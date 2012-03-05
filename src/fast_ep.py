#!/usr/bin/env python
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# fast_ep - main program.

import os
import sys
import time
import shutil
import math
from multiprocessing import Pool

from iotbx import mtz
from libtbx.phil import parse
from cctbx.sgtbx import space_group, space_group_symbols

if not 'FAST_EP_ROOT' in os.environ:
    raise RuntimeError, 'FAST_EP_ROOT not set'

fast_ep_lib = os.path.join(os.environ['FAST_EP_ROOT'], 'lib')

if not fast_ep_lib in sys.path:
    sys.path.append(fast_ep_lib)

from generate_possible_spacegroups import generate_chiral_spacegroups_unique, \
     spacegroup_enantiomorph, spacegroup_full, sanitize_spacegroup
from number_sites_estimate import number_sites_estimate, \
     number_residues_estimate
from guess_the_atom import guess_the_atom
from run_job import run_job, run_job_cluster, is_cluster_job_finished
from fast_ep_helpers import autosharp
from fast_ep_shelxd import run_shelxd_cluster, run_shelxd_local, analyse_res
from fast_ep_shelxe import run_shelxe_cluster, run_shelxe_local

def useful_number_sites(cell, pointgroup):
    nha = number_sites_estimate(cell, pointgroup)

    result = []

    for f in [0.25, 0.5, 1.0, 2.0, 4.0]:
        nha_test = int(round(f * nha))
        if nha_test and not nha_test in result:
            result.append(nha_test)

    return result

def modify_ins_text(ins_text, spacegroup, nsites):
    '''Update the text in a SHELXD .ins file to handle the correct number
    of sites and spacegroup symmetry operations.'''

    new_text = []

    symm = [op.as_xyz().upper() for op in
            space_group(space_group_symbols(spacegroup).hall()).smx()]

    for record in ins_text:
        if 'SYMM' in record:
            if not symm:
                continue
            for op in symm:
                if op == 'X,Y,Z':
                    continue
                new_text.append(('SYMM %s' % op))
            symm = None
        elif 'FIND' in record:
            new_text.append(('FIND %d' % nsites))
        else:
            new_text.append(record.strip())

    return new_text

class logger:
    def __init__(self):
        self._fout = open('fast_ep.log', 'w')
        return

    def __del__(self):
        self._fout.close()
        self._cout = None
        return

    def __call__(self, line):
        print line
        self._fout.write('%s\n' % line)
        return

class Fast_ep_parameters:
    '''A class to wrap up the parameters for fast_ep e.g. the number of machines
    to use, the number of cpus on each machine, the input reflection file.'''

    def __init__(self):
        self._phil = parse("""
fast_ep {
  machines = 1
    .type = int
  cpu = 8
    .type = int
  input = 'fast_dp.mtz'
    .type = str
}
""")
        argument_interpreter = self._phil.command_line_argument_interpreter(
            home_scope = 'fast_ep')
        for argv in sys.argv[1:]:
            command_line_phil = argument_interpreter.process(arg = argv)
            self._phil = self._phil.fetch(command_line_phil)
        self._parameters = self._phil.extract()
        
        return

    def get_machines(self):
        return self._parameters.fast_ep.machines

    def get_cpu(self):
        return self._parameters.fast_ep.cpu

    def get_input(self):
        return self._parameters.fast_ep.input

class Fast_ep:
    '''A class to run shelxc / d / e to very quickly establish (i) whether
    experimental phasing is likely to be successful and (ii) what the
    correct parameteters and number of heavy atom sites.'''

    def __init__(self, parameters):
        '''Instantiate class and perform initial processing needed before the
        real work is done.'''
        
        self._hklin = parameters.get_input()
        self._cpu = parameters.get_cpu()
        self._machines = parameters.get_machines()

        if self._machines == 1:
            self._cluster = False
        else:
            self._cluster = True

        self._wd = os.getcwd()
        self._log = logger()

        self._log('Using %d cpus / %d machines' % (self._cpu, self._machines))

        # pull information we'll need from the input MTZ file - the unit cell,
        # the pointgroup and the number of reflections in the file

        m = mtz.object(self._hklin)

        self._pointgroup = m.space_group().type().number()
        for crystal in m.crystals():
            if crystal.name() != 'HKL_base':
                self._unit_cell = crystal.unit_cell().parameters()

        self._nrefl = m.n_reflections()

        self._log('Input:     %s' % self._hklin)
        self._log('Unit cell: %.2f %.2f %.2f %.2f %.2f %.2f' % \
                  self._unit_cell)
        self._log('Nrefl:     %d' % self._nrefl)

        # Now set up the job - run shelxc, assess anomalous signal, compute
        # possible spacegroup options, generate scalepack format reflection
        # file etc.

        # FIXME check status of output here

        run_job('mtz2sca', [self._hklin, 'sad.sca'])

        # in here run shelxc to generate the ins file (which will need to be
        # modified) and the hkl files, which will need to be copied.

        self._spacegroups = generate_chiral_spacegroups_unique(self._pointgroup)
        self._nsites = useful_number_sites(self._unit_cell, self._pointgroup)

        spacegroup = self._spacegroups[0]
        nsite = self._nsites[0]

        shelxc_output = run_job(
            'shelxc', ['sad'],
            ['sad sad.sca',
             'cell %.3f %.3f %.3f %.3f %.3f %.3f' % self._unit_cell,
             'spag %s' % sanitize_spacegroup(spacegroup),
             'find %d' % nsite,
             'mind -3.5',
             'ntry 200'])

        # FIXME in here perform some analysis of the shelxc output - how much
        # anomalous signal was reported?
 
        open('shelxc.log', 'w').write(''.join(shelxc_output))

        # store the ins file text - will need to modify this when we come to
        # run shelxd...
        
        self._ins_text = open('sad_fa.ins', 'r').readlines()

        return

    def find_sites(self):
        '''Actually perform the substructure calculation.'''

        t0 = time.time()

        cluster = self._cluster
        njobs = self._machines
        ncpu = self._cpu

        # set up N x M shelxd jobs

        jobs = [ ]

        nrefl = 1 + int(math.floor(self._cpu * self._nrefl / 100000.0))

        for spacegroup in self._spacegroups:
            for nsite in self._nsites:
                wd = os.path.join(self._wd, spacegroup, str(nsite))
                if not os.path.exists(wd):
                    os.makedirs(wd)

                new_text = modify_ins_text(self._ins_text, spacegroup, nsite)

                shutil.copyfile(os.path.join(self._wd, 'sad_fa.hkl'),
                                os.path.join(wd, 'sad_fa.hkl'))

                open(os.path.join(wd, 'sad_fa.ins'), 'w').write(
                    '\n'.join(new_text))

                jobs.append({'nrefl':nrefl, 'ncpu':ncpu, 'wd':wd})

        # actually execute the tasks - either locally or on a cluster, allowing
        # for potential for fewer available machines than jobs

        self._log('Running %d x shelxd_mp jobs' % len(jobs))
        
        pool = Pool(min(njobs, len(jobs)))

        if cluster:
            pool.map(run_shelxd_cluster, jobs)
        else:
            pool.map(run_shelxd_local, jobs)

        # now gather up all of the results, find the one with best cfom

        best_cfom = 0.0

        results = { }

        for spacegroup in self._spacegroups:
            for nsite in self._nsites:
                wd = os.path.join(self._wd, spacegroup, str(nsite))
                res = open(os.path.join(wd, 'sad_fa.res')).readlines()

                cc, cc_weak, cfom, nsite_real = analyse_res(res)

                results[(spacegroup, nsite)] = (cc, cc_weak, cfom, nsite_real)

                if cfom > best_cfom:
                    best_cfom = cfom
                    best_spacegroup = spacegroup
                    best_nsite = nsite
                    best_nsite_real = nsite_real

        for spacegroup in self._spacegroups:
            if spacegroup == best_spacegroup:
                self._log('Spacegroup: %s (best)' % spacegroup)
            else:
                self._log('Spacegroup: %s' % spacegroup)
                
            for nsite in self._nsites:
                (cc, cc_weak, cfom, nsite_real) = results[(spacegroup, nsite)]
                if (spacegroup, nsite) == (best_spacegroup, best_nsite):
                    self._log('%3d %6.2f %6.2f %6.2f %3d (best)' %
                              (nsite, cc, cc_weak, cfom, nsite_real))
                else:
                    self._log('%3d %6.2f %6.2f %6.2f %3d' %
                              (nsite, cc, cc_weak, cfom, nsite_real))

        t1 = time.time()
        self._log('Time: %.2f' % (t1 - t0))

        self._log('Best spacegroup: %s' % best_spacegroup)
        self._log('Best nsites:     %d' % best_nsite_real)

        self._best_spacegroup = best_spacegroup
        self._best_nsite = best_nsite_real
        
        # copy back result files

        best = os.path.join(self._wd, best_spacegroup, str(best_nsite))

        for ending in 'lst', 'pdb', 'res':
            shutil.copyfile(os.path.join(best, 'sad_fa.%s' % ending),
            os.path.join(self._wd, 'sad_fa.%s' % ending))

        return

    def phase(self):
        '''Perform the phasing following from the substructure determination,
        using the best solution found. This will test a range of sensible
        solvent fractions'''

        t0 = time.time()
        
        cluster = self._cluster
        njobs = self._machines
        ncpu = self._cpu

        solvent_fractions = [0.25 + 0.05 * j for j in range(11)]

        jobs = [ ]

        for solvent_fraction in solvent_fractions:
            wd = os.path.join(self._wd, '%.2f' % solvent_fraction)
            if not os.path.exists(wd):
                os.makedirs(wd)
            shutil.copyfile(os.path.join(self._wd, 'sad.hkl'),
                            os.path.join(wd, 'sad.hkl'))
            for ending in 'lst', 'pdb', 'res', 'hkl':
                shutil.copyfile(os.path.join(self._wd, 'sad_fa.%s' % ending),
                                os.path.join(wd, 'sad_fa.%s' % ending))

            jobs.append({'nsite':self._best_nsite, 'solv':solvent_fraction,
                         'hand':'original', 'wd':wd})
            jobs.append({'nsite':self._best_nsite, 'solv':solvent_fraction,
                         'hand':'inverted', 'wd':wd})

        self._log('Running %d x shelxe jobs' % len(jobs))

        pool = Pool(min(njobs * ncpu, len(jobs)))

        if cluster:
            pool.map(run_shelxe_cluster, jobs)
        else:
            pool.map(run_shelxe_local, jobs)

        results = { }

        best_fom = 0.0
        best_solvent = None
        best_hand = None

        for solvent_fraction in solvent_fractions:
            wd = os.path.join(self._wd, '%.2f' % solvent_fraction)
            for record in open(os.path.join(wd, 'sad.lst')):
                if 'Estimated mean FOM =' in record:
                    fom_orig = float(record.split()[4])
            for record in open(os.path.join(wd, 'sad_i.lst')):
                if 'Estimated mean FOM =' in record:
                    fom_inv = float(record.split()[4])
            results[solvent_fraction] = (fom_orig, fom_inv)

            if fom_orig > best_fom:
                best_fom = fom_orig
                best_solvent = solvent_fraction
                best_hand = 'original'

            if fom_inv > best_fom:
                best_fom = fom_inv
                best_solvent = solvent_fraction
                best_hand = 'inverted'

        self._log('Solv. Orig. Inv.')
        for solvent_fraction in solvent_fractions:
            fom_orig, fom_inv = results[solvent_fraction]
            if solvent_fraction == best_solvent:
                self._log(
                    '%.2f %.3f %.3f (best)' % (solvent_fraction, fom_orig,
                                               fom_inv))
            else:
                self._log('%.2f %.3f %.3f' % (solvent_fraction, fom_orig,
                                              fom_inv))

        self._log('Best solvent: %.2f' % best_solvent)
        self._log('Best hand:    %s' % best_hand)

        wd = os.path.join(self._wd, '%.2f' % best_solvent)

        if best_hand == 'original':
            for ending in 'phs', 'pha', 'lst', 'hat':
                shutil.copyfile(os.path.join(wd, 'sad.%s' % ending),
                                os.path.join(self._wd, 'sad.%s' % ending))
        else:
            for ending in 'phs', 'pha', 'lst', 'hat':
                shutil.copyfile(os.path.join(wd, 'sad_i.%s' % ending),
                                os.path.join(self._wd, 'sad.%s' % ending))
            self._best_spacegroup = spacegroup_enantiomorph(
                self._best_spacegroup)

        self._log('Best spacegroup: %s' % self._best_spacegroup)
                
        run_job('convert2mtz', ['-hklin', 'sad.phs', '-mtzout', 'sad.mtz',
                                '-colin', 'F FOM PHI SIGF',
                                '-cell', '%f %f %f %f %f %f' % self._unit_cell,
                                '-spacegroup',
                                spacegroup_full(self._best_spacegroup)],
            [], self._wd)

        t1 = time.time()
        self._log('Time: %.2f' % (t1 - t0))

        return
    
if __name__ == '__main__':
    fast_ep = Fast_ep(Fast_ep_parameters())
    fast_ep.find_sites()
    fast_ep.phase()
