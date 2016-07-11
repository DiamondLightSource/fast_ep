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
import traceback
import string
import json
from os.path import basename, splitext
from multiprocessing import Pool

from iotbx import mtz
from iotbx import pdb
from libtbx.phil import parse
from cctbx.sgtbx import space_group, space_group_symbols
from iotbx.scalepack import merge as merge_scalepack
from libtbx import introspection

if not 'FAST_EP_ROOT' in os.environ:
    raise RuntimeError, 'FAST_EP_ROOT not set'

fast_ep_lib = os.path.join(os.environ['FAST_EP_ROOT'], 'lib')

if not fast_ep_lib in sys.path:
    sys.path.append(fast_ep_lib)

from xml_output import write_ispyb_xml, xmlfile2json

from generate_possible_spacegroups import generate_chiral_spacegroups_unique, \
     spacegroup_enantiomorph, spacegroup_full, sanitize_spacegroup
from number_sites_estimate import number_sites_estimate, \
     number_residues_estimate
from guess_the_atom import guess_the_atom
from run_job import run_job, run_job_cluster, is_cluster_job_finished
from fast_ep_shelxd import run_shelxd_drmaa_array, run_shelxd_local, analyse_res, \
     happy_shelxd_log
from fast_ep_shelxe import run_shelxe_drmaa_array, run_shelxe_local

def useful_number_sites(_cell, _pointgroup):
    nha = number_sites_estimate(_cell, _pointgroup)

    result = []

    for f in [0.25, 0.5, 1.0, 2.0, 4.0]:
        nha_test = int(round(f * nha))
        if nha_test and not nha_test in result:
            result.append(nha_test)

    if not result:
        result = [1]

    return result

def modify_ins_text(_ins_text, _spacegroup, _nsites, _rlimit):
    '''Update the text in a SHELXD .ins file to handle the correct number
    of sites and spacegroup symmetry operations.'''

    new_text = []

    symm = [op.as_xyz().upper() for op in
            space_group(space_group_symbols(_spacegroup).hall()).smx()]

    for record in _ins_text:
        if 'SYMM' in record:
            if not symm:
                continue
            for op in symm:
                if op == 'X,Y,Z':
                    continue
                new_text.append(('SYMM %s' % op))
            symm = None
        elif 'FIND' in record:
            new_text.append(('FIND %d' % _nsites))
        elif 'SHEL' in record:
            new_text.append(('SHEL 999 %.1f' % _rlimit))
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

    def __call__(self, _line):
        sys.stdout.write('%s\n' % _line)
        self._fout.write('%s\n' % _line)
        return

class Fast_ep_parameters:
    '''A class to wrap up the parameters for fast_ep e.g. the number of machines
    to use, the number of cpus on each machine, the input reflection file.'''

    def __init__(self):
        self._phil = parse("""
fast_ep {
  machines = 1
    .type = int
  cpu = %d
    .type = int
  sad = 'fast_dp.mtz'
    .type = str
  native = None
    .type = str
  atom = Se
    .type = str
  spg = None
    .type = str
  nsites = None
    .type = int
  ntry = 200
    .type = int
  xml = ''
    .type = str
  plot = False
    .type = bool
  trace = False
    .type = bool
}
""" % introspection.number_of_processors(return_value_if_unknown = 1))
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

    def get_sad(self):
        return self._parameters.fast_ep.sad

    def get_native(self):
        return self._parameters.fast_ep.native

    def get_atom(self):
        return self._parameters.fast_ep.atom

    def get_spg(self):
        return self._parameters.fast_ep.spg

    def get_nsites(self):
        return self._parameters.fast_ep.nsites

    def get_ntry(self):
        return self._parameters.fast_ep.ntry

    def get_xml(self):
        return self._parameters.fast_ep.xml

    def get_plot(self):
        return self._parameters.fast_ep.plot

    def get_trace(self):
        return self._parameters.fast_ep.trace

class Fast_ep:
    '''A class to run shelxc / d / e to very quickly establish (i) whether
    experimental phasing is likely to be successful and (ii) what the
    correct parameteters and number of heavy atom sites.'''

    def __init__(self, _parameters):
        '''Instantiate class and perform initial processing needed before the
        real work is done. This includes assessment of the anomalous signal,
        the signal to noise and conversion to pseudo-scalepack format of the
        data for input into shelxc, which is used to compute FA values.'''

        self._hklin = _parameters.get_sad()
        self._native_hklin = _parameters.get_native()
        self._cpu = _parameters.get_cpu()
        self._machines = _parameters.get_machines()
        self._atom = _parameters.get_atom()
        self._ntry = _parameters.get_ntry()
        self._plot = _parameters.get_plot()
        self._trace = _parameters.get_trace()
        self._data = None

        if self._machines == 1:
            self._cluster = False
        else:
            self._cluster = True

        self._wd = os.getcwd()
        self._log = logger()

        self._log('Using %d cpus / %d machines' % (self._cpu, self._machines))

        self._full_command_line = ' '.join(sys.argv)

        # pull information we'll need from the input MTZ file - the unit cell,
        # the pointgroup and the number of reflections in the file. select
        # first Miller array in file which has anomalous data

        # --- SAD DATA ---

        m = mtz.object(self._hklin)
        mas = m.as_miller_arrays()

        for ma in mas:
            if not ma.anomalous_flag():
                continue
            if str(ma.observation_type()) != 'xray.intensity':
                continue
            self._data = ma
            break

        if not self._data:
            raise RuntimeError, 'no anomalous intensity data found in %s' % \
                self._hklin

        # --- NATIVE DATA ---

        if self._native_hklin:
            mn = mtz.object(self._native_hklin)
            masn = m.as_miller_arrays()

            for ma in masn:
                if ma.anomalous_flag():
                    continue
                if str(ma.observation_type()) != 'xray.intensity':
                    continue
                self._native = ma
                break

        else:
            self._native = None

        if not self._data:
            raise RuntimeError, 'no intensity data found in %s' % self._hklin

        self._pointgroup = self._data.space_group().type().number()
        self._unit_cell = self._data.unit_cell().parameters()

        self._nrefl = m.n_reflections()

        # write out a nice summary of the data set properties and what columns
        # were selected for analysis

        self._log('Input:       %s' % self._hklin)
        self._log('Columns:     %s' % self._data.info().label_string())
        if self._native_hklin:
            self._log('Native:      %s' % self._native_hklin)
            self._log('Columns:     %s' % self._native.info().label_string())
        self._log('Unit cell:   %.2f %.2f %.2f %.2f %.2f %.2f' % \
                  self._unit_cell)
        self._log('N try:       %d' % self._ntry)
        self._log('Pointgroup:  %s' % m.space_group().type().lookup_symbol())
        self._log('Resolution:  %.2f - %.2f' % self._data.resolution_range())

        self._xml_name = _parameters.get_xml()
        self._xml_results= {}
        self._xml_results['LOWRES'] = self._data.resolution_range()[0]
        self._xml_results['HIGHRES'] = self._data.resolution_range()[1]

        self._log('Nrefl:       %d / %d' % (self._nrefl,
                                            self._data.n_bijvoet_pairs()))
        self._log('DF/F:        %.3f' % self._data.anomalous_signal())

        differences = self._data.anomalous_differences()

        self._log('dI/sig(dI):  %.3f' % (sum(abs(differences.data())) /
                                         sum(differences.sigmas())))

        # Now set up the job - run shelxc, assess anomalous signal, compute
        # possible spacegroup options, generate scalepack format reflection
        # file etc.

        self._data = self._data.apply_scaling(target_max = 1.0e5)

        merge_scalepack.write(file_name = 'sad.sca',
                              miller_array = self._data)

        if self._native_hklin:
            merge_scalepack.write(file_name = 'native.sca',
                                  miller_array = self._native)

        # in here run shelxc to generate the ins file (which will need to be
        # modified) and the hkl files, which will need to be copied.

        if _parameters.get_spg():
            self._spacegroups = [_parameters.get_spg(),]
        else:
            self._spacegroups = generate_chiral_spacegroups_unique(self._pointgroup)

        self._log('Spacegroups: %s' % ' '.join(self._spacegroups))

        if _parameters.get_nsites():
            self._nsites = [_parameters.get_nsites(),]
        else:
            self._nsites = useful_number_sites(self._unit_cell, self._pointgroup)

        spacegroup = self._spacegroups[0]
        nsite = self._nsites[0]
        ntry = self._ntry

        self._xml_results['SHELXC_SPACEGROUP_ID'] = space_group_symbols(
            spacegroup).number()

        if self._native_hklin:
            shelxc_output = run_job(
                'shelxc', ['sad'],
                ['sad sad.sca',
                 'nat native.sca',
                'cell %.3f %.3f %.3f %.3f %.3f %.3f' % self._unit_cell,
                'spag %s' % sanitize_spacegroup(spacegroup),
                'sfac %s' % string.upper(_parameters.get_atom()),
                'find %d' % nsite,
                'mind -3.5',
                'ntry %d' % ntry])

        else:
            shelxc_output = run_job(
                'shelxc', ['sad'],
                ['sad sad.sca',
                 'cell %.3f %.3f %.3f %.3f %.3f %.3f' % self._unit_cell,
                 'spag %s' % sanitize_spacegroup(spacegroup),
                 'sfac %s' % string.upper(_parameters.get_atom()),
                 'find %d' % nsite,
                 'mind -3.5',
                 'ntry %d' % ntry])

        # FIXME in here perform some analysis of the shelxc output - how much
        # anomalous signal was reported?

        open('shelxc.log', 'w').write(''.join(shelxc_output))

        table = { }

        for record in shelxc_output:
            if record.strip().startswith('Resl.'):
                resolutions = map(float, record.replace(' - ', ' ').split()[2:])
                table['dmin'] = resolutions
            if record.strip().startswith('<I/sig>'):
                table['isig'] = map(float, record.split()[1:])
            if record.strip().startswith('%Complete'):
                table['comp'] = map(float, record.split()[1:])
            if record.strip().startswith('<d"/sig>'):
                table['dsig'] = map(float, record.split()[1:])

        for row in ['isig', 'comp', 'dsig']:
            pad = len(table['dmin']) - len(table[row])
            if pad > 0:
                table[row] += [float('nan')] * pad

        shells = len(table['dmin'])

        self._log('SHELXC summary:')
        self._log('Dmin  <I/sig>  %comp  <d"/sig>')
        for j in range(shells):
            self._log('%5.2f  %6.2f  %6.2f  %5.2f' %
                      (table['dmin'][j], table['isig'][j],
                       table['comp'][j], table['dsig'][j]))

        dmin = self._data.resolution_range()[1]
        self._ano_rlimits =  [dmin, dmin + 0.25, dmin + 0.5]

        self._log('Anomalous limits: %s' %  ' '.join(["%.1f" % v for v in self._ano_rlimits]))

        # store the ins file text - will need to modify this when we come to
        # run shelxd...

        self._ins_text = open('sad_fa.ins', 'r').readlines()

        return

    def find_sites(self):
        '''Actually perform the substructure calculation, using many runs of
        shelxd_mp (one per spacegroup per nsites to test). This requires
        copying into place the input files generated by shelxc, making
        modifications to set the spacegroup (as symmetry operations) and the
        number of sites to look for.'''

        t0 = time.time()

        cluster = self._cluster
        njobs = self._machines
        ncpu = self._cpu

        # set up N x M shelxd jobs

        jobs = [ ]

        # the shelx programs are fortran so we need to tell them how much space
        # to allocate on the command-line - this is done by passing -LN on the
        # command line where N is calculated as follows:

        nrefl = max(10, 1 + 2 * int(1 + math.floor(self._nrefl / 100000.0)))

        # modify the instruction file (.ins) for the number of sites and
        # symmetry operations for each run

        for spacegroup in self._spacegroups:
            for nsite in self._nsites:
                for rlimit in self._ano_rlimits:
                    wd = os.path.join(self._wd, spacegroup.replace(':', '-'), str(nsite), "%.1f" % rlimit)
                    if not os.path.exists(wd):
                        os.makedirs(wd)

                    new_text = modify_ins_text(self._ins_text, spacegroup, nsite, rlimit)

                    shutil.copyfile(os.path.join(self._wd, 'sad_fa.hkl'),
                                    os.path.join(wd, 'sad_fa.hkl'))

                    open(os.path.join(wd, 'sad_fa.ins'), 'w').write(
                        '\n'.join(new_text))

                    jobs.append({'nrefl':nrefl, 'ncpu':ncpu, 'wd':wd})

        # actually execute the tasks - either locally or on a cluster, allowing
        # for potential for fewer available machines than jobs

        self._log('Running %d x shelxd_mp jobs' % len(jobs))

        if cluster:
            run_shelxd_drmaa_array(self._wd, nrefl, ncpu, njobs, jobs)
        else:
            pool = Pool(min(njobs, len(jobs)))
            pool.map(run_shelxd_local, jobs)

        # now gather up all of the results, find the one with best cfom

        best_cfom = 0.0
        best_cc = 0.0
        best_ccweak = 0.0
        best_spacegroup = None
        best_nsite = 0
        best_nsite_real = 0
        best_ano_rlimit = 0

        results = { }

        for spacegroup in self._spacegroups:
            for nsite in self._nsites:
                for rlimit in self._ano_rlimits:
                    wd = os.path.join(self._wd, spacegroup.replace(':', '-'), str(nsite), "%.1f" % rlimit)

                    shelxd_log = os.path.join(wd, 'sad_fa.lst')

                    if self._plot:
                        from fast_ep_helpers import plot_shelxd_cc
                        shelxd_plot = os.path.join(wd, 'sad_fa.png')

                        plot_shelxd_cc(shelxd_log, shelxd_plot, spacegroup, nsite, rlimit)

                    if happy_shelxd_log(shelxd_log):
                        res = open(os.path.join(wd, 'sad_fa.res')).readlines()
                        cc, cc_weak, cfom, nsite_real = analyse_res(res)

                        results[(spacegroup, nsite, rlimit)] = (cc, cc_weak, cfom,
                                                        nsite_real)

                        if cfom > best_cfom and nsite_real > 0:
                            best_cfom = cfom
                            best_cc = cc
                            best_ccweak = cc_weak
                            best_spacegroup = spacegroup
                            best_nsite = nsite
                            best_nsite_real = nsite_real
                            best_ano_rlimit = rlimit

                    else:
                        results[(spacegroup, nsite, rlimit)] = (0.0, 0.0, 0.0, 0)

        if not best_spacegroup:
            raise RuntimeError, 'All shelxd jobs failed'

        for spacegroup in self._spacegroups:
            if spacegroup == best_spacegroup:
                self._log('Spacegroup: %s (best)' % spacegroup)
                self._xml_results['SPACEGROUP'] = space_group_symbols(
                    spacegroup).number()
            else:
                self._log('Spacegroup: %s' % spacegroup)

            self._log('No.  Res.  CCall  CCweak CFOM  No. found')
            for nsite in self._nsites:
                for rlimit in self._ano_rlimits:
                    (cc, cc_weak, cfom, nsite_real) = results[(spacegroup, nsite, rlimit)]
                    if (spacegroup, nsite, rlimit) == (best_spacegroup, best_nsite, best_ano_rlimit):
                        self._log('%3d  %.1f  %6.2f %6.2f %6.2f %3d (best)' %
                            (nsite, rlimit, cc, cc_weak, cfom, nsite_real))
                    else:
                        self._log('%3d  %.1f  %6.2f %6.2f %6.2f %3d' %
                            (nsite, rlimit, cc, cc_weak, cfom, nsite_real))

        t1 = time.time()
        self._log('Time: %.2f' % (t1 - t0))

        self._log('Best spacegroup: %s' % best_spacegroup)
        self._log('Best nsites:     %d' % best_nsite_real)
        self._log('Best resolution: %.1f A' % best_ano_rlimit)
        self._log('Best CC / weak:  %.2f / %.2f' % (best_cc, best_ccweak))

        self._best_spacegroup = best_spacegroup
        self._best_nsite = best_nsite_real
        self._best_ano_rlimit = best_ano_rlimit
        self._best_cfom = best_cfom
        self._best_cc = best_cc
        self._best_ccweak = best_ccweak

        # copy back result files

        best = os.path.join(self._wd, best_spacegroup.replace(':', '-'), str(best_nsite), "%.1f" % best_ano_rlimit)

        endings = ['lst', 'pdb', 'res']
        if self._plot:
            endings.append('png')

        for ending in endings:
            shutil.copyfile(os.path.join(best, 'sad_fa.%s' % ending),
                            os.path.join(self._wd, 'sad_fa.%s' % ending))

        return

    def phase(self):
        '''Perform the phasing following from the substructure determination,
        using the best solution found, using shelxe. This will be run for a
        range of sensible solvent fractions between 25% and 75% and for
        both hands of the substructure. N.B. for a chiral screw axis (e.g. P41)
        this will also invert the spacegroup (to e.g. P43) which needs to
        be remembered in transforming the output.'''

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


        if cluster:
            run_shelxe_drmaa_array(self._wd, njobs, jobs)
        else:
            pool = Pool(min(njobs * ncpu, len(jobs)))
            pool.map(run_shelxe_local, jobs)

        results = { }

        best_fom = 0.0
        best_solvent = None
        best_hand = None

        for solvent_fraction in solvent_fractions:
            wd = os.path.join(self._wd, '%.2f' % solvent_fraction)

            if self._plot:
                from fast_ep_helpers import plot_shelxe_contrast
                plot_shelxe_contrast(os.path.join(wd, 'sad.lst'),
                                     os.path.join(wd, 'sad_i.lst'),
                                     os.path.join(wd, 'sad.png'),
                                     solvent_fraction)

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

        self._best_fom = best_fom
        self._best_solvent = best_solvent
        self._best_hand = best_hand

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
            filename_to_read = 'sad.lst'
        elif best_hand == 'inverted':
            filename_to_read = 'sad_i.lst'
        else:
            raise RuntimeError, 'unknown hand'
        file_to_read = os.path.join(wd, filename_to_read)

        self.get_fom_mapCC(file_to_read)

        self._xml_results['FOM'] = best_fom
        self._xml_results['SOLVENTCONTENT'] = best_solvent
        self._xml_results['ENANTIOMORPH'] = (best_hand=='inverted')

        # copy the result files from the most successful shelxe run into the
        # working directory, before converting to mtz format for inspection with
        # e.g. coot.

        # FIXME in here map correct site file to ASU

        from fast_ep_helpers import map_sites_to_asu
        if best_hand == 'original':
            map_sites_to_asu(self._best_spacegroup,
                             os.path.join(wd, 'sad_fa.pdb'),
                             os.path.join(self._wd, 'sites.pdb'))
        else:
            map_sites_to_asu(self._best_spacegroup,
                             os.path.join(wd, 'sad_fa.pdb'),
                             os.path.join(self._wd, 'sites.pdb'),
                             invert=True)

        if best_hand == 'original':
            for ending in ['phs', 'pha', 'lst', 'hat']:
                shutil.copyfile(os.path.join(wd, 'sad.%s' % ending),
                                os.path.join(self._wd, 'sad.%s' % ending))
        else:
            for ending in ['phs', 'pha', 'lst', 'hat']:
                shutil.copyfile(os.path.join(wd, 'sad_i.%s' % ending),
                                os.path.join(self._wd, 'sad.%s' % ending))
            self._best_spacegroup = spacegroup_enantiomorph(
                self._best_spacegroup)

        if self._plot:
            for ending in ['png']:
                shutil.copyfile(os.path.join(wd, 'sad.%s' % ending),
                                os.path.join(self._wd, 'sad.%s' % ending))

        self._log('Best spacegroup: %s' % self._best_spacegroup)

        if self._trace:
            # rerun shelxe to trace the chain
            self._nres_trace = 0
            arguments = ['sad', 'sad_fa', '-h%d' % self._best_nsite,
                         '-s%.2f' % best_solvent, '-a3', '-m20']
            if not best_hand == 'original':
                arguments.append('-i')
            output = run_job('shelxe', arguments, [], self._wd)
            for record in output:
                if 'residues left after pruning' in record:
                    self._nres_trace = int(record.split()[0])
            self._log('Traced:       %d' % self._nres_trace)
            shutil.copyfile(os.path.join(self._wd, 'sad.pdb'),
                            os.path.join(self._wd, 'sad_trace.pdb'))

        # convert sites to pdb, inverting if needed

        xs = pdb.input(os.path.join(
            self._wd, 'sad_fa.pdb')).xray_structure_simple()
        if best_hand == 'inverted':
            open('sad.pdb', 'w').write(xs.change_hand().as_pdb_file())
        else:
            open('sad.pdb', 'w').write(xs.as_pdb_file())

        o = run_job('convert2mtz', ['-hklin', 'sad.phs', '-mtzout', 'sad.mtz',
                                   '-colin', 'F FOM PHI SIGF',
                                   '-cell', '%f %f %f %f %f %f' % self._unit_cell,
                                   '-spacegroup',
                                   spacegroup_full(self._best_spacegroup)],
            [], self._wd)

        open('convert2mtz.log', 'w').write('\n'.join(o))

        t1 = time.time()
        self._log('Time: %.2f' % (t1 - t0))

        return

    def write_results(self):
        '''Write a little data file which can be used for subsequent analysis
        with other phasing programs, based on what we have learned in the
        analysis above.'''

        open(os.path.join(self._wd, 'fast_ep.dat'), 'w').write('\n'.join([
            'SPACEGROUP: %s' % self._best_spacegroup,
            'NSITE: %d' % self._best_nsite,
            'SOLVENT: %.2f' % self._best_solvent, '']))

        json_dict = {'_hklin'          : self._hklin,
                     '_cpu'            : self._cpu,
                     '_machines'       : self._machines,
                     '_use_cluster'    : True,
                     '_spacegroup'     : [self._best_spacegroup],
                     '_nsites'         : [self._best_nsite],
                     'nsite_real'      : self._best_nsite,
                     'solv'            : self._best_solvent,
                     'cc'              : self._best_cc,
                     'cc_weak'         : self._best_ccweak,
                     'cfom'            : self._best_cfom,
                     'fom'             : self._best_fom,
                     'inverted'        : self._best_hand=='inverted',
                     '_trace'          : False,
                     '_xml_name'       : None
                    }

        json_data = json.dumps(json_dict, indent=4, separators=(',', ':'))
        with open(os.path.join(self._wd, 'fast_ep_data.json'), 'w') as json_file:
            json_file.write(json_data)

        return

    def get_fom_mapCC(self, file_to_read):
        for record in open(file_to_read):
            check_string = 'd    inf'
            #special due to initial case and separated resolutions
            if check_string in record:
                resolution_ranges = record.split(check_string)[1].split(' - ')
                for resolution_number in xrange(0, len(resolution_ranges) - 1):
                    resolution_number_name = str(resolution_number).zfill(2)
                    k = 'RESOLUTION_LOW' + resolution_number_name
                    if resolution_number == 0:
                        self._xml_results[k] = self._data.resolution_range()[0]
                    else:
                        self._xml_results[k] = float(resolution_ranges[
                            resolution_number])
                    k = 'RESOLUTION_HIGH' + resolution_number_name
                    self._xml_results[k] = float(resolution_ranges[
                        resolution_number + 1])

            parse_pairs = [['<FOM>', 'FOM'], ['<mapCC>', 'MAPCC'],
                           ['N   ', 'NREFLECTIONS']]
            for check_string, field_name in parse_pairs:
                self._find_line_parse_string(check_string, record, field_name)

    # look for a line starting with check_string, then split the rest of the
    # string into values and store them

    def _find_line_parse_string(self, check_string, record, field_name):
        if check_string in record:
            # remove the check_string, then split the remainder
            field_values = record.split(check_string)[1].split()
            for field_number in xrange(0, len(field_values)):
                field_number_name = str(field_number).zfill(2)
                k = field_name + field_number_name
                self._xml_results[k] = field_values[field_number]

    def write_xml(self):
        if self._xml_name == '':
            return
        filename = os.path.join(self._wd, self._xml_name)
        write_ispyb_xml(filename, self._full_command_line, self._wd,
                        self._xml_results)

        try:
            json_file_name = '.'.join([splitext(basename(self._xml_name))[0],
                                       'json'])
            with open(os.path.join(
                self._wd, json_file_name), 'w') as json_file:
                json_data = xmlfile2json(filename)
                json_file.write(json_data)
        except ImportError, e:
            pass

if __name__ == '__main__':
    fast_ep = Fast_ep(Fast_ep_parameters())
    try:
        fast_ep.find_sites()
    except RuntimeError, e:
        fast_ep._log('*** FIND: %s ***' % str(e))
        traceback.print_exc(file = open('fast_ep.error', 'w'))
        sys.exit(1)

    try:
        fast_ep.phase()
    except RuntimeError, e:
        fast_ep._log('*** PHASE %s ***' % str(e))
        traceback.print_exc(file = open('fast_ep.error', 'w'))
        sys.exit(1)

    try:
        fast_ep.write_results()
    except RuntimeError, e:
        fast_ep._log('*** WRITE_RESULTS %s ***' % str(e))
        traceback.print_exc(file = open('fast_ep.error', 'w'))
        sys.exit(1)

    try:
        fast_ep.write_xml()
    except RuntimeError, e:
        fast_ep._log('*** FINISH %s ***' % str(e))
        traceback.print_exc(file = open('fast_ep.error', 'w'))
        sys.exit(1)
