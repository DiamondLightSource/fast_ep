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
import json
from os.path import basename, splitext
from multiprocessing import Pool
from operator import and_
from itertools import product
from math import isnan

from iotbx import pdb
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.phil import parse
from cctbx.sgtbx import space_group_symbols
from cctbx.xray import observation_types
from iotbx.scalepack import merge as merge_scalepack
from libtbx import introspection

if 'FAST_EP_ROOT' in os.environ:
    sys.path.append(os.environ['FAST_EP_ROOT'])

from lib.report import render_html_report
from lib.xml_output import write_ispyb_xml, xmlfile2json, store_string_xml
from lib.generate_possible_spacegroups import generate_chiral_spacegroups_unique, \
     spacegroup_enantiomorph, spacegroup_full, sanitize_spacegroup
from lib.run_job import run_job
from src.fast_ep_helpers import modify_ins_text, useful_number_sites
from src.fast_ep_shelxd import get_shelxd_results, get_average_ranks, \
     log_rank_table, log_shelxd_results, run_shelxd_drmaa_array,\
    run_shelxd_local, get_shelxd_result_ranks, log_shelxd_results_advanced
from src.fast_ep_shelxe import run_shelxe_drmaa_array, run_shelxe_local,\
    read_shelxe_log
from src.fast_ep_plots import plot_shelxd_cc, plot_shelxe_contrast,\
    hist_shelxd_cc, plot_shelxe_fom_mapcc, plot_shelxe_mean_fom_cc,\
    plot_anom_shelxc
from datetime import datetime

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
  data = 'fast_dp.mtz'
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
  ncycle = 20
    .type = int
  rlims = None
    .type = floats(value_min=0)
  xml = ''
    .type = str
  mode = *basic advanced
    .help = "fast_ep operation setting"
    .type = choice
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

    def get_data(self):
        return self._parameters.fast_ep.data

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

    def get_ncycle(self):
        return self._parameters.fast_ep.ncycle

    def get_rlims(self):
        return self._parameters.fast_ep.rlims

    def get_xml(self):
        return self._parameters.fast_ep.xml

    def get_trace(self):
        return self._parameters.fast_ep.trace
    
    def get_mode(self):
        return self._parameters.fast_ep.mode

class Fast_ep:
    '''A class to run shelxc / d / e to very quickly establish (i) whether
    experimental phasing is likely to be successful and (ii) what the
    correct parameteters and number of heavy atom sites.'''

    def __init__(self, _parameters):
        '''Instantiate class and perform initial processing needed before the
        real work is done. This includes assessment of the anomalous signal,
        the signal to noise and conversion to pseudo-scalepack format of the
        data for input into shelxc, which is used to compute FA values.'''

        self._hklin = _parameters.get_data()
        self._native_hklin = _parameters.get_native()
        self._native = None
        self._cpu = _parameters.get_cpu()
        self._machines = _parameters.get_machines()
        self._atom = _parameters.get_atom()
        self._ntry = _parameters.get_ntry()
        self._ncycle = _parameters.get_ncycle()
        self._ano_rlimits = _parameters.get_rlims()
        self._trace = _parameters.get_trace()
        self._mode = _parameters.get_mode()
        self._data = None
        self._spacegroups = [_parameters.get_spg(),] if _parameters.get_spg() else None
        self._nsites = [_parameters.get_nsites(),] if _parameters.get_nsites() else None
        self._xml_name = _parameters.get_xml()
        self._start_time = datetime.now().strftime("%c")

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

        reader = any_reflection_file(self._hklin)
        self._file_type = reader.file_type()
    
        if self._file_type == 'ccp4_mtz':
            self._file_content = reader.file_content()
            self._is_merged = False if self._file_content.n_batches() > 0 else True
            self._all_data = [m for m in reader.as_miller_arrays(merge_equivalents=self._is_merged)
                    if type(m.observation_type()) is observation_types.intensity and (m.anomalous_flag() if self._is_merged else True)]
            if not self._all_data:
                raise RuntimeError, 'no intensity data found in %s' % self._hklin
            
            if self._native_hklin:
                native_reader = any_reflection_file(self._native_hklin)
                try:
                    self._native = next(m for m in native_reader.as_miller_arrays(merge_equivalents=True)
                                        if (type(m.observation_type()) is observation_types.intensity and
                                            not m.anomalous_flag()))
                except StopIteration:
                    self._native = None
        else:
            raise RuntimeError, 'Unsupported input file type: %s' % self._file_type

        self._nrefl = self._file_content.n_reflections()
        self._pointgroup = self._file_content.space_group_number()
        self._dmin, self._dmax = self._file_content.max_min_resolution()

        self._xml_results= {}
        self._xml_results['LOWRES'] = self._dmin
        self._xml_results['HIGHRES'] = self._dmax

        return


    def fa_values(self):

        self._log('Input:       %s' % self._hklin)
        self._log('N try:       %d' % self._ntry)
        dataset_names = ['sad',] if len(self._all_data) == 1 else ['peak', 'infl', 'hrem', 'lrem']
        if 'sad' in dataset_names:
            self._xml_results['SUBSTRUCTURE_METHOD'] = 'SAD'
        else:
            self._xml_results['SUBSTRUCTURE_METHOD'] = 'MAD'
        zip_dataset_names = zip(dataset_names, self._all_data)
        if self._native:
            zip_dataset_names.append(('nat', self._native))

        # write out a nice summary of the data set properties and what columns
        # were selected for analysis
        self._dataset_table = []
        for dtname, data in zip_dataset_names:
            self._log('Dataset:     %s' % dtname)
            self._log('Columns:     %s' % data.info().label_string())
            self._unit_cell = data.unit_cell().parameters()
            self._log('Unit cell:   %.2f %.2f %.2f %.2f %.2f %.2f' % \
                      self._unit_cell)
            self._log('Pointgroup:  %s' % data.crystal_symmetry().space_group().type().lookup_symbol())
            self._log('Resolution:  %.2f - %.2f' % data.resolution_range())

            self._log('Nrefl:       %d / %d' % (data.size(),
                                                data.n_bijvoet_pairs() if data.anomalous_flag() else 0))
            if data.anomalous_flag():
                self._log('DF/F:        %.3f' % data.anomalous_signal())

                differences = data.anomalous_differences()

                self._log('dI/sig(dI):  %.3f' % (sum(abs(differences.data())) /
                                                 sum(differences.sigmas())))

            table_vals = {'dtname': dtname,
                          'col_labels': data.info().label_string(),
                          'unit_cell': self._unit_cell,
                          'pg': data.crystal_symmetry().space_group().type().lookup_symbol(),
                          'resol_range': data.resolution_range(),
                          'nrefl': data.size(),
                          'n_pairs': data.n_bijvoet_pairs() if data.anomalous_flag() else 0,
                          'anom_flg': data.anomalous_flag()
                         }
            if data.anomalous_flag():
                table_vals.update({'anom_signal': data.anomalous_signal(),
                                   'anom_diff': (sum(abs(data.anomalous_differences().data())) / 
                                                 sum(data.anomalous_differences().sigmas())),
                                  })

            self._dataset_table.append(table_vals)

            
            # Now set up the job - run shelxc, assess anomalous signal, compute
            # possible spacegroup options, generate scalepack format reflection
            # file etc.

            if self._is_merged:
                intensities = data
            else:
                indices = self._file_content.extract_original_index_miller_indices()
                intensities = data.customized_copy(indices=indices, info=data.info())
            
            merge_scalepack.write(file_name = '.'.join([dtname, 'sca']),
                                  miller_array = intensities)

        # in here run shelxc to generate the ins file (which will need to be
        # modified) and the hkl files, which will need to be copied.

        if not self._spacegroups:
            self._spacegroups = generate_chiral_spacegroups_unique(self._pointgroup)

        self._log('Spacegroups: %s' % ' '.join(self._spacegroups))

        if not self._nsites:
            self._nsites = useful_number_sites(self._unit_cell, self._pointgroup)

        spacegroup = self._spacegroups[0]
        nsite = self._nsites[0]
        ntry = self._ntry

        self._xml_results['SHELXC_SPACEGROUP_ID'] = space_group_symbols(spacegroup).number()

        shelxc_input_files = ['%s %s.sca' % (v[0],v[0]) for v in zip_dataset_names]
        shelxc_stdin = ['cell %.3f %.3f %.3f %.3f %.3f %.3f' % self._unit_cell,
                 'spag %s' % sanitize_spacegroup(spacegroup),
                 'sfac %s' % self._atom.upper(),
                 'find %d' % nsite,
                 'mind -3.5',
                 'ntry %d' % ntry]
        shelxc_output = run_job('shelxc', ['sad'],
                                shelxc_input_files +
                                shelxc_stdin)

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

        plot_anom_shelxc(table['dmin'], table['isig'], table['dsig'], 'shelxc_anom.png')

        if self._ano_rlimits == [0]:
            self._ano_rlimits = [self._dmax]
        elif not self._ano_rlimits:
            self._ano_rlimits = [self._dmax] if self._mode == 'basic' \
                                else [self._dmax, self._dmax + 0.25, self._dmax + 0.5]

        self._log('Anomalous limits: %s' %  ' '.join(["%.2f" % v for v in self._ano_rlimits]))

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
        timeout = max(600, self._ntry)

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
                    wd = os.path.join(self._wd, spacegroup.replace(':', '-'), str(nsite), "%.2f" % rlimit)
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
            run_shelxd_drmaa_array(self._wd, nrefl, ncpu, njobs, jobs, timeout)
        else:
            pool = Pool(min(njobs, len(jobs)))
            pool.map(run_shelxd_local, jobs)

        # now gather up all of the results, find the one with best cfom

        results = get_shelxd_results(self._wd,
                                     self._spacegroups,
                                     self._nsites,
                                     self._ano_rlimits,
                                     self._mode == 'advanced')
        if self._mode == 'basic':
            best_keys, best_stats = max(filter(lambda (k, v): (v['nsites'] > 0), results.iteritems()),
                                        key=lambda (k, v): v['CFOM'])
        else:
            result_ranks = get_shelxd_result_ranks(results,
                                                   self._spacegroups,
                                                   self._nsites,
                                                   self._ano_rlimits)
            aver_ranks = get_average_ranks(self._spacegroups, self._nsites, self._ano_rlimits, results, result_ranks)
            best_sg, _ = min(filter(lambda (_, v): reduce(and_, (not isnan(x) for x in v.values())),
                                    aver_ranks.iteritems()),
                             key=lambda (_, v): min(v.values()))
        
            best_keys, best_stats = max(filter(lambda (k, v): (v['nsites'] > 0) and (k[0] == best_sg),
                                           results.iteritems()),
                                    key=lambda (k, v): v['CCres'])
            log_rank_table(self._log, aver_ranks, self._spacegroups, best_sg)
        
        (best_spacegroup,
         best_nsite,
         best_ano_rlimit) = best_keys

        (best_cc,
         best_ccweak,
         best_cfom,
         best_nsite_real) = [best_stats[col] for col in ['CCall', 'CCweak', 'CFOM', 'nsites']]
        
        if self._mode == 'basic':
            log_shelxd_results(self._log, results, self._spacegroups, best_keys, self._xml_results)
        else:
            log_shelxd_results_advanced(self._log, results, result_ranks, self._spacegroups, best_keys, self._xml_results)
        
        try:
            plot_shelxd_cc(self._wd, results, self._spacegroups, 'shelxd_cc.png')
            plot_shelxd_cc(self._wd,
                           dict([(k,v) for k,v in results.iteritems() if k[1] == best_nsite and k[2] == best_ano_rlimit]),
                           self._spacegroups, 'shelxd_cc_best.png')
            hist_shelxd_cc(self._wd, results, self._spacegroups)
        except:
            self._log("Exception thrown while plotting SHELXD results.")

        t1 = time.time()
        self._log('Time: %.2f' % (t1 - t0))

        if not best_spacegroup:
            raise RuntimeError, 'All shelxd jobs failed'

        self._log('Best spacegroup: %s' % best_spacegroup)
        self._log('Best nsites:     %d' % best_nsite_real)
        self._log('Best resolution: %.2f A' % best_ano_rlimit)
        self._log('Best CC / weak:  %.2f / %.2f' % (best_cc, best_ccweak))

        self._best_spacegroup = best_spacegroup
        self._best_nsite = best_nsite_real
        self._best_ano_rlimit = best_ano_rlimit
        self._best_cfom = best_cfom
        self._best_cc = best_cc
        self._best_ccweak = best_ccweak
        
        self._results = results

        # copy back result files

        best = os.path.join(self._wd, best_spacegroup.replace(':', '-'), str(best_nsite), "%.2f" % best_ano_rlimit)

        endings = ['lst', 'pdb', 'res']

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
        timeout = 600 + self._ncycle

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
                         'ncycle':self._ncycle,
                         'hand':'original', 'wd':wd})
            jobs.append({'nsite':self._best_nsite, 'solv':solvent_fraction,
                         'ncycle':self._ncycle,
                         'hand':'inverted', 'wd':wd})

        self._log('Running %d x shelxe jobs' % len(jobs))


        if cluster:
            run_shelxe_drmaa_array(self._wd, njobs, jobs, timeout)
        else:
            pool = Pool(min(njobs * ncpu, len(jobs)))
            pool.map(run_shelxe_local, jobs)

        shelxe_stats = read_shelxe_log(self._wd, solvent_fractions)
        skey = lambda s: '%.2f' % s
        
        best_solvent, best_hand, best_fom = max(((solv, hand, shelxe_stats['mean_fom_cc'][skey(solv)][hand]['mean_fom'])
                                                for solv, hand in product(solvent_fractions, ['original', 'inverted'])),
                                                key=lambda v: v[-1])
        
        try:
            plot_shelxe_contrast({best_solvent: shelxe_stats['contrast'][skey(best_solvent)]},
                                 os.path.join(self._wd, 'sad_best.png'), True)
            plot_shelxe_contrast(shelxe_stats['contrast'],
                                 os.path.join(self._wd, 'sad.png'))
            plot_shelxe_fom_mapcc(shelxe_stats['fom_mapcc'],
                                 os.path.join(self._wd, 'fom_mapcc.png'))
            plot_shelxe_mean_fom_cc(shelxe_stats['mean_fom_cc'],
                                 os.path.join(self._wd, 'mean_fom_cc.png'))
        except:
            self._log("Exception thrown while plotting SHELXE results.")

        self._best_fom = best_fom
        self._best_solvent = best_solvent
        self._best_hand = best_hand

        self._log('Solv. Orig. Inv.')
        for solvent_fraction in solvent_fractions:
            fom_orig, fom_inv = [shelxe_stats['mean_fom_cc'][skey(solvent_fraction)][hand]['pseudo_cc']
                                 for hand in ['original', 'inverted']]
            if solvent_fraction == best_solvent:
                self._log(
                    '%.2f %.3f %.3f (best)' % (solvent_fraction, fom_orig,
                                               fom_inv))
            else:
                self._log('%.2f %.3f %.3f' % (solvent_fraction, fom_orig,
                                              fom_inv))

        self._log('Best solvent: %.2f' % best_solvent)
        self._log('Best hand:    %s' % best_hand)

        wd = os.path.join(self._wd, skey(best_solvent))

        best_fom_mapcc = shelxe_stats['fom_mapcc'][skey(best_solvent)][best_hand]
        parse_pairs = [([self._dmin,] + best_fom_mapcc['resol'][:-1], 'RESOLUTION_LOW'),
                       (best_fom_mapcc['resol'], 'RESOLUTION_HIGH'),
                       (best_fom_mapcc['fom'], 'FOM'),
                       (best_fom_mapcc['mapcc'], 'MAPCC'),
                       (best_fom_mapcc['nrefl'], 'NREFLECTIONS')]
        for field_values, field_name in parse_pairs:
            store_string_xml(self._xml_results, field_values, field_name)
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

        self._log('Best spacegroup: %s' % self._best_spacegroup)

        if self._trace:
            # rerun shelxe to trace the chain
            self._nres_trace = 0
            arguments = ['sad', 'sad_fa', '-h%d' % self._best_nsite,
                         '-s%.2f' % best_solvent, '-d%.2f' % self._best_ano_rlimit, '-a3', '-m20']
            if not best_hand == 'original':
                arguments.append('-i')
            output = run_job('shelxe', arguments, [], self._wd)
            for record in output:
                if 'residues left after pruning' in record:
                    self._nres_trace = int(record.split()[0])
            pdb_org = os.path.join(self._wd, 'sad.pdb')
            pdb_inv = os.path.join(self._wd, 'sad_i.pdb')
            pdb_final = os.path.join(self._wd, 'sad_trace.pdb')
            try:
                if best_hand == 'inverted':
                    shutil.copyfile(pdb_inv, pdb_final)
                else:
                    shutil.copyfile(pdb_org, pdb_final)
                self._log('Traced:       %d residues' % self._nres_trace)
            except IOError:
                self._log('Chain tracing was unsuccessful.')

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


        template_dict = {'cpu'          : self._cpu,
                         'machines'     : self._machines,
                         'hklin'        : self._hklin,
                         'wd'           : self._wd,
                         'start_time'   : self._start_time,
                         'ntry'         : self._ntry,
                         'dataset_table': self._dataset_table,
                         'best_sg'      : self._best_spacegroup,
                         'nsite_real'   : self._best_nsite,
                         'best_rlim'    : '%.2f' % self._best_ano_rlimit,
                         'solv'         : '%2d' % (self._best_solvent * 100,),
                         'cc'           : self._best_cc,
                         'cc_weak'      : self._best_ccweak,
                         'cfom'         : self._best_cfom,
                         'fom'          : self._best_fom,
                         'hand'         : self._best_hand,
                         }
        try:
            render_html_report(template_dict)
        except:
            self._log("Exception thrown while generating HTML summary.")
        
        return


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
        fast_ep.fa_values()
    except RuntimeError, e:
        fast_ep._log('*** FA: %s ***' % str(e))
        traceback.print_exc(file = open('fast_ep.error', 'w'))
        sys.exit(1)

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
