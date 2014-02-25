#!/usr/bin/env python
#
# fast_mr ->
#
# Fast molecular replacement in the spirit of fast_dp, starting from coordinate
# files and using brute force (and educated guesses) to get everything going.
#
# fast_mr - main program.

import os
import sys
import time
import shutil
import math
import traceback
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

from xml_output import write_ispyb_xml

from generate_possible_spacegroups import generate_chiral_spacegroups, \
     spacegroup_enantiomorph, spacegroup_full, sanitize_spacegroup
from run_job import run_job, run_job_cluster, is_cluster_job_finished
from fast_mr_phaser import run_phaser_cluster
from parse_pdb import pdb_file_nres

class logger:
    def __init__(self):
        self._fout = open('fast_mr.log', 'w')
        return

    def __del__(self):
        self._fout.close()
        self._cout = None
        return

    def __call__(self, _line):
        sys.stdout.write('%s\n' % _line)
        self._fout.write('%s\n' % _line)
        return

class Fast_mr:

    def __init__(self, hklin, xyzin_and_ids):
        
        self._hklin = os.path.abspath(hklin)
        self._xyzins = [os.path.abspath(xyzin_and_id[0])
                        for xyzin_and_id in xyzin_and_ids]
        self._ids = [xyzin_and_id[1] for xyzin_and_id in xyzin_and_ids]
        self._cpu = 2
        self._machines = 10
        self._wd = os.getcwd()
        self._log = logger()

        self._log('Using %d cpus / %d machines' % (self._cpu, self._machines))

        self._full_command_line = ' '.join(sys.argv)

        # pull information we'll need from the input MTZ file - the unit cell,
        # the pointgroup and the number of reflections in the file. select
        # first Miller array in file which has native data

        # --- SAD DATA ---

        m = mtz.object(self._hklin)
        mas = m.as_miller_arrays()

        self._data = None

        for ma in mas:
            if str(ma.observation_type()) != 'xray.amplitude':
                continue
            self._data = ma
            break

        if not self._data:
            raise RuntimeError, 'no intensity data found in %s' % \
                self._hklin

        self._pointgroup = self._data.space_group().type().number()
        self._unit_cell = self._data.unit_cell().parameters()

        self._nrefl = m.n_reflections()
        self._spacegroups = generate_chiral_spacegroups(self._pointgroup)

        # write out a nice summary of the data set properties and what columns
        # were selected for analysis

        self._log('Input:       %s' % self._hklin)
        self._log('Columns:     %s' % self._data.info().label_string())
        self._log('Unit cell:   %.2f %.2f %.2f %.2f %.2f %.2f' % \
                  self._unit_cell)
        self._log('Pointgroup:  %s' % m.space_group().type().lookup_symbol())
        self._log('Resolution:  %.2f - %.2f' % self._data.resolution_range())
        self._log('Nrefl:       %d' % self._nrefl)
        self._log('Spacegroups: %s' % ' '.join(self._spacegroups))

        self._log('Input coordinate files:')
        self._nres = []
        for xyzin, _id in zip(self._xyzins, self._ids):
            nres = pdb_file_nres(xyzin)
            self._nres.append(nres)
            self._log('%40s %8d %.3f' % (os.path.split(xyzin)[1], nres, _id))

        total_nres = sum(self._nres)
        # FIXME calculate probable number of complexes in here

        self._copies = 1
        
        return

    def do_mr(self):
        t0 = time.time()

        cluster = True
        njobs = self._machines
        ncpu = self._cpu

        # set up N phaser jobs

        jobs = [ ]

        for spacegroup in self._spacegroups:
            wd = os.path.join(self._wd, spacegroup)
            if not os.path.exists(wd):
                os.makedirs(wd)
            commands = ['mode mr_auto',
                        'spacegroup %s' % spacegroup,
                        'hklin %s' % self._hklin,
                        'labin F=F SIGF=SIGF',
                        'root mr%s' % spacegroup]
            for j, (xyzin, _id, nres) in enumerate(
                zip(self._xyzins, self._ids, self._nres)):
                commands.append('ensemble m%d pdb %s identity %f' %
                                (j, xyzin, 100 * _id))
            for j, (xyzin, _id, nres) in enumerate(
                zip(self._xyzins, self._ids, self._nres)):
                commands.append('composition protein nres %d num %d' %
                                (nres, self._copies))
            for j, (xyzin, _id, nres) in enumerate(
                zip(self._xyzins, self._ids, self._nres)):
                commands.append('search ensemble m%d num %d' %
                                (j, self._copies))
            
            jobs.append((wd, commands))

        # actually execute the tasks - either locally or on a cluster, allowing
        # for potential for fewer available machines than jobs

        self._log('Running %d x phaser jobs' % len(jobs))
        
        pool = Pool(min(njobs, len(jobs)))

        if cluster:
            pool.map(run_phaser_cluster, jobs)
        else:
            print 1/0

        # now look for the results
        worked = []
        for job in jobs:
            wd = job[0]
            spacegroup = os.path.split(wd)[-1]
            if os.path.exists(os.path.join(wd, 'mr%s.sol' % spacegroup)):
                worked.append(os.path.join(wd, 'mr%s.sol' % spacegroup))

        for w in worked:
            sol = open(w).read()
            for record in sol.split('\n'):
                if 'SOLU SPAC' in record:
                    spacegroup = record.replace(
                        'SOLU SPAC', '').replace(' ', '')
                if 'SOLU SET' in record:
                    tfz = float(record.replace('=', ' ').split()[5])
            print 'Solution: %s %.2f' % (spacegroup, tfz)

        t1 = time.time()
        self._log('Time: %.2f' % (t1 - t0))

if __name__ == '__main__':
    xyzin_and_ids = []
    for arg in sys.argv[2:]:
        if ':' in arg:
            xyzin = arg.split(':')[0]
            _id = float(arg.split(':')[1])
            if _id > 1.0:
                _id /= 100.0
            xyzin_and_ids.append((xyzin, _id))
        else:
            xyzin_and_ids.append((arg, 1.0))
            
    fast_mr = Fast_mr(sys.argv[1], xyzin_and_ids)
    try:
        fast_mr.do_mr()
    except RuntimeError, e:
        fast_mr._log('*** MR: %s ***' % str(e))
        traceback.print_exc(file = open('fast_mr.error', 'w'))
        sys.exit(1)

