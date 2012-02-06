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

from iotbx import mtz
from cctbx.sgtbx import space_group, space_group_symbols

if not 'FAST_EP_ROOT' in os.environ:
    raise RuntimeError, 'FAST_EP_ROOT not set'

fast_ep_lib = os.path.join(os.environ['FAST_EP_ROOT'], 'lib')

if not fast_ep_lib in sys.path:
    sys.path.append(fast_ep_lib)

from generate_possible_spacegroups import generate_chiral_spacegroups_unique
from number_sites_estimate import number_sites_estimate
from run_job import run_job, run_job_cluster, is_cluster_job_finished

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
            space_group(space_group_symbols(spacegroup).hall()).all_ops()]

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
    
def fast_ep(hklin):
    '''Quickly (i.e. with the aid of lots of CPUs) try to experimentally
    phase the data in here.'''

    working_directory = os.getcwd()

    m = mtz.object(hklin)
    
    pointgroup = m.space_group().type().number()

    for crystal in m.crystals():
        if crystal.name() != 'HKL_base':
            unit_cell = crystal.unit_cell().parameters()

    # first run mtz2sca on this input => sad.sca

    run_job('mtz2sca', [hklin, 'sad.sca'])

    # in here run shelxc to generate the ins file (which will need to be
    # modified) and the hkl files, which will need to be copied.

    spacegroup_0 = generate_chiral_spacegroups_unique(pointgroup)[0]
    nsites_0 = useful_number_sites(unit_cell, pointgroup)[0]

    run_job('shelxc', ['sad'],
            ['sad sad.sca',
             'cell %.3f %.3f %.3f %.3f %.3f %.3f' % unit_cell,
             'spag %s' % spacegroup_0,
             'find %d' % nsites_0,
             'mind -3.5',
             'ntry 200'])
    
    # then for all possible spacegroups and all possible numbers of
    # sites modify the ins file and run shelxd on the cluster.

    ins_text = open('sad_fa.ins', 'r').readlines()

    job_ids = []

    for spacegroup in generate_chiral_spacegroups_unique(pointgroup):
        for nsites in useful_number_sites(unit_cell, pointgroup):
            
            wd = os.path.join(working_directory, spacegroup, str(nsites))
            if not os.path.exists(wd):
                os.makedirs(wd)

            new_text = modify_ins_text(ins_text, spacegroup, nsites)

            shutil.copyfile(os.path.join(working_directory, 'sad_fa.hkl'),
                            os.path.join(wd, 'sad_fa.hkl'))

            open(os.path.join(wd, 'sad_fa.ins'), 'w').write(
                '\n'.join(new_text))

            job_id = run_job_cluster('shelxd_mp', ['-L4', 'sad_fa'], [], wd)
            job_ids.append(job_id)

    for job_id in job_ids:
        while not is_cluster_job_finished(job_id):
            time.sleep(1)

    results = { }

    best_cfom = 0.0
    best_spacegroup = None
    best_nsites = None

    for spacegroup in generate_chiral_spacegroups_unique(pointgroup):
        for nsites in useful_number_sites(unit_cell, pointgroup):
            res = open(os.path.join(working_directory, spacegroup, str(nsites),
                                    'sad_fa.res'), 'r').readlines()

            cc = float(res[0].split()[5])
            cc_weak = float(res[0].split()[7])
            cfom = float(res[0].split()[9])

            results[(spacegroup, nsites)] = (cc, cc_weak, cfom)

            if cfom > best_cfom:
                best_cfom = cfom
                best_spacegroup = spacegroup
                best_nsites = nsites

    for spacegroup in generate_chiral_spacegroups_unique(pointgroup):
        for nsites in useful_number_sites(unit_cell, pointgroup):

            (cc, cc_weak, cfom) = results[(spacegroup, nsites)]

            if (spacegroup, nsites) == (best_spacegroup, best_nsites):
                print '%10s %2d %6.2f %6.2f %6.2f ***' % \
                      (spacegroup, nsites, cc, cc_weak, cfom)
            else:
                print '%10s %2d %6.2f %6.2f %6.2f' % \
                      (spacegroup, nsites, cc, cc_weak, cfom)
                
if __name__ == '__main__':
    fast_ep(sys.argv[1])
    
