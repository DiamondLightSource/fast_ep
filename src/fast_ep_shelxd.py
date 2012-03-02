#!/usr/bin/env python
#
# fast_ep_shelxd ->
# 
# code to run shelxd and manage the jobs - this will be called from within
# a multiprocess task so life is easier if the input is provided in the form
# of a dictionary with the name of the "problem" we're working on and the 
# number of reflections to make room for. 

import os
import sys
import time

if not 'FAST_EP_ROOT' in os.environ:
    raise RuntimeError, 'FAST_EP_ROOT not set'

fast_ep_lib = os.path.join(os.environ['FAST_EP_ROOT'], 'lib')

if not fast_ep_lib in sys.path:
    sys.path.append(fast_ep_lib)

from run_job import run_job, run_job_cluster, is_cluster_job_finished

def run_shelxd_cluster(settings)
    '''Run shelxd_mp on cluster with settings given in dictionary, containing:
    
    nrefl = 1 + floor(nref / 100000) - space to allocate
    ncpu - number of cpus to use
    wd - working directory'''

    nrefl = settings['nrefl']
    ncpu = settings['ncpu']
    wd = settings['wd']

    job_id = run_job_cluster(
        'shelxd_mp', ['-L%d' nrefl, 'sad_fa', '-t%d' % ncpu],
        [], wd, ncpu)

    while not is_cluster_job_finished(job_id):
        time.sleep(1)

    return

def run_shelxd_local(settings)
    '''Run shelxd_mp locally settings given in dictionary, containing:
    
    nrefl = 1 + floor(nref / 100000) - space to allocate
    ncpu - number of cpus to use
    wd - working directory'''

    nrefl = settings['nrefl']
    ncpu = settings['ncpu']
    wd = settings['wd']

    job_output = run_job(
        'shelxd_mp', ['-L%d' nrefl, 'sad_fa', '-t%d' % ncpu], [], wd)

    return

def analyse_res(res):
    cc = float(res[0].split()[5])
    cc_weak = float(res[0].split()[7])
    cfom = float(res[0].split()[9])
    
    # estimate real # sites - as drop below 30% relative occupancy
    
    nsites_real = 0

    for record in res:
        if not 'SE' in record[:2]:
            continue
        if float(record.split()[5]) > 0.3:
            nsites_real += 1

    return cc, cc_weak, cfom, nsites_real
