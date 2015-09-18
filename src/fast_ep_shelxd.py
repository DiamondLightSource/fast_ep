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

from run_job import run_job, run_job_cluster, is_cluster_job_finished, setup_job_drmaa

def run_shelxd_cluster(_settings):
    '''Run shelxd on cluster with settings given in dictionary, containing:

    nrefl = 1 + floor(nref / 100000) - space to allocate
    ncpu - number of cpus to use
    wd - working directory'''

    nrefl = _settings['nrefl']
    ncpu = _settings['ncpu']
    wd = _settings['wd']

    job_id = run_job_cluster(
        'shelxd', ['-L%d' % nrefl, 'sad_fa', '-t%d' % ncpu],
        [], wd, ncpu, timeout = 600)

    while not is_cluster_job_finished(job_id):
        time.sleep(1)

    return

def run_shelxd_drmaa(njobs, job_settings):
    '''Run shelxd on cluster with settings given in dictionary, containing:

    nrefl = 1 + floor(nref / 100000) - space to allocate
    ncpu - number of cpus to use
    wd - working directory'''

    import drmaa
    with drmaa.Session() as session:

        job = session.createJobTemplate()

        batches = range(0, len(job_settings), njobs)
        for idx in batches:
            jobs = []
            for _settings in job_settings[idx:idx+njobs]:

                nrefl = _settings['nrefl']
                ncpu = _settings['ncpu']
                wd = _settings['wd']

                setup_job_drmaa(job,
                                'shelxd', ['-L%d' % nrefl, 'sad_fa', '-t%d' % ncpu],
                                [], wd, ncpu, timeout = 600)
                jobs.append(session.runJob(job))
            session.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
        session.deleteJobTemplate(job)
    return

def run_shelxd_local(_settings):
    '''Run shelxd locally settings given in dictionary, containing:

    nrefl = 1 + floor(nref / 100000) - space to allocate
    ncpu - number of cpus to use
    wd - working directory'''

    nrefl = _settings['nrefl']
    ncpu = _settings['ncpu']
    wd = _settings['wd']

    job_output = run_job(
        'shelxd', ['-L%d' % nrefl, 'sad_fa', '-t%d' % ncpu], [], wd)

    open(os.path.join(wd, 'shelxd.log'), 'w').write(''.join(job_output))

    return

def analyse_res(_res):

    try:
        cc = float(_res[0].split()[5])
    except (ValueError, IndexError):
        cc = float('nan')

    try:
        cc_weak = float(_res[0].split()[7])
    except (ValueError, IndexError):
        cc_weak = float('nan')

    try:
        cfom = float(_res[0].split()[9])
    except (ValueError, IndexError):
        cfom = float('nan')

    # estimate real # sites - as drop below 30% relative occupancy

    nsites_real = 0

    try:
        for record in _res:
            if not 'SE' in record[:2]:
                continue
            if float(record.split()[5]) > 0.3:
                nsites_real += 1
    except (ValueError, IndexError):
        nsites_real = float('nan')

    return cc, cc_weak, cfom, nsites_real

def happy_shelxd_log(_shelxd_lst_file):
    for record in open(_shelxd_lst_file):
        if '** NO SUITABLE PATTERSON VECTORS FOUND **' in record:
            return False
        if '** CANNOT ALLOCATE ENOUGH MEMORY **' in record:
            return False

    best_cfom = 0.0

    # columns can get merged in output, so watch for that

    for record in open(_shelxd_lst_file):
        if '*****' in record:
            continue
        if record.startswith(' Try'):
            best_cfom_token = record.replace(',', ' ').split()[-3]
            best_cfom = float(best_cfom_token.replace('best', ''))

    if best_cfom == 0.0:
        return False

    return True
