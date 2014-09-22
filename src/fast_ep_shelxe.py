#!/usr/bin/env python
#
# fast_ep_shelxe ->
#
# code to run shelxe and manage the jobs - this will be called from within
# a multiprocess task so life is easier if the input is provided in the form
# of a dictionary with the information we need.

import os
import sys
import time

if not 'FAST_EP_ROOT' in os.environ:
    raise RuntimeError, 'FAST_EP_ROOT not set'

fast_ep_lib = os.path.join(os.environ['FAST_EP_ROOT'], 'lib')

if not fast_ep_lib in sys.path:
    sys.path.append(fast_ep_lib)

from run_job import run_job, run_job_cluster, is_cluster_job_finished

def run_shelxe_cluster(_settings):
    '''Run shelxe on cluster with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    hand - original or inverted
    wd - working directory'''

    nsite = _settings['nsite']
    solv = _settings['solv']
    hand = _settings['hand']
    wd = _settings['wd']

    if hand == 'original':
        job_id = run_job_cluster(
            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                       '-s%f' % solv, '-m20'], [], wd, 1, timeout = 600)
    else:
        job_id = run_job_cluster(
            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                       '-s%f' % solv, '-m20', '-i'], [], wd, 1, timeout = 600)

    while not is_cluster_job_finished(job_id):
        time.sleep(1)

    return

def run_shelxe_local(_settings):
    '''Run shelxe locally with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    hand - original or inverted
    wd - working directory'''

    nsite = _settings['nsite']
    solv = _settings['solv']
    hand = _settings['hand']
    wd = _settings['wd']

    if hand == 'original':
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                        '-s%f' % solv, '-m20'], [], wd)
    else:
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                        '-s%f' % solv, '-m20', '-i'], [], wd)

    return
