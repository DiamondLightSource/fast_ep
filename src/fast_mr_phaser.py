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

def run_phaser_cluster(wd_commands):
    wd, commands = wd_commands
    job_id = run_job_cluster(
        'phaser', [], commands, wd, 1, timeout = 600)

    while not is_cluster_job_finished(job_id):
        time.sleep(1)
            
    return

