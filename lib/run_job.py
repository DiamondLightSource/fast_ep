#!/usr/bin/env python
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# run_job - a tool for running a job

import subprocess
import os
import random
import string

def random_string():
    return ''.join(random.sample(string.lowercase, 6))

def run_job(executable, arguments = [], stdin = [], working_directory = None):
    '''Run a program with some command-line arguments and some input,
    then return the standard output when it is finished.'''

    if working_directory is None:
        working_directory = os.getcwd()

    command_line = '%s' % executable
    for arg in arguments:
        command_line += ' "%s"' % arg

    popen = subprocess.Popen(command_line,
                             bufsize = 1,
                             stdin = subprocess.PIPE,
                             stdout = subprocess.PIPE,
                             stderr = subprocess.STDOUT,
                             cwd = working_directory,
                             universal_newlines = True,
                             shell = True)

    for record in stdin:
        popen.stdin.write('%s\n' % record)

    popen.stdin.close()

    output = []

    while True:
        record = popen.stdout.readline()
        if not record:
            break

        output.append(record)

    return output

def run_job_cluster(executable, arguments = [], stdin = [],
                    working_directory = None, ncpu = 1):
    '''Run a program with some command-line arguments and some input,
    then return the standard output when it is finished.'''

    if working_directory is None:
        working_directory = os.getcwd()

    rs = random_string()

    script = open(os.path.join(working_directory, 'FEP_%s.sh' % rs), 'w')

    script.write('#!/bin/bash\n')

    command_line = '%s' % executable
    for arg in arguments:
        command_line += ' "%s"' % arg

    if stdin:
        script.write('%s << eof\n' % command_line)
        for record in stdin:
            script.write('%s\n' % record)
        script.write('eof\n')
    else:
        script.write('%s\n' % command_line)        

    script.close()

    if ncpu > 1:
        qsub_output = run_job(
            'qsub', ['-V', '-pe', 'smp', str(ncpu), '-cwd', '-q', 'medium.q',
                     'FEP_%s.sh' % rs], [], working_directory)
    else:
        qsub_output = run_job(
            'qsub', ['-V', '-cwd', '-q', 'medium.q',
                     'FEP_%s.sh' % rs], [], working_directory)
        
    job_id = int(qsub_output[0].split()[2])

    return job_id

def is_cluster_job_finished(job_id):

    qstat_output = run_job('qstat')

    for record in qstat_output:
        if str(job_id) in record:
            return False

    return True

