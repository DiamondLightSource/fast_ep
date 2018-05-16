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
                             shell = True,
                             env = os.environ)

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
                    working_directory = None, ncpu = 1, timeout = None, sge_project=None):
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

    if timeout:
        timeout_tokens = ['-l', 'h_rt=%d' % timeout]
    else:
        timeout_tokens = []

    if sge_project:
        project_tokens = ['-P %s' % sge_project]
    else:
        project_tokens = []

    queue = 'medium.q'

    if ncpu > 1:
        qsub_output = run_job(
            'qsub', timeout_tokens + project_tokens + ['-V', '-pe', 'smp', str(ncpu),
                                      '-l', 'release="*"',
                                      '-cwd', '-q', queue,
                                      'FEP_%s.sh' % rs], [], working_directory)
    else:
        qsub_output = run_job(
            'qsub', timeout_tokens + project_tokens + ['-V', '-cwd', '-q', queue,
                                      'FEP_%s.sh' % rs], [], working_directory)

    if 'Unable to run job' in qsub_output[0]:
        raise RuntimeError, 'error submitting job to queue'


    job_id = None
    for record in qsub_output:
        if 'Your job' in record:
            job_id = int(record.split()[2])

    return job_id

def is_cluster_job_finished(job_id):

    qstat_output = run_job('qstat')

    for record in qstat_output:
        if str(job_id) in record:
            return False

    return True

def setup_job_drmaa(job, executable, arguments = [], stdin = [],
                    working_directory = None, ncpu = 1, timeout = None):
    '''Generate a script to run a program with some command-line arguments and
    setup cluster job for submission using DRMAA API.'''

    if working_directory is None:
        working_directory = os.getcwd()

    rs = random_string()
    script_path = os.path.join(working_directory, 'FEP_%s.sh' % rs)
    with open(script_path, 'w') as script:

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

    job.jobName = 'FEP_%s.sh' % rs
    job.remoteCommand = 'sh'
    job.workingDirectory = working_directory
    job.args = [script_path]

    qsub_args = ['-V',]
    if timeout:
        qsub_args += ['-l', 'h_rt=%s' % timeout]
    if ncpu > 1:
        qsub_args += ['-pe', 'smp', str(ncpu)]

    job.nativeSpecification = ' '.join(qsub_args)
    job.jobCategory = 'medium'
