#!/usr/bin/env python
#
# fast_ep_shelxe ->
#
# code to run shelxe and manage the jobs - this will be called from within
# a multiprocess task so life is easier if the input is provided in the form
# of a dictionary with the information we need.

import os
import time

from lib.run_job import run_job, run_job_cluster, is_cluster_job_finished, setup_job_drmaa


def run_shelxe_cluster(_settings):
    '''Run shelxe on cluster with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    resol - high resolution limit
    hand - original or inverted
    wd - working directory'''

    nsite = _settings['nsite']
    solv = _settings['solv']
    hand = _settings['hand']
    resol = _settings['resol']
    wd = _settings['wd']

    if hand == 'original':
        job_id = run_job_cluster(
            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite, '-d%f' % resol,
                       '-s%f' % solv, '-m20'], [], wd, 1, timeout = 600)
    else:
        job_id = run_job_cluster(
            'shelxe', ['sad', 'sad_fa', '-h%d' % nsite, '-d%f' % resol,
                       '-s%f' % solv, '-m20', '-i'], [], wd, 1, timeout = 600)

    while not is_cluster_job_finished(job_id):
        time.sleep(1)

    return


def run_shelxe_drmaa(njobs, job_settings):
    '''Run shelxe on cluster with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    resol - high resolution limit
    hand - original or inverted
    wd - working directory'''

    import drmaa
    with drmaa.Session() as session:

        job = session.createJobTemplate()

        batches = range(0, len(job_settings), njobs)
        for idx in batches:
            jobs = []
            for _settings in job_settings[idx:idx+njobs]:


                nsite = _settings['nsite']
                solv = _settings['solv']
                hand = _settings['hand']
                resol = _settings['resol']
                wd = _settings['wd']

                if hand == 'original':
                    setup_job_drmaa(job,
                                    'shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                    '-s%f' % solv, '-d%f' % resol, '-m20'], [], wd, 1, timeout = 600)
                else:
                    setup_job_drmaa(job,
                                    'shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                    '-s%f' % solv, '-d%f' % resol, '-m20', '-i'], [], wd, 1, timeout = 600)

                jobs.append(session.runJob(job))
            session.synchronize(jobs, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
        session.deleteJobTemplate(job)
    return


def run_shelxe_drmaa_array(wd, njobs, job_settings, timeout):
    '''Run shelxe on cluster with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    resol - high resolution limit
    hand - original or inverted
    wd - working directory'''

    script_path = os.path.join(wd, 'shelxe_batch.sh')
    with open(script_path, 'w') as script:

        script.write('#!/bin/bash\n')

        for idx, _settings in enumerate(job_settings, start=1):

            hand = '-i' if _settings['hand'] == 'inverted' else ''
            script.write('WORKING_DIR_{idx}={wd}\n'.format(idx=idx,
                                                           wd= _settings['wd']))
            script.write('COMMAND_{idx}="shelxe sad sad_fa -h{nsite} -s{solv} -m{ncycle} {hand}"\n'.format(idx=idx,
                                                                                                     nsite=_settings['nsite'],
                                                                                                     solv=_settings['solv'],
                                                                                                     ncycle=_settings['ncycle'],
                                                                                                     hand=hand))

        script.write('TASK_WORKING_DIR=WORKING_DIR_${SGE_TASK_ID}\n')
        script.write('TASK_COMMAND=COMMAND_${SGE_TASK_ID}\n')
        script.write('cd ${!TASK_WORKING_DIR}\n')
        script.write('${!TASK_COMMAND} > ${!TASK_WORKING_DIR}/FEP_shelxe_${SGE_TASK_ID}.out  2> ${!TASK_WORKING_DIR}/FEP_shelxe_${SGE_TASK_ID}.err')

    import drmaa
    with drmaa.Session() as session:
        job = session.createJobTemplate()
        job.jobName = 'FEP_shelxe'
        job.workingDirectory = wd
        job.remoteCommand = 'sh'
        args = [script_path,]
        job.args = args

        if os.environ.get('USER', '') == 'gda2':
            job.jobCategory = 'high'
        else:
            job.jobCategory = 'medium'

        job.nativeSpecification = '-V -l h_rt={timeout} -tc {njobs}  -o /dev/null -e /dev/null'.format(timeout=timeout,
                                                                                                       njobs=njobs)

        job_ids = session.runBulkJobs(job, 1, len(job_settings), 1)
        session.synchronize(job_ids, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
        session.deleteJobTemplate(job)


def run_shelxe_local(_settings):
    '''Run shelxe locally with settings given in dictionary, containing:

    nsite - number of sites
    solv - solvent fraction
    resol - high resolution limit
    hand - original or inverted
    wd - working directory'''

    nsite = _settings['nsite']
    solv = _settings['solv']
    hand = _settings['hand']
    resol = _settings['resol']
    wd = _settings['wd']

    if hand == 'original':
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                        '-s%f' % solv, '-d%f' % resol, '-m20'], [], wd)
    else:
        job_output = run_job('shelxe', ['sad', 'sad_fa', '-h%d' % nsite,
                                        '-s%f' % solv, '-d%f' % resol, '-m20', '-i'], [], wd)

    return


def read_shelxe_log(pth, solvent_fractions):

    res = {'contrast': {},
           'fom_mapcc': {},
           'mean_fom_cc': {}}
    
    for solvent_fraction in solvent_fractions:
        solv_key = '%.2f' % solvent_fraction
        for stats in res.itervalues():
            stats[solv_key] = {}
        
        wd = os.path.join(pth, solv_key)
        for hand, lst in [('original', os.path.join(wd, 'sad.lst'),),
                          ('inverted', os.path.join(wd, 'sad_i.lst'))]:
            
            contrast = []
            fom_mapcc = {}
            
            for record in open(lst):
                if 'Contrast' in record and 'Connect' in record:
                    tokens = record.replace(',', ' ').split()
                    contrast.append(float(tokens[5]))
                elif 'Estimated mean FOM =' in record:
                    mean_fom = float(record.split()[4])
                    pseudo_cc = float(record.split()[-2]) / 100.
                    res['mean_fom_cc'][solv_key][hand] = {'mean_fom': mean_fom,
                                                          'pseudo_cc': pseudo_cc}
                else:
                    for lbl, k, lst in [('d    inf', ' - ', 'resol'),
                                        ('<FOM>', None, 'fom'),
                                        ('<mapCC>', None, 'mapcc'),
                                        ('N      ', None, 'nrefl')]:
                        if lbl in record:
                            fom_mapcc[lst] = [float(s) for s in record.split(k)[1:]]
            res['contrast'][solv_key][hand] = contrast
            res['fom_mapcc'][solv_key][hand] = fom_mapcc
    
    return res

