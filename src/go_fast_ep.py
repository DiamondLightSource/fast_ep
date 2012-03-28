#!/usr/bin/env python
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# a little program to decide to go (or no)
#
# criteria -
#
# low resolution Rmerge < 5%
# overall completeness > 95%
# then mid slope > 1.3 or mid slope > 1.15 and resolution better than 1.5

import sys
from iotbx import mtz
from cctbx.miller import build_set
from cctbx.crystal import symmetry as crystal_symmetry

def parse_fast_dp_log(fast_dp_log):

    resolution = None
    completeness = None
    mid_slope = None
    rmerge_low = None

    for record in open(fast_dp_log):
        if 'High resolution' in record:
            resolution = float(record.split()[2])
        if 'Rmerge' in record:
            rmerge_low = float(record.split()[2])
        if 'Anom. Completeness' in record:
            completeness = float(record.split()[2])
        if 'Mid-slope' in record:
            mid_slope = float(record.split()[1])

    return resolution, completeness, rmerge_low, mid_slope

def go_fast_ep(fast_dp_log):
    resolution, completeness, rmerge_low, mid_slope = parse_fast_dp_log(
        fast_dp_log)

    if rmerge_low > 0.05:
        return 'No'
    if completeness < 95:
        return 'No'
    if mid_slope > 1.3:
        return 'Go'
    if mid_slope > 1.15 and resolution < 1.5:
        return 'Go'
    return 'No'

def go_fast_ep_from_data(_mtz_file):
    '''Decide whether to run fast_ep or no based on the actual data based on
    the following criteria:

    completeness > 80% to dmin or 2.0 whichever the lower
    dI / s(dI) > 1.0 if resolution lower than 2.0, > 0.8 if better than
    1.5, smoothly varying in between.'''

    m = mtz.object(_mtz_file)

    mas = m.as_miller_arrays()
    data = None

    for ma in mas:
        if not ma.anomalous_flag():
            continue
        data = ma
        break

    if not data:
        return False

    d_min, d_max = sorted(data.resolution_range())

    if d_min > 2.0:
        differences = data.anomalous_differences()
        signal_to_noise = sum(abs(differences.data())) / \
            sum(differences.sigmas())
        completeness = data.completeness()

        if completeness < 0.8:
            return False
        if signal_to_noise < 1.0:
            return False
        return True
    
    else:
        data2 = data.resolution_filter(d_min = 2.0)
        differences = data2.anomalous_differences()
        signal_to_noise = sum(abs(differences.data())) / \
            sum(differences.sigmas())
        completeness = data2.completeness()

        if completeness < 0.8:
            return False
        if signal_to_noise < 0.8:
            return False

        return True
    
if __name__ == '__main__':

    status = False
    
    try:
        status = go_fast_ep_from_data(sys.argv[1])
    except:
        status = False

    if status:
        print 'Go'
    else:
        print 'No'
