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
    if completeness < 0.95:
        return 'No'
    if mid_slope > 1.3:
        return 'Go'
    if mid_slope > 1.15 and resolution < 1.5:
        return 'Go'
    return 'No'

if __name__ == '__main__':
    print go_fast_ep(sys.argv[1])
        
