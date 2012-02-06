#!/usr/bin/env python
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# number_sites_estimate - estimate the number of sites likely for a given
# unit cell and spacegroup.

from cctbx.sgtbx import space_group, space_group_symbols
from cctbx.uctbx import unit_cell

def number_residues_estimate(cell, pointgroup):
    '''Guess the number of residues in the ASU, assuming 50% solvent etc.'''

    sg = space_group(space_group_symbols(pointgroup).hall())
    uc = unit_cell(cell)

    n_ops = len(sg.all_ops())

    v_asu = uc.volume() / n_ops

    return int(round(v_asu / (2.7 * 128)))

def number_sites_estimate(cell, pointgroup):
    '''Guess # heavy atoms likely to be in here (as a floating point number)
    based on Matthews coefficient, average proportion of methionine in
    protein sequences and typical mass of an amino acid.'''

    sg = space_group(space_group_symbols(pointgroup).hall())
    uc = unit_cell(cell)

    n_ops = len(sg.all_ops())

    v_asu = uc.volume() / n_ops

    return 0.023 * v_asu / (2.7 * 128)

def test():
    cell = [42.18, 42.18, 39.24, 90.00, 90.00, 90.00]
    spacegroup = 'P422'
    number_sites = number_sites_estimate(cell, spacegroup)
    assert(number_sites < 1)
    assert(number_sites > 0)
    print 'OK'

if __name__ == '__main__':

    import sys

    if len(sys.argv) == 3:
        cell = map(float, sys.argv[1].split(','))
        spacegroup = sys.argv[2]
        print '%.1f' % number_sites_estimate(cell, spacegroup)
    else:
        test()
