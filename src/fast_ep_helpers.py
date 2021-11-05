#!/usr/bin/env python
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# format_for_autosharp - format an input file for autosharp.

import os
from math import isinf, isnan
from cctbx.sgtbx import space_group, space_group_symbols

from lib.number_sites_estimate import number_sites_estimate
from lib.run_job import run_job

_autosharp_file = '''
# ----------------------------------------
# General information about project
# ----------------------------------------
  autoSHARP_proj="random"
  autoSHARP_jobi="0"
  autoSHARP_titl="random"
  autoSHARP_type="MAD"
  autoSHARP_rate="5"
  autoSHARP_molw="0.0"
  autoSHARP_nres="{nres:d}"
  autoSHARP_pirf=""
  autoSHARP_spgr="None"
  autoSHARP_resl="9999.0"
  autoSHARP_resh="0.01"
  autoSHARP_nset="1"
  autoSHARP_user="{user:s}"
  autoSHARP_ulvl="3"
  autoSHARP_chtm="no"
  autoSHARP_csum="no"
# ----------------------------------------
# autoSHARP protocol information
# ----------------------------------------
  autoSHARP_EntryPoint="2"
  autoSHARP_EntryPoint3_Path="3"
  autoSHARP_EntryPoint3_Path3_Opt="7"
  autoSHARP_DetectPgm="shelx"
# ----------------------------------------
# Ersatz CCP4i stuff
# ----------------------------------------
  autoSHARP_RunningType="ccp4i"
  autoSHARP_ccp4i_workdir="{working_dir:s}"
# ----------------------------------------
# Dataset 1
# ----------------------------------------
  autoSHARP_iden_1="peak"
  autoSHARP_wave_1="{wavelength:.4f}"
  autoSHARP_hatm_1="{atom:s}"
  autoSHARP_nsit_1="{nsites:d}"
  autoSHARP_sitf_1=""
  autoSHARP_fone_1="0.0"
  autoSHARP_ftwo_1="0.0"
  autoSHARP_fmid_1="F"
  autoSHARP_smid_1="SIGF"
  autoSHARP_dano_1="DANO"
  autoSHARP_sano_1="SIGDANO"
  autoSHARP_isym_1="ISYM"
  autoSHARP_dtyp_1="MTZ"
  autoSHARP_data_1="{hklin:s}"
'''

def ctruncate_anomalous_signal(hklin):
    '''

    Estimated limits of anomalous signal
      Wang limit (deltaI/I) > 0.6% : 1.3 A
      anomalous limit (deltaI/sig) > 1.3 : 2.1 A
      measurability limit (Nanon/Nov) > 5% : 1.8 A

Use

      measurability limit (Nanon/Nov) > 5% : 1.8 A

      + / - 0.2 A


from

ctruncate -mtzin ../AUTOMATIC_DEFAULT_scaled.mtz -mtzout truncated.mtz -colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'

'''
    ctruncate_output = run_job(' '.join(['ctruncate',
                                         '-mtzin', hklin,
                                         '-mtzout', 'truncated.mtz',
                                         '-colano', '"/*/*/[I(+),SIGI(+),I(-),SIGI(-)]"']),
                                         [], [])
    open('ctruncate.log', 'w').write(''.join(ctruncate_output))

    maximum_resolution = float('nan')
    rlimit = float('nan')

    for record in ctruncate_output:
        if 'Maximum resolution =' in record:
            maximum_resolution = float(record.split()[-2])
        if "measurability limit (Nanon/Nov)" in record:
            rlimit = float(record.split()[-2])
            break

    if not isinf(rlimit) and not isnan(rlimit) and \
        not isinf(maximum_resolution) and not isnan(maximum_resolution):
        return [max(maximum_resolution, rlimit - 0.2), rlimit, rlimit + 0.2]

    return None

def autosharp(nres, user, wavelength, atom, nsites, hklin):
    return _autosharp_file.format(
        nres = nres,
        user = user,
        working_dir = os.getcwd(),
        wavelength = wavelength,
        atom = atom,
        nsites = nsites,
        hklin = hklin)


def map_sites_to_asu(spacegroup,
                     pdb_in,
                     pdb_out,
                     invert=False):
    '''Map sites to asu of input spacegroup (as sites from shelxd claim
    P1 in CRYST1 record) inverting if necessary. N.B. if inverting sites
    also need to invert spacegroup.'''

    from cctbx.crystal import symmetry, direct_space_asu
    from iotbx.pdb import hierarchy
    from scitbx.array_family import flex

    sg = space_group(space_group_symbols(spacegroup).hall())
    coords = hierarchy.input(file_name=pdb_in)
    cs = coords.input.crystal_symmetry()
    uc = cs.unit_cell()
    cs2 = symmetry(unit_cell=uc, space_group=sg)
    xs = coords.xray_structure_simple().customized_copy(crystal_symmetry=cs2)

    if invert:
        xs = xs.change_hand()

    am = xs.crystal_symmetry().asu_mappings(0.0)
    xyz = xs.sites_cart()
    am.process_sites_cart(xyz)
    xyz = flex.vec3_double()
    for m in am.mappings():
        xyz.append(m[0].mapped_site())
    xs.set_sites_cart(xyz)

    open(pdb_out, 'w').write(xs.as_pdb_file())

    return

def useful_number_sites(_cell, _pointgroup):
    nha = number_sites_estimate(_cell, _pointgroup)

    result = []

    for f in [0.25, 0.5, 1.0, 2.0, 4.0]:
        nha_test = int(round(f * nha))
        if nha_test and not nha_test in result:
            result.append(nha_test)

    if len(result) < 3:
        result = [1, 2, 3]

    return result

def modify_ins_text(_ins_text, _spacegroup, _nsites, _rlimit):
    '''Update the text in a SHELXD .ins file to handle the correct number
    of sites and spacegroup symmetry operations.'''

    new_text = []

    symm = [op.as_xyz().upper() for op in
            space_group(space_group_symbols(_spacegroup).hall()).smx()]

    for record in _ins_text:
        if 'SYMM' in record:
            if not symm:
                continue
            for op in symm:
                if op == 'X,Y,Z':
                    continue
                new_text.append(('SYMM %s' % op))
            symm = None
        elif 'FIND' in record:
            new_text.append(('FIND %d' % _nsites))
        elif 'SHEL' in record:
            new_text.append(('SHEL 999 %.1f' % _rlimit))
        else:
            new_text.append(record.strip())

    return new_text

def get_scipy():
  import os, sys
  # make sure we can get scipy, if not try failing over to version in CCP4
  try:
    import scipy.cluster
    found = True
  except ImportError, e:
    found = False

  if not found and 'CCP4' in os.environ:
    sys.path.append(os.path.join(os.environ['CCP4'], 'lib', 'python2.7',
                                 'site-packages'))
    try:
      import scipy.cluster
      found = True
    except ImportError, e:
      found = False

  if not found:
      raise RuntimeError, 'fast_ep needs scipy'


if __name__ == '__main__':

    print(autosharp(180, 'gw56', 0.97966, 'Se', 8, 'xia2.txt'))
