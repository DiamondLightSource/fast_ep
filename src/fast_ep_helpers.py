#!/usr/bin/env python
#
# fast_ep ->
#
# Fast experimental phasing in the spirit of fast_dp, starting from nothing
# and using brute force (and educated guesses) to get everything going.
#
# format_for_autosharp - format an input file for autosharp.

import os

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

def autosharp(nres, user, wavelength, atom, nsites, hklin):
    return _autosharp_file.format(
        nres = nres,
        user = user,
        working_dir = os.getcwd(),
        wavelength = wavelength,
        atom = atom,
        nsites = nsites,
        hklin = hklin)

def plot_shelxd_cc(fa_lst_file, png_file, spacegroup, sites):
    '''Plot cc weak vs. cc from shelxd run.'''

    # first scrape out the cc values

    cc_all = []
    cc_weak = []
    
    for record in open(fa_lst_file):
        if 'Try' in record and 'CC All/Weak' in record:
            tokens = record.replace(',', ' ').replace('CPU', 'CPU ').split()
            cc_all.append(float(tokens[6]))
            cc_weak.append(float(tokens[8]))

    # now generate plot 
            
    from matplotlib import pyplot

    pyplot.xlabel('CC (all)')
    pyplot.ylabel('CC (weak)')
    pyplot.title('Substructure search %d sites in %s' % (sites, spacegroup))
    pyplot.scatter(cc_all, cc_weak, label = 'CC')
    pyplot.axis([-10, 100, -10, 100])
    pyplot.legend()
    pyplot.savefig(png_file)
    pyplot.close()
    
    return

def plot_shelxe_contrast(original_lst, other_lst, png_file, solvent):
    '''Plot contrast vs. cycle number from shelxe.'''

    # first scrape out the contrast values

    contrast_orig = []
    contrast_other = []
    
    for record in open(original_lst):
        if 'Contrast' in record and 'Connect' in record:
            tokens = record.replace(',', ' ').split()
            contrast_orig.append(float(tokens[5]))
        
    for record in open(other_lst):
        if 'Contrast' in record and 'Connect' in record:
            tokens = record.replace(',', ' ').split()
            contrast_other.append(float(tokens[5]))

    cycles = [j + 1 for j in range(len(contrast_orig))]
            
    # now generate plot 
            
    from matplotlib import pyplot

    pyplot.xlabel('Cycle')
    pyplot.ylabel('Contrast')
    pyplot.title('Phasing contrast for solvent fraction %.2f' % solvent)
    pyplot.plot(cycles, contrast_orig, label = 'Original')
    pyplot.plot(cycles, contrast_other, label = 'Inverse')
    pyplot.axis([0, len(contrast_orig), 0, 1.5])
    pyplot.legend()
    pyplot.savefig(png_file)
    pyplot.close()

    return

if __name__ == '__main__':

    print autosharp(180, 'gw56', 0.97966, 'Se', 8, 'xia2.txt')
