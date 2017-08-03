#!/usr/bin/env python
#
# fast_ep_plots ->
#
# Plotting routines for SHELX results
#

import os.path
from itertools import product
from math import ceil

import matplotlib
matplotlib.use('Agg')
from matplotlib import ticker
import matplotlib.pylab as plt
import matplotlib.cm as cm

from src.fast_ep_shelxd import read_shelxd_log

params = {
   'axes.labelsize': 12,
   'font.size': 12,
   'legend.fontsize': 12,
   'legend.columnspacing': 0.5,
   'legend.scatterpoints' : 5,
   'legend.handletextpad' : 0.3,
   'axes.labelpad': 15.0,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.subplot.hspace'  : .3,
   }

matplotlib.rcParams.update(params)


def plot_anom_shelxc(resol, isig, dsig, png_file):
    '''Plot <I/sig> and <d"/sig> vs. resolution from SHELXC'''

    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax2 = ax1.twinx()

    x = range(len(resol))
    ax1.plot(x, isig, lw=1, label='<I/sig>', c='r')
    ax2.plot(x, dsig, lw=1, label='<d"/sig>', c='b')

    ax1.set_xlabel('Resolution / $\mathsf{\AA}$', fontsize=14)
    ax1.set_ylabel('<I/sig>', fontsize=14)
    ax2.set_ylabel('<d"/sig>', fontsize=14)

    plt.xticks(x, resol)
    ax1.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    lgd1 = ax1.legend(h1 + h2, l1 + l2, bbox_to_anchor=[.25, 1.02], loc=3, ncol=2, fontsize=14)
    plt.savefig(png_file, bbox_extra_artists=(lgd1, ), bbox_inches='tight')
    plt.close()


def plot_shelxd_cc(pth, results, spacegroups, png_file):
    '''Summary plot of CCall / CCweak values from SHELXD for all
    space group, heavy atom number and resolution combinations'''

    _, nsites, ano_rlimits = map(sorted, map(set, zip(*results.keys())))

    fig, axarr = plt.subplots(len(ano_rlimits), len(nsites), squeeze=False,
                              figsize=(4*len(nsites), 3*len(ano_rlimits)))
    shelxd_plot = os.path.join(pth, png_file)
    for (i, rlimit), (j, nsite) in product(enumerate(ano_rlimits), enumerate(nsites)):
        tmp_ax = axarr[i,j]
        cc_all_sg, cc_weak_sg = [], []
        for clr, spacegroup in enumerate(spacegroups):
            color = cm.jet(float(clr)/len(spacegroups))
            wd = os.path.join(pth, spacegroup.replace(':', '-'), str(nsite), "%.2f" % rlimit)
            shelxd_log = os.path.join(wd, 'sad_fa.lst')
            cc_all, cc_weak, _ = read_shelxd_log(shelxd_log)
            cc_all_sg.extend(cc_all)
            cc_weak_sg.extend(cc_weak)

            tmp_ax.scatter(cc_all, cc_weak, c=color, s=5, label=spacegroup, lw=0, alpha=0.75)
            textstr = 'HA : %s  Resol: %.2f$\mathsf{\AA}$' % (nsite, rlimit)
            tmp_ax.set_title(textstr, fontsize=10)
            if i == len(ano_rlimits) - 1:
                tmp_ax.set_xlabel('CCall', fontsize=10)
            if j == 0:
                tmp_ax.set_ylabel('CCweak', fontsize=10)
            if i == 0 and j == len(nsites) - 1:
                lgd = tmp_ax.legend(bbox_to_anchor=(1.02, 1), loc=2, fontsize=10)
    plt.savefig(shelxd_plot, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def hist_shelxd_cc(pth, results, spacegroups, rows=2):
    '''Histogram plot of SHELXD results for best solutions found for all
    combinations of space group, heavy atom number and resolution parameters'''

    _, nsite_set, resol_set = map(sorted, map(set, zip(*results.keys())))

    cols = ['CCall', 'CCweak', 'CFOM', 'nsites']
    plt_labels = dict(zip(cols, ['CCall', 'CCweak', 'CFOM', 'No. Found HA sites']))
    ncols = int(ceil(float(len(cols))/rows))
    fig, axarr = plt.subplots(rows, ncols, figsize=(5*ncols, 6))
    for i, sg in enumerate(spacegroups):
        x_vals = []
        y_vals = dict([(col, []) for col in cols])
        width = 0.8 / len(spacegroups)
        offset = i * width
        color = cm.jet(float(i)/len(spacegroups))
        for nsite in nsite_set:
            for resol in resol_set:
                x_vals.append(' | '.join([str(nsite), '%.2f'%resol]))
                try:
                    vals = results[(sg, nsite, resol)]
                except:
                    continue
                for col, y_val in y_vals.iteritems():
                    y_val.append(vals[col])

        x_pos = [x + offset for x in range(len(x_vals))]
        for idx, col in enumerate(cols):
            x, y = divmod(idx, rows)
            axarr[y,x].bar(x_pos, y_vals[col], width=width, color=color, label=sg, alpha=0.75)
            axarr[y,x].set_title(plt_labels[col], fontsize=12)
            if y == 1 and i == len(spacegroups)/2:
                axarr[y,x].set_xticks(x_pos)
                axarr[y,x].set_xticklabels(x_vals, rotation='vertical')
                axarr[y,x].set_xlabel('HA Sites | Resolution')
            elif y == 0:
                axarr[y,x].set_xticks([])
    lgd =  axarr[0,1].legend(bbox_to_anchor=(1.02, 1), loc=2)
    hist_plot = os.path.join(pth, 'hist_shelxd_cc.png')
    plt.savefig(hist_plot, bbox_extra_artists=(lgd,), bbox_inches='tight')


def plot_shelxe_contrast(shelxe_contrast, png_file, add_legend=False):
    '''Plot contrast vs. cycle number from SHELXE.'''

    fig, ax = plt.subplots(figsize=(6, 4))
    for l, solvent_fraction in enumerate(sorted(shelxe_contrast.keys()), 1):
        contrast_vals = shelxe_contrast[solvent_fraction]
        contrast_orig, contrast_other = contrast_vals['original'], contrast_vals['inverted']
        cycles = [j + 1 for j in range(len(contrast_orig))]

        lb_orig = 'Orig. {}'.format(solvent_fraction)
        lb_inv = 'Inv. {}'.format(solvent_fraction)
        color = cm.Paired(float(l)/12)
        ax.plot(cycles, contrast_orig, lw=1, c=color, label=lb_orig)
        ax.plot(cycles, contrast_other, ls='dashed', label=lb_inv, lw=1, c=color)

    plt.xlabel('Cycle', fontsize=14)
    plt.ylabel('Contrast', fontsize=14)
    plt.tick_params(labelsize=14)
    if add_legend:
        lgd = ax.legend(bbox_to_anchor=[.2, 1.02], loc=3, ncol=2)
        plt.savefig(png_file, bbox_extra_artists=(lgd,), bbox_inches='tight')
    else:
        plt.savefig(png_file, bbox_inches='tight')

    plt.close()


def plot_shelxe_fom_mapcc(fom_mapcc, png_file):
    '''Plot contrast vs. cycle number from SHELXE.'''

    fig, (ax1, ax2) = plt.subplots(2, figsize=(6, 4), sharex=True)
    for l, solvent_fraction in enumerate(sorted(fom_mapcc.keys()), 1):
        vals = fom_mapcc[solvent_fraction]
        orig, other = vals['original'], vals['inverted']

        lb_orig = 'Orig. {}'.format(solvent_fraction)
        lb_inv = 'Inv.'
        color = cm.Paired(float(l)/12)
        x = range(len(orig['resol']))
        ax1.plot(x, orig['fom'], lw=1, label=lb_orig, c=color)
        ax1.plot(x, other['fom'], ls='dashed', label=lb_inv, lw=1, c=color)
        ax2.plot(x, orig['mapcc'], lw=1, label=lb_orig, c=color)
        ax2.plot(x, other['mapcc'], ls='dashed', label=lb_inv, lw=1, c=color)

    plt.xlabel('Resolution / $\mathsf{\AA}$')
    ax1.set_ylabel('<FOM>')
    ax2.set_ylabel('<mapCC>')
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p: orig['resol'][p]))
    lgd = ax1.legend(bbox_to_anchor=[1.02, 1.], loc=2, ncol=2, fontsize=10)

    plt.savefig(png_file, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def plot_shelxe_mean_fom_cc(mean_fom_cc, png_file):
    '''Plot <FOM> and <mapCC>  vs. Resolution from SHELXE'''

    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax2 = ax1.twinx()
    solv_axis = sorted(mean_fom_cc.keys())
    fom_orig, fom_inv, mapcc_orig, mapcc_inv = [[mean_fom_cc[solv][hand][stat] for solv in solv_axis]
                                                for stat, hand in product(['mean_fom', 'pseudo_cc'],
                                                                          ['original', 'inverted'])]
    x = range(len(solv_axis))
    ax1.plot(x, fom_orig, lw=1, label='Est.<FOM> Orig.', c='r')
    ax1.plot(x, fom_inv, ls='dashed', label='Inv.', lw=1, c='r')
    ax2.plot(x, mapcc_orig, lw=1, label='pseudo-free CC Orig.', c='b')
    ax2.plot(x, mapcc_inv, ls='dashed', label='Inv.', lw=1, c='b')

    ax1.set_xlabel('Solvent content', fontsize=14)
    ax1.set_ylabel('Est.<FOM>', fontsize=14)
    ax2.set_ylabel('pseudo-free CC', fontsize=14)

    plt.xticks(x, solv_axis)
    ax1.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    lgd1 = ax1.legend(h1 + h2, l1 + l2, bbox_to_anchor=[.15, 1.02], loc=3, ncol=2)

    plt.savefig(png_file, bbox_extra_artists=(lgd1, ), bbox_inches='tight')
    plt.close()
