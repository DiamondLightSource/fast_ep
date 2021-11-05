#!/usr/bin/env python
#
# fast_ep_plots ->
#
# Plotting routines for SHELX results
#

import codecs
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


def plot_anom_shelxc(resol, isig, dsig, chi2, cc12, png_file):
    '''Plot <I/sig> and <d"/sig> vs. resolution from SHELXC'''

    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax2 = ax1.twinx()
    if chi2:
        ax3 = ax1.twinx()
        ax3.spines["right"].set_position(("axes", 1.12))
    if cc12:
        ax4 = ax1.twinx()
        ax4.spines["right"].set_position(("axes", 1.25))

    x = range(len(resol))
    plt1, = ax1.plot(x, isig, lw=1, label='$\mathregular{<I/\sigma>}$', c='r')
    plt2, = ax2.plot(x, dsig, lw=1, label='$\mathregular{<d^{\prime\prime}/\sigma>}$', c='b')
    ax1.set_xlabel('Resolution / $\mathregular{\AA}$', fontsize=14)
    ax1.set_ylabel('$\mathregular{<I/\sigma>}$', fontsize=14, color=plt1.get_color(), labelpad=2)
    ax2.set_ylabel('$\mathregular{<d^{\prime\prime}/\sigma>}$', fontsize=14, color=plt2.get_color(), labelpad=2)

    plt.xticks(x, resol)
    ax1.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ha = h1 + h2
    la = l1 + l2
    bbox_to_anchor, ncol = [.25, 1.02], 2

    try:
        plt3, = ax3.plot(x, chi2, lw=1, label='$\mathregular{\chi^2}$', c='c')
        plt4, = ax4.plot(x, cc12, lw=1, label='$\mathregular{CC_{anom}}$', c='g')

        ax3.tick_params(labelsize=14)
        ax4.tick_params(labelsize=14)
        ax3.set_ylabel('$\mathregular{\chi^2}$', fontsize=14, color=plt3.get_color(), labelpad=3)
        ax4.set_ylabel('$\mathregular{CC_{anom}}$', fontsize=14, color=plt4.get_color(), labelpad=3)

        h3, l3 = ax3.get_legend_handles_labels()
        h4, l4 = ax4.get_legend_handles_labels()
        ha += h3 + h4
        la += l3 + l4
        bbox_to_anchor, ncol = [.12, 1.02], 4
    except:
        pass

    lgd1 = ax1.legend(ha, la, bbox_to_anchor=bbox_to_anchor, loc=3, ncol=ncol, fontsize=14)
    plt.savefig(png_file, bbox_extra_artists=(lgd1, ), bbox_inches='tight')
    plt.close()


def plot_shelxd_cc(pth, results, spacegroups, png_file):
    '''Summary plot of CCall / CCweak values from SHELXD for all
    space group, heavy atom number and resolution combinations'''

    _, nsites, ano_rlimits = map(sorted, map(set, zip(*list(results.keys()))))

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
            textstr = 'HA : %s  Resol: %.2f$\mathregular{\AA}$' % (nsite, rlimit)
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

    _, nsite_set, resol_set = map(sorted, map(set, zip(*list(results.keys()))))

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
                for col, y_val in y_vals.items():
                    y_val.append(vals[col])

        x_pos = [x + offset for x in range(len(x_vals))]
        for idx, col in enumerate(cols):
            x, y = divmod(idx, rows)
            axarr[y,x].bar(x_pos, y_vals[col], width=width, color=color, edgecolor='k', label=sg, alpha=0.75)
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
        cycles_orig, contrast_orig = contrast_vals['original']
        cycles_other, contrast_other = contrast_vals['inverted']

        lb_orig = 'Orig. {}'.format(solvent_fraction)
        lb_inv = 'Inv. {}'.format(solvent_fraction)
        color = cm.Paired(float(l)/12)
        ax.plot(cycles_orig, contrast_orig, lw=1, c=color, label=lb_orig)
        ax.plot(cycles_other, contrast_other, ls='dashed', label=lb_inv, lw=1, c=color)

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
        try:
            orig, other = vals['original'], vals['inverted']
            x_vals = range(len(orig['resol']))
            x_labels = orig['resol']
        except KeyError:
            continue

        lb_orig = 'Orig. {}'.format(solvent_fraction)
        lb_inv = 'Inv.'
        color = cm.Paired(float(l)/12)
        ax1.plot(x_vals, orig['fom'], lw=1, label=lb_orig, c=color)
        ax1.plot(x_vals, other['fom'], ls='dashed', label=lb_inv, lw=1, c=color)
        ax2.plot(x_vals, orig['mapcc'], lw=1, label=lb_orig, c=color)
        ax2.plot(x_vals, other['mapcc'], ls='dashed', label=lb_inv, lw=1, c=color)

    plt.xlabel('Resolution / $\mathregular{\AA}$')
    ax1.set_ylabel('<FOM>')
    ax2.set_ylabel('<mapCC>')
    ax1.xaxis.set_ticks(x_vals)
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p: x_labels[p]))
    lgd = ax1.legend(bbox_to_anchor=[1.02, 1.], loc=2, ncol=2, fontsize=10)

    plt.savefig(png_file, bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()


def plot_shelxe_mean_fom_cc(mean_fom_cc, png_file):
    '''Plot <FOM> and <mapCC>  vs. Resolution from SHELXE'''

    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax2 = ax1.twinx()
    solv_axis = sorted(mean_fom_cc.keys())
    fom_orig, fom_inv, mapcc_orig, mapcc_inv = [[(idx, mean_fom_cc[solv][hand][stat]) for idx, solv in enumerate(solv_axis)]
                                                for stat, hand in product(['mean_fom', 'pseudo_cc'],
                                                                          ['original', 'inverted'])]
    ax1.plot([x for x,_ in fom_orig],
             [y for _,y in fom_orig],
             lw=1, label='Est.<FOM> Orig.', c='r')
    ax1.plot([x for x,_ in fom_inv],
             [y for _,y in fom_inv],
             ls='dashed', label='Inv.', lw=1, c='r')
    ax2.plot([x for x,_ in mapcc_orig],
             [y for _,y in mapcc_orig],
             lw=1, label='Pseudo-free CC Orig.', c='b')
    ax2.plot([x for x,_ in mapcc_inv],
             [y for _,y in mapcc_inv],
             ls='dashed', label='Inv.', lw=1, c='b')

    ax1.set_xlabel('Solvent content', fontsize=14)
    ax1.set_ylabel('Est.<FOM>', fontsize=14)
    ax2.set_ylabel('Pseudo-free CC', fontsize=14)

    plt.xticks(range(len(solv_axis)), solv_axis)
    ax1.tick_params(labelsize=14)
    ax2.tick_params(labelsize=14)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    lgd1 = ax1.legend(h1 + h2, l1 + l2, bbox_to_anchor=[.15, 1.02], loc=3, ncol=2)

    plt.savefig(png_file, bbox_extra_artists=(lgd1, ), bbox_inches='tight')
    plt.close()


def plot_b64encoder(plot_list):
    '''Encode image to embed into HTML summary'''

    res = {}
    for plt in plot_list:
        with open(plt, 'rb') as fh:
            enc_data = codecs.encode(fh.read(), encoding="base64").decode("ascii")
            tag,_ = os.path.splitext(plt)
            res[tag] = enc_data.replace("\n", "")

    return res
