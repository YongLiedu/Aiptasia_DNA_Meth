#!/usr/bin/env python3

"""
> plot_meth_pct_over_norm_gene.py <

Uses seaborn to plot a 2D-kde-histogram of the methylation percentage across 
a normalized genic region.
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import seaborn as sns

WINDOW = 2000

parser = argparse.ArgumentParser(description="""
Uses seaborn to plot a 2D-kde-histogram of the methylation percentage across 
a normalized genic region, TSS and TTS.""")

parser.add_argument('augmented_bismark_cov', metavar="tsv_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="Augmented Bismark file.")
parser.add_argument('--savefig', action='store_true',
                    help="Saves a pdf of the diagram.")

args = parser.parse_args()

# read data
meth_pct = []
genic_pos = []

tsv_reader = csv.reader(args.augmented_bismark_cov, delimiter='\t')
for row in tsv_reader:
    if not row: continue
    
    # check that the meth positions are in genes
    if row[6] not in ['intergenic', 'no_info']:
        meth_pct.append(float(row[3]))
        genic_pos.append(float(row[7]))

    ## TSS - part I: get +1 to +WINDOW
    ## check that the meth positions are in genes
    #if row[6] not in ['intergenic', 'no_info']:
    #    if int(row[8]) < WINDOW:
    #        meth_pct.append(float(row[3]))
    #        genic_pos.append(int(row[8]) + 1)
    
    # TSS - part II: get -1 to -WINDOW
    # check that the meth positions are intergenic regions
    elif row[6] == 'intergenic':
        if row[10] == 'upstream_crick' and int(row[8]) < WINDOW:
            meth_pct.append(float(row[3]))
            genic_pos.append(0.5*(-int(row[8]) - 1)/WINDOW)
        elif row[11] == 'downstream_watson' and abs(int(row[9])) <= WINDOW:
            meth_pct.append(float(row[3]))
            genic_pos.append(0.5*int(row[9])/WINDOW)

    ## TTS - part I: get -WINDOW to -1 (genic region)
    ## check that the meth positions are in genes
    #if row[6] not in ['intergenic', 'no_info']:
    #    if abs(int(row[9])) <= WINDOW:
    #        meth_pct.append(float(row[3]))
    #        genic_pos.append(int(row[9]))
    
    # TTS - part II: get +1 to +WINDOW (intergenic region)
    # check that the meth positions are intergenic regions
    #if row[6] == 'intergenic':
        elif row[10] == 'upstream_watson' and int(row[8]) < WINDOW:
            meth_pct.append(float(row[3]))
            genic_pos.append(0.5*((int(row[8]) + 1)/WINDOW)+1)
        elif row[11] == 'downstream_crick' and abs(int(row[9])) <= WINDOW:
            meth_pct.append(float(row[3]))
            genic_pos.append(0.5*(abs(int(row[9]))/WINDOW)+1)

# plotting with seaborn
# seaborn's default colour palette is 'deep':  
#   ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]
#      blue       green       red      purple     yellow    lightblue
#sns.set(style="white")
#sns.set(style="ticks")
rc={'font.size': 18, 'axes.labelsize': 18, 'legend.fontsize': 15, 
    'axes.titlesize': 18, 'xtick.labelsize': 15, 'ytick.labelsize': 15}
sns.set(rc=rc)
sns.set_style("white")
sns.set_style("ticks")
p = sns.jointplot(np.array(genic_pos), np.array(meth_pct), kind="kde",  
    stat_func=None, color='#4C72B0', size=5, space=0, xlim=(-0.5,1.5), ylim=(0,100))
# for hex plots
# p = sns.jointplot(np.array(genic_pos), np.array(meth_pct), kind="hex", size=10, 
    # xlim=(0,1), ylim=(0,100), space=0, #joint_kws={'gridsize': 50}, 
    # marginal_kws={'bins': 50})        # for tweaking bin numbers
p.set_axis_labels('Normalized gene length', 'Methylation %')
p.ax_joint.set_xticks([-0.5, 0, 1, 1.5])
p.ax_joint.set_xticklabels(['-2000', 'TSS', 'TTS', '2000'])


if args.savefig:
    fig = plt.gcf()
    
    # without bbox_inches, the saved figure has truncated axes.
    output_filename = args.augmented_bismark_cov.name[:-4] + '.wholeGene.pdf'
    fig.savefig(output_filename, bbox_inches='tight', format='pdf')
