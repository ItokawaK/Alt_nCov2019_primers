#!/usr/bin/env python

try:
    import matplotlib
    matplotlib.use('Agg')
finally:
    import matplotlib.pyplot as plt
    import matplotlib.lines as lines
    import matplotlib.patches as patches

import subprocess
import sys
import os
import pandas as pd
import numpy as np
import re

# samtools depth analysis runner
def samtools_depth(bam_file):
    p1 = subprocess.Popen(['samtools','depth','-a', bam_file],
                       stdout = subprocess.PIPE)
    out = p1.communicate()[0]
    out = [l.split('\t') for l in out.decode().rstrip().split('\n')]
    out_tbl = pd.DataFrame({'POS': [int(i[1]) for i in out],
                            'DEPTH':  [int(i[2]) if int(i[2]) > 0 else  + 0.9
                                      for i in out]
                           })

    return out_tbl

# samtools stats analysis runner
def samtools_stats(bam_file):
    p1 = subprocess.Popen(['samtools','stats', bam_file],
                          stdout = subprocess.PIPE)

    out= p1.communicate()[0]
    total_l = 0

    for l in out.decode().rstrip().split("\n"):
        l = l.split('\t')
        if l[0] != 'SN':
            continue
        if l[1] == 'total length:':
            total_l = int(l[2])
        elif l[1] == 'percentage of properly paired reads (%):':
            per_paired = float(l[2]) / 100
    return((total_l, per_paired))

# samtools mpileup runner
def samtools_mpileup(bam_file, ref_fa):
    def readbase_parser(read_base_str):
        i = 0
        split_bases = []
        while i < len(read_base_str):
            if read_base_str[i] in 'ATGCNatgcn.,*':
                split_bases.append(read_base_str[i])
                i += 1
            elif read_base_str[i] == '^':
                i += 2
                split_bases.append(read_base_str[i])
                i += 1
            elif read_base_str[i] == '$':
                i += 1
            elif read_base_str[i] in '+-':
                indel_len_str = re.findall('^[\+-](\d+)', read_base_str[i:])[0]
                i += 1
                i += len(indel_len_str)
                i += int(indel_len_str)
            else:
                print('Unknown char: ' + read_base_str[i], file = sys.stderr)
                print('Unknown char: ' + read_base_str, file = sys.stderr)
                return split_bases

        return split_bases

    p1 = subprocess.Popen(['samtools','mpileup', '-f', ref_fa, '-ax', bam_file],
                       stdout = subprocess.PIPE)
    out = p1.communicate()[0]
    out = [l.split('\t') for l in out.decode().rstrip().split('\n')]
    #return out
    out_tbl = pd.DataFrame({'POS': [int(i[1]) for i in out],
                            'REF_BASE': [i[2] for i in out],
                            'DEPTH':  [int(i[3]) for i in out],
                            'READ_BASES':  [readbase_parser(i[4]) for i in out]
                           })
    def count_mismtaches(read_bases):
        cnt = 0
        for base in read_bases:
            if base in 'ATGCatgc':
                cnt += 1
        return cnt

    out_tbl['MISMATCHES'] = [count_mismtaches(tmp) for tmp in out_tbl['READ_BASES']]

    return out_tbl

# Adding lines for mismatches in plot
def add_mismatch(tbl, ax, threashold = 0.8):
    xy = []
    for index, row in tbl.iterrows():
        if row['MISMATCHES'] > row['DEPTH'] * threashold:
            x = row['POS']
            y = row['MISMATCHES']
            ax.plot([x,x],[1,y],
                    color = 'r',
                    linewidth=0.5)

# Adding gene boxes
def add_genes(ax):
    gene_tbl = pd.DataFrame(
        {
         'Name':["ORF1ab", "S", "ORF3a" , "E", "M", "ORF6", "ORF7a" , "ORF8", "N", "ORF10"],
         'Start':[265,21562,25392,26244,26522,27201,27393,27893,28273,29557],
         'End':[21555,25384,26220,26472,27191,27387,27759,28259,29533,29674]
        }
                            )

    for index, row in gene_tbl.iterrows():
        g_start = row['Start']
        g_end = row['End']
        g_name = row['Name']
        y1 = 0.2
        y2 = 0.3
        if index % 2 == 1:
            y1 += 0.05
            y2 += 0.05
            text_y = 0.4
        else:
            text_y =0.15

        r = patches.Rectangle(xy=(g_start, y1),
                        width= g_end - g_start,
                        height = y2-y1,
                        fc='yellow',
                        ec='k',
                        zorder=101)
        ax.add_patch(r)
        ax.text(x=(g_start + g_end)/2,
          y=text_y,
          s=g_name,
          ha='center',
          va='center',
          fontsize=3,
          weight = 'bold',
          zorder=102
          )

def plot_depths(tbl, primer_region, ax, meta_data=None, hline=10):

    # Setting x and y labels
    ax.set_ylabel('Depth')
    ax.set_xlabel('Genome position nt')

    # Primer regions
    primer_hash = {}
    for i in range(196):
        tmp = primer_region.loc[i, 3].split('_')
        if tmp[1] not in primer_hash:
            primer_hash[tmp[1]] = {}
            primer_hash[tmp[1]]['start'] = 0
            primer_hash[tmp[1]]['end'] = 0
        if 'LEFT' in tmp[2]:
            primer_hash[tmp[1]]['start'] = primer_region.loc[i, 1]
        else:
            primer_hash[tmp[1]]['end'] = primer_region.loc[i, 2]

    # Adding vertical strips for highlight
    # for i in range(1,99):
    #   if i%2 == 1:
    #       my_col = 'red'
    #   else:
    #       my_col = 'green'
    #
    #   ax.axvspan(primer_hash[str(i)]['start'],
    #               primer_hash[str(i)]['end'],
    #               alpha=0.2,
    #               color=my_col
    #             )

    # Plotting depths

    ax.fill([np.min(tbl['POS'])] + list(tbl['POS']) + [np.max(tbl['POS'])],
            [0.9] + list(tbl['DEPTH'] ) + [0.9],
            'b')

    # Adding targets on the bottom of the plotting area
    shade = patches.Rectangle(xy=(-0, 0.51),
                              width=30000,
                              height=0.48,
                              fc='w',
                              zorder=100)
    ax.add_patch(shade)

    for i in range(1,99):
        base1 = 0.5
        base2 = base1 * np.sqrt(2)
        h1 = (base2 - base1) * 0.85
        h2 = (1- base2) * 0.85
        if i%2 == 1:
            base = base1
            h = h1
        else:
            base = base2
            h = h2

        r = patches.Rectangle(xy=(primer_hash[str(i)]['start'], base),
                              width=primer_hash[str(i)]['end']-primer_hash[str(i)]['start'],
                              height = h,
                              fc='#42bff5',
                              ec='k',
                              zorder=101)
        ax.add_patch(r)
        ax.text(x=(primer_hash[str(i)]['end'] + primer_hash[str(i)]['start'])/2,
                y=base + 0.05,
                s=str(i),
                ha='center',
                fontsize=3,
                weight = 'bold',
                zorder=102
                )
    ax.axhline(1, color='k',zorder=102)

    # Adding a horizontal line at hline
    ax.axhline(hline,
               color='k',
               linestyle=':',
               linewidth=0.8,
               zorder=103)

    # Ticks
    #ymax = np.ceil(tbl['DEPTH'].max() / 10000) * 10000
    ax.set_xticks([i * 1000 for i in range(1,30)], minor=True)
    ax.set_xticks([i * 5000 for i in range(7)])
    ax.set_xticklabels([str(i * 5000) for i in range(7)], fontsize='8')

    ax.set_yscale('log')
    ax.set_ylim(0.12, 10000)

    ax.set_yticks([i * ii  for ii in (1,10,100,1000) for i in range(1,11)],
                  minor=True)
    ax.set_yticks([1, 10, 100, 1000, 10000])
    ax.set_yticklabels([str(i) for i in (1,10,100,1000,10000)], fontsize='8')


    # Some information at the right side of figure

    if meta_data != None:
        x,y = np.array([[1, 1.04], [0.85, 0.85]])
        h_offset = 0.01
        for i in range(len(meta_data)):
          ax.text(x[0] + h_offset,
                  y[0] -0.04- (0.04 * i),
                  '{}: {}'.format(meta_data[i][0], meta_data[i][1]),
                  fontsize='6',
                  clip_on=False,
                  transform=ax.transAxes)

def main(bam_files, outpdf, primer_bed, add_mismatches=False, fa_file = None):

    if add_mismatches and fa_file is None:
        sys.exit('Requires the reference fasta')

    n_sample = len(bam_files)
    primer_region = pd.read_csv(primer_bed, sep='\t', header = None)

    fig_hi = 3 * n_sample
    fig = plt.figure(figsize=(8, fig_hi))
    plt.subplots_adjust(right=0.85,
                        hspace=0.5,
                        bottom=0.5/fig_hi,
                        top=1-0.5/fig_hi)

    for i in range(n_sample):
        if add_mismatches:
            tbl = samtools_mpileup(bam_files[i], fa_file)
        else:
            tbl = samtools_depth(bam_files[i])
        align_stats = samtools_stats(bam_files[i])
        meta_data = [('Total Seq.', '{:.1f} Mb'.format(align_stats[0]/1e6)),
                     ('Paired properly', '{:.1%} '.format(align_stats[1]))]
        title = os.path.basename(bam_files[i])
        ax = fig.add_subplot(n_sample, 1, i+1)
        ax.set_title(title)
        plot_depths(tbl,
                    primer_region,
                    ax,
                    meta_data = meta_data,
                    hline=10)
        add_genes(ax)
        if add_mismatches:
            add_mismatch(tbl, ax, threashold = 0.8)

    plt.savefig(outpdf, format='pdf')

if __name__=='__main__':
    import argparse
    import sys
    import os

    _version = 0.6

    parser = argparse.ArgumentParser(description='Output depth plot in PDF. Ver: {}'.format(_version))
    parser.add_argument('-i',
                        '--bams',
                        nargs='*',
                        help='Paths for input BAMs')
    parser.add_argument('-p',
                        '--primer',
                        help='primer_region in BED format')
    parser.add_argument('-o',
                        '--out',
                        help='Output PDF file name')
    parser.add_argument('-r',
                        '--ref_fa',
                        help='Reference fasta file [optional]')

    args = parser.parse_args()

    if not args.out:
        sys.exit('-o (--out) option is mandate')
    if not args.bams:
        sys.exit('-i (--bams) option is mandate')
    if not args.primer:
        sys.exit('-p (--primer) option is mandate')

    if args.ref_fa:
        add_mismatches = True
    else:
        add_mismatches = False

    for file in args.bams:
        if not os.path.isfile(file):
            sys.exit('{} was not found'.format(file))

    if not os.path.isfile(args.primer):
        sys.exit('{} was not found'.format(args.primer))

    main(args.bams,
         args.out,
         args.primer,
         add_mismatches=add_mismatches,
         fa_file=args.ref_fa)
