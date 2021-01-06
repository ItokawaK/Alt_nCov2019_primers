#!/usr/bin/env python

try:
    import matplotlib
    matplotlib.use('Agg')
finally:
    import matplotlib.pyplot as plt
    import matplotlib.lines as lines
    import matplotlib.patches as patches
    from matplotlib.ticker import NullFormatter

import subprocess
import sys
import os
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import signal
import matplotlib.transforms as transforms

signal.signal(signal.SIGINT, signal.SIG_DFL)

color_scheme = {'plot_fc': 'gray',
                'primerP1_fc':'#ff5353',
                'primerP2_fc':'#5fefc7',
                'gene_fc':'#ffe17f',
                'mismatch_normal':'#ffe17f',
                'mismatch_primer':'red',
                }

# samtools depth analysis runner

def filter_softclipped(bam_file):
    p0 = subprocess.Popen(['samtools','view','-h', bam_file],
                          stdout = subprocess.PIPE)
    p1 = subprocess.Popen(['awk', '/^@/ || $6!~/S/'],
                          stdin = p0.stdout,
                          stdout = subprocess.PIPE)

    return p1

def samtools_depth(bam_file, remove_softclipped=False, min_readlen=0):

    if remove_softclipped:
        p1 = filter_softclipped(bam_file)
        p2 = subprocess.Popen(['samtools','depth','-a', '-l', str(min_readlen),'-'],
                              stdin = p1.stdout,
                              stdout = subprocess.PIPE)
    else:
        p2 = subprocess.Popen(['samtools','depth','-a','-l', str(min_readlen), bam_file],
                              stdout = subprocess.PIPE)

    out = p2.communicate()[0]
    out = [l.split('\t') for l in out.decode().rstrip().split('\n')]
    out_tbl = pd.DataFrame({'POS': [int(i[1]) for i in out],
                            'DEPTH':  [int(i[2]) if int(i[2]) > 0 else  0.9 for i in out]
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

class Mismatch:
    def __init__(self, pos, refbase, altbase, count, type):
        self.pos = pos
        self.refbase = refbase
        self.altbase = altbase
        self.count = count
        self.type = type # snp, del. ins

    @property
    def show_mismatch(self):

        MAX_CHAR = 20
        if self.type == 'snp':
            return f'{self.refbase}{self.pos}{self.altbase}'
        elif self.type == 'del':
            if len(self.altbase) < MAX_CHAR:
                _alt = self.altbase
            else:
                _alt = self.altbase[0:MAX_CHAR] + '...'
            return f'del_{self.pos}{_alt}'
        elif self.type == 'ins':
            if len(self.altbase) < MAX_CHAR:
                _alt = self.altbase
            else:
                _alt = self.altbase[0:MAX_CHAR] + '...'
            return f'ins_{self.pos}{_alt}'
    @property
    def length(self):
        return len(self.altbase)

# samtools mpileup runner
def samtools_mpileup(bam_file, ref_fa, threashold=0.8):

    def readbase_parser(mpileup_line):
        read_base_str = mpileup_line[4]
        pos = int(mpileup_line[1])
        refbase = mpileup_line[2]
        depth = int(mpileup_line[3])
        MIN_NUM_MISMATH = 2
        MIN_MISMATH_RATIO = threashold
        MIN_INDEL_RATIO = threashold * 0.9

        i = 0
        # split_bases = []
        mismatches = defaultdict(int)
        indels = defaultdict(int)
        while i < len(read_base_str):
            first_char = read_base_str[i]
            if first_char in 'n.,*$':
                i += 1
            elif first_char in 'ATGCatgc':
                # split_bases.append(first_char.upper())
                mismatches[first_char.upper()] += 1
                i += 1
            elif first_char == '^':
                i += 2
                # split_bases.append(read_base_str[i])
            # elif read_base_str[i] == '$':
            #     i += 1
            elif read_base_str[i] in '+-':
                indel_type, indel_len_str = re.findall('^([\+-])(\d+)', read_base_str[i:])[0]
                i += 1
                i += len(indel_len_str)
                indel_seq = read_base_str[i:(i + int(indel_len_str))]
                if indel_type == '-':
                    indel_str = 'del_' +  indel_seq.upper()
                elif indel_type == '+':
                    indel_str = 'ins_' +  indel_seq.upper()
                indels[indel_str] += 1
                i += int(indel_len_str)
            else:
                print('Unknown char: ' + read_base_str[i], file = sys.stderr)
                print('Unknown char: ' + read_base_str, file = sys.stderr)
                # return split_bases

        # return split_bases

        out_mismatch = None
        out_indel = None
        if mismatches:
            max_mismatch =  max([(cnt, base) for base, cnt in mismatches.items()])
            mismatch_cnt = max_mismatch[0]
            altbase = max_mismatch[1]
            if (mismatch_cnt >= MIN_NUM_MISMATH) and (mismatch_cnt > (depth * MIN_MISMATH_RATIO)):
                out_mismatch = Mismatch(pos, refbase, altbase, mismatch_cnt, 'snp')
        # else:
        #     out_mismatch = None

        if indels:
            max_indel =  max([(cnt, base) for base, cnt in indels.items()])
            indel_cnt = max_indel[0]
            indel_str = max_indel[1]
            if (indel_cnt >= MIN_NUM_MISMATH) and (indel_cnt > (depth * MIN_INDEL_RATIO)):
                indel_type, indel_seq = indel_str.split("_")
                out_indel = Mismatch(pos+1, refbase, indel_seq, indel_cnt, indel_type)

        #
        # else:
        #     out_indel = None

        return (out_mismatch, out_indel)

    # def count_mismtaches(read_bases):
    #     cnt = 0
    #     for base in read_bases:
    #         if base in 'ATGCatgc':
    #             cnt += 1
    #     return cnt
    #
    # def count_mismtaches2(read_bases):
    #     alt_bases = {'A':0, 'T':0, 'G':0, 'C':0}
    #     indels = defaultdict(int)
    #     for base in read_bases:
    #         if base in 'Aa':
    #             alt_bases['A'] += 1
    #         elif base in 'Tt':
    #             alt_bases['T'] += 1
    #         elif base in 'Gg':
    #             alt_bases['G'] += 1
    #         elif base in 'Cc':
    #             alt_bases['C'] += 1
    #
    #     inverse = [(cnt, base) for base, cnt in alt_bases.items()]
    #     max_cnt, max_base = max(inverse)
    #     return (max_cnt, max_base)

    p1 = subprocess.Popen(['samtools','mpileup', '-f', ref_fa, '-ax', bam_file],
                       stdout = subprocess.PIPE)
    out = p1.communicate()[0]
    out = [l.split('\t') for l in out.decode().rstrip().split('\n')]
    #return out
    out_tbl = pd.DataFrame({'POS': [int(row[1]) for row in out],
                            'REF_BASE': [row[2] for row in out],
                            'DEPTH':  [int(row[3]) if int(row[3]) > 0 else  0.9 for row in out],
                            'READ_BASE':  [row[4] for row in out],
                            'MISMATCHES': [readbase_parser(row) for row in out]
                           })

    # for i, row in out_tbl.iterrows():
    #     max_mistmatch, max_indel = readbase_parser(row.READ_BASE)
    #
    #     if max_mistmatch:
    #         if max_mistmatch[0] > row.DEPTH *

    return out_tbl

# Adding lines for mismatches in plot
# Detecting mismatches contained in primer region
def add_mismatch(tbl,
                 ax,
                 threashold = 0.8,
                 primer_bed=None,
                 seq_name=None,
                 refseq_vector=None):

    # MIN_NUM_MISMATH = 2

    if primer_bed:
        df = pd.read_csv(primer_bed, sep='\t', header = None)
        def is_contained(pos, df):
            for index, row in df.iterrows():
                start = row[1]
                end = row[2]
                # if pos -1 < start:
                #     return False
                if pos -1 >= start and pos <= end:
                    return True

            return False

    xy = []
    count = 0
    for index, row in tbl.iterrows():
        mismatch_obj = row['MISMATCHES'][0]
        indel_obj = row['MISMATCHES'][1]

        if mismatch_obj:
            x = mismatch_obj.pos
            y = mismatch_obj.count
            col = color_scheme['mismatch_normal']
            if primer_bed and is_contained(x, df):
                col = color_scheme['mismatch_primer']
            ax.plot([x,x],[1,y],
                    color = col,
                    linewidth=0.5,
                    zorder= 120)
            # Adding a label for the mismatch base and position
            ax.text(x=x,
                    y=np.sqrt(y) / ((count % 2) + 1),
                    s=mismatch_obj.show_mismatch,
                    fontsize=2,
                    zorder=120,
                    color='#363636')
            count += 1

        if indel_obj:
            x = indel_obj.pos
            y = indel_obj.count
            col = 'blue'
            if primer_bed and is_contained(x, df):
                col = color_scheme['mismatch_primer']

            r = patches.Rectangle(xy=(x, 1),
                                  width=indel_obj.length,
                                  height = y,
                                  fc=col,
                                  ec=col,
                                  linewidth=0.5,
                                  zorder=120)
            ax.add_patch(r)
            # Adding a label for the mismatch base and position
            ax.text(x=x,
                    y=np.sqrt(y) / ((count % 2) + 1),
                    s=indel_obj.show_mismatch,
                    fontsize=2,
                    zorder=120,
                    color='#363636')
            count += 1

    if refseq_vector:
        MIN_DEPTH_CONSENSUS = 10

        _tbl = tbl.sort_values(by='POS', ascending=False)
        for index, row in _tbl.iterrows():
            mismatch_obj = row['MISMATCHES'][0]
            indel_obj = row['MISMATCHES'][1]

            if row.DEPTH < MIN_DEPTH_CONSENSUS:
                try:
                    refseq_vector[row.POS - 1] = 'N'
                except:
                    pass

                continue

            if mismatch_obj:
                idx = (mismatch_obj.pos - 1)
                refseq_vector[idx] = mismatch_obj.altbase
            if indel_obj:
                idx = (indel_obj.pos - 1)
                if indel_obj.type == 'del':
                    del refseq_vector[idx:(idx + len(indel_obj.altbase))]
                if indel_obj.type == 'ins':
                    refseq_vector[idx:idx] = indel_obj.altbase

        fasta_writer(seq_name,''.join(refseq_vector))

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
        ymin = 1/6 * 0.2
        height = 1/6 * 0.15
        if index % 2 == 1:
            ymin += 1/6 *0.05
            text_y = ymin + height
            text_va = 'bottom'
        else:
            text_y = ymin - 1/6 * 0.02
            text_va = 'top'

        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

        r = patches.Rectangle(xy=(g_start, ymin),
                        width= g_end - g_start,
                        height = height,
                        fc=color_scheme['gene_fc'],
                        ec='k',
                        linewidth=0.5,
                        zorder=101,
                        transform=trans)
        ax.add_patch(r)
        ax.text(x=(g_start + g_end)/2,
          y=text_y,
          s=g_name,
          ha='center',
          va=text_va,
          fontsize=3,
          weight='bold',
          zorder=102,
          transform=trans
          )

def add_amplicons(primer_bed, ax, highlights=[], ymax=1/6 *0.92, ymin=1/6*0.55):

    df = pd.read_csv(primer_bed, sep='\t', header = None)

    n_primer = len(df)

    # Primer regions
    amplicon_dict = {}
    for index, row in df.iterrows():
        start = row[1]
        end = row[2]
        primer_name = row[3]
        prefix, amplicom_id, primer_direction, *_ = primer_name.split('_')
        if amplicom_id not in amplicon_dict:
            amplicon_dict[amplicom_id] = {}
            amplicon_dict[amplicom_id]['start'] = 0
            amplicon_dict[amplicom_id]['end'] = 0
        if 'LEFT' in primer_direction:
            amplicon_dict[amplicom_id]['start'] = start
        else:
            amplicon_dict[amplicom_id]['end'] = end

    amplicon_ids = amplicon_dict.keys()

    '''
        ______      __ top2
        | 2  |
        |____|      __ base2
            ______  __ top1
            | 3  |
            |____|  __ base1
    '''

    base1 = ymin
    top1 = (ymax + ymin)/2
    base2 = (ymax + ymin)/2
    top2 = ymax

    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

    for id in amplicon_ids:
        id_num = int(id)

        if id_num % 2 == 1:
            base = base1
            top = top1
            h = top1 - base1
            fc = color_scheme['primerP1_fc']
        else:
            base = base2
            top = top2
            h = top2 - base2
            fc = color_scheme['primerP2_fc']

        x1 = amplicon_dict[id]['start']
        x2 = amplicon_dict[id]['end']

        r = patches.Rectangle(xy=(x1, base),
                              width=x2-x1,
                              height=h,
                              linewidth = 0.5,
                              fc=fc,
                              ec='k',
                              zorder=101,
                              transform = trans)
        ax.add_patch(r)
        ax.text(x=(x1+x2)/2,
                y=(base+top)/2,
                s=id,
                ha='center',
                va='center',
                fontsize=3.5,
                weight='bold',
                zorder=102,
                transform= trans)

        if id in highlights:

            r = patches.Rectangle(xy=(x1, 1/6),
                                  width=x2-x1,
                                  height=5/6,
                                  linewidth = 0,
                                  fc=fc,
                                  alpha=0.5,
                                  zorder=50,
                                  transform = trans)
            ax.add_patch(r)

            if id_num%2==1:
                y = 0.97
            else:
                y = 1.0

            ax.text(
                x=(x1+x2)/2,
                y=y,
                s=id,
                ha='center',
                va='top',
                fontsize=5,
                weight='bold',
                zorder=55,
                transform=trans
            )

def set_plot_area(ax, max_hight=10000):

    # Setting x and y labels
    ax.set_ylabel('Depth')
    ax.set_xlabel('Genome position nt')

    ax.axhline(1, color='k', linewidth=0.5, zorder=102) # line at depth=1

    # Ticks
    #ymax = np.ceil(tbl['DEPTH'].max() / 10000) * 10000

    ax.set_xlim(0, 30000)
    ax.set_xticks([i * 1000 for i in range(1,30)], minor=True)
    ax.set_xticks([i * 5000 for i in range(7)])
    ax.set_xticklabels([str(i * 5000) for i in range(7)], fontsize='8')

    ax.set_yscale('log')
    ymin = 10**-(np.log10(max_hight)/5)
    # For linear scale
    # ymin = -(max_hight)/5

    ax.set_ylim(ymin, max_hight)
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    y_major_ticks = [el for el in ax.get_yticks() if el >=1 and el <= max_hight]
    y_minor_ticks = [el for el in ax.get_yticks(minor=True) if el >1 and el < max_hight]
    ax.set_yticks(y_major_ticks)
    ax.set_yticks(y_minor_ticks, minor=True)
    ax.set_yticklabels([str(int(el)) for el in y_major_ticks])

    # ax.set_ylim(ymin, max_hight)
    #
    # ax.set_yticks([i * ii  for ii in (1,10,100,1000) for i in range(1,11)],
    #               minor=True)
    # ax.set_yticks([1, 10, 100, 1000, 10000])
    # ax.set_yticklabels([str(i) for i in (1,10,100,1000,10000)], fontsize='8')

def fasta_parser(file_path):
    seq = ''
    with open (file_path) as f:
        header = next(f)
        if not header.startswith('>'):
            exit(f'Error: {file_path} does not have a heder line')

        for l in f:
            if not l.startswith('>'):
                seq += l.rstrip()
            else:
                break
    return seq

def fasta_writer(name, inseq):
    print(f'>{name}')
    while inseq:
        print(inseq[0:80])
        inseq = inseq[80:]

def main(bam_files,
         outpdf,
         primer_bed=None,
         highlight_arg=None,
         fa_file=None,
         num_cpu=1,
         mismatches_thresh=0.8,
         remove_softclipped=False,
         min_readlen=0,
         out_consensus=False):

    if fa_file == None:
        with ProcessPoolExecutor(max_workers = num_cpu) as executor:
            executed1 = [executor.submit(samtools_depth, bam, remove_softclipped, min_readlen) for bam in bam_files]
    else:
        with ProcessPoolExecutor(max_workers = num_cpu) as executor:
            executed1 = [executor.submit(samtools_mpileup, bam, fa_file, threashold=mismatches_thresh) for bam in bam_files]

    depth_tbls = [ex.result() for ex in executed1]

    with ProcessPoolExecutor(max_workers = num_cpu) as executor:
        executed2 = [executor.submit(samtools_stats, bam) for bam in bam_files]

    stats = [ex.result() for ex in executed2]

    if highlight_arg:
        highlights=highlight_arg.split(',')
    else:
        highlights=[]

    n_sample = len(bam_files)

    fig_hi = 3 * n_sample
    fig = plt.figure(figsize=(8, fig_hi))
    plt.subplots_adjust(right=0.85,
                        hspace=0.5,
                        bottom=0.5/fig_hi,
                        top=1-0.5/fig_hi)

    if fa_file:
        refseq_str = fasta_parser(fa_file)

    for i in range(n_sample):
        # if add_mismatches:
        #     tbl = samtools_mpileup(bam_files[i], fa_file)
        # else:
        #     tbl = samtools_depth(bam_files[i])

        align_stats = stats[i]
        meta_data = ['Total Seq: {:.1f} Mb'.format(align_stats[0]/1e6),
                     'Paired properly: {:.1%}'.format(align_stats[1])]
        title = os.path.basename(bam_files[i]).rstrip('.bam')
        ax = fig.add_subplot(n_sample, 1, i+1)
        ax.set_title(title)
        set_plot_area(ax, max_hight=10000)
        tbl = depth_tbls[i]

        ax.fill([np.min(tbl['POS'])] + list(tbl['POS']) + [np.max(tbl['POS'])],
                [0.9] + list(tbl['DEPTH']) + [0.9],
                color_scheme['plot_fc'],
                zorder=52)

        # Horizontal line at 10
        ax.axhline(10,
               color='k',
               linestyle=':',
               linewidth=0.8,
               zorder=103)

        # Adding supportive data
        ax.text(1.01, 0.7,
                '\n'.join(meta_data),
                fontsize=5,
                ha='left',
                transform=ax.transAxes)

        if primer_bed != None:
            add_amplicons(primer_bed, ax, highlights=highlights)
        add_genes(ax)
        if fa_file != None:
            refseq_vector = None
            if out_consensus:
                refseq_vector = list(refseq_str)

            add_mismatch(tbl,
                         ax,
                         primer_bed=primer_bed,
                         seq_name = title,
                         refseq_vector=refseq_vector,)

        labels = [item.get_text() for item in ax.get_yticklabels()]

    plt.savefig(outpdf, format='pdf')

if __name__=='__main__':
    import argparse
    import sys
    import os

    _version = 0.9

    parser = argparse.ArgumentParser(description='Output depth plot in PDF. Ver: {}'.format(_version))
    parser.add_argument('-i',
                        '--bams',
                        nargs='*',
                        help='Paths for input BAMs')
    parser.add_argument('-o',
                        '--out',
                        help='Output PDF file name')
    parser.add_argument('-p',
                        '--primer', default=None,
                        help='Primer regions in BED format [optional]')
    parser.add_argument('-l',
                        '--highlights', default=None,
                        help='Add highlights on selected amplicons. '
                             'Give amplicon numbers delimited by comma (e.g. 18,76,...) '
                             'Can only be used with the -p --primer option. [optional]')
    parser.add_argument('-r',
                        '--ref_fa', default=None,
                        help='Reference fasta file [optional]')
    parser.add_argument('-t',
                        '--threads', default=1, type=int,
                        help='Num tasks to process concurrently [optional]')
    parser.add_argument('-m',
                        '--mismatches_thresh', default=0.8, type=float,
                        help='Show mismatches higher than this ratio (default=0.8). '
                             'Only effective with the -r option [optional]')
    parser.add_argument('-s',
                        '--ignore_softclipped', action='store_true',
                        help='Ignore softclipped reads (default=False). [optional]')
    parser.add_argument('--min_readlen', default=0, type=int,
                        help='Minumum length of read (default=0). [optional]')
    parser.add_argument('--out_consensus', action='store_true',
                        help='Output consensus to STDOUT. Experimental.')

    args = parser.parse_args()

    if not args.out:
        sys.exit('-o (--out) option is mandate')

    if not args.bams:
        sys.exit('-i (--bams) option is mandate')

    if args.highlights and not args.primer:
        sys.exit('-l can be used only with the -p (--primer) option')

    for file in args.bams:
        if not os.path.isfile(file):
            sys.exit('{} was not found'.format(file))

    if args.primer and not os.path.isfile(args.primer):
        sys.exit('{} was not found'.format(args.primer))

    if args.ref_fa and not os.path.isfile(args.ref_fa):
        sys.exit('{} was not found'.format(args.ref_fa))

    if args.out_consensus:
        print('WARNIG!!!!: Consensus generation is an experimental function.',
                file=sys.stderr)
        if not args.ref_fa:
            sys.exit('--out_consensus can be used only with the -r (--ref_fa) option')

    main(args.bams,
         args.out,
         primer_bed=args.primer,
         highlight_arg=args.highlights,
         fa_file=args.ref_fa,
         num_cpu=args.threads,
         mismatches_thresh=args.mismatches_thresh,
         remove_softclipped=args.ignore_softclipped,
         min_readlen=args.min_readlen,
         out_consensus=args.out_consensus)
