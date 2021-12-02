#!/usr/bin/env python3

import sys
import argparse
import re
import gzip

from Alignment import Alignment
from Fragment import Fragment
from Primer_range import Primer_range

parser = argparse.ArgumentParser(description='This tool trims portion of primer from aligned reads. '
                                              'Input name sorted sam or output of the bwa mem directry by PIPE.')
parser.add_argument('primer_bed', help='bed file describing primer coorinates')
parser.add_argument('fastq_1', help='name of output fastq file 1')
parser.add_argument('fastq_2', help='name of output fastq file 2')
parser.add_argument('--gzip', action='store_true', help='gzip output fastq files. Adds .gz extention automatically.')
parser.add_argument('--no_merge', action='store_true', help='Does not merge primer ranges. Not recommended.')
parser.add_argument('--verbose', action='store_true', help='output detail infomation for debugging')


args = parser.parse_args()

BED_FILE = args.primer_bed

primer_range = Primer_range(BED_FILE)
if not args.no_merge:
    primer_range.primer_merge()

# [+strand alignment, -strand alignment]
alignment_bucket = [None, None]

out_buffer1 = []
out_buffer2 = []
cnt = 0
current_read = 'default' # Current read name
for sam_line in sys.stdin:
    # Skip header
    if sam_line.startswith('@'):
        continue

    alignment = Alignment(sam_line.rstrip())

    # Initialize alignment backet
    if current_read != alignment.read_name:
        alignment_bucket = [None, None]
        current_read = alignment.read_name

    # Skip non-primary and supplemental alignments
    if alignment.flag & (256 + 2048):
        continue

    if alignment.strand == '+':
        alignment_bucket[0] = alignment
    elif alignment.strand == '-':
        alignment_bucket[1] = alignment

    # process if the bucket is full
    if alignment_bucket[0] and alignment_bucket[1]:

        fragment = Fragment(alignment_bucket[0], alignment_bucket[1])

        # Check if fragment ends are conteined in a primer region
        # Return range if True, None otherwise.
        range_left = primer_range.is_contained(fragment.ref_start, 'left')
        range_right = primer_range.is_contained(fragment.ref_end, 'right')

        # Chop ends of the fragment overlapping to primer
        fragment.slice(range_left)
        fragment.slice(range_right)

        if args.verbose:
            sys.stderr.write(current_read + ':')
            sys.stderr.write('  Fragment interval: {}-{}\n'.format(fragment.ref_start, fragment.ref_end))
            sys.stderr.write('  Left part overlapped with: {}\n'.format(range_left))
            sys.stderr.write('  Right part overlapped with: {}\n'.format(range_right))
            sys.stderr.write('    Left clipped: {}\n'.format(fragment.left_trimmed))
            sys.stderr.write('    Right clipped: {}\n'.format(fragment.right_trimmed))

        # Get fasta string list
        fastq_lines = fragment.get_fastqlines()

        # Create fasta string
        out_read1 = "\n".join(fastq_lines[0]) + "\n"
        out_read2 = "\n".join(fastq_lines[1]) + "\n"

        # Push fasta string to output buffer
        if args.gzip:
            out_read1 = out_read1.encode()
            out_read2 = out_read2.encode()

        out_buffer1.append(out_read1)
        out_buffer2.append(out_read2)

# Write fasta
if args.gzip:
    f1 = gzip.open(args.fastq_1 + '.gz', 'wb', compresslevel=3)
    f2 = gzip.open(args.fastq_2 + '.gz', 'wb', compresslevel=3)
    f1.write(b"".join(out_buffer1))
    f2.write(b"".join(out_buffer2))
else:
    f1 = open(args.fastq_1, 'w')
    f2 = open(args.fastq_2, 'w')
    f1.write("".join(out_buffer1))
    f2.write("".join(out_buffer2))

f1.close()
f2.close()
