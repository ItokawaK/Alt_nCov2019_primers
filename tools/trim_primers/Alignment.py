import re

class Alignment():
    def __init__(self, sam_line):
        sam_line_split = sam_line.split('\t')
        self.read_name = sam_line_split[0]
        self.flag = int(sam_line_split[1])
        self.ref_name = sam_line_split[2]
        self.ref_start = int(sam_line_split[3]) - 1 # SAM has 1-based coordinate
        self.mapq = int(sam_line_split[4])
        self.cigar = sam_line_split[5]
        self.next_ref_name = sam_line_split[6]
        self.next_ref_start = int(sam_line_split[7]) - 1
        self.fragment_len = int(sam_line_split[8])
        self.read_seq = sam_line_split[9]
        self.read_q = sam_line_split[10]
        self.opt_tags = sam_line_split[11:]

        cigar_sub = [
                (int(tmp[:-1]), tmp[-1])
                  for tmp in re.findall(r'\d+\w', self.cigar)
                              ]
        self.cigar_expand = "".join([
                                     tmp[1] * tmp[0]
                                      for tmp in cigar_sub
                                     ])

        if self.flag & 64:
            self.read_order = 0
        if self.flag & 128:
            self.read_order = 1
        if self.flag & 16:
            self.strand = '-'
        else:
            self.strand = '+'

        self.read2refcorr = []

        walk = 0
        for cigletter in self.cigar_expand:
            if cigletter in "S":
                self.read2refcorr.append("S")
            if cigletter in 'MDN=X':
                self.read2refcorr.append(self.ref_start + walk)
                walk += 1

        self.ref_end = self.ref_start + walk

    def get_samline():
        pass
