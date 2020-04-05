import Alignment

class Fragment():
    '''
    This object represents a molecule inserted between the adapters.
    '''
    def __init__(self, alignment_l, alignment_r):

        self.alignments = (alignment_l, alignment_r)
        self.read_name = alignment_l.read_name

        self.ref_start = alignment_l.ref_start
        self.ref_end = alignment_r.ref_end
        self.size = self.ref_end - self.ref_start
        self.left_trimmed = 0
        self.right_trimmed = 0

    def get_fastqlines(self):
        '''
        Returns a list of sequences and
        qualities of both reads.
        '''

        def revcomp_dna(DNA_SEQ):
            comp_base_dict = {"A":"T", "T":"A",
                              "G":"C", "C":"G",
                              "N":"N"}
            return "".join([comp_base_dict[letter] for letter in DNA_SEQ[::-1]])

        out_list = [None, None]

        for ali in self.alignments:
            if ali.strand == "+":
               seq = ali.read_seq
               qual = ali.read_q
            else:
               seq = revcomp_dna(ali.read_seq)
               qual = ali.read_q[::-1]

            out_list[ali.read_order] = ["@" + self.read_name,
                                        seq,
                                        "+",
                                        qual]
        return out_list

    def slice(self, range):
        '''
        Slices the fragment at the pos.
        cut_pos value represents the 0-based index of the left most base of the
        right part.
        Returns tuple for clipped lengths on read1 and 2.
        '''

        if range == None:
            return 0

        direction = range[2]

        if direction == "+":
            cut_pos = range[1]
        else:
            cut_pos = range[0]

        if not (self.ref_start < cut_pos and cut_pos <= self.ref_end):
            return 0

        def walk_along(alignment, cut_pos):
            walk = 0

            for pos in alignment.read2refcorr:
                if pos == "S":
                    walk += 1
                elif pos >= cut_pos:
                    break
                else:
                    walk += 1
            return walk

        walk0 = walk_along(self.alignments[0], cut_pos)
        walk1 = walk_along(self.alignments[1], cut_pos)

        if direction == "+":
            self.left_trimmed = (walk0, walk1)
        if direction == "-":
            self.right_trimmed = (len(self.alignments[0].read_seq) - walk0,
                                  len(self.alignments[1].read_seq) - walk1)

        if direction == "+":
            self.alignments[0].read_seq = self.alignments[0].read_seq[walk0:]
            self.alignments[0].read_q = self.alignments[0].read_q[walk0:]
            self.alignments[0].read2refcorr = self.alignments[0].read2refcorr[walk0:]
            self.alignments[1].read_seq = self.alignments[1].read_seq[walk1:]
            self.alignments[1].read_q = self.alignments[1].read_q[walk1:]
            self.alignments[1].read2refcorr = self.alignments[1].read2refcorr[walk1:]

        if direction == "-":
            self.alignments[0].read_seq = self.alignments[0].read_seq[:walk0]
            self.alignments[0].read_q = self.alignments[0].read_q[:walk0]
            self.alignments[0].read2refcorr = self.alignments[0].read2refcorr[:walk0]
            self.alignments[1].read_seq = self.alignments[1].read_seq[:walk1]
            self.alignments[1].read_q = self.alignments[1].read_q[:walk1]
            self.alignments[1].read2refcorr = self.alignments[1].read2refcorr[:walk1]
