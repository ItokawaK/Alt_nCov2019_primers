class Primer_range():
    def __init__(self, BED_FILE):
        self.ranges_fwd = []
        self.ranges_rev = []
        with open(BED_FILE, 'r') as bed_h:
            for bed_line in bed_h:
                bed_line = bed_line.rstrip()
                bed_fields = bed_line.split('\t')
                start = int(bed_fields[1])
                end = int(bed_fields[2])
                strand = bed_fields[5]
                if strand == "+":
                    self.ranges_fwd.append((start, end, strand))
                if strand == "-":
                    self.ranges_rev.append((start, end, strand))

    def is_contained(self, end_pos, end_type):
        '''
        This function return range of primer region containing a given position.
        Returns None if there is no containing interval.
        '''
        if end_type == "left":
            for range in self.ranges_fwd:

                if range[0] > end_pos:
                    return None
                if range[0] <= end_pos and end_pos < range[1]:
                    return range
        if end_type == "right":
            for range in self.ranges_rev:
                if range[0] > end_pos:
                    return None
                if range[0] < end_pos and end_pos <= range[1]:
                    return range
