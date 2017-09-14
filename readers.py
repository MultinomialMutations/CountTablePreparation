import bisect
import math


class BinDiscretizer:
    """
    Discretize with given breaks
    """

    def __init__(self, breaks):
        assert all(breaks[i] < breaks[i + 1] for i in xrange(len(breaks) - 1))
        self.breaks = breaks

    def __getitem__(self, x):
        return bisect.bisect(self.breaks, x) - 1


class SortedBedGraph:
    def __init__(self, path, bedname, anno_type, cat, f_start=0):
        self.name = bedname
        self.f = open(path)
        self.type = anno_type
        if self.type == "binary":
            self.num = 2
            self.cat = [0, 1]
        if self.type == "categorical":
            self.cat = cat
            self.num = len(cat) + 1
        if self.type == "bin":
            self.cat = cat
            self.num = len(cat)
            self.bin = BinDiscretizer(cat)
        if self.type == "log":
            self.cat = cat
            self.num = len(cat) + 1
            self.bin = BinDiscretizer(cat)
        self.f_start = f_start
        self.f.seek(self.f_start)
        t = self.f.readline().split()
        self.bed_chr = t[0]
        self.bed_start = int(t[1])
        self.bed_end = int(t[2])
        self.bed_value = t[-1]

    def __getitem__(self, x):
        if self.type == "binary":
            return x
        if x == self.num - 1:
            return "NA"
        else:
            return self.cat[x]

    def query(self, chr, pos):
        while self.bed_chr < chr or (
                        self.bed_chr == chr and self.bed_end <= pos):
            t = self.f.readline().split()
            if len(t) == 0:
                #print "end of annotation file !"
                if self.type == "binary":
                    return 0
                return self.num - 1
            self.bed_chr = t[0]
            self.bed_start = int(eval(t[1]))
            self.bed_end = int(eval(t[2]))
            self.bed_value = t[-1]

        if (chr == self.bed_chr) and (self.bed_start <= pos < self.bed_end):
            if self.type == "value":
                return float(self.bed_value)
            if self.type == "bin":
                return self.bin[float(self.bed_value)]
            if self.type == "categorical":
                return self.cat.index(self.bed_value)
            if self.type == "binary":
                return 1
            if self.type == "log":
                v = float(self.bed_value)
                if v < 0.0000001:
                    return self.bin[-7]
                else:
                    return self.bin[round(math.log10(v))]
        else:
            if self.type == "binary":
                return 0
            return self.num - 1
