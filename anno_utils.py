from readers import *
import re
num_format = re.compile("^[1-9][0-9]*$")

def get_readers(f):
    cov_num = int(f.readline().strip())
    names = []
    readers = []
    bins = []
    for i in xrange(cov_num):
        path = f.readline().strip()
        name = f.readline().strip()
        anno_type = f.readline().strip()
        cat = f.readline().split()
        if anno_type == "bin" or anno_type == "log":
            cat = [float(x) for x in cat]
        t = SortedBedGraph(path, name, anno_type, cat)
        readers.append(t)
        names.append(name)
        bins.append(t.num)
    return names, readers, bins

class Subs:
    def __init__(self, chr, pos, mut_from, mut_to, sample, cancer):
        self.chr = chr

        #only for mutation files that starts from 1
        self.pos = int(pos) - 1 

        self.mut_from = mut_from
        self.mut_to = mut_to
        self.sample = sample
        self.cancer = cancer


class SubsReader:
    """
    Read mutation file as csv format.
    The mutation file should contain at least: chr start from to sample.id cancer.id
    """

    def __init__(self, f):
        features = ["chr", "start", "from", "to", "sample", "cancer"]
        self.mut_file = f
        self.feature = f.readline().split()
        self.feature_index = [self.feature.index(i) for i in features]


    def next(self):
        t = self.mut_file.readline().split()
        if len(t):
            tp = (t[i] for i in self.feature_index)
            return Subs(*tp)
        else:
            return Subs('z', 0, 0, 0, 0, 0)
