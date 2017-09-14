import argparse
import pickle
from sets import Set
from anno_utils import *
from bx.seq.twobit import *

base = "ACTG"
parser = argparse.ArgumentParser(
    description='''Counts number of covered positions and the number of
     mutated positions for each combination of different genome annotations
     and sequence contexts''')
parser.add_argument('annotations',
                    help='''File with path and discretization method for each
                         annotation track''',
                    type=argparse.FileType('r'))
parser.add_argument('covered_region',
                    help='''Sorted bed-file with those regions where mutations
                     could be called''', type=argparse.FileType('r'))
parser.add_argument('mutations',
                    help='''Sorted file with for each mutation''',
                    type=argparse.FileType('r'))
parser.add_argument('sample',
                    help='''Sample cancer combine''',
                    type=argparse.FileType('r'))
parser.add_argument('ref',
                    help=''' reference genome''',
                    type=argparse.FileType('r'))
parser.add_argument('output',
                    help='output file',
                    type=argparse.FileType('w'))
args = parser.parse_args()

names, readers, bins = get_readers(args.annotations)
head = names + ["left", "right", "sample_id", "cancer_type",
                "from", "to", "count"]
anno_num = len(readers)


mutation = SubsReader(args.mutations)
genome = TwoBitFile(args.ref)
mut = mutation.next()

sample_set = Set([])
sample_cancer = {}
sample_array = []
for line in args.sample:
    sample, cancer = line.split()
    sample_set.add(sample)
    sample_cancer[sample] = cancer
    sample_array.append(sample)

summary = {} 

for line in args.covered_region:
    chr, start, end = line.split()
    start = int(start)
    end = int(end)
    if chr == "chrX" or chr == "chrY":
        continue
    for pos in xrange(start, end+1):
        # genome reference
        context = genome[chr][pos-1: pos+2].upper()
        if 'N' in context or len(context) != 3:
            continue

        covs = []
        for reader in readers:
            covs.append(reader.query(chr, pos))

        covs = covs + [base.index(context[0]), base.index(context[2])]

        sample_with_mut = Set([])
        while mut.chr < chr or (mut.chr == chr and mut.pos < pos):
            mut = mutation.next()

        while mut.chr == chr and mut.pos == pos:
            sample_with_mut.add(mut.sample)

            covs_with_mut = covs + [sample_array.index(mut.sample),
                                    base.index(mut.mut_from),
                                    base.index(mut.mut_to)]
            if tuple(covs_with_mut) not in summary:
                summary[tuple(covs_with_mut)] = 0
            summary[tuple(covs_with_mut)] += 1
            mut = mutation.next()

        for i, sample in enumerate(sample_array):
            if sample in sample_with_mut:
                continue
            ref = base.index(context[1])
            covs_without_mut = covs + [i, ref, ref]
            if tuple(covs_without_mut) not in summary:
                summary[tuple(covs_without_mut)] = 0
            summary[tuple(covs_without_mut)] += 1

pickle.dump(summary, args.output)
