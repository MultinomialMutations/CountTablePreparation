import argparse
import pickle
from sets import Set
from anno_utils import *

base = "ACTG"

parser = argparse.ArgumentParser()
parser.add_argument('annotations',
                    help='''File with path and discretization method for each
                         annotation track''',
                    type=argparse.FileType('r'))
parser.add_argument('sample',
                    help='''Sample cancer combine''',
                    type=argparse.FileType('r'))
parser.add_argument('subcount',
                    help='''subset count path ''',
                    type=argparse.FileType('r'))
parser.add_argument('output',
                    help='output file',
                    type=argparse.FileType('w'))
args = parser.parse_args()

names, readers, bins = get_readers(args.annotations)
head = names + ["left", "right", "sample_id", "cancer_type",
                "from", "to", "count"]
anno_num = len(readers)


sample_set = Set([])
sample_cancer = {}
sample_array = []
for line in args.sample:
    sample, cancer = line.split()
    sample_set.add(sample)
    sample_cancer[sample] = cancer
    sample_array.append(sample)

summary = {}

for file_path in args.subcount:
    f = open(file_path.strip(), 'r')
    record = pickle.load(f)
    for tp in record:
        if tp not in summary:
            summary[tp] = 0
        summary[tp] += record[tp] 

args.output.write('\t'.join(head) + "\n")
for tp in summary:
    tmp = []
    for j in xrange(anno_num):
        tmp.append(str(readers[j][tp[j]]))

    tmp += [base[tp[anno_num]], base[tp[anno_num+1]]]
    sample_id = sample_array[tp[anno_num + 2]]
    tmp += [sample_id, sample_cancer[sample_id]]
    tmp.append(base[tp[anno_num + 3]])
    tmp.append(base[tp[anno_num + 4]])
    args.output.write("\t".join(tmp + [str(summary[tp])]) + "\n")
