# CountTablePreparation

## Introduction
The python scripts 'CountTablePreparation' was written to Prepare count table for a given genomic region and a set of annotations in cancer genome. In order to generate count tables fast enough, we separate regions into small size and run the scripts individually and then combine the counts together. 

In the 'multi_covs_summary', counts are summarized for each possible combination of the levels of genomic covariates. The value of these covariates are discritized in certain methods to prevent countless levels (see 'Discritization methods' for more information). The results for each region are stored as dictionary in a python pickle file.

The function 'combine_counts' takes in a filelist of pickle files and combine the dictionaries. It then generates a tab separated count file summaries the information of all the regions.

## Usage
```
multi_cov_summary.py
usage: multi_cov_summary.py [-h]
                            annotations covered_region mutations sample ref
                            output
```

```
combine_counts.py
usage: combine_counts.py [-h] annotations sample subcount output
```

To use the scripts, you need to first edit an annotation file beginning with the number of annotations and then for each annotation 4 lines indicating the path to the annotation file, the name/abbreviation of the annotation, the discretization methods and the levels. (see 'sample_data/anno.txt').

The covered_region is a sorted .bed file of genomic regions. The regions should be sorted by chromosome and positions. (see 'sample_data/region1.txt')

The mutation file records the mutation calls. Each line represents one mutation for a specific sample. Each record should contain features as 'sample', 'cancer', 'chr', 'start', 'from', 'to'. We now only consider SNVs. (see 'sample_data/mutation.tsv')

The sample file is the map from sample to cancer type. Each line represents one sample. (see 'sample_data/sample_cancer.txt')

The ref file is the human genome annotation in .2bit format. We use hg19.2bit. 

The subcounts file is a list of pickle files. (see sample_output/filelist)

## Discretization methods
In order to prevent huge count tables for each possible value of a given annotation, we use different discretization methods for continuous variables. These can be edited in the annotation files.

'bin': manually set breaks and bin continuous variables to the closest break.

'log': take log value and bin to intergers (used for local mutation rate)

'binary': mark whether a given site is annotated or not

Categorical variables keep there levels in the count tables.

## Scripts used for the analysis in Bertl, J., Guo, Q. et al (2017)
https://figshare.com/projects/MultinomialMutations/24685



### The package is used in: 

Bertl, J.; Guo, Q.; Rasmussen, M. J.; Besenbacher, S; Nielsen, M. M.; Hornshøj, H.; Pedersen, J. S. & Hobolth, A. A Site Specific Model And Analysis Of The Neutral Somatic Mutation Rate In Whole-Genome Cancer Data. bioRxiv, 2017. doi: https://doi.org/10.1101/122879 

Juul, M.; Bertl, J.; Guo, Q.; Nielsen, M. M.; Świtnicki, M.; Hornshøj, H.; Madsen, T.; Hobolth, A. & Pedersen, J. S. Non-coding cancer driver candidates identified with a sample- and position-specific model of the somatic mutation rate. eLife, 2017, 6, e21778. https://doi.org/10.7554/eLife.21778


