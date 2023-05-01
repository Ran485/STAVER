#!/bin/bash


# python  ./staver_pipeline.py \
#         --thread_numbers < The CPU worker numbers, Default to nmax-2 > \
#         --input < The DIA data input directory> \
#         --output_peptide < The processed DIA peptide data output directory > \
#         --output_protein < The processed DIA protein data output directory > \
#         --count_cutoff_same_libs < Default to 1 > \
#         --count_cutoff_diff_libs < Default to 2 > \
#         --proteins_cv_thresh < Default to 0.3 > \
#         --na_threshold < Default to 0.3 > \
#         --top_precursor_ions < Default to 3 > \
#         --file_suffix < Default to "_F1_R1" >  \


python  ./staver_pipeline.py \
        --thread_numbers 16 \
        --input ../test_data/ \
        --output_peptide ../test_dataresults/DIA_shicuan/peptides/ \
        --output_protein /Volumes/T7_Shield/staver/results/DIA_shicuan/proteins/ \
        --count_cutoff_same_libs 1 \
        --count_cutoff_diff_libs 2 \
        --peptides_cv_thresh 0.3 \
        --proteins_cv_thresh 0.3 \
        --na_threshold 0.3 \
        --top_precursor_ions 3 \
        --file_suffix _F1_R2 \
