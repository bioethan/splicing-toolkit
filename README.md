# Splicing Toolkit

## Purpose
The goal of the `splicing_toolkit` package is to make performing splicing classification on long read sequencing data simpler and easily repreducible. The current package consists of two functions stored in `utils.py`, but I am working on creating a CLI that allows for easy access to those utility functions, among others that will assist with UTR exploration and polyadenylation.

## Functions
`gtf_to_bed12(path_gtf, out_dir)`
This function takes a path to a GTF file and returns three separate bed files, one for genes, one for exons, and one for introns. It processes the GTF to generate these files for use in the splicing classifier.

`splicing_classifier(LRS_url, output_save)`
This function takes the three bed files generated by `gtf_to_bed12` and the path to the long read sequencing data (`LRS_url`) and generates a file that classifies the splicing status for each of the reads in the LRS data. 

This repo is under construction!!


