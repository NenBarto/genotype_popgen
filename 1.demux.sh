#!/bin/bash

#A de novo analysis in Stacks proceeds in six major stages. First, reads are demultiplexed and cleaned by the process_radtags program. The next three stages comprise the main Stacks pipeline: building loci (ustacks), creating the catalog of loci (cstacks), and matching against the catalog (sstacks). In the fifth stage, the gstacks program is executed to assemble and merge paired-end contigs, call variant sites in the population and genotypes in each sample. In the final stage, the populations program is executed, which can filter data, calculate population genetics statistics, and export a variety of data formats.

#A reference-based analysis in Stacks proceeds in three major stages. The process_radtags program is executed as in a de novo analysis, then the gstacks program is executed, which will read in aligned reads to assemble loci (and genotype them as in a de novo analysis), followed by the populations program.



