# lumpr_paper
This repository contains a collection of Latex files to reproduce the paper of Pilz et al. (2017) (directory __Script/__), that has been submitted to GMD, including R scripts to reproduce (most of) the graphics therein (__analysis/__). The Paper describes the R package _lumpR_ for the hillslope-based discretisation of landscapes to be used in hydrological models.

# lumpR
The R package lumpR is available at https://github.com/tpilz/lumpR

See the README and Wiki therein as well as the documentation in R for more information.

# Review contribution
As the Paper has been submitted to an open access journal, everybody is welcome to contribute to the review process. To do so, follow the standard procedure of paper review or try out the following:

1. Clone the repository
2. Go into the Script/ directory, create a personal copy of the script.tex file
3. Adjust that copy as you wish
4. To make comments I suggest the _todonotes_ package
5. To highlight your contributions from the original text, I suggest the tool _latexdiff_ which creates a new *.tex file which, in turn, you can compile to produce a pdf highlighting the differences in a word-like style
6. Publish your contributions/comments directly at the GMD discussion page

A tutorial on how to collaboratively write a paper with Latex (i.e., how to use _todonotes_ and _latexdiff_ etc.) can be found in directory __howto_review_with_latex/__ (courtesy of Aline Murawski).

# Notes on analyses scripts
For the paper, an ensemble analysis has been performed with 12,250 parameter realisations. This took several days on a high performance cluster employing 96 CPUs. The template R script used for these experiments has been included (__analysis/model_runs/apply_lumpR_template.R__) but it is not possible to reproduce the experiments on the fly. In addtion, due to unclear legal status, not all of the used input data can be made publicly available.

Moreover, several hundred GB of data has been produced. The raw results are, thus, not included. The results relevant for the analyses, however, have been collected and compressed and are available: __analysis/collect_results/*.Rdata__ (to be opened within R via `load()`).

Scripts to reproduce the reservoir plots can be found in __analysis/comparison_results/__.

Files for the calculation of streamflow indices and for the sensitivity analysis are in __analysis/sensitivity/__.
