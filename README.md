# lumpR Paper
This repository contains a collection of Latex files to reproduce the paper of Pilz et al. (2017) (directory __Script/__), that was published in GMD, including R scripts to reproduce (most of) the graphics therein (__analysis/__). The Paper describes the R package _lumpR_ for the hillslope-based discretisation of landscapes to be used in hydrological models.

Reference (and link to GMD):

[Pilz, T., Francke, T., and Bronstert, A.: lumpR 2.0.0: an R package facilitating landscape discretisation for hillslope-based hydrological models, Geosci. Model Dev., 10, 3001-3023, doi: 10.5194/gmd-10-3001-2017, 2017.](https://www.geosci-model-dev.net/10/3001/2017/gmd-10-3001-2017.html)


# lumpR
The R package lumpR is available at https://github.com/tpilz/lumpR

See the README and Wiki therein as well as the documentation in R for more information.


# Notes on analyses scripts
For the paper, an ensemble analysis has been performed with 12,250 parameter realisations. This took several days on a high performance cluster employing 96 CPUs. The template R script used for these experiments has been included (__analysis/model_runs/apply_lumpR_template.R__) but it is not possible to reproduce the experiments on the fly. In addtion, due to unclear legal status, not all of the used input data can be made publicly available.

Moreover, several hundred GB of data has been produced. The raw results are, thus, not included. The results relevant for the analyses, however, have been collected and compressed and are available: __analysis/collect_results/*.Rdata__ (to be opened within R via `load()`).

Scripts to reproduce the reservoir plots can be found in __analysis/comparison_results/__.

Files for the calculation of streamflow indices and for the sensitivity analysis are in __analysis/sensitivity/__.
