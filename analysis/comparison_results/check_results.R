# This file check_results.R is part of the analysis scripts for the paper Pilz et al. (2017), GMD
# Copyright (C) 2017 Tobias Pilz
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



# This script is for (simplistic) checking the Bengue reservoir inflow data for completeness and plausible values.
# See command line or check_results.Rout for the results.

library(xts)

# working dir
setwd("/mnt/scratch/users/stud06/tpilz/LUMP_paper/analysis/comparison_results/")

# R object with collected simulation data
dat <- "../collect_results/alldat_Bengue_inflow_m3s_xts.Rdata"

# expected length of data vectors
dat_len <- length(seq.Date(as.Date("2001-01-01"), as.Date("2013-12-31"), by="day"))

# range of accepted data values (values smaller or larger will be reported)
val_range <- c(0,500)

# runs (i.e. list element names) expected to be within the data object
params <- list(
  thresh_sub =  c(1000, # 57 subbasins
                  2000, # 20 subbasins
                  5000, # 10 subbasins
                  10000, # 7 subbasins
                  30000), # 1 subbasin
  eha_thres =  c(25, # 826 EHAs (many are discarded for being too small, having not enough sample points)
                 50, # 1132 EHAs
                 100, # 906 EHAs
                 200, # 530 EHAs
                 500, # 220 EHAs
                 750, # 147 EHAs
                 1000), # 115 EHAs
  no_LUs = c(5,10,20,50,75,100,150,200,250,300),
  no_atts = c(1:7),
  cTC1 = c(1:5)
)
combinations <- expand.grid(params)
run_names <-   dbname <- paste0("lump_paper", "_s", combinations$thresh_sub, "_eha", combinations$eha_thres, 
                                "_lunum", combinations$no_LUs, "_luatts", combinations$no_atts,
                                "_tc", combinations$cTC1) 

### CALCULATIONS ###

# load data
load(dat)

# check run_names within data and store those that could not be found for report
names_na <- run_names[!(run_names %in% names(dat_all))]

# check lengths of data vectors
len_err <- sapply(dat_all, function(x) length(x) != dat_len)
len_err <- which(len_err)

# identify non finite (NA etc.) values
nas <- sapply(dat_all, function(x) which(!is.finite(x)))
nas <- nas[lapply(nas, length) > 0]

# out-of-range values
outofrange_min <- sapply(dat_all, function(x) which(x < val_range[1]))
outofrange_min <- outofrange_min[lapply(outofrange_min, length) > 0]

outofrange_max <- sapply(dat_all, function(x) which(x > val_range[2]))
outofrange_max <- outofrange_max[lapply(outofrange_max, length) > 0]


# REPORT #
message("SHORT REPORT ON DATA CHECKS (analyse R objects for more information):\n", 
        "\n",
        "Runs expected but not within dataset: ", length(names_na), " of ", length(run_names),
        "\n",
        "Runs with length of data vectors not as expected: ", length(len_err), " of ", ncol(dat_all),
        "\n",
        "Runs with non finite values: ", length(nas), " of ", ncol(dat_all),
        "\n",
        "Runs with values smaller expected value range: ", length(outofrange_min), " of ", ncol(dat_all),
        "\n",
        "Runs with values larger expected value range: ", length(outofrange_max), " of ", ncol(dat_all))

if (length(names_na)>0) {
	write.table(data.frame(runs_na=names_na), "runs_missing.log", row.names=F, quote=F, sep="\t")
}
	