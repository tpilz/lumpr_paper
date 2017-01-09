# This file collect_results.R is part of the analysis scripts for the paper Pilz et al. (2017), GMD
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



# RAW DATA NOT ON GITHUB DUE TO SPACE LIMITATIONS! SCRIPT JUST FOR DEMONSTRATIONS / CONFIRMABILITY!
#
# Script collects all model run data and writes them into a compressed R file for later re-use
# Adapted to collect Bengue inflow data in m3/s and store them as xts object; column names are model run names
# The same also for daily rainfall sums over the whole catchment in hm3 and number of actual generated spatial units

library(xts)
library(plyr)
library(dplyr)

# working dir
setwd("/mnt/scratch/users/stud06/tpilz/LUMP_paper/analysis/collect_results/")

# directory with model runs (should contain <model_rundir> -> "output" and "output" with usual WASA input/output sub-directories and files)
model_runs <- "/mnt/scratch/users/stud06/tpilz/LUMP_paper/model_runs/"

# directory with discretisations
discretisations <- "/mnt/scratch/users/stud06/tpilz/LUMP_paper/discretisations/"

# name of R file containing all collected data
name_out_inflow <- "alldat_Bengue_inflow_m3s_xts.Rdata"
name_out_level <- "alldat_Bengue_reslevel_hm3_xts.Rdata"
name_out_precip <- "alldat_precip_sum_hm3_xts.Rdata"

### collect data ###

# reservoir inflow data #
dat_all <- dir(path = model_runs, pattern = "res_1_watbal.out", recursive = T, full.names = T) %>%
  # apply for each model output
  lapply(function(x) {
    # read data
    dat <- read.table(x, header = T) %>%
      # select values and build xts object
      mutate(res = xts(inflow.m..3.s.., as.Date(paste(year., day., sep="-"), format = "%Y-%j"))) %>%
      #mutate(res = xts(volume.m..3./10^6, as.Date(paste(year., day., sep="-"), format = "%Y-%j"))) %>%
      select(res) %>%
      xts(index(.$res))
    # name columns (name of model run)
    names(dat) <- head(tail(unlist(strsplit(x, split = "/")),3), 1)
    return(dat)
  }
  ) %>%
  # xts object with multiple columns instead of list with multiple xts objects
  do.call(cbind, .)

# save as R object
save(dat_all, file=name_out_inflow)

rm(dat_all)
gc()




# reservoir level data #
dat_all <- dir(path = model_runs, pattern = "res_1_watbal.out", recursive = T, full.names = T) %>%
  # apply for each model output
  lapply(function(x) {
    # read data
    dat <- read.table(x, header = T) %>%
      # select values and build xts object
      mutate(res = xts(volume.m..3./10^6, as.Date(paste(year., day., sep="-"), format = "%Y-%j"))) %>%
      select(res) %>%
      xts(index(.$res))
    # name columns (name of model run)
    names(dat) <- head(tail(unlist(strsplit(x, split = "/")),3), 1)
    return(dat)
  }
  ) %>%
  # xts object with multiple columns instead of list with multiple xts objects
  do.call(cbind, .)

# save as R object
save(dat_all, file=name_out_level)

rm(dat_all)
gc()




# Daily rainfall sums over the whole catchment in hm3 #
dat_all <- list.dirs(model_runs, recursive=F) %>%
  # apply for each model directory
  lapply(function(x) {
    # read subbasin information
    dat_sub <- readLines(dir(x, pattern="hymo.dat", recursive=T, full.names=T))[-c(1:2)] %>%
      sapply(function(y) {
        as.numeric(unlist(strsplit(y, "\t"))[c(1:2)])
      })
    dat_sub <- data.frame(subbasin=dat_sub[1,], area=dat_sub[2,])
    row.names(dat_sub) <- NULL
    dat_sub <- arrange(dat_sub, subbasin)
    # read precipitation time series
    dat_pr <- read.table(dir(x, pattern="rain_daily.dat", recursive=T, full.names=T), skip = 2, header=T)
    dat_pr <- xts(dat_pr[,-c(1:2)], as.Date(sprintf("%08d", dat_pr[,1]), format="%d%m%Y"))
    # calculate actual precipitation amount in hm3 for the whole catchment area
    dat_out <- apply(coredata(dat_pr), 1, function(z) sum(z * dat_sub$area * 1e-3))
    dat_out <- xts(dat_out, index(dat_pr))
    colnames(dat_out) <- tail(unlist(strsplit(x, split = "/")),1)
 
    return(dat_out)
  }
  ) %>%
  # xts object with multiple columns instead of list with multiple xts objects
  do.call(cbind, .)

# save as R object
save(dat_all, file=name_out_precip)
