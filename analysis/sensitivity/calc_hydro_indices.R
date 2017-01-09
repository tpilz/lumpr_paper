# This file calc_hydro_indices.R is part of the analysis scripts for the paper Pilz et al. (2017), GMD
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



# Script calculates selected streamflow indices (and stores them in a compressed R file) and gives an overview plot.
# For indices see declaration of function hydInd() below.

library(hydrostats)
library(hydroTSM)
library(xts)
library(ggplot2)
library(plyr)
library(dplyr)

#setwd("/mnt/scratch/users/stud06/tpilz/LUMP_paper/analysis/sensitivity/")
setwd("/home/tobias/Promotion/Paper/lumpR/analysis/sensitivity/")

# R object with streamflow data
dat <- "../collect_results/alldat_Bengue_inflow_m3s_xts.Rdata"

# R object with precipitation time series
dat_pr <- "../collect_results/alldat_precip_sum_hm3_xts.Rdata"

# time frame that shall be used for analysis (e.g., warm up period can be excluded)
period <- c(as.Date("2006-01-01"), as.Date("2013-12-31"))

# filename of histogram plots for an overview over hydrological index realisation
out_hist_file <- "hist_hydIndices_overview.png"

# filename for R object containing hydrological index values
out_r_file <- "hydIndices_realisations.Rdata"


# PARAMETERS #
# all flow values <= this value [m3/s] will be treated as zero flow
thresh_zero <- 0.01

# the quantile of all non-zero flow values (i.e. all values from all realisations) used to define high flood events
# i.e. the derived value of this quantile is used as threshold for defining high flow events (parameter threshold of hydrostats::high.spells)
quant_flood_thresh <- 0.9


# PLOT OPTIONS #

# which indices to use and in which order facets should appear (for original names see function declaration hydInd() below)
ind_sel <- c("rc", "p_z", "m_h", "fdc_s", "f_l", "f_h", "r_r", "r_f")

# how indices shall be named in the plot; formatting as with expression() in R (or change label_parsed, see R help)
ind_labs <- as_labeller(c(rc="RR~('%')", m_h="Q[avmax]~(m^3~s^-1)", p_z="P[flow]~('%')", f_h="f[high]~(year^-1)", f_l="f[low]~(year^-1)",
                          r_r="atop(RC[rise],(m^3~s^-1~day^-1))", r_f="atop(RC[fall],(m^3~s^-1~day^-1))", fdc_s="SFDC~('-')"), 
                        label_parsed)



# EXTERNAL FUNCTION #
# function for calculation of hydrological indices to be applied over each time series (takes time series as xts objects with one column)
hydInd <- function(dat.q.ts, dat.pr.ts, na.rm = F, ignore.zeros = T, thresh.zero, flood.thresh) {

  # remove NAs (if TRUE)
  if(na.rm && any(is.na(dat.q.ts)))
    dat.q.ts <- dat.q.ts[!is.na(dat.q.ts)]
  
  # remove zeros if desired (only for certain quantile-based calculations)
  if(ignore.zeros)
    dat.q.ts.nozero <- dat.q.ts[!(dat.q.ts <= thresh.zero)]
  else
    dat.q.ts.nozero <- dat.q.ts
  
  # data.frame object required for external hydrostats-functions
  df.hydrostats <- data.frame(Q=as.numeric(dat.q.ts), Date=index(dat.q.ts))
  
  # calculate hydrostats-based statistics (relevant information assigned later on)
  low.hydrostats <- low.spells(df.hydrostats, plot = F)
  high.hydrostats <- high.spells(df.hydrostats, threshold=flood.thresh, ignore.zeros = ignore.zeros, ctf.threshold = thresh.zero, plot = F)
  #colwells.hydrostats <- Colwells(data.frame(Q=as.numeric(dat.q.ts.nozero), Date=as.Date(names(dat.q.ts.nozero))))
  
  
  
  # MAGNITUDE
  # runoff ratio over the whole given period [%]
  # calculate sums of flows in hm3 per day
  dat.q.sum <- dat.q.ts * 86400 / 1e6
  # calculate runoff ratio over the whole period
  rc <- sum(dat.q.sum) / sum(dat.pr.ts) * 100
  
  # exceedance probability of zero flows [%]
  p_z <- length(dat.q.ts.nozero)/length(dat.q.ts) * 100
  
  # high flow: average annual maximum flow [V/T]
  m_h <- high.hydrostats$avg.max.ann
  
  
  # FLOW REGIME
  # average slope of flow duration curve for medium range (33% to 66%) of non-zero flows
  # the higher the value the more variable the flow regime; see Sawicz et al. (2011)
  quants <- quantile(dat.q.ts.nozero*100, probs=c(1-0.33,1-0.66))
  fdc_s <- (log(quants[1]) - log(quants[2])) / (0.66-0.33)
  

  # FREQUENCY
  # low flow: average number of low flow events per year [1/year]
  f_l <- low.hydrostats$low.spell.freq
  # high flow: average number of high spell events per year [1/year]
  f_h <- high.hydrostats$spell.freq
  
  
  # CONCENTRATION
  # Rate of change in flow events
  # averge absolute daily change during rise periods within high flow events [V/T/day]
  r_r <- high.hydrostats$avg.rise
  r_f <- high.hydrostats$avg.fall
  
  
#   # NOT USED IN PAPER
#
#   # streamflow elasticity
#   # sensitivity of a catchment's streamflow response to changes in precipitation at annual scale
#   #   - median( dQ/dP * P/Q), where dQ (dP) is difference between previous year's streamflow (precip) to current year's streamflow (precip)
#   #   - value of one: 1% change in precip leads to 1% change in streamflow; greater (less) than one: catchment is elastic (inelastic), i.e. (in)sensitive to a change of precipitation
#   #   - see Sawicz et al. (2011)
#   dat.q.sumy <- apply.yearly(dat.q.sum, sum)
#   dat.pr.sumy <- apply.yearly(dat.pr.ts, sum)
#   elast <- median( (diff(dat.q.sumy) / diff(dat.pr.sumy)) * 1/rc, na.rm = T )
#   
#   # Magnitude of flow events
#   # average flow: ratio of 25th/75th percentiles of (non-zero) daily flows over all years [-]
#   m_a <- quantile(dat.q.ts.nozero, 0.25) / quantile(dat.q.ts.nozero, 0.75)
#   # low flow: average range of monthly minimum flows in a year [V/T]
#   tmp <- aggregate(as.numeric(dat.q.ts), list(year=as.numeric(format(index(dat.q.ts), "%Y")),
#                                               month=as.numeric(format(index(dat.q.ts), "%m"))), FUN = min)
#   tmp <- tapply(tmp$x, list(tmp$year), function(x) abs(diff(range(x))))
#   m_l <- mean(tmp)
#   
#   # Duration of flow events
#   # low flow: average duration of low flow spells [days]
#   d_l <- low.hydrostats$avg.low.spell.duration
#   # high flow: average duration of high flow spells [days]
#   d_h <- high.hydrostats$avg.high.spell.duration
#   
#   # Timing of flow events
#   test <- Colwells(df.hydrostats, fn="median")
  
  out <- c(rc, m_h, p_z, f_h, f_l, r_r, r_f, fdc_s)
  names(out) <- c("rc", "m_h", "p_z", "f_h", "f_l", "r_r", "r_f", "fdc_s")
  return(out)
}

### CALCULATIONS ###

# load data
load(dat)
dat_flow <- dat_all
load(dat_pr)
dat_prec <- dat_all

rm(dat_all)
gc()

# make sure the same time frame is used
dat_flow <- dat_flow[seq.Date(period[1], period[2], by="day"),]
dat_prec <- dat_prec[seq.Date(period[1], period[2], by="day"),]

# derive flow value for defining flood events
dat_v <- as.numeric(dat_flow)
dat_v <- dat_v[dat_v > thresh_zero]
flood_thresh <- quantile(dat_v, probs=quant_flood_thresh)
rm(dat_v)
gc()

# calculation of indices
ind <- sapply(colnames(dat_flow), simplify = F, USE.NAMES = T,
              FUN = function(x) hydInd(dat_flow[,x], dat_prec[,x], thresh.zero = thresh_zero, flood.thresh=flood_thresh))

# create data.frame for ggplot
ind_df <- data.frame(value=unlist(ind),
                     realisation=rep(names(ind), each=length(ind[[1]])),
                     index=names(ind[[1]]))
rownames(ind_df) <- NULL

# store as R object for re-use
save(ind_df, file=out_r_file)



# PLOT #

load(out_r_file) # in case only the plot shall be reproduced ...

# produce histogramms of all indices
ind_df <- filter(ind_df, index %in% ind_sel)

gp <- ggplot(data = ind_df, aes(x = value)) +
  geom_histogram() +
  ylab("Number of experiments") +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_blank(), strip.text.x = element_text(size = 18)) +
  facet_grid(. ~ factor(index, levels=ind_sel), scales = "free", labeller = ind_labs)

ggsave(out_hist_file, plot = gp, width=18, height=6, units="in", dpi=300)

