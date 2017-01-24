# This file plot_compare_results.R is part of the analysis scripts for the paper Pilz et al. (2017), GMD
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




# This script produces overview plots of simulated discharge time series:
# - <file_out_inflow_overview>: daily, monthly, and yearly averaged discharge plots; all simulated time series and a observed one (to give an overview; NOT in the paper)
# - <file_out_reslevel_paper>: Simulated reservoir levels with precipitation input and comparison with observed levels
# - <file_out_prec_overview>: Rainfall input for all experiments (daily, monthly, yearly; NOT in the paper)
# - <file_out_prec_diffs>: Absolute maximum rainfall difference over all experiments for a certain time step (daily, monthly, yearly)
# - <file_out_prec_rc_compare>: Yearly precipitation sums in comparison to values reported by Medeiros & de Araújo (2014) and yearly runoff coefficients in comparison to values reported by Medeiros & de Araújo (2014)
#
# ATTENTION: Script has been originally run on a high performance cluster. Calculations consume a lot of memory!

library(xts)
library(ggplot2)
library(broom)
library(plyr)
library(dplyr)
library(tidyr)
library(gridExtra)
#library(doMC)
#registerDoMC(cores=32)

# working dir
setwd("/mnt/scratch/users/stud06/tpilz/LUMP_paper/analysis/comparison_results/")

# R object with collected simulation data
dat <- "../collect_results/alldat_Bengue_inflow_m3s_xts.Rdata"

# R object with precipitation data
dat_pr <- "../collect_results/alldat_precip_sum_hm3_xts.Rdata"

# R object with collected simulated reservoir filling volume data
dat_level <- "../collect_results/alldat_Bengue_reslevel_hm3_xts.Rdata"


# output files #
# Bengue inflows: overview
file_out_inflow_overview <- "Bengue_inflow_overview.png"
file_out_reslevel_paper <- "Bengue_resvol_overview.png"
file_out_prec_overview <- "rainfall_area-weighted_overview.png"
file_out_prec_diffs <- "rainfall_maxabsdiff_experiments_overview.png"
file_out_prec_rc_compare <- "Bengue_prec_rc_compare.pdf"

file_out <- "Bengue_inflow"

file_out_vol <- "Bengue_resvol"

# observed data file
file_obs <- "./lista_aportes_diario_177.csv"

# time frame that shall be used for analysis (e.g., warm up period can be excluded); only applied to certain plots
period <- c(as.Date("2006-01-01"), as.Date("2013-12-31"))

### CALCULATIONS ###

OVERVIEW OF BENGUE INFLOWS #
load data
load(dat)

dat_obs <- read.table(file_obs, header=T, sep=";", skip = 4)
dat_obs_xts <- xts(dat_obs$Aporte.Total..hm..*10^6/86400, as.Date(dat_obs$Data))
colnames(dat_obs_xts) <- "obs"
rm(dat_obs)

# merge observed and simulated data
dat_all <- merge(dat_all, dat_obs_xts, all=F)
rm(dat_obs_xts)

# time series of monthly and yearly means
dat_all_m <- apply.monthly(dat_all, mean, na.rm=F)
dat_all_y <- apply.yearly(dat_all, mean, na.rm=F)

# re-shape xts object for plotting with ggplot2
dat_tidy_d <- cbind(tidy(dat_all), aggregation="daily")
dat_tidy_m <- cbind(tidy(dat_all_m), aggregation="monthly")
dat_tidy_y <- cbind(tidy(dat_all_y), aggregation="yearly")
rm(dat_all, dat_all_m, dat_all_y)
dat_tidy <- rbind(dat_tidy_d, dat_tidy_m, dat_tidy_y)
rm(dat_tidy_d, dat_tidy_m, dat_tidy_y)

# plot general overview of daily, monthly, and yearly discharge time series with observed data
png(file_out_inflow_overview, width=12, height=8, units="in", res=300)

ggplot(dat_tidy, aes(x=index, y=value, group=series)) +
  geom_line(data = subset(dat_tidy, series != "obs"), colour="gray", lwd=0.5) +
  geom_line(data = subset(dat_tidy, series == "obs"), colour="black", lwd=0.5) +
  facet_grid(aggregation ~ ., scale="free_y") +
  scale_x_date("Year", limits = c(as.Date("2001-01-01"), as.Date("2013-12-31")), date_breaks = "2 years", date_minor_breaks = "1 year") +
  ylab(expression(paste("Discharge / ", m^3, s^-1))) +
  theme_bw()

dev.off()

rm(dat_tidy)




# RESERVOIR LEVEL #
load(dat_level)
dat_level <- dat_all[seq.Date(period[1], period[2], by="day"),]

load(dat_pr)
dat_prec <- dat_all[seq.Date(period[1], period[2], by="day"),1] # only one column as area-weighted precipitation sum over the catchment is almost the same (has been checked)
dat_prec <- tidy(dat_prec)
dat_prec$series <- "prec"
rm(dat_all)

dat_obs <- read.table(file_obs, header=T, sep=";", skip = 4)
dat_obs_level_xts <- xts(dat_obs$Volume..hm.., as.Date(dat_obs$Data))
colnames(dat_obs_level_xts) <- "obs"
rm(dat_obs)

dat_all <- merge(dat_level,dat_obs_level_xts, all=F)
rm(dat_level, dat_obs_level_xts)

dat_all <- tidy(dat_all)

dat_all <- rbind(dat_all, dat_prec)

dat_all <- rbind(dat_all, data.frame(index=NA, series="NA", value=as.numeric(NA))) # quick and dirty solution to get the legend as desired

gp <- ggplot(dat_all, aes(x=index, y=value, group=series)) +
  geom_line(data = subset(dat_all, series == "NA"), lwd=0.5, mapping = aes(colour="Catchment precipitation")) +
  geom_bar(data = subset(dat_all, series == "prec"), mapping = aes(colour="Catchment precipitation"), stat = "identity", show.legend = F) +
  geom_line(data = subset(dat_all, !(series %in% c("obs", "prec"))), lwd=0.5, mapping = aes(colour="Simulated reservoir levels")) +
  geom_line(data = subset(dat_all, series == "obs"), lwd=0.5, mapping = aes(colour="Observed reservoir level")) +
  geom_hline(mapping = aes(colour="Maximum reservoir level", yintercept = 19.6), linetype="dashed") +
  scale_x_date("Date", date_breaks = "1 year", date_minor_breaks = "1 year") +
  scale_color_manual(values=c("lightblue", "#F8766D", "black", "gray50")) +
  ylab(expression(paste("Water Volume (", hm^3, ")"))) +
  theme_bw(base_size = 18) +
  theme(legend.position=c(.99,.99), legend.justification = c(.99,.99), legend.key.width = unit(25,"pt"), legend.key = element_blank()) +
  guides(colour = guide_legend(title=NULL, override.aes = list(size=1.2, bg="white")))

ggsave(file_out_reslevel_paper, plot = gp, width=14, height=6, units="in", dpi=300)

rm(dat_all)




# PRECIPITATION #
load(dat_pr)
dat_all <- dat_all * 1.08 # get values in mm by incorporating catchment area via factor: 1 hm3 / 926 km2 = 1e9 l / 0.926e9 m2 = 1.08 mm

# time series of monthly and yearly sums
dat_all_m <- apply.monthly(dat_all, colSums, na.rm=F)
index(dat_all_m) <- index(dat_all_m) - 15 # yearly values shall be plotted at the middle of the respective month (approximately)
dat_all_y <- apply.yearly(dat_all, colSums, na.rm=F)
index(dat_all_y) <- index(dat_all_y) - 365/2 # yearly values shall be plotted at the middle of the respective year (approximately)

# re-shape xts object for plotting with ggplot2
dat_tidy_d <- cbind(tidy(dat_all), aggregation="daily")
dat_tidy_m <- cbind(tidy(dat_all_m), aggregation="monthly")
dat_tidy_y <- cbind(tidy(dat_all_y), aggregation="yearly")
rm(dat_all, dat_all_m, dat_all_y)
dat_tidy <- rbind(dat_tidy_d, dat_tidy_m, dat_tidy_y)
rm(dat_tidy_d, dat_tidy_m, dat_tidy_y)

# plot
gp <- ggplot(dat_tidy, aes(x=index, y=value, group=series)) +
  geom_line(colour="gray", lwd=0.5) +
  facet_grid(aggregation ~ ., scale="free_y") +
  scale_x_date("Date", limits = c(as.Date("2001-01-01"), as.Date("2013-12-31")), date_breaks = "2 years", date_minor_breaks = "1 year") +
  ylab("Rainfall (mm)") +
  theme_bw(base_size = 18)
ggsave(file_out_prec_overview, plot = gp, width=14, height=8, units="in", dpi=300)

# calculate and plot absolute differences between experiments
dat_diff <- dat_tidy %>%
  ddply(c("index", "aggregation"), function(x) { # had some problems with .parallel=T
    x %>%
      summarise(diffs = max(abs(diff(value)), na.rm=T))
  })

gp <- ggplot(dat_diff, aes(x=index, y=diffs, group=aggregation)) +
  geom_bar(colour="gray", stat = "identity") +
  facet_grid(aggregation ~ ., scale="free_y") +
  scale_x_date("Date", limits = c(as.Date("2001-01-01"), as.Date("2013-12-31")), date_breaks = "2 years", date_minor_breaks = "1 year") +
  ylab("Maximum absolute difference in rainfall input (mm)") +
  theme_bw(base_size = 18)
ggsave(file_out_prec_diffs, plot = gp, width=14, height=7, units="in", dpi=300)




# RUNOFF RATIOS #
# for comparative study in Bengue catchment see Medeiros & de Araujo (2014)
load(dat_pr)
dat_prec <- dat_all
load(dat)
dat_flow <- dat_all
rm(dat_all)

# calculate sums of flows in hm3 per year
dat_flow_sumy <- apply.yearly(dat_flow*86400/1e6, colSums)
index(dat_flow_sumy) <- as.Date(format(index(dat_flow_sumy), "%Y-01-01"))
dat_flow_tidy <- tidy(dat_flow_sumy)
rm(dat_flow,dat_flow_sumy)

# calculate yearly precip sum in hm3 (daily values already in hm3)
dat_prec_sumy <- apply.yearly(dat_prec, colSums)
index(dat_prec_sumy) <- as.Date(format(index(dat_prec_sumy), "%Y-01-01"))
dat_prec_tidy <- tidy(dat_prec_sumy)
rm(dat_prec,dat_prec_sumy)

# join dfs and calculate runoff ratio for every year and realisation
dat_rc <- inner_join(dat_prec_tidy, dat_flow_tidy, by=c("index", "series"), suffix=c(".prec", ".flow"), all=F) %>%
  mutate(rc=value.flow/value.prec,
         value.prec.mm = value.prec * 1.08)
rm(dat_prec_tidy,dat_flow_tidy)

# fill in values from Medeiros & de Araujo (2014) manually
comparisons <- data.frame(index = as.Date(c("2004-01-01", "2005-01-01", "2006-01-01", "2007-01-01", "2008-01-01", "2009-01-01", "2010-01-01", "2011-01-01")),
						  series = "MedeirosAraujo",
						  rc = c(2.3, 0.1, 0.2, 1.0, 1.9, 1.2, 0.7, 2.3)/100,
						  value.prec.mm = c(1002, 518, 561, 739, 733, 743, 449, 598)
						  )
dat_rc <- rbind.fill(dat_rc, comparisons)

# plot
gp1 <- ggplot(dat_rc, aes(x=index, y=rc*100, group=series)) +
  geom_point(data = subset(dat_rc, series != "MedeirosAraujo"), colour="black") +
  geom_point(data = subset(dat_rc, series == "MedeirosAraujo"), colour="red") +
  scale_y_continuous("Runoff ratio (%)", breaks=seq(0, 10, by=2)) +
  scale_x_date("Year", limits = c(as.Date("2006-01-01"), as.Date("2011-01-01")), date_breaks = "1 year", date_minor_breaks = "1 year", date_labels = "%Y") +
  theme_bw(base_size = 18)
#ggsave(file_out_prec_compare, plot = gp, width=7, height=6, units="in", dpi=300)

gp2 <- ggplot(dat_rc, aes(x=index, y=value.prec.mm, group=series)) +
  geom_point(data = subset(dat_rc, series != "MedeirosAraujo"), colour="black") +
  geom_point(data = subset(dat_rc, series == "MedeirosAraujo"), colour="red") +
  scale_y_continuous("Precipitation sum (mm)") +
  scale_x_date("Year", limits = c(as.Date("2006-01-01"), as.Date("2011-01-01")), date_breaks = "1 year", date_minor_breaks = "1 year", date_labels = "%Y") +
  theme_bw(base_size = 18)
#ggsave(file_out_rc_compare, plot = gp, width=7, height=6, units="in", dpi=300)

pdf(file_out_prec_rc_compare, width=14, height=4)
grid.arrange(gp1, gp2, widths=c(.5,.5), ncol=2, nrow=1)
dev.off()
