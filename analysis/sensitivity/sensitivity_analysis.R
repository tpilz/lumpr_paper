# This file sensitivity_analysis.R is part of the analysis scripts for the paper Pilz et al. (2017), GMD
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



# This script does the sensitivity analysis based on streamflow index data produced with calc_hydro_indices.R.
# As results it produces a number of plots; see below.
# 
# This script includes two approaches:
# - density-based approach: independent of statistical moments and better suited in case of skewed value
#   distributions; PAWN method by Pianosi & Wagner (2015); compares unconditional and conditional
#   cumulated density functions by Kolmogorov-Smirnov statistic
# - variance-based approach (NOT IN PAPER): output distributions for some index/input factor combinations multi-modal
#   and/or highly skewed (see plots) which violates assumption that output variance is sufficient to
#   fully characterise output uncertainty

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(doMC)
registerDoMC(cores=8)

#setwd("/mnt/scratch/users/stud06/tpilz/LUMP_paper/analysis/sensitivity/") 
setwd("/home/tobias/Promotion/Paper/lumpR/analysis/sensitivity/")

# R object with data.frame of hydrological indices per realisation
dat_file <- "./hydIndices_realisations.Rdata"

# how discretisation parameters shall be named in the plot
par_labs <- c(s="SUB_thresh", eha="EHA_thresh", lunum="LU_no", tc="TC_no", luatts="LU_atts")

# which indices/parameters to use and in which order facets should appear (for original names see raw data)
ind_sel <- c("rc", "p_z", "m_h", "fdc_s", "f_l", "f_h", "r_r", "r_f")
par_sel <- c("SUB_thresh", "EHA_thresh", "LU_no", "LU_atts", "TC_no")

# how indices shall be named in the plot; formatting as with expression() in R (or change label_parsed, see R help)
ind_labs <- as_labeller(c(rc="RR", m_h="Q[avmax]", p_z="P[flow]", f_h="f[high]", f_l="f[low]", r_r="RC[rise]", r_f="RC[fall]", fdc_s="SFDC"), 
                        label_parsed)
# same as ind_labs but with units
ind_labs_units <- as_labeller(c(rc="RR~('%')", m_h="Q[avmax]~(m^3~s^-1)", p_z="P[flow]~('%')", f_h="f[high]~(year^-1)", f_l="f[low]~(year^-1)", r_r="atop(RC[rise],(m^3~s^-1~day^-1))", r_f="atop(RC[fall],(m^3~s^-1~day^-1))", fdc_s="SFDC~('-')",
                                SUB_thresh="SUB_thresh", EHA_thresh="EHA_thresh", LU_no="LU_no", TC_no="TC_no", LU_atts="LU_atts"), 
                        label_parsed)

# PAWN plots
plot_PAWN_KS_pars <- "sensitivity_PAWN_KS_pars.png"
plot_PAWN_bars <- "sensitivity_PAWN_bars.png"
plot_PAWN_ecdf <- "sensitivity_PAWN_ecdf.pdf"

# Variance-based SA plots (NOT IN PAPER)
plot_variance_bar <- "sensitivity_variance_barplots.png"
plot_variance_pattern <- "sensitivity_variance_tiles.png"



### UTILITY FUNCTION ###
# function to extract a legend from a ggplot for separate plotting
# obtained from: http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
} 


### CALCULATIONS ###

# get data
load(dat_file)

# explicitly list input factors (i.e. parameters) and realisations
dat_ind <- ind_df %>% 
  filter(index %in% ind_sel) %>%
  mutate(index = factor(index, levels=ind_sel)) %>%
  separate(realisation, c("p1", "p2", "s", "eha", "lunum", "luatts", "tc"), sep="_", extra="merge") %>%
  select(-p1,-p2)
dat_ind[,c("s", "eha", "lunum", "luatts", "tc")] <- lapply(dat_ind[,c("s", "eha", "lunum", "luatts", "tc")], function(x) as.numeric(gsub("[a-z]", "", x)))

params <- names(select(dat_ind, -value))


### PAWN density-based method (Pianosi & Wagener, 2015) ###
# calculation of conditional and unconditional empirical cumulated density values
dat_sens <- dat_ind %>%
  # apply over all indices
  ddply("index", .parallel=T, function(x){
    # apply over parameters to calculate condition cdfs
    lapply(params[-which(params == "index")], function(y){
      x %>%
        mutate(ecd_uncon = ecdf(value)(value), # calculate unconditional empirical cumulative densitiy values
               Nu = length(value)) %>% # number of values used to estimate ecdf
          # apply over parameter realisations of current parameter y
          ddply(y, function(z){
            z %>%
              select(value, ecd_uncon, Nu) %>%
              mutate(ecd = ecdf(value)(value)) %>% # calculate conditional empirical cumulative densitiy values for each parameter realisation
              mutate(ecd_diff = abs(ecd - ecd_uncon)) # calculate absolute difference between conditional and unconditional cdf
          })
    }) %>%
      rbind.fill()
  }) %>%
  gather("param", "param_val", one_of(params[-which(params == "index")]), na.rm=T) %>%
  mutate(param = as.factor(revalue(param, replace = par_labs))) %>% # replace parameter labels
  mutate(param = factor(param, levels=par_sel)) # ensure specified order of labels

# calculate Kolmogorov-Smirnoff statistic (i.e. maximum difference between conditional and unconditional ecdf)
dat_sens_ks <- dat_sens %>%
  ddply(c("index", "param", "param_val"), .parallel=T, function(x){
    x %>%
      summarise(KS = max(ecd_diff),
                Nu = mean(Nu),
                Nc = length(ecd_diff))
  })

# calculate mean and median of KS values and critical KS value by two-sample test with confidence level alpha = 0.05 (i.e. c(alpha) = 1.36)
dat_sens_stats <- dat_sens_ks %>%
  ddply(c("index", "param"), .parallel=T, function(x) {
    x %>%
      summarise(KS_mean = mean(KS),
                KS_med = median(KS),
                Nu = mean(Nu),
                Nc_mean = mean(Nc)
      ) %>%
      mutate(KS_crit = 1.36 * sqrt( (Nu+Nc_mean) / (Nu*Nc_mean) ) ) # use Nc on average used to compute conditional CDFs)
  })



#PLOTS #

# plot KS values vs. parameters with KS statistics
png(plot_PAWN_KS_pars, width=18, height=8, unit="in", res=300)
ggplot(dat_sens_ks, aes(x=param_val, y=KS)) +
  geom_line() +
  geom_point() +
  geom_hline(data = dat_sens_stats, aes(yintercept=KS_mean)) +
  geom_hline(data = dat_sens_stats, aes(yintercept=KS_med), linetype="dashed") +
  geom_hline(data = dat_sens_stats, aes(yintercept=KS_crit), colour="red") +
  coord_cartesian(ylim=c(0,1)) +
  facet_grid(index ~ param, scale="free_x")
dev.off()


# barplots of parameter vs. sensitivity measure for each streamflow index
gp <- ggplot(dat_sens_stats, aes(x=param, y=KS_med, group=index)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=KS_crit,ymax=KS_crit), colour=hcl(h = 15, l = 65, c = 100)) +
  labs(y="PAWN sensitivity index (-)", x="Discretisation parameter") +
  #scale_x_discrete(labels = parse(text = levels(dat_sens_stats$param))) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5), strip.text.x = element_text(size = 20)) +
  facet_grid(. ~ factor(index, levels=ind_sel), scales = "free", labeller = ind_labs)
ggsave(plot_PAWN_bars, plot = gp, width=18, height=8, units="in", dpi=300)


# plot ecdfs
dat_plot_ecd <- dat_sens %>%
  gather("ecd_type", "ecd_val", ecd_uncon, ecd)

# create plot for each parameter (to get individual color scales)
gp <- dat_plot_ecd %>%
  dlply("param", function(x){
    gp_t <- ggplot(x, aes(x=value, y=ecd_val)) +
      geom_line(data = subset(x, ecd_type!="ecd_uncon"), aes(group=param_val, colour=factor(param_val))) +
      geom_line(data = subset(x, ecd_type=="ecd_uncon"), colour="black") +
      scale_y_continuous(breaks = round(seq(0, 1, by = 0.5),1)) +
      labs(x="Streamflow index value", y="Empirical cumulated density (-)") +
      theme_bw(base_size = 18) +
      theme(legend.position = "none",  strip.text.x = element_text(size = 16)) +
      facet_grid(param ~ index, scale="free_x", labeller = ind_labs_units)
    ggplotGrob(gp_t)
  })

# combine all plots
gp.all <- do.call("rbind.gtable", gp)

# alterations (remove axes, reduces spaces, add legends, etc. for proper arrangement of plots)
# - push ylab of first sub plot across width of united plot
ylabs <- which(grepl("axis.title.y", sapply(gp.all$grobs, function(x) x$name)))
gp.all$layout[first(ylabs),"b"] <- gp.all$layout[last(ylabs),"b"]
# - remove ylab from all but first sub plot
gp.all$grobs[ylabs[-1]] <- lapply(gp.all$grobs[ylabs[-1]], function(x) x <- zeroGrob())

# remove bottom axes from all but the last plot (one bottom axis for every column, i.e. every index)
no_index <- length(levels(dat_plot_ecd$index))
bot_axes <- head(grep("axis-b", gp.all$layout$name), -no_index)
gp.all$grobs[bot_axes] <- lapply(gp.all$grobs[bot_axes], function(x) x <- zeroGrob())

# remove xlab from all but the last plot
xlabs <- which(grepl("axis.title.x", sapply(gp.all$grobs, function(x) x$name)))
gp.all$grobs[head(xlabs, -1)] <- lapply(gp.all$grobs[head(xlabs, -1)], function(x) x <- zeroGrob())

# remove panels from all but the first plot
panels <- grep("panel", gp.all$layout$name)
top <- unique(gp.all$layout$t[panels])[-1]
gp.all <- gp.all[-(top-1),]

# adjust heights between facet rows
#gp.all$heights[c(7,8,13,14,19,20,25,26,5,11,17,23)] <- unit(0, "pt")
gp.all$heights[c(11,12,21,22,31,32,41,42,8,18,28,38)] <- unit(0, "pt")
gp.all$heights[c(9,19,29,39)] <- unit(5.5, "pt")

# legend (generalized, manual)
df <- data.frame(
  x = runif(100),
  y = runif(100),
  z = runif(100)
)
leg1 <- ggplot(df, aes(x, y, colour=z)) +
  geom_point(aes(colour=z)) +
  scale_colour_gradientn("Parameter value for conditional density functions   ", breaks=c(.01,.5,.99), labels=c("small", "medium", "large"), colours = rainbow(10), limits=c(0,1)) +
  guides(colour = guide_colourbar(title.theme = element_text(size=18, angle = 0), title.hjust = 0, title.vjust = .8,
                                  label.theme = element_text(size=15, angle=0),
                                  barwidth = 15, barheight = 1.5, ticks = F, direction="horizontal")) +
  theme(legend.position = "bottom", legend.justification=c(.1,0), legend.margin = margin(t=15, b=5))

leg2 <- ggplot(df, aes(x, y, colour=1)) +
  geom_line() +
  scale_colour_continuous("Unconditional density function") +
  guides(colour = guide_legend(label.theme = element_blank(), title.theme = element_text(size=18, angle = 0), 
                               override.aes = list(linetype=1, colour="black", size=1.2))) +
  theme(legend.position = "bottom", legend.justification=c(0,0), legend.margin = margin(t=15, b=20), 
        legend.key = element_rect(fill="white"), legend.key.width = unit(30, "pt"))
  
# draw plot
pdf(plot_PAWN_ecdf, width=18, height=9)
lay <- rbind(c(1,1),
             c(2,3))
grid.arrange(gp.all, g_legend(leg1), g_legend(leg2), layout_matrix=lay, heights=c(.9,.1), widths=c(.5,.5))
dev.off()






 
### VARIANCE-BASED: NOT IN PAPER###
# variance-based sensitivities
dat_sens <- dat_ind %>%
  # apply over all indices
  ddply("index", .parallel=T, function(x){
    params_t <- params[-which(params == "index")]
    
    # calculate total variance
    var_tot <- var(x$value)
    
    # apply over parameters
    lapply(params_t, function(y){
      # calculate first-order-index
      first <- x %>%
        # apply over parameter realisations of current parameter y
        ddply(y, function(z){
          z %>%
            summarise(mean_i = mean(value)) # calculate mean over current realisations of y
        }) %>%
          summarise(first = var(mean_i) / var_tot) %>% # first-order index: variance of means divided by total variance
            mutate(parameter = y)
      
      # calculate total-order-index
      param_nr <- which(params_t == y) # get current parameter
      params_remain <- params_t[-param_nr] # extract current parameter from dataset
      tmp <- x # dataset as tmp
      names(tmp)[which(names(tmp) == y)] <- "y" # name of current parameter is "y"
      total <- tmp %>%
        ddply(params_remain, function(z) { # apply over remaining parameter combinations
          z %>%
            summarise(var_noni = var(value, na.rm=T)) # calculate variance over current realisation of 'non-y'
        }) %>%
          summarise(total = mean(var_noni, na.rm=T) / var_tot) %>% # total-order index: mean of variances divided by total variance
            mutate(parameter = y)
      
      return(full_join(total, first, by="parameter"))
      }) %>%
      rbind.fill()
  }) %>%
  gather("sens_ind", "sensit", total, first)

# barplots of parameter vs. sensitivity measure for each streamflow index
png(plot_variance_bar, width=16, height=6, unit="in", res=300)
ggplot(dat_sens, aes(x=reorder(parameter, sensit), y=sensit, group=sens_ind)) +
  geom_bar(aes(fill=sens_ind), stat="identity", position="dodge") +
  facet_grid(. ~ index) + 
  labs(y="Sensitivity (-)", x="Parameter", fill="Sensitivity\nindex") +
  theme(axis.text.x = element_text(angle = 45))
dev.off()

# tile plot with for parameter vs. streamflow index with tiles coloured according to sensitivity measure
png(plot_variance_pattern, width=8, height=6, unit="in", res=300)
ggplot(dat_sens, aes(x=reorder(parameter, sensit), y=index)) +
  geom_tile(aes(fill=sensit, width=.8, height=.8), colour="gray") +
  coord_fixed(ratio=1) +
  scale_fill_gradient2(low = "white",high = "darkblue", mid="lightblue") +
  labs(y="Streamflow index", x="Parameter", fill="Sensitivity") +
  facet_grid(. ~ sens_ind) + 
  theme_bw()
dev.off()


