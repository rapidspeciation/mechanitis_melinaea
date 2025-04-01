require("HDInterval")
require(tidyverse)
options(scipen = 999)

# functions
mcmcSummary = function(x) {
  mcmcSum = data.frame(stat = colnames(x))
  mcmcSum$median = -1
  mcmcSum$hpd_025 = -1
  mcmcSum$hpd_975 = -1
  
  for (cl in 2:ncol(x)) {
    cl
    hdi = hdi(x[,cl])
    mcmcSum$median[cl] = median(x[,cl])
    mcmcSum$hpd_025[cl] = hdi[1]
    mcmcSum$hpd_975[cl] = hdi[2]
  }
  
  mSum = mcmcSum[2:nrow(mcmcSum),]
  mSum_t = setNames(data.frame(t(mSum[,-1])), mSum[,1])
  return(mSum_t)
}


# parameters
dir = "path/to/mcmc/dir/"
popset = "mmp_msci"
model = "modelA"

# read mcmc
mcmc1 = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep1/",popset,".mcmc.txt"), header=T)
mcmc2 = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep2/",popset,".mcmc.txt"), header=T)
mcmc3 = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep3/",popset,".mcmc.txt"), header=T)
mcmc4 = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep4/",popset,".mcmc.txt"), header=T)
mcmc5 = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep5/",popset,".mcmc.txt"), header=T)
mcmc6 = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep6/",popset,".mcmc.txt"), header=T)

# summary of mcmc runs
m1 = mcmcSummary(mcmc1); m1
m2 = mcmcSummary(mcmc2); m2
m3 = mcmcSummary(mcmc3); m3
m4 = mcmcSummary(mcmc4); m4
m5 = mcmcSummary(mcmc5); m5
m6 = mcmcSummary(mcmc6); m6



## Traces ****************************************************************************************
# read mcmc
mcmcAll = data.frame()
for (i in 1:6) {
  inp = read.table(paste0(dir,popset,"/bpp-msci-noncod.l2000_g10000/",model,"/rep",i,"/",popset,".mcmc.txt"), header=T)
  inp$rep = as.character(i)
  mcmcAll = rbind(mcmcAll, inp) 
}
# plot trace
mcmcAll %>%
  filter(row_number() %% 100 == 0) %>%
  pivot_longer(-c(Gen,rep), names_to = "parameter", values_to = "value" ) %>%
  # filter(str_detect(parameter, "^phi")) %>%
  filter(!str_detect(parameter, "lnL")) %>%
  ggplot() +
  geom_line(aes(x=Gen/1000000, y=value, col=rep)) +
  scale_x_continuous(breaks = seq(0,1,.1), name = expression("MCMC Iteration (x"~10^6~")")) +
  scale_color_brewer() +
  facet_wrap(~parameter, ncol=4, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none") +
  ggtitle(paste0(popset," - ",model))

# plot trace (density)
tet = mcmcAll %>% 
  pivot_longer(-c(Gen,rep), names_to = "parameter", values_to = "value" ) %>%
  filter(str_detect(parameter, "^theta")) %>%
  ggplot() +
  geom_density(aes(value, group = rep, col=rep)) +
  facet_wrap(~parameter, ncol=4) +
  theme_classic() +
  theme(legend.position = "none")
tau = mcmcAll %>% 
  pivot_longer(-c(Gen,rep), names_to = "parameter", values_to = "value" ) %>%
  filter(str_detect(parameter, "^tau")) %>%
  ggplot() +
  geom_density(aes(value, group = rep, col=rep)) +
  facet_wrap(~parameter, ncol=4) +
  theme_classic() +
  theme(legend.position = "none")
phi = mcmcAll %>% 
  pivot_longer(-c(Gen,rep), names_to = "parameter", values_to = "value" ) %>%
  filter(str_detect(parameter, "^phi")) %>%
  ggplot() +
  geom_density(aes(value, group = rep, col=rep)) +
  facet_wrap(~parameter, ncol=4) +
  theme_classic() +
  theme(legend.position = "none")

# combine plots
require(cowplot)
t = plot_grid(tet)
b = plot_grid(tau, phi, nrow=1, rel_widths = c(4,1))
plot_grid(t,b, nrow=2, rel_heights = c(2,1))
