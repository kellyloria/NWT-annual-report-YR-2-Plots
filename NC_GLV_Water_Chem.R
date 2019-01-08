### Compling Nel Caine's Stream Water Chem Data ####
## by Kelly Loria ##
# 2018-12-20

library(tidyverse)

nav <- read.csv("navasolu.nc.data.csv", header=T, 
                na = c("NaN", "DNS",  "EQCL", "N/A", "NP", "NSS", "NV", "u", "QNS", NA, " ", ""))
ari <- read.csv("ariksolu.nc.data.csv", header=T, 
                na = c("NaN", "DNS",  "EQCL", "N/A", "NP", "NSS", "NV", "u", "QNS", NA, " ", ""))
gl5 <- read.csv("gre5solu.nc.data.csv", header=T,
                na = c("NaN", "DNS",  "EQCL", "N/A", "NP", "NSS", "NV", "u", "QNS", NA, " ", ""))
alb <- read.csv("albisolu.nc.data.csv", header=T, 
                na = c("NaN", "DNS",  "EQCL", "N/A", "NP", "NSS", "NV", "u", "QNS", NA, " ", ""))
gl4 <- read.csv("gre4solu.nc.data.csv", header=T, 
                na = c("NaN", "DNS",  "EQCL", "N/A", "NP", "NSS", "NV", "u", "QNS", NA, " ", ""))

newdat <- rbind(nav, ari, gl5, gl4, alb)
summary(newdat)

#select columns of interest
DOC_dat <- select(newdat, c("local_site", "year", "date", "DOC"))
colnames(DOC_dat)

summary(DOC_dat$DOC)

ggplot(data = DOC_dat, # data
       aes(x = year, # aesthetics
           y = DOC,
           color = local_site)) +
  geom_point() + geom_smooth(method = "lm") +
  # geom 
  facet_grid(~local_site) + # facet
  scale_color_brewer(palette = "Spectral")

### nav ###
nav.DOC.mod <- lm(DOC ~ year, data = nav)
summary(nav.DOC.mod)
range(nav$year)

### ari ###
ari.DOC.mod <- lm(DOC ~ year, data = ari)
summary(ari.DOC.mod)
range(ari$year)

### gl4 ###
gl4.DOC.mod <- lm(DOC ~ year, data = gl4)
summary(gl4.DOC.mod)
range(gl4$year)

### gl5 ###
gl5.DOC.mod <- lm(DOC ~ year, data = gl5)
summary(gl5.DOC.mod)
range(gl5$year)

### alb ###
alb.DOC.mod <- lm(DOC ~ year, data = alb)
summary(alb.DOC.mod)
range(alb$year)


