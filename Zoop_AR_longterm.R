### Long-Term Zooplankton ###

# Created by Kelly Loria #
# 11/09/2018
##########################

# Load in packages
library(ggplot2)
library(scales)
library(dplyr)
library(tidyverse)
library(lme4)
library(lmerTest)

# load in extended summer PCA values 
# choose NWT_sumallPCclimate_19822017.csv
dt_summer <- read.csv(file.choose())
head(dt_summer)

# load zooplankton data 
infile1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-nwt/161/1/9b37e40db2134f4f163cdac258d54560" 
infile1 <- sub("^https","http",infile1) 
dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "LTER_site",     
                 "local_site",     
                 "year",     
                 "date",     
                 "depth",     
                 "net_um",     
                 "water_sampled",     
                 "taxon_name",     
                 "num_of_taxon",     
                 "taxon_dens",     
                 "avg_taxon_length",     
                 "tot_indiv_meas",     
                 "tot_counted",     
                 "tot_density"    ), check.names=TRUE)


# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$LTER_site)!="factor") dt1$LTER_site<- as.factor(dt1$LTER_site)
if (class(dt1$local_site)!="factor") dt1$local_site<- as.factor(dt1$local_site)                                   
# attempting to convert dt1$date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1date<-as.Date(dt1$date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1date) == length(tmp1date[!is.na(tmp1date)])){dt1$date <- tmp1date } else {print("Date conversion failed for dt1$date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1date) 
if (class(dt1$depth)=="factor") dt1$depth <-as.numeric(levels(dt1$depth))[as.integer(dt1$depth) ]
if (class(dt1$net_um)=="factor") dt1$net_um <-as.numeric(levels(dt1$net_um))[as.integer(dt1$net_um) ]
if (class(dt1$water_sampled)=="factor") dt1$water_sampled <-as.numeric(levels(dt1$water_sampled))[as.integer(dt1$water_sampled) ]
if (class(dt1$taxon_name)!="factor") dt1$taxon_name<- as.factor(dt1$taxon_name)
if (class(dt1$num_of_taxon)=="factor") dt1$num_of_taxon <-as.numeric(levels(dt1$num_of_taxon))[as.integer(dt1$num_of_taxon) ]
if (class(dt1$taxon_dens)=="factor") dt1$taxon_dens <-as.numeric(levels(dt1$taxon_dens))[as.integer(dt1$taxon_dens) ]
if (class(dt1$avg_taxon_length)=="factor") dt1$avg_taxon_length <-as.numeric(levels(dt1$avg_taxon_length))[as.integer(dt1$avg_taxon_length) ]
if (class(dt1$tot_indiv_meas)=="factor") dt1$tot_indiv_meas <-as.numeric(levels(dt1$tot_indiv_meas))[as.integer(dt1$tot_indiv_meas) ]
if (class(dt1$tot_counted)=="factor") dt1$tot_counted <-as.numeric(levels(dt1$tot_counted))[as.integer(dt1$tot_counted) ]
if (class(dt1$tot_density)=="factor") dt1$tot_density <-as.numeric(levels(dt1$tot_density))[as.integer(dt1$tot_density) ]

# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1) 

#Subset for all cladoceran or copepod densities 
summary(dt1$taxon_name)
names(dt1)

clad_den <- subset(dt1, taxon_name=="Total cladoceran",
                     select=LTER_site:tot_density)
summary(clad_den$taxon_name) 

cop_den <- subset(dt1, taxon_name=="Total copepod",
                   select=LTER_site:tot_density)
summary(cop_den$taxon_name) 


# function to calculate standard error around the mean 
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

## annual average cladoceran densities GL4
clad_den_ave <- aggregate(taxon_dens ~ year, data=clad_den, FUN=mean) 
clad_den_ave$SE <- aggregate(taxon_dens ~ year, data=clad_den, FUN=se)[,2]
clad_den_ave$label <- "clad"

## annual average copepod densities GL4
cop_den_ave <- aggregate(taxon_dens ~ year, data=cop_den, FUN=mean) 
cop_den_ave$SE <- aggregate(taxon_dens ~ year, data=cop_den, FUN=se)[,2]
cop_den_ave$label <- "cop"

## annual average all zoop densities GL4
zoop_den_ave <- aggregate(tot_density ~ year, data=dt1, FUN=mean) 
zoop_den_ave$SE <- aggregate(tot_density ~ year, data=dt1, FUN=se)[,2]
#change column name to be able to compare
names(zoop_den_ave)[2]<-"taxon_dens"
zoop_den_ave$label <- "all"

###### combind averages for groups of zoop
zoop_ave_den_group <- rbind(zoop_den_ave, cop_den_ave, clad_den_ave)
###### combind averages for 2 groups of zoop
zoop_ave_den_Dgroup <- rbind(cop_den_ave, clad_den_ave)


## annual max cladoceran densities GL4
clad_den_max <- aggregate(taxon_dens ~ year, data=clad_den, FUN=max) 
clad_den_max$SE <- aggregate(taxon_dens ~ year, data=clad_den, FUN=se)[,2]
clad_den_max$label <- "clad"

## annual total  copepod densities GL4
cop_den_max <- aggregate(taxon_dens ~ year, data=cop_den, FUN=max) 
cop_den_max$SE <- aggregate(taxon_dens ~ year, data=cop_den, FUN=se)[,2]
cop_den_max$label <- "cop"

## annual max all zoop densities GL4
zoop_den_max <- aggregate(tot_density ~ year, data=dt1, FUN=max) 
zoop_den_max$SE <- aggregate(tot_density ~ year, data=dt1, FUN=se)[,2]
names(zoop_den_max)[2]<-"taxon_dens"
zoop_den_max$label <- "all"

###### combind averages for groups of zoop
zoop_max_den_group <- rbind(zoop_den_max, cop_den_max, clad_den_max)
###### combind max for 2 groups of zoop
zoop_max_den_Dgroup <- rbind(cop_den_max, clad_den_max)


## Average annual zoop densities for copepods and cladocerans 


png("Zoop_ave_den_overtime.png",  
    width = 7.5,
    height = 4.2,
    units = "in",
    res = 1200,
    pointsize = 8)

plot(NA, NA, xlim = c(2012, 2017), ylim = c(0, 8.5), xlab = "Year",
     ylab = "Taxon average density in GL4 (perL)", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2012:2017), cex.axis=1.2)

lines(y = clad_den_ave$taxon_dens, x= clad_den_ave$year, lty =1, lwd = 2, col= "sienna1")
lines(y = cop_den_ave$taxon_dens,  x= cop_den_ave$year, lty =1, lwd = 2, col= "tomato")
lines(y = zoop_den_ave$taxon_dens,  x= zoop_den_ave$year, lty =1, lwd = 2, col= "tomato3")


points(y = clad_den_ave$taxon_dens, x= clad_den_ave$year, col= "sienna1", pch=19) 
arrows(clad_den_ave$year, (clad_den_ave$taxon_dens + clad_den_ave$SE), 
       clad_den_ave$year, (clad_den_ave$taxon_dens - clad_den_ave$SE), 
       length=0.0, code=3, lwd=1, col= "sienna1")


points(y = cop_den_ave$taxon_dens, x= cop_den_ave$year, col= "tomato", pch=19) 
arrows(cop_den_ave$year, (cop_den_ave$taxon_dens + cop_den_ave$SE), 
       cop_den_ave$year, (cop_den_ave$taxon_dens - cop_den_ave$SE), 
       length=0.0, code=3, lwd=1, col= "tomato")

points(y = zoop_den_ave$taxon_dens, x= zoop_den_ave$year, col= "tomato3", pch=19) 
arrows(zoop_den_ave$year, (zoop_den_ave$taxon_dens + zoop_den_ave$SE), 
       zoop_den_ave$year, (zoop_den_ave$taxon_dens - zoop_den_ave$SE), 
       length=0.0, code=3, lwd=1, col= "tomato3")

legend("topright", legend=c("Cladoceran", "Copepod", "All taxa"),
       col=c("sienna1", "tomato", "tomato3"), lty=1, lwd=2, cex=1.2)

dev.off()

png("Zoop_max_den_overtime.png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

plot(NA, NA, xlim = c(2012, 2017), ylim = c(0, 12), xlab = "Year",
     ylab = "Taxon max density in GL4 (perL)", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2012:2017), cex.axis=1.5)

lines(y = clad_den_max$taxon_dens, x= clad_den_max$year, lty =1, lwd = 1, col= "sienna1")
lines(y = cop_den_max$taxon_dens,  x= cop_den_max$year, lty =1, lwd = 1, col= "tomato")
lines(y = zoop_den_max$taxon_dens,  x= zoop_den_max$year, lty =1, lwd = 1, col= "tomato3")


points(y = clad_den_max$taxon_dens, x= clad_den_max$year, col= "sienna1", pch=19) 
arrows(clad_den_max$year, (clad_den_max$taxon_dens + clad_den_max$SE), 
       clad_den_max$year, (clad_den_max$taxon_dens - clad_den_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "sienna1")


points(y = cop_den_max$taxon_dens, x= cop_den_max$year, col= "tomato", pch=19) 
arrows(cop_den_max$year, (cop_den_max$taxon_dens + cop_den_max$SE), 
       cop_den_max$year, (cop_den_max$taxon_dens - cop_den_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "tomato")

points(y = zoop_den_max$taxon_dens, x= zoop_den_max$year, col= "tomato3", pch=19) 
arrows(zoop_den_max$year, (zoop_den_max$taxon_dens + zoop_den_max$SE), 
       zoop_den_max$year, (zoop_den_max$taxon_dens - zoop_den_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "tomato3")

legend("topright", legend=c("Cladoceran", "Copepod", "All taxa"),
       col=c("sienna1", "tomato", "tomato3"), lty=1, cex=1.2)

dev.off()

# zoop den all averages and pca for extened summer
##zooplankton and data and PCA
comb_dat_zoop_ave <- left_join(zoop_ave_den_group, dt_summer[c("eco_year", "sumallPC1")],
                         by = c("year" = "eco_year"))

comb_dat_Dzoop_ave <- left_join(zoop_ave_den_Dgroup, dt_summer[c("eco_year", "sumallPC1")],
                               by = c("year" = "eco_year"))

comb_dat_allzoop_ave <- left_join(zoop_den_ave, dt_summer[c("eco_year", "sumallPC1")],
                                by = c("year" = "eco_year"))

# zoop den MAX averages and pca for extened summer
comb_dat_zoop_max <- left_join(zoop_max_den_group, dt_summer[c("eco_year", "sumallPC1")],
                               by = c("year" = "eco_year"))

comb_dat_Dzoop_max <- left_join(zoop_max_den_Dgroup, dt_summer[c("eco_year", "sumallPC1")],
                                by = c("year" = "eco_year"))


all_den <- subset(dt1, taxon_name=="Chydoridae sp." | taxon_name=="Daphnia pulicaria" 
                  | taxon_name=="Hesperodiaptomus shoshone" | taxon_name=="Nauplii "
                  | taxon_name=="Neonate"| taxon_name=="Bosminidae sp.",
                  select=LTER_site:tot_density)
summary(all_den$taxon_name) 



comb_dat_all_raw <- left_join(all_den, dt_summer[c("eco_year", "sumallPC1")],
                                by = c("year" = "eco_year"))

##### PCA Plots######

ave_zoop_plot <- qplot(sumallPC1, taxon_dens, data = comb_dat_zoop_ave, geom="point", ylab = "Average annual taxa density", color=label) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

ave_zoop_plot + scale_colour_manual(values = c("tomato3", "sienna1", "tomato"))

ggsave("PCA_ave_zoop_den.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

zoop.mod_ave_all <- glm(taxon_dens ~ scale(sumallPC1) + label, data = comb_dat_zoop_ave)
summary(zoop.mod_ave_all)

zoop.mod_ave_all.1 <- lmer(taxon_dens ~ scale(sumallPC1) + (1|label), data = comb_dat_zoop_ave)
summary(zoop.mod_ave_all.1)

### just clad and cop
ave_zoop_plot <- qplot(sumallPC1, taxon_dens, data = comb_dat_Dzoop_ave, geom="point", ylab = "Average annual taxa density", color=label) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

ave_zoop_plot + scale_colour_manual(values = c("sienna1", "tomato"))

ggsave("PCA_ave_zoop_den.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

Dzoop.mod_ave_all <- glm(taxon_dens ~ scale(sumallPC1) + label, data = comb_dat_Dzoop_ave)
summary(Dzoop.mod_ave_all)

Dzoop.mod_ave_all.1 <- lmer(taxon_dens ~ scale(sumallPC1) + (1|label), data = comb_dat_Dzoop_ave)
summary(Dzoop.mod_ave_all.1)

Azoop.mod_ave_all <- glm(taxon_dens ~ scale(sumallPC1), data = comb_dat_allzoop_ave)
summary(Azoop.mod_ave_all)

### PCA Plots for max 
max_zoop_plot <- qplot(sumallPC1, taxon_dens, data = comb_dat_zoop_max, geom="point", ylab = "Average annual taxa density", color=label) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

max_zoop_plot + scale_colour_manual(values = c("tomato3", "sienna1", "tomato"))

ggsave("PCA_max_zoop_den.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


Dzoop.mod_max_all <- glm(taxon_dens ~ scale(sumallPC1) + label, data = comb_dat_Dzoop_max)
summary(Dzoop.mod_max_all)

Dzoop.mod_max_all.1 <- lmer(taxon_dens ~ scale(sumallPC1) + (1|label), data = comb_dat_Dzoop_max)
summary(Dzoop.mod_max_all.1)

names(dt1)
head(dt1)
dt1$taxon_dens

zoop.mod_all <- glm(taxon_dens ~ scale(sumallPC1), data = comb_dat_all_raw)
summary(zoop.mod_all)