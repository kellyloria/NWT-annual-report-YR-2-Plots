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

## annual average copepod densities GL4
cop_den_ave <- aggregate(taxon_dens ~ year, data=cop_den, FUN=mean) 
cop_den_ave$SE <- aggregate(taxon_dens ~ year, data=cop_den, FUN=se)[,2]

## annual average all zoop densities GL4
zoop_den_ave <- aggregate(tot_density ~ year, data=dt1, FUN=mean) 
zoop_den_ave$SE <- aggregate(tot_density ~ year, data=dt1, FUN=se)[,2]


## annual max cladoceran densities GL4
clad_den_max <- aggregate(taxon_dens ~ year, data=clad_den, FUN=max) 
clad_den_max$SE <- aggregate(taxon_dens ~ year, data=clad_den, FUN=se)[,2]

## annual max copepod densities GL4
cop_den_max <- aggregate(taxon_dens ~ year, data=cop_den, FUN=max) 
cop_den_max$SE <- aggregate(taxon_dens ~ year, data=cop_den, FUN=se)[,2]

## annual total  copepod densities GL4
cop_den_max <- aggregate(taxon_dens ~ year, data=cop_den, FUN=max) 
cop_den_max$SE <- aggregate(taxon_dens ~ year, data=cop_den, FUN=se)[,2]

## annual max all zoop densities GL4
zoop_den_max <- aggregate(tot_density ~ year, data=dt1, FUN=max) 
zoop_den_max$SE <- aggregate(tot_density ~ year, data=dt1, FUN=se)[,2]

