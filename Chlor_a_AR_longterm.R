### Long-Term Chlorophyll-a plots for all locations within GL4 ###

# Created by Kelly Loria #
# 11/08/2018
##########################

# Load in packages
library(ggplot2)
library(scales)
library(dplyr)
library(tidyverse)


# load in extended summer PCA values 
# choose NWT_sumallPCclimate_19822017.csv
dt_summer <- read.csv(file.choose())
head(dt_summer)
# x axis with PC1: earlier <---> later in summer

# load on water quality data with chla values:
# water_quality_GLV.dm.data_KL_11_07_2018.csv
WQdata <- read.csv(file.choose())
names(WQdata)

#Subset the data for GL4 observations with chla observations (>0)
WQdat1 <- subset(WQdata, local_site=="GL4" & chl_a > -0.5,
                 select=local_site:chl_a)
names(WQdat1)

# Modify format of dates so it is not seen as a factor 
WQdat1$date
range(WQdat1$date)
WQdat1$date1 <- as.Date(WQdat1$date, format="%Y-%m-%d")
range(WQdat1$date1)


### subset for locations in GL4###

#Subset for all within lake observations (0-9m)
WQdat_lake <- subset(WQdat1, location=="LAKE" & chl_a > -0.5,
                     select=local_site:date1)

# subset for with in lake observations at Surface
WQdat_lake_0m <- subset(WQdat_lake, location=="LAKE" & depth  < 1,
                        select=local_site:date1)

# subset for with in lake observations at metalimnion
WQdat_lake_3m <- subset(WQdat_lake, location=="LAKE" & depth  > 1 & depth < 7,
                        select=local_site:date1)

# subset for with in lake observations at hypolimnion
WQdat_lake_9m <- subset(WQdat_lake, location=="LAKE" & depth  > 7,
                        select=local_site:date1)

# subset for with in lake observations at lake inlet
WQdat_in_lake<- subset(WQdat1, location== "INLET",
                       select=local_site:date1)

# subset for with in lake observations at lake outlet
WQdat_out_lake<- subset(WQdat1, location== "OUTLET",
                        select=local_site:date1)

### r binds to merge subsetted data inlet and outlet
in_out_wc <- rbind(WQdat_in_lake, WQdat_out_lake)


### Calculate average values for lake locations ###

# function to calculate standard error around the mean 
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

## average chlor-a for all inlake depths by year 
WQdat_lake_ave <- aggregate(chl_a ~ year, data=WQdat_lake, FUN=mean) 
WQdat_lake_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake, FUN=se)[,2]

# average annual chla at inlet
WQdat_IN_ave <- aggregate(chl_a ~ year, data= WQdat_in_lake, FUN=mean) 
WQdat_IN_ave$SE <- aggregate(chl_a ~ year, data=WQdat_in_lake, FUN=se)[,2]

# average annual chla at outlet
WQdat_out_ave <- aggregate(chl_a ~ year, data= WQdat_out_lake, FUN=mean) 
WQdat_out_ave$SE <- aggregate(chl_a ~ year, data=WQdat_out_lake, FUN=se)[,2]

# average annual chla at lake surface
WQdat_0m_ave <- aggregate(chl_a ~ year, data= WQdat_lake_0m, FUN=mean) 
WQdat_0m_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake_0m, FUN=se)[,2]

# average annual chla at metalimnion
WQdat_3m_ave <- aggregate(chl_a ~ year, data= WQdat_lake_3m, FUN=mean) 
WQdat_3m_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake_3m, FUN=se)[,2]

# average annual chla at hypolimnion 
WQdat_9m_ave <- aggregate(chl_a ~ year, data= WQdat_lake_9m, FUN=mean) 
WQdat_9m_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake_9m, FUN=se)[,2]


### Combine extended summer PCA with subsetted data ###

##chlora for all GL4 (Surface, Meta, Hypo, IN and OUT) data and PCA
comb_dat_WQ <- left_join(WQdat1, dt_summer[c("eco_year", "sumallPC1")],
                         by = c("year" = "eco_year"))

## chlora for all in lake data and PCA
comb_dat_lake_WQ <- left_join(WQdat_lake, dt_summer[c("eco_year", "sumallPC1")],
                              by = c("year" = "eco_year"))

## Annual Average chlora for all inlake data and PCA
comb_dat_L_WQ <- left_join(WQdat_lake_ave, dt_summer[c("eco_year", "sumallPC1")],
                           by = c("year" = "eco_year"))

## chlora for all inlet lake data and PCA
comb_dat_in_WQ <- left_join(WQdat_in_lake, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## chlora for inlet and outlet data and PCA
comb_dat_in_out_WQ <- left_join(in_out_wc, dt_summer[c("eco_year", "sumallPC1")],
                                by = c("year" = "eco_year"))

## Annual Average chlora for hypo and PCA
comb_dat_9m_WQ <- left_join(WQdat_9m_ave, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

#########################################################################
### Visualize relationships of Chlor overtime by observation location ###


#png("Chlor_A overtime(all locations).png",  
#    width = 5.25,
#    height = 3.25,
#    units = "in",
#    res = 1200,
#    pointsize = 4 )

## Average annual chla for all lake locations

plot(NA, NA, xlim = c(2000, 2017), ylim = c(0, 10), xlab = "Year",
     ylab = "Average annual chlor-A in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2000:2017), cex.axis=1.5)

lines(y = WQdat_IN_ave$chl_a, x= WQdat_IN_ave$year, lty =1, lwd = 1, col= "lightseagreen")
lines(y = WQdat_out_ave$chl_a,  x= WQdat_out_ave$year, lty =1, lwd = 1, col= "cadetblue3")
lines(y = WQdat_0m_ave$chl_a,  x= WQdat_0m_ave$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WQdat_3m_ave$chl_a,  x= WQdat_3m_ave$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WQdat_9m_ave$chl_a,  x= WQdat_9m_ave$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WQdat_IN_ave$chl_a, x= WQdat_IN_ave$year, col= "lightseagreen", pch=19) 
arrows(WQdat_IN_ave$year, (WQdat_IN_ave$chl_a + WQdat_IN_ave$SE), 
       WQdat_IN_ave$year, (WQdat_IN_ave$chl_a - WQdat_IN_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WQdat_out_ave$chl_a, x= WQdat_out_ave$year, col= "cadetblue3", pch=19) 
arrows(WQdat_out_ave$year, (WQdat_out_ave$chl_a + WQdat_out_ave$SE), 
       WQdat_out_ave$year, (WQdat_out_ave$chl_a - WQdat_out_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")

points(y = WQdat_0m_ave$chl_a, x= WQdat_0m_ave$year, col= "darkolivegreen3", pch=19) 
arrows(WQdat_0m_ave$year, (WQdat_0m_ave$chl_a + WQdat_0m_ave$SE), 
       WQdat_0m_ave$year, (WQdat_0m_ave$chl_a - WQdat_0m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WQdat_3m_ave$chl_a, x= WQdat_3m_ave$year, col= "darkolivegreen4", pch=19) 
arrows(WQdat_3m_ave$year, (WQdat_3m_ave$chl_a + WQdat_3m_ave$SE), 
       WQdat_3m_ave$year, (WQdat_3m_ave$chl_a - WQdat_3m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WQdat_9m_ave$chl_a, x= WQdat_9m_ave$year, col= "darkgreen", pch=19) 
arrows(WQdat_9m_ave$year, (WQdat_9m_ave$chl_a + WQdat_9m_ave$SE), 
       WQdat_9m_ave$year, (WQdat_9m_ave$chl_a - WQdat_9m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topleft", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

#dev.off()



#png("Chlor_A overtime(surface,meta,hypo).png",  
#    width = 5.25,
#    height = 3.25,
#    units = "in",
#    res = 1200,
#    pointsize = 4 )

## Average annual chla for within lake locations

plot(NA, NA, xlim = c(2000, 2017), ylim = c(0, 10), xlab = "Year",
     ylab = "Average annual chlor-A in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2000:2017), cex.axis=1.5)

lines(y = WQdat_0m_ave$chl_a,  x= WQdat_0m_ave$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WQdat_3m_ave$chl_a,  x= WQdat_3m_ave$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WQdat_9m_ave$chl_a,  x= WQdat_9m_ave$year, lty =1, lwd = 1, col= "darkgreen")

points(y = WQdat_0m_ave$chl_a, x= WQdat_0m_ave$year, col= "darkolivegreen3", pch=19) 
arrows(WQdat_0m_ave$year, (WQdat_0m_ave$chl_a + WQdat_0m_ave$SE), 
       WQdat_0m_ave$year, (WQdat_0m_ave$chl_a - WQdat_0m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WQdat_3m_ave$chl_a, x= WQdat_3m_ave$year, col= "darkolivegreen4", pch=19) 
arrows(WQdat_3m_ave$year, (WQdat_3m_ave$chl_a + WQdat_3m_ave$SE), 
       WQdat_3m_ave$year, (WQdat_3m_ave$chl_a - WQdat_3m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")


points(y = WQdat_9m_ave$chl_a, x= WQdat_9m_ave$year, col= "darkgreen", pch=19) 
arrows(WQdat_9m_ave$year, (WQdat_9m_ave$chl_a + WQdat_9m_ave$SE), 
       WQdat_9m_ave$year, (WQdat_9m_ave$chl_a - WQdat_9m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topleft", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

#dev.off()


#png("Chlor_A overtime(hypo).png",  
#    width = 5.25,
#    height = 3.25,
#    units = "in",
#    res = 1200,
#    pointsize = 4 )

## Average annual chla for hypo locations

plot(NA, NA, xlim = c(2000, 2017), ylim = c(0, 10), xlab = "Year",
     ylab = "Average annual chlor-A in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2000:2017), cex.axis=1.5)


lines(y = WQdat_9m_ave$chl_a,  x= WQdat_9m_ave$year, lty =1, lwd = 1, col= "darkgreen")

points(y = WQdat_9m_ave$chl_a, x= WQdat_9m_ave$year, col= "darkgreen", pch=19) 
arrows(WQdat_9m_ave$year, (WQdat_9m_ave$chl_a + WQdat_9m_ave$SE), 
       WQdat_9m_ave$year, (WQdat_9m_ave$chl_a - WQdat_9m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topleft", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

#dev.off()




# plot chlor-a for in lake averages and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_L_WQ, geom="point", ylab = "Average annual in lake chlor-a") + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

chlora.mod_1yr <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_L_WQ)
summary(chlora.mod_1yr) # year average not significant

# plot chlor-a for in lake averages and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_WQ, geom="point", ylab = "Average chlor-a ", color = factor(location), shape=factor(depth)) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

chlora.mod1 <- lmer(chl_a ~ scale(sumallPC1) + (1 | location), data = comb_dat_WQ)
summary(chlora.mod1) 
#inlet outlet might be more responsive 

# plot chlor-a for in lake and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_lake_WQ, geom="point", ylab = "Average chlor-a ", color = factor(depth)) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

chlora.mod2 <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_lake_WQ)
summary(chlora.mod2) # not significant 


# plot chlor-a for inlet and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_in_WQ, geom="point", ylab = "Average chlor-a inlet", color = factor(depth)) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

chlora.mod2 <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_in_WQ)
summary(chlora.mod2) # inlet is significant 

# plot chlor-a for inlet and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = , geom="point", ylab = "Average chlor-a inlet", color = factor(depth)) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

chlora.mod2 <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_in_WQ)
summary(chlora.mod2) # inlet is significant 


#comb_dat_in_out_WQ

# plot chlor-a for 9m and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_9m_WQ, geom="point", ylab = "Average chlor-a ") + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

chlora.mod_9m <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_9m_WQ)
summary(chlora.mod_9m) # inlet and outlet is more significant 

#Inlet and outlet chlorophyll-a seem to more sensitive that inlake values

## check out chlor-a max 





