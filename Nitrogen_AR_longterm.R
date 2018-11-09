### Long-Term Nitrogen plots for all locations within GL4 ###

# Created by Kelly Loria #
# 11/08/2018
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
# x axis with PC1: earlier <---> later in summer

# load on water chem data with chla values:
# glvwatsolu.dm.data_KL_11_07_2018.csv
WCdata <- read.csv(file.choose(),
                   na = c("NaN", "DNS",  "EQCL", "N/A", "NP", "NSS", "NV", "u", "QNS", NA, " ", ""))
names(WCdata)





#####################
### Nitrate plot ###

# Modify format of dates so it is not seen as a factor 
WCdata$date
range(WCdata$date)
WCdata$date1 <- as.Date(WCdata$date, format="%Y-%m-%d")
range(WCdata$date1)

### subset for locations in GL4###

#Subset the data for GL4 observations with chla observations (>0)
WCdat1 <- subset(WCdata, local_site=="GL4",
                 select=local_site:date1)
summary(WCdat1$local_site) #all set


#Subset for all within lake observations (0-9m)
WCdat_lake <- subset(WCdat1, location=="LAKE",
                     select=local_site:date1)
summary(WCdat_lake$location)

#Subset for all 0m in lake samples: 
#Subset for all within lake observations (0m)
WCdat_0m <- subset(WCdat_lake, location=="LAKE" & depth  < 1,
                     select=local_site:date1)
summary(WCdat_0m)

#Subset for all within lake observations (3m)
WCdat_3m <- subset(WCdat_lake, location=="LAKE" & depth  > 1 & depth < 7,
                   select=local_site:date1)
summary(WCdat_3m)

#Subset for all within lake observations (9m)
WCdat_9m <- subset(WCdat_lake, location=="LAKE" & depth > 7,
                   select=local_site:date1)
summary(WCdat_9m)

#Subset for all within lake observations INLET
WCdat_IN <- subset(WCdat1, location=="INLET",
                   select=local_site:date1)
summary(WCdat_IN)

#Subset for all within lake observations OUTLET
WCdat_OUT <- subset(WCdat1, location=="OUTLET",
                   select=local_site:date1)
summary(WCdat_OUT)


### Calculate average values for lake locations ###

# function to calculate standard error around the mean 
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}


# average annual nitrate at inlet
WCdat_IN$TD
WCdat_IN_ave <- aggregate(TDN ~ year, data= WCdat_IN, FUN=mean) 
WCdat_IN_ave$SE <- aggregate(TDN ~ year, data=WCdat_IN, FUN=se)[,2]

WCdat_IN_ave$depth <- "IN"

# average annual chla at outlet
WCdat_out_ave <- aggregate(TDN ~ year, data= WCdat_OUT, FUN=mean) 
WCdat_out_ave$SE <- aggregate(TDN ~ year, data =WCdat_OUT, FUN=se)[,2]

WCdat_out_ave$depth <- "OUT"

# average annual chla at lake surface
WCdat_0m_ave <- aggregate(TDN ~ year, data= WCdat_0m, FUN=mean) 
WCdat_0m_ave$SE <- aggregate(TDN ~ year, data=WCdat_0m, FUN=se)[,2]

WCdat_0m_ave$depth <- "0m"

# average annual chla at metalimnion
WCdat_3m_ave <- aggregate(TDN ~ year, data= WCdat_3m, FUN=mean) 
WCdat_3m_ave$SE <- aggregate(TDN ~ year, data= WCdat_3m, FUN=se)[,2]

WCdat_3m_ave$depth <- "3m"

# average annual chla at hypolimnion 
WCdat_9m_ave <- aggregate(TDN ~ year, data= WCdat_9m, FUN=mean) 
WCdat_9m_ave$SE <- aggregate(TDN ~ year, data=WCdat_9m, FUN=se)[,2]

WCdat_9m_ave$depth <- "9m"

ave_annual_TDN <- rbind(WCdat_IN_ave, WCdat_out_ave, 
                        WCdat_0m_ave, WCdat_3m_ave, WCdat_9m_ave)

## join nwt climate dataset (year and sumallPC1) with water quality dataset
all_TDN_combined_dat <- left_join(ave_annual_TDN, dt_summer[c("eco_year", "sumallPC1")],
                          by = c("year" = "eco_year"))


png("TDN overtime(all locations).png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

## Average annual TDN for all lake locations

plot(NA, NA, xlim = c(1998, 2016), ylim = c(0, 32), xlab = "Year",
     ylab = "Average annual TDN in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(1998:2016), cex.axis=1.5)

lines(y = WCdat_IN_ave$TDN, x= WCdat_IN_ave$year, lty =1, lwd = 1, col= "lightseagreen")
lines(y = WCdat_out_ave$TDN,  x= WCdat_out_ave$year, lty =1, lwd = 1, col= "cadetblue3")
lines(y = WCdat_0m_ave$TDN,  x= WCdat_0m_ave$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WCdat_3m_ave$TDN,  x= WCdat_3m_ave$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WCdat_9m_ave$TDN,  x= WCdat_9m_ave$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WCdat_IN_ave$TDN, x= WCdat_IN_ave$year, col= "lightseagreen", pch=19) 
arrows(WCdat_IN_ave$year, (WCdat_IN_ave$TDN + WCdat_IN_ave$SE), 
       WCdat_IN_ave$year, (WCdat_IN_ave$TDN - WCdat_IN_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WCdat_out_ave$TDN, x= WCdat_out_ave$year, col= "cadetblue3", pch=19) 
arrows(WCdat_out_ave$year, (WCdat_out_ave$TDN + WCdat_out_ave$SE), 
       WCdat_out_ave$year, (WCdat_out_ave$TDN - WCdat_out_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")

points(y = WCdat_0m_ave$TDN, x= WCdat_0m_ave$year, col= "darkolivegreen3", pch=19) 
arrows(WCdat_0m_ave$year, (WCdat_0m_ave$TDN + WCdat_0m_ave$SE), 
       WCdat_0m_ave$year, (WCdat_0m_ave$TDN - WCdat_0m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WCdat_3m_ave$TDN, x= WCdat_3m_ave$year, col= "darkolivegreen4", pch=19) 
arrows(WCdat_3m_ave$year, (WCdat_3m_ave$TDN + WCdat_3m_ave$SE), 
       WCdat_3m_ave$year, (WCdat_3m_ave$TDN - WCdat_3m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WCdat_9m_ave$TDN, x= WCdat_9m_ave$year, col= "darkgreen", pch=19) 
arrows(WCdat_9m_ave$year, (WCdat_9m_ave$TDN + WCdat_9m_ave$SE), 
       WCdat_9m_ave$year, (WCdat_9m_ave$TDN - WCdat_9m_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topright", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

dev.off()


# plot TDN for in lake (surface, meta, hypo) averages and pca for extened summer add in linear model 
PCA_allinlake_ave_TDN_plot <- qplot(sumallPC1, TDN, data = all_TDN_combined_dat, geom="point", ylab = "Average annual TDN", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

PCA_allinlake_ave_TDN_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))


ggsave("PCA_allinlake_aveTDN.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

TDN.mod_inlake <- lmer(TDN ~ scale(sumallPC1) + (1|depth), data = all_TDN_combined_dat)
summary(TDN.mod_inlake) # year average not significant

######### Max TDN not Average 
# average annual nitrate at inlet
WCdat_IN$TD
WCdat_IN_max <- aggregate(TDN ~ year, data= WCdat_IN, FUN=max) 
WCdat_IN_max$SE <- aggregate(TDN ~ year, data=WCdat_IN, FUN=se)[,2]

WCdat_IN_max$depth <- "IN"

# average annual chla at outlet
WCdat_out_max <- aggregate(TDN ~ year, data= WCdat_OUT, FUN=max) 
WCdat_out_max$SE <- aggregate(TDN ~ year, data =WCdat_OUT, FUN=se)[,2]

WCdat_out_max$depth <- "OUT"

# average annual chla at lake surface
WCdat_0m_max <- aggregate(TDN ~ year, data= WCdat_0m, FUN=max) 
WCdat_0m_max$SE <- aggregate(TDN ~ year, data=WCdat_0m, FUN=se)[,2]

WCdat_0m_max$depth <- "0m"

# average annual chla at metalimnion
WCdat_3m_max <- aggregate(TDN ~ year, data= WCdat_3m, FUN=max) 
WCdat_3m_max$SE <- aggregate(TDN ~ year, data= WCdat_3m, FUN=se)[,2]

WCdat_3m_max$depth <- "3m"

# average annual chla at hypolimnion 
WCdat_9m_max <- aggregate(TDN ~ year, data= WCdat_9m, FUN=max) 
WCdat_9m_max$SE <- aggregate(TDN ~ year, data=WCdat_9m, FUN=se)[,2]

WCdat_9m_max$depth <- "9m"

max_annual_TDN <- rbind(WCdat_IN_max, WCdat_out_max, 
                        WCdat_0m_max, WCdat_3m_max, WCdat_9m_max)

## join nwt climate dataset (year and sumallPC1) with water quality dataset
all_TDN_combined_dat_max <- left_join(max_annual_TDN, dt_summer[c("eco_year", "sumallPC1")],
                                  by = c("year" = "eco_year"))

## Annual max TDN for all lake locations

png("TDN max overtime(all locations).png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

plot(NA, NA, xlim = c(1998, 2016), ylim = c(0, 45), xlab = "Year",
     ylab = "Max annual TDN in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(1998:2016), cex.axis=1.5)

lines(y = WCdat_IN_max$TDN, x= WCdat_IN_max$year, lty =1, lwd = 1, col= "lightseagreen")
lines(y = WCdat_out_max$TDN,  x= WCdat_out_max$year, lty =1, lwd = 1, col= "cadetblue3")
lines(y = WCdat_0m_max$TDN,  x= WCdat_0m_max$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WCdat_3m_max$TDN,  x= WCdat_3m_max$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WCdat_9m_max$TDN,  x= WCdat_9m_max$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WCdat_IN_max$TDN, x= WCdat_IN_max$year, col= "lightseagreen", pch=19) 
arrows(WCdat_IN_max$year, (WCdat_IN_max$TDN + WCdat_IN_max$SE), 
       WCdat_IN_max$year, (WCdat_IN_max$TDN - WCdat_IN_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WCdat_out_max$TDN, x= WCdat_out_max$year, col= "cadetblue3", pch=19) 
arrows(WCdat_out_max$year, (WCdat_out_max$TDN + WCdat_out_max$SE), 
       WCdat_out_max$year, (WCdat_out_max$TDN - WCdat_out_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")

points(y = WCdat_0m_max$TDN, x= WCdat_0m_max$year, col= "darkolivegreen3", pch=19) 
arrows(WCdat_0m_max$year, (WCdat_0m_max$TDN + WCdat_0m_max$SE), 
       WCdat_0m_max$year, (WCdat_0m_max$TDN - WCdat_0m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WCdat_3m_max$TDN, x= WCdat_3m_max$year, col= "darkolivegreen4", pch=19) 
arrows(WCdat_3m_max$year, (WCdat_3m_max$TDN + WCdat_3m_max$SE), 
       WCdat_3m_max$year, (WCdat_3m_max$TDN - WCdat_3m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WCdat_9m_max$TDN, x= WCdat_9m_max$year, col= "darkgreen", pch=19) 
arrows(WCdat_9m_max$year, (WCdat_9m_max$TDN + WCdat_9m_max$SE), 
       WCdat_9m_max$year, (WCdat_9m_max$TDN - WCdat_9m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topright", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

dev.off()

# plot TDN for in lake (surface, meta, hypo) averages and pca for extened summer add in linear model 
PCA_allinlake_max_TDN_plot <- qplot(sumallPC1, TDN, data = all_TDN_combined_dat_max, geom="point", ylab = "Max annual TDN", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

PCA_allinlake_max_TDN_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))


ggsave("PCA_allinlake_maxTDN.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

TDN.mod_inlake_max <- lmer(TDN ~ scale(sumallPC1) + (1|depth), data = all_TDN_combined_dat_max)
summary(TDN.mod_inlake_max) # year average not significant

##### For Nitrate ######
# average annual nitrate at inlet
WCdat_IN$NO3.
WCdat_IN_NO3_ave <- aggregate(NO3. ~ year, data= WCdat_IN, FUN=mean) 
WCdat_IN_NO3_ave$SE <- aggregate(NO3. ~ year, data=WCdat_IN, FUN=se)[,2]

WCdat_IN_NO3_ave$depth <- "IN"

# average annual chla at outlet
WCdat_out_NO3_ave <- aggregate(NO3. ~ year, data= WCdat_OUT, FUN=mean) 
WCdat_out_NO3_ave$SE <- aggregate(NO3. ~ year, data =WCdat_OUT, FUN=se)[,2]

WCdat_out_NO3_ave$depth <- "OUT"

# average annual chla at lake surface
WCdat_0m_NO3_ave <- aggregate(NO3. ~ year, data= WCdat_0m, FUN=mean) 
WCdat_0m_NO3_ave$SE <- aggregate(NO3. ~ year, data=WCdat_0m, FUN=se)[,2]

WCdat_0m_NO3_ave$depth <- "0m"

# average annual chla at metalimnion
WCdat_3m_NO3_ave <- aggregate(NO3. ~ year, data= WCdat_3m, FUN=mean) 
WCdat_3m_NO3_ave$SE <- aggregate(NO3. ~ year, data= WCdat_3m, FUN=se)[,2]

WCdat_3m_NO3_ave$depth <- "3m"

# average annual chla at hypolimnion 
WCdat_9m_NO3_ave <- aggregate(NO3. ~ year, data= WCdat_9m, FUN=mean) 
WCdat_9m_NO3_ave$SE <- aggregate(NO3. ~ year, data=WCdat_9m, FUN=se)[,2]

WCdat_9m_NO3_ave$depth <- "9m"

ave_annual_NO3 <- rbind(WCdat_IN_NO3_ave, WCdat_out_NO3_ave, 
                        WCdat_0m_NO3_ave, WCdat_3m_NO3_ave, WCdat_9m_NO3_ave)

## join nwt climate dataset (year and sumallPC1) with water quality dataset
all_TDN_combined_dat <- left_join(ave_annual_TDN, dt_summer[c("eco_year", "sumallPC1")],
                                  by = c("year" = "eco_year"))

###### Nitrate plot
png("NO3 ave overtime(all locations).png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

plot(NA, NA, xlim = c(1998, 2016), ylim = c(0, 27), xlab = "Year",
     ylab = "Average annual NO3 in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(1998:2016), cex.axis=1.5)

lines(y = WCdat_IN_NO3_ave$NO3., x= WCdat_IN_NO3_ave$year, lty =1, lwd = 1, col= "lightseagreen")
lines(y = WCdat_out_NO3_ave$NO3.,  x= WCdat_out_NO3_ave$year, lty =1, lwd = 1, col= "cadetblue3")
lines(y = WCdat_0m_NO3_ave$NO3.,  x= WCdat_0m_NO3_ave$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WCdat_3m_NO3_ave$NO3.,  x= WCdat_3m_NO3_ave$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WCdat_9m_NO3_ave$NO3.,  x= WCdat_9m_NO3_ave$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WCdat_IN_NO3_ave$NO3., x= WCdat_IN_NO3_ave$year, col= "lightseagreen", pch=19) 
arrows(WCdat_IN_NO3_ave$year, (WCdat_IN_NO3_ave$NO3. + WCdat_IN_NO3_ave$SE), 
       WCdat_IN_NO3_ave$year, (WCdat_IN_NO3_ave$NO3. - WCdat_IN_NO3_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WCdat_out_NO3_ave$NO3., x= WCdat_out_NO3_ave$year, col= "cadetblue3", pch=19) 
arrows(WCdat_out_NO3_ave$year, (WCdat_out_NO3_ave$NO3. + WCdat_out_NO3_ave$SE), 
       WCdat_out_NO3_ave$year, (WCdat_out_NO3_ave$NO3. - WCdat_out_NO3_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")

points(y = WCdat_0m_NO3_ave$NO3., x= WCdat_0m_NO3_ave$year, col= "darkolivegreen3", pch=19) 
arrows(WCdat_0m_NO3_ave$year, (WCdat_0m_NO3_ave$NO3. + WCdat_0m_NO3_ave$SE), 
       WCdat_0m_NO3_ave$year, (WCdat_0m_NO3_ave$NO3. - WCdat_0m_NO3_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WCdat_3m_NO3_ave$NO3., x= WCdat_3m_NO3_ave$year, col= "darkolivegreen4", pch=19) 
arrows(WCdat_3m_NO3_ave$year, (WCdat_3m_NO3_ave$NO3. + WCdat_3m_NO3_ave$SE), 
       WCdat_3m_NO3_ave$year, (WCdat_3m_NO3_ave$NO3. - WCdat_3m_NO3_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WCdat_9m_NO3_ave$NO3., x= WCdat_9m_NO3_ave$year, col= "darkgreen", pch=19) 
arrows(WCdat_9m_NO3_ave$year, (WCdat_9m_NO3_ave$NO3. + WCdat_9m_NO3_ave$SE), 
       WCdat_9m_NO3_ave$year, (WCdat_9m_NO3_ave$NO3. - WCdat_9m_NO3_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topright", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

dev.off()

ave_annual_TDN$NO3.
## join nwt climate dataset (year and sumallPC1) with water quality dataset
all_NO3_combined_dat <- left_join(ave_annual_TDN, dt_summer[c("eco_year", "sumallPC1")],
                                  by = c("year" = "eco_year"))

# plot TDN for in lake (surface, meta, hypo) averages and pca for extened summer add in linear model 
PCA_allinlake_aveNO3_plot <- qplot(sumallPC1, NO3., data = all_NO3_combined_dat, geom="point", ylab = "Average annual NO3", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

PCA_allinlake_aveNO3_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))


##### max nitrate #######
##### For Nitrate ######
# average annual nitrate at inlet
WCdat_IN$NO3.
WCdat_IN_NO3_max <- aggregate(NO3. ~ year, data= WCdat_IN, FUN=max) 
WCdat_IN_NO3_max$SE <- aggregate(NO3. ~ year, data=WCdat_IN, FUN=se)[,2]

WCdat_IN_NO3_max$depth <- "IN"

# average annual chla at outlet
WCdat_out_NO3_max <- aggregate(NO3. ~ year, data= WCdat_OUT, FUN=max) 
WCdat_out_NO3_max$SE <- aggregate(NO3. ~ year, data =WCdat_OUT, FUN=se)[,2]

WCdat_out_NO3_max$depth <- "OUT"

# average annual chla at lake surface
WCdat_0m_NO3_max <- aggregate(NO3. ~ year, data= WCdat_0m, FUN=max) 
WCdat_0m_NO3_max$SE <- aggregate(NO3. ~ year, data=WCdat_0m, FUN=se)[,2]

WCdat_0m_NO3_max$depth <- "0m"

# average annual chla at metalimnion
WCdat_3m_NO3_max <- aggregate(NO3. ~ year, data= WCdat_3m, FUN=max) 
WCdat_3m_NO3_max$SE <- aggregate(NO3. ~ year, data= WCdat_3m, FUN=se)[,2]

WCdat_3m_NO3_max$depth <- "3m"

# average annual chla at hypolimnion 
WCdat_9m_NO3_max <- aggregate(NO3. ~ year, data= WCdat_9m, FUN=max) 
WCdat_9m_NO3_max$SE <- aggregate(NO3. ~ year, data=WCdat_9m, FUN=se)[,2]

WCdat_9m_NO3_max$depth <- "9m"

ave_annual_NO3 <- rbind(WCdat_IN_NO3_max, WCdat_out_NO3_max, 
                        WCdat_0m_NO3_max, WCdat_3m_NO3_max, WCdat_9m_NO3_max)

## join nwt climate dataset (year and sumallPC1) with water quality dataset
all_NO3_max_combined_dat <- left_join(ave_annual_NO3, dt_summer[c("eco_year", "sumallPC1")],
                                  by = c("year" = "eco_year"))

###### Nitrate plot
png("NO3 max overtime(all locations).png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

plot(NA, NA, xlim = c(1998, 2016), ylim = c(0, 40), xlab = "Year",
     ylab = "Max annual NO3 in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(1998:2016), cex.axis=1.5)

lines(y = WCdat_IN_NO3_max$NO3., x= WCdat_IN_NO3_max$year, lty =1, lwd = 1, col= "lightseagreen")
lines(y = WCdat_out_NO3_max$NO3.,  x= WCdat_out_NO3_max$year, lty =1, lwd = 1, col= "cadetblue3")
lines(y = WCdat_0m_NO3_max$NO3.,  x= WCdat_0m_NO3_max$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WCdat_3m_NO3_max$NO3.,  x= WCdat_3m_NO3_max$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WCdat_9m_NO3_max$NO3.,  x= WCdat_9m_NO3_max$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WCdat_IN_NO3_max$NO3., x= WCdat_IN_NO3_max$year, col= "lightseagreen", pch=19) 
arrows(WCdat_IN_NO3_max$year, (WCdat_IN_NO3_max$NO3. + WCdat_IN_NO3_max$SE), 
       WCdat_IN_NO3_max$year, (WCdat_IN_NO3_max$NO3. - WCdat_IN_NO3_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WCdat_out_NO3_max$NO3., x= WCdat_out_NO3_max$year, col= "cadetblue3", pch=19) 
arrows(WCdat_out_NO3_max$year, (WCdat_out_NO3_max$NO3. + WCdat_out_NO3_max$SE), 
       WCdat_out_NO3_max$year, (WCdat_out_NO3_max$NO3. - WCdat_out_NO3_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")

points(y = WCdat_0m_NO3_max$NO3., x= WCdat_0m_NO3_max$year, col= "darkolivegreen3", pch=19) 
arrows(WCdat_0m_NO3_max$year, (WCdat_0m_NO3_max$NO3. + WCdat_0m_NO3_max$SE), 
       WCdat_0m_NO3_max$year, (WCdat_0m_NO3_max$NO3. - WCdat_0m_NO3_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WCdat_3m_NO3_max$NO3., x= WCdat_3m_NO3_max$year, col= "darkolivegreen4", pch=19) 
arrows(WCdat_3m_NO3_max$year, (WCdat_3m_NO3_max$NO3. + WCdat_3m_NO3_max$SE), 
       WCdat_3m_NO3_max$year, (WCdat_3m_NO3_max$NO3. - WCdat_3m_NO3_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WCdat_9m_NO3_max$NO3., x= WCdat_9m_NO3_max$year, col= "darkgreen", pch=19) 
arrows(WCdat_9m_NO3_max$year, (WCdat_9m_NO3_max$NO3. + WCdat_9m_NO3_max$SE), 
       WCdat_9m_NO3_max$year, (WCdat_9m_NO3_max$NO3. - WCdat_9m_NO3_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topright", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

dev.off()

ave_annual_TDN$NO3.
## join nwt climate dataset (year and sumallPC1) with water quality dataset
all_NO3_combined_dat_max <- left_join(all_NO3_max_combined_dat, dt_summer[c("eco_year", "sumallPC1")],
                                  by = c("year" = "eco_year"))

# plot TDN for in lake (surface, meta, hypo) averages and pca for extened summer add in linear model 
PCA_allinlake_maxNO3_plot <- qplot(sumallPC1, NO3., data = all_NO3_max_combined_dat, geom="point", ylab = "Max annual NO3", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

PCA_allinlake_maxNO3_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))


ggsave("PCA_allinlake_maxNO3.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

NO3.mod_inlake_max <- lmer(NO3. ~ scale(sumallPC1) + (1|depth), data = all_NO3_max_combined_dat)
summary(NO3.mod_inlake_max) # year average not significant


