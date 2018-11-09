### Long-Term Chlorophyll-a plots for all locations within GL4 ###

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

## annual average all chlor-a data aty GL4
WQdat_all_ave <- aggregate(chl_a ~ year, data=WQdat1, FUN=mean) 
WQdat_all_ave$SE <- aggregate(chl_a ~ year, data=WQdat1, FUN=se)[,2]


# average annual chla at inlet
WQdat_IN_ave <- aggregate(chl_a ~ year, data= WQdat_in_lake, FUN=mean) 
WQdat_IN_ave$SE <- aggregate(chl_a ~ year, data=WQdat_in_lake, FUN=se)[,2]

WQdat_IN_ave$depth <- "IN"

# average annual chla at outlet
WQdat_out_ave <- aggregate(chl_a ~ year, data= WQdat_out_lake, FUN=mean) 
WQdat_out_ave$SE <- aggregate(chl_a ~ year, data=WQdat_out_lake, FUN=se)[,2]

WQdat_out_ave$depth <- "OUT"

# average annual chla at lake surface
WQdat_0m_ave <- aggregate(chl_a ~ year, data= WQdat_lake_0m, FUN=mean) 
WQdat_0m_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake_0m, FUN=se)[,2]

WQdat_0m_ave$depth <- "0m"

# average annual chla at metalimnion
WQdat_3m_ave <- aggregate(chl_a ~ year, data= WQdat_lake_3m, FUN=mean) 
WQdat_3m_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake_3m, FUN=se)[,2]

WQdat_3m_ave$depth <- "3m"

# average annual chla at hypolimnion 
WQdat_9m_ave <- aggregate(chl_a ~ year, data= WQdat_lake_9m, FUN=mean) 
WQdat_9m_ave$SE <- aggregate(chl_a ~ year, data=WQdat_lake_9m, FUN=se)[,2]

WQdat_9m_ave$depth <- "9m"


chlor_ave_all <- rbind(WQdat_IN_ave, WQdat_out_ave, WQdat_0m_ave, WQdat_3m_ave, WQdat_9m_ave)
chlor_ave_inlake <- rbind(WQdat_0m_ave, WQdat_3m_ave, WQdat_9m_ave)
in_out_ave <- rbind(WQdat_IN_ave, WQdat_out_ave)

### Combine extended summer PCA with subsetted data ###

##chlora for all GL4 (Surface, Meta, Hypo, IN and OUT) data and PCA
comb_dat_WQ <- left_join(WQdat1, dt_summer[c("eco_year", "sumallPC1")],
                         by = c("year" = "eco_year"))


## chlora for all inlet lake data and PCA
comb_dat_in_WQ <- left_join(WQdat_in_lake, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## chlora for inlet and outlet data and PCA
comb_dat_in_out_WQ <- left_join(in_out_wc, dt_summer[c("eco_year", "sumallPC1")],
                                by = c("year" = "eco_year"))




############################################

## Annual Average all chlora at GL4 and  PCA
comb_dat_all_WQ <- left_join(WQdat_all_ave, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## Annual all inlake chlora  with depth at GL4 and  PCA
comb_dat_inlake_WQ <- left_join(chlor_ave_inlake, dt_summer[c("eco_year", "sumallPC1")],
                             by = c("year" = "eco_year"))

## Annual Average chlora for hypo and PCA
comb_dat_9m_WQ <- left_join(WQdat_9m_ave, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## Annual Average chlora for meta and PCA
comb_dat_3m_WQ <- left_join(WQdat_3m_ave, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## Annual Average chlora for surface and PCA
comb_dat_0m_WQ <- left_join(WQdat_0m_ave, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## Annual Average chlora for IN and OUT and PCA
comb_dat_INOUT_WQ <- left_join(in_out_ave, dt_summer[c("eco_year", "sumallPC1")],
                            by = c("year" = "eco_year"))

## chlora for all average in lake data and PCA
comb_dat_lake_WQ <- left_join(chlor_ave_all, dt_summer[c("eco_year", "sumallPC1")],
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

#png("Chlor_A overtime(IN and OUT).png",  
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

points(y = WQdat_IN_ave$chl_a, x= WQdat_IN_ave$year, col= "lightseagreen", pch=19) 
arrows(WQdat_IN_ave$year, (WQdat_IN_ave$chl_a + WQdat_IN_ave$SE), 
       WQdat_IN_ave$year, (WQdat_IN_ave$chl_a - WQdat_IN_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WQdat_out_ave$chl_a, x= WQdat_out_ave$year, col= "cadetblue3", pch=19) 
arrows(WQdat_out_ave$year, (WQdat_out_ave$chl_a + WQdat_out_ave$SE), 
       WQdat_out_ave$year, (WQdat_out_ave$chl_a - WQdat_out_ave$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")


legend("topleft", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

#dev.off()




# plot chlor-a for in lake averages and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_all_WQ, geom="point", ylab = "Average annual chlor-a") + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

#ggsave("PCA_all_ave_chla.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_1yr <- lmer(chl_a ~ scale(sumallPC1) + (1|depth), data = comb_dat_all_WQ)
summary(chlora.mod_1yr) # year average not significant


# plot chlor-a for in lake hypo averages and pca for extened summer add in linear model 
PCA_hypo_ave_chla<- qplot(sumallPC1, chl_a, data = comb_dat_9m_WQ, geom="point", ylab = "Average annual chlor-a (Hypo)", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())
PCA_hypo_ave_chla + scale_colour_manual(values = c("darkgreen"))


ggsave("PCA_hypo_ave_chla.pdf", scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_9myr <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_9m_WQ)
summary(chlora.mod_9myr) # year average not significant



# plot chlor-a for in lake (surface, meta, hypo) averages and pca for extened summer add in linear model 
PCA_allinlake_ave_chla_plot <- qplot(sumallPC1, chl_a, data = comb_dat_inlake_WQ, geom="point", ylab = "Average annual chlor-a (Surface, Meta, Hypo)", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

PCA_allinlake_ave_chla_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))


ggsave("PCA_allinlake_ave_chla.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_inlake <- lmer(chl_a ~ scale(sumallPC1) + (1|depth), data = comb_dat_inlake_WQ)
summary(chlora.mod_inlake) # year average not significant

# plot chlor-a for in lake surface averages and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_0m_WQ, geom="point", ylab = "Average annual chlor-a (Surface)") + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

#ggsave("PCA_surface_ave_chla.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_0m <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_0m_WQ)
summary(chlora.mod_0m)


# plot chlor-a for in lake meta averages and pca for extened summer add in linear model 
qplot(sumallPC1, chl_a, data = comb_dat_3m_WQ, geom="point", ylab = "Average annual chlor-a (Meta)") + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())

#ggsave("PCA_meta_ave_chla.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_3m <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_3m_WQ)
summary(chlora.mod_3m)

# plot chlor-a for IN and OUT averages and pca for extened summer add in linear model 
PCA_INOUT_ave_chla <- qplot(sumallPC1, chl_a, data = comb_dat_INOUT_WQ, geom="point", ylab = "Average annual chlor-a (Meta)", color=depth) + 
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text())
PCA_INOUT_ave_chla+scale_colour_manual(values = c("lightseagreen", "cadetblue3"))


ggsave("PCA_INOUT_ave_chla.pdf", scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_IN_OUT <- lmer(chl_a ~ scale(sumallPC1) + (1|depth), data = comb_dat_INOUT_WQ)
summary(chlora.mod_IN_OUT)




# plot chlor-a for all averages and pca for extened summer add in linear model 
ave_chl_plot <- qplot(sumallPC1, chl_a, data = comb_dat_lake_WQ, geom="point", ylab = "Average annual chlor-a (Meta)", color=depth) + 
  geom_smooth(method = "lm", color=depth) +
  theme_classic() + theme(axis.text.x = element_text())

ave_chl_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))

ggsave("PCA_all_ave_chla.pdf",
      scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_ave_all <- glm(chl_a ~ scale(sumallPC1) + depth, data = comb_dat_lake_WQ)
summary(chlora.mod_ave_all)

chlora.mod_ave_all <- lmer(chl_a ~ scale(sumallPC1) + (1|depth), data = comb_dat_lake_WQ)
summary(chlora.mod_ave_all)
















##### Averages not giving conclusive resluts. 
#### Look at max values 

## annual average all chlor-a data aty GL4
WQdat_all_ave <- aggregate(chl_a ~ year, data=WQdat1, FUN=max) 
WQdat_all_ave$SE <- aggregate(chl_a ~ year, data=WQdat1, FUN=se)[,2]


## average chlor-a for all inlake depths by year 
WQdat_lake_max <- aggregate(chl_a ~ year, data=WQdat_lake, FUN=max) 
WQdat_lake_max$SE <- aggregate(chl_a ~ year, data=WQdat_lake, FUN=se)[,2]


# average annual chla at inlet
WQdat_IN_max <- aggregate(chl_a ~ year, data= WQdat_in_lake, FUN=max) 
WQdat_IN_max$SE <- aggregate(chl_a ~ year, data=WQdat_in_lake, FUN=se)[,2]

WQdat_IN_max$depth <- "IN"

# average annual chla at outlet
WQdat_out_max <- aggregate(chl_a ~ year, data= WQdat_out_lake, FUN=max) 
WQdat_out_max$SE <- aggregate(chl_a ~ year, data=WQdat_out_lake, FUN=se)[,2]

WQdat_out_max$depth <- "OUT"

# average annual chla at lake surface
WQdat_0m_max <- aggregate(chl_a ~ year, data= WQdat_lake_0m, FUN=max) 
WQdat_0m_max$SE <- aggregate(chl_a ~ year, data=WQdat_lake_0m, FUN=se)[,2]

WQdat_0m_max$depth <- "0m"

# average annual chla at metalimnion
WQdat_3m_max <- aggregate(chl_a ~ year, data= WQdat_lake_3m, FUN=max) 
WQdat_3m_max$SE <- aggregate(chl_a ~ year, data=WQdat_lake_3m, FUN=se)[,2]

WQdat_3m_max$depth <- "3m"

# average annual chla at hypolimnion 
WQdat_9m_max<- aggregate(chl_a ~ year, data= WQdat_lake_9m, FUN=max) 
WQdat_9m_max$SE <- aggregate(chl_a ~ year, data=WQdat_lake_9m, FUN=se)[,2]

WQdat_9m_max$depth <- "9m"

# rbind all aggregated max values
WQdat_max_alldepth <- rbind(WQdat_IN_max, WQdat_out_max, WQdat_0m_max, WQdat_3m_max, WQdat_9m_max)

## Annual max all chlora at GL4 and  PCA
comb_dat_all_max_WQ <- left_join(WQdat_max_alldepth, dt_summer[c("eco_year", "sumallPC1")],
                             by = c("year" = "eco_year"))


#########################################
####### Visulalize Max Chla##############

png("Chlor_A-max_overtime(all locations).png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

## Average annual chla for all lake locations

plot(NA, NA, xlim = c(2000, 2017), ylim = c(0, 21), xlab = "Year",
     ylab = "Max annual chlor-A in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2000:2017), cex.axis=1.5)

lines(y = WQdat_IN_max$chl_a, x= WQdat_IN_max$year, lty =1, lwd = 1, col= "lightseagreen")
lines(y = WQdat_out_max$chl_a,  x= WQdat_out_max$year, lty =1, lwd = 1, col= "cadetblue3")
lines(y = WQdat_0m_max$chl_a,  x= WQdat_0m_max$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WQdat_3m_max$chl_a,  x= WQdat_3m_max$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WQdat_9m_max$chl_a,  x= WQdat_9m_max$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WQdat_IN_max$chl_a, x= WQdat_IN_max$year, col= "lightseagreen", pch=19) 
arrows(WQdat_IN_max$year, (WQdat_IN_max$chl_a + WQdat_IN_max$SE), 
       WQdat_IN_max$year, (WQdat_IN_max$chl_a - WQdat_IN_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "lightseagreen")


points(y = WQdat_out_max$chl_a, x= WQdat_out_max$year, col= "cadetblue3", pch=19) 
arrows(WQdat_out_max$year, (WQdat_out_max$chl_a + WQdat_out_max$SE), 
       WQdat_out_max$year, (WQdat_out_max$chl_a - WQdat_out_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "cadetblue3")

points(y = WQdat_0m_max$chl_a, x= WQdat_0m_max$year, col= "darkolivegreen3", pch=19) 
arrows(WQdat_0m_max$year, (WQdat_0m_max$chl_a + WQdat_0m_max$SE), 
       WQdat_0m_max$year, (WQdat_0m_max$chl_a - WQdat_0m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WQdat_3m_max$chl_a, x= WQdat_3m_max$year, col= "darkolivegreen4", pch=19) 
arrows(WQdat_3m_max$year, (WQdat_3m_max$chl_a + WQdat_3m_max$SE), 
       WQdat_3m_max$year, (WQdat_3m_max$chl_a - WQdat_3m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WQdat_9m_max$chl_a, x= WQdat_9m_max$year, col= "darkgreen", pch=19) 
arrows(WQdat_9m_max$year, (WQdat_9m_max$chl_a + WQdat_9m_max$SE), 
       WQdat_9m_max$year, (WQdat_9m_max$chl_a - WQdat_9m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


legend("topleft", legend=c("Inlet", "Outlet", "0m", "3m", "9m"),
       col=c("lightseagreen", "cadetblue3", "darkolivegreen3", "darkolivegreen4", "darkgreen" ), lty=1, cex=1.2)

dev.off()


# plot chlor-a for max chlora and pca for extened summer add in linear model 
pac_max_plot <- qplot(sumallPC1, chl_a, data = comb_dat_all_max_WQ, geom="point", ylab = "Maximum annual chlor-a", color=depth)  +
  geom_smooth(method = "lm") +
  theme_classic() + theme(axis.text.x = element_text()) 

pac_max_plot + scale_colour_manual(values = c("darkolivegreen3", "darkolivegreen4", "darkgreen", "lightseagreen", "cadetblue3"))

ggsave("PCA_max_chla.pdf",
       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)

chlora.mod_all_max <- glm(chl_a ~ scale(sumallPC1), data = comb_dat_all_max_WQ)
summary(chlora.mod_all_max)




##########################
## For Zoop comparison ##
#########################


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

# function to calculate standard error around the mean 
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

## annual max all zoop densities GL4
zoop_den_max <- aggregate(tot_density ~ year, data=dt1, FUN=max) 
zoop_den_max$SE <- aggregate(tot_density ~ year, data=dt1, FUN=se)[,2]
names(zoop_den_max)[2]<-"taxon_dens"
zoop_den_max$label <- "all"


png("Zoop_mden_chlom_0vertime.png",  
    width = 5.25,
    height = 3.25,
    units = "in",
    res = 1200,
    pointsize = 4 )


plot(NA, NA, xlim = c(2000, 2017), ylim = c(0, 21), xlab = "Year",
     ylab = "Max annual chlor-A in GL4", cex.lab=1.5, cex.axis=1.2, 
     xaxt="n")
axis(1, at=c(2000:2017), cex.axis=1.5)

lines(y = WQdat_0m_max$chl_a,  x= WQdat_0m_max$year, lty =1, lwd = 1, col= "darkolivegreen3")
lines(y = WQdat_3m_max$chl_a,  x= WQdat_3m_max$year, lty =1, lwd = 1, col= "darkolivegreen4")
lines(y = WQdat_9m_max$chl_a,  x= WQdat_9m_max$year, lty =1, lwd = 1, col= "darkgreen")


points(y = WQdat_0m_max$chl_a, x= WQdat_0m_max$year, col= "darkolivegreen3", pch=19) 
arrows(WQdat_0m_max$year, (WQdat_0m_max$chl_a + WQdat_0m_max$SE), 
       WQdat_0m_max$year, (WQdat_0m_max$chl_a - WQdat_0m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen3")

points(y = WQdat_3m_max$chl_a, x= WQdat_3m_max$year, col= "darkolivegreen4", pch=19) 
arrows(WQdat_3m_max$year, (WQdat_3m_max$chl_a + WQdat_3m_max$SE), 
       WQdat_3m_max$year, (WQdat_3m_max$chl_a - WQdat_3m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkolivegreen4")

points(y = WQdat_9m_max$chl_a, x= WQdat_9m_max$year, col= "darkgreen", pch=19) 
arrows(WQdat_9m_max$year, (WQdat_9m_max$chl_a + WQdat_9m_max$SE), 
       WQdat_9m_max$year, (WQdat_9m_max$chl_a - WQdat_9m_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "darkgreen")


# add in zoop
## Allow a second plot on the same graph
par(new=TRUE)
plot(zoop_den_max$year, zoop_den_max$taxon_dens,type = "l", xaxt = "n", yaxt = "n",
     ylab = "", xlab = "", col = "tomato3", lty = 1)
axis(side = 4, ylim= c(0, 20))


lines(y = zoop_den_max$taxon_dens,  x= zoop_den_max$year, lty =1, lwd = 1, col= "tomato3")
points(y = zoop_den_max$taxon_dens, x= zoop_den_max$year, col= "tomato3", pch=19) 
arrows(zoop_den_max$year, (zoop_den_max$taxon_dens + zoop_den_max$SE), 
       zoop_den_max$year, (zoop_den_max$taxon_dens - zoop_den_max$SE), 
       length=0.0, code=3, lwd=0.5, col= "tomato3")

legend("topright", legend=c("0m", "3m", "9m", "Zoop density"),
       col=c("darkolivegreen3", "darkolivegreen4", "darkgreen", "tomato3" ), lty=1, cex=1.2)

dev.off()
