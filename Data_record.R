wsol <- read.csv(file.choose())
head(wsol)
names(wsol)
y1998 <- subset(wsol, year==1998,
              select=local_site:POC)
summary(y1998$year)
summary(y1998$date)
#1998 GL1= 1, GL3=22, GL4=18

y1999 <- subset(wsol, year==1999,
                select=local_site:POC)
summary(y1999$year)
summary(y1999$local_site)
#1999 GL3= 3, GL4=11


y2000 <- subset(wsol, year==2000,
                select=local_site:POC)
summary(y2000$year)
summary(y2000$local_site)
#2000 GL4=20

y2001 <- subset(wsol, year==2001,
                select=local_site:POC)
summary(y2001$year)
summary(y2001$local_site)

y2002 <- subset(wsol, year==2002,
                select=local_site:POC)
summary(y2002$year)
summary(y2002$local_site)

y2003 <- subset(wsol, year==2003,
                select=local_site:POC)
summary(y2003$year)
summary(y2003$local_site)

y2004 <- subset(wsol, year==2004,
                select=local_site:POC)
summary(y2004$year)
summary(y2004$local_site)

y2005 <- subset(wsol, year==2005,
                select=local_site:POC)
summary(y2005$year)
summary(y2005$local_site)

y2006 <- subset(wsol, year==2006,
                select=local_site:POC)
summary(y2006$year)
summary(y2006$local_site)

y2007 <- subset(wsol, year==2007,
                select=local_site:POC)
summary(y2007$year)
summary(y2007$local_site)

y2008 <- subset(wsol, year==2008,
                select=local_site:POC)
summary(y2008$year)
summary(y2008$local_site)

y2009 <- subset(wsol, year==2009,
                select=local_site:POC)
summary(y2009$year)
summary(y2009$local_site)

y2010 <- subset(wsol, year==2010,
                select=local_site:POC)
summary(y2010$year)
summary(y2010$local_site)

y2011 <- subset(wsol, year==2011,
                select=local_site:POC)
summary(y2011$year)
summary(y2011$local_site)

y2012 <- subset(wsol, year==2012,
                select=local_site:POC)
summary(y2012$year)
summary(y2012$local_site)

y2013 <- subset(wsol, year==2013,
                select=local_site:POC)
summary(y2013$year)
summary(y2013$local_site)

y2014 <- subset(wsol, year==2014,
                select=local_site:POC)
summary(y2014$year)
summary(y2014$local_site)

y2016 <- subset(wsol, year==2018,
                select=local_site:POC)
summary(y2016$year)
summary(y2016$local_site)

















