## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load, echo=FALSE---------------------------------------------------------
library(isoWater)

## ----values-------------------------------------------------------------------
wiDB_values(c("types", "countries"))

## ----query1-------------------------------------------------------------------
ls = wiDB_sites(countries = "US", types = "Lake_or_pond")
ps = wiDB_sites(countries = "US", types = "Precipitation")

## ----siteData, include=FALSE--------------------------------------------------
ls = readRDS("lsl.rds")
ps = readRDS("psl.rds")

## ----sites, fig.width=5, fig.asp=0.6------------------------------------------
omar = par("mar")
par(mar = c(4, 4, 1, 1))
plot(ls[, 2:1], xlim = c(-125, -68), ylim = c(25, 50))
points(ps[, 2:1], col = "red")
par(mar = omar)

## ----query2-------------------------------------------------------------------
ld = wiDB_data(minLat = 41, maxLat = 42, minLong = -94, maxLong = -93, types = "Lake")
pd = wiDB_data(minLat = 41, maxLat = 42, minLong = -94, maxLong = -93, types = "Precipitation")

## ----dataData, include = FALSE------------------------------------------------
ld = readRDS("ldl.rds")
pd = readRDS("pdl.rds")

## ----data---------------------------------------------------------------------
ld$data
pd$projects

## ----plot, fig.width=5, fig.asp=0.8-------------------------------------------
omar = par("mar")
par(mar = c(5, 5, 1, 1))
plot(pd$data$d18O, pd$data$d2H, xlab = expression(delta^{18}*"O"),
     ylab = expression(delta^{2}*"H"))
abline(10, 8)
points(ld$data$d18O, ld$data$d2H, col = "red")
par(mar = omar)

## ----mwl, fig.width=5, fig.asp=0.8--------------------------------------------
#extract the precipitation H and O isotope values
HO = data.frame(pd$data$d2H, pd$data$d18O)
MWL = mwl(HO)

## ----iso----------------------------------------------------------------------
#we will analyze one of our lake water samples, and assume realistic
#values of analytical uncertainty and no error covariance
obs = iso(ld$data$d2H[1], ld$data$d18O[1], 0.5, 0.1)

## ----ELslope------------------------------------------------------------------
#assumed value that is reasonable for lake water evaporation in this
#location/climate
slope = c(5.2, 0.3)

## ----mwlSource1, results='hide'-----------------------------------------------
mwls1 = mwlSource(obs, MWL, slope, ngens = 5e3)

## ----mwlSource1out------------------------------------------------------------
mwls1$summary

## ----traceplot, fig.width=5, fig.asp=0.8--------------------------------------
plot(mwls1$results$source_d2H[1:2500], type = "l", ylim = range(mwls1$results$source_d2H))
lines(mwls1$results$source_d2H[2501:5000], col=2)
lines(mwls1$results$source_d2H[5001:7500], col=3)

## ----mwlSource2, results='hide'-----------------------------------------------
mwls2 = mwlSource(obs, MWL, slope, stype = 2, ngens = 5e3)

## ----mulSource2out, fig.width=5, fig.asp=0.8----------------------------------
omar = par("mar")
par(mar = c(5, 5, 1, 1))
plot(HO[, 2:1], xlim = c(-20, 0), ylim = c(-140, 0),
     xlab = expression(delta^{18}*"O"),
     ylab = expression(delta^{2}*"H"))
abline(MWL[2], MWL[1])
points(obs[2:1], col = "red")
points(mwls1$results$source_d18O, mwls1$results$source_d2H, col = "blue")
points(mwls2$results$source_d18O, mwls2$results$source_d2H, col = "green")
par(mar = omar)

## ----sources------------------------------------------------------------------
#prep our data - we'll average the precipitation values by quarter
q = quarters(as.POSIXct(pd$data$Collection_Date))
qu = sort(unique(q))
ql = length(qu)

pd_q = data.frame("H" = numeric(ql), "O" = numeric(ql),
                  "Hsd" = numeric(ql), "Osd" = numeric(ql),
                  "HOc" = numeric(ql))
for(i in seq_along(qu)){
  pd_q$H[i] = mean(pd$data$d2H[q == qu[i]], na.rm = TRUE)
  pd_q$O[i] = mean(pd$data$d18O[q == qu[i]], na.rm = TRUE)
  pd_q$Hsd[i] = sd(pd$data$d2H[q == qu[i]], na.rm = TRUE)
  pd_q$Osd[i] = sd(pd$data$d18O[q == qu[i]], na.rm = TRUE)
  pd_q$HOc[i] = cov(pd$data$d18O[q == qu[i]], pd$data$d2H[q == qu[i]],
                    use = "pairwise.complete.obs")
}
pd_q

#make the iso object, providing the stats calculated above
sources = iso(pd_q$H, pd_q$O, pd_q$Hsd, pd_q$Osd, pd_q$HOc)

## ----mixSource1---------------------------------------------------------------
mixs1 = mixSource(obs, sources, slope, ngens = 2e4, ncores = 2)
mixs1$summary

## ----mixSource2, fig.width=6, fig.asp=0.7-------------------------------------
mixs2 = mixSource(obs, sources, slope, prior = c(2, 1, 1, 1), ngens = 2e4, ncores = 2)
boxplot(mixs1$results[, 3:6], outline = FALSE)
boxplot(mixs2$results[, 3:6], add = TRUE, col = rgb(1, 0, 0, 0.5),
        outline = FALSE)

