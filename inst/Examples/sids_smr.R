## This example deals with the use of mixture model to obtain empirical Bayes estimates
## for Poisson rates in rare disease mapping.
## The data consists of sudden infant death syndrome counts in 100 counties in north Caroline and are
## available at http://sal.agecon.uiuc.edu/datasets/sids.zip and also in the examples of
## Bugs (http://www.mrc-bsu.cam.ac.uk/bugs/). The version we use here is
## nc.sids in package=spdep since it is very complete with North Carolina maps and centroide.
## This example only use the 1974 data.                   

library(moc)                   ## first load the required moc library
data(nc.sids,package="spdep")  ## and the data.

## Poisson distribution is appropriate for counts of rare events.

poiss <- function(x,la,...) {dpois(x,la)}

## The mean of the poisson process is taken as the base rate (exp(p) > 0) times
## the expected number of cases which is proportional to the size of the county.

smr.ecases <- function(p,expected.cases) {exp(p)*cbind(expected.cases)}

nc.sids$ecases <- nc.sids$BIR74*sum(nc.sids$SID74)/sum(nc.sids$BIR74)

gmu.smr.sids <- list(G1=function(p) smr.ecases(p[1],nc.sids$ecases),
                     G2=function(p) smr.ecases(p[2],nc.sids$ecases),
                     G3=function(p) smr.ecases(p[3],nc.sids$ecases),
                     G4=function(p) smr.ecases(p[4],nc.sids$ecases))

## In fact it is coded such that the analysis corresponds to that of the
## standardized mortality rate (SMR = O/E : observed cases/expected cases)
## which we know is higly sensitive to the size of the county extreme SMR
## having much chances to appear in smaller counties. As illustrated by the
## following plot:
plot(nc.sids$ecases,nc.sids$SID74/nc.sids$ecases)

## Fitting mixture model would allow to account for over(under)-dispersion
## and to capture the range of regressions on the size of the counties.

nc.smr <- moc(nc.sids$SID74, density = poiss, gmu = gmu.smr.sids[1], 
    pgmu = log(c(1)))


nc.smr2 <- moc(nc.sids$SID74, density = poiss, groups=2, gmu = gmu.smr.sids[1:2], 
    pgmu = log(c(1,3)), pgmix=0)

nc.smr3 <- moc(nc.sids$SID74, density = poiss, groups=3, gmu = gmu.smr.sids[1:3], 
     pgmu = log(c(0.9,2.5,5)), pgmix=c(0,0))

nc.smr4 <- moc(nc.sids$SID74, density = poiss, groups=4, gmu = gmu.smr.sids, 
    pgmu = log(c(0.2,1,6,10)), pgmix=c(0,0,0))

AIC(nc.smr,nc.smr2,nc.smr3,nc.smr4,k="BIC")
nc.smr2

## As we can see the 2 groups models is sufficient to capture all essential heterogeneity.
plot(nc.sids$ecases,nc.sids$SID74)
lines(c(0,max(nc.sids$ecases)),exp(nc.smr2$coef[1])*c(0,max(nc.sids$ecases)))
lines(c(0,max(nc.sids$ecases)),exp(nc.smr2$coef[2])*c(0,max(nc.sids$ecases)))

## As ce can see the empirical posterior mean provides a good shrinkage of the
## SMR estimates which is less dependent on county size
plot(nc.sids$ecases,nc.sids$SID74/nc.sids$ecases)
points(nc.sids$ecases,apply(post(nc.smr2),1,function(x) sum(x*exp(nc.smr2$coef[1:2]))),col="red",pch=5)
lines(c(0,max(nc.sids$ecases)),exp(nc.smr2$coef[1])*c(1,1))
lines(c(0,max(nc.sids$ecases)),exp(nc.smr2$coef[2])*c(1,1))

## Using the posterior empirical mean estimates is equivalent to the use
## of the posterior probabilities that we can use to mix basic colors
## to map on the counties.
smr2.col <- mix.colors.moc(nc.smr2,group.colors=c("blue","red"))

library(maptools)   
plot(sidspolys, main="SMR_EB of SIDS in North Carolina",border="grey",smr2.col,forcefill=FALSE)
## If you don't use the library maptools the following commands are equivalent.
## plot(attr(sidspolys, "maplim"),type="n",asp=1)
## for(i in 1:100) polygon(sidspolys[[i]], border="grey",col=smr2.col[i])
legend(-84,34,format(apply(cbind(c(0,0.25,0.5,0.75,1),c(1,0.75,0.5,0.25,0)),1,function(x) sum(x*exp(nc.smr2$coef[1:2]))),
                     digits=4),fill=c(rgb(1,0,0),rgb(0.75,0,0.25),rgb(0.5,0,0.5),rgb(0.25,0,0.75),rgb(0,0,1)))

## As we can see there are some spatial heterogeneity in the spatial distribution of sids.
## Part of it can be due to some remaining effect of county size, however we can
## use the spatial structure to account for possible spatial dependency and achieve
## better smoothing.

## A simple way to do this is to make the mixture probabilities of each counties depend on the
## SMR of the neighborhood counties.

nc.sids.smr.nb <- sapply(ncCR85.nb,function(nb) sum(nc.sids[nb,"SID74"])/sum(nc.sids[nb,"BIR74"]))/sum(nc.sids$SID74)*sum(nc.sids$BIR74)

smr.mix.nb <- function(pmix) t(apply(cbind(pmix[1]+pmix[2]*nc.sids.smr.nb),1,inv.glogit))


nc.smr2.nb <- moc(nc.sids$SID74, density = poiss, groups=2, gmu = gmu.smr.sids[1:2], 
    gmixture=smr.mix.nb,pgmu =  nc.smr2$coef[1:2], pgmix=c(nc.smr2$coef[3],0.1))

## The positive coefficient of regression on neighborhood counties has a small p-value
## and the decrease in -likelihood is important showing the presence spatial dependance.
nc.smr2.nb
AIC(nc.smr2,nc.smr2.nb,k="BIC")

## Mapping the empirical posterior SMR estimates again
smr2.nb.col <- mix.colors.moc(nc.smr2.nb,group.colors=c("blue","red"))

plot(sidspolys, main="SMR_EB_NB of SIDS in North Carolina",border="grey",smr2.nb.col,forcefill=FALSE)
legend(-84,34,format(apply(cbind(c(0,0.25,0.5,0.75,1),c(1,0.75,0.5,0.25,0)),1,function(x) sum(x*exp(nc.smr2.nb$coef[1:2]))),
                     digits=4),fill=c(rgb(1,0,0),rgb(0.75,0,0.25),rgb(0.5,0,0.5),rgb(0.25,0,0.75),rgb(0,0,1)))

## now shows more spatial homogeneity and the values shows that we have achived slightly better shrinkage.
 plot(nc.sids$ecases,nc.sids$SID74/nc.sids$ecases)
 points(nc.sids$ecases,apply(post(nc.smr2),1,function(x) sum(x*exp(nc.smr2$coef[1:2]))),col="red",pch=5)
 points(nc.sids$ecases,apply(post(nc.smr2.nb),1,function(x) sum(x*exp(nc.smr2.nb$coef[1:2]))),col="blue",pch=8)

## The preceeding model require the identification of the neighborhood of each county
## and restrict spatial dependance to adjacent counties. An alternative approach
## would be to consider the dependance of each county on other counties
## to be inversely proportional to their distance. A simple way to
## achieve this would be to compute for each county, the weighted SMR of the other counties
## with the weigths being the inverse of the squared distance between the counties.
sids.dist.nb <- t(apply(sidscents,1,function(x) apply((x-t(sidscents))^2,2,sum)))
for(i in 1:100) sids.dist.nb[i,-i] <-  1/sids.dist.nb[i,-i]/sum(1/sids.dist.nb[i,-i])

nc.smr.dist.nb <- apply(nc.sids$SID74*t(sids.dist.nb),2,sum)/apply(nc.sids$BIR74*t(sids.dist.nb),2,sum)/sum(nc.sids$SID74)*sum(nc.sids$BIR74)

smr.mix.dist.nb <- function(pmix) t(apply(cbind(pmix[1]+pmix[2]*nc.smr.dist.nb),1,inv.glogit))

nc.smr2.dist.nb <- moc(nc.sids$SID74, density = poiss, groups=2, gmu = gmu.smr.sids[1:2], 
    gmixture=smr.mix.dist.nb,pgmu =  nc.smr2$coef[1:2], pgmix=c(nc.smr2$coef[3],0.1))

nc.smr2.dist.nb
AIC(nc.smr2,nc.smr2.nb,nc.smr2.dist.nb,k="BIC")
entropy(nc.smr2,nc.smr2.nb,nc.smr2.dist.nb)

## As we can see the two spatial models are quite comparable, the second model showing
## a slightly smaller prior entropy denoting a stronger spatial prediction resulting
## in a little more homogeneity in the SMR mapping.
smr2.dist.nb.col <- mix.colors.moc(nc.smr2.dist.nb,group.colors=c("blue","red"))
plot(sidspolys, main="SMR_EB_NB of SIDS in North Carolina",border="grey",smr2.dist.nb.col,forcefill=FALSE)
legend(-84,34,format(apply(cbind(c(0,0.25,0.5,0.75,1),c(1,0.75,0.5,0.25,0)),1,function(x) sum(x*exp(nc.smr2.dist.nb$coef[1:2]))),
                     digits=4),fill=c(rgb(1,0,0),rgb(0.75,0,0.25),rgb(0.5,0,0.5),rgb(0.25,0,0.75),rgb(0,0,1)))

## The second model however uses a specific function of the distance between counties which may not
## always be appropriate.
