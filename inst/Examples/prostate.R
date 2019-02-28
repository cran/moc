## This example concerns the prostate clinical trial data of
## Byar and Green and is available at Statlib in 
## http://lib.stat.cmu.edu/datasets/Andrews/T46.1
## This is a large data set with a lot of variables of mixed
## types (continuous, categorical) and it was used by
## Jorgensen and Hunt to illustrate the MULTIMIX
##  (http://www.stats.waikato.ac.nz/Staff/maj/multimix)
## approach to cluster pre-trial variables.
## Here we use the more complete first version of this data set
## to illustrate the MULTIMIX approach in MOC. Note that the data set
## contain an AP (serum prostatic acid phosphatase) measurement of 9999
## which look like an extreme outlier and can be taught to be a miscoding
## of missing value which are -9999. However this value was retained in the
## analysis of the previous authors, AP is log transformed and this
## value then appears to be less influential. The previous authors also
## removed all subjects with missing values, which is highly questionable,
## here we are going to use all subjects.


#If you have an internet connection
#prostate <- read.table(file="http://lib.stat.cmu.edu/datasets/Andrews/T46.1",
#col.names=c("T46","One","id1","id2","stage", "rx", "u2", "u3", "u4", "dtime", "status", 
#	"age", "wti", "pf", "hx", "sbp", "dbp", "ekg", "hg", "sz", "sg", "ap", "bm"),
#na.strings=c("-9999",""),
#colClasses=c(rep("NULL",4),"numeric","factor",rep("numeric",4),"factor",rep("numeric",2),rep("factor",2),
#	rep("numeric",2),"factor",rep("numeric",4),"factor"))
#
#levels(prostate$rx) <- c("placebo","0.2 mg estrogen","1.0 mg estrogen","5.0 mg estrogen")
#
#levels(prostate$status) <- c("alive","dead - prostatic cancer","dead - heart or vascular",
#	"dead - cerebrovascular","dead - pulmonary embolus","dead - other cancer",
#	"dead - respiratory disease","dead - other non-cancer","dead - unspecified non-cancer","dead - unknown cause")
#
#levels(prostate$pf) <- c("normal activity","in bed < 50% daytime","in bed > 50% daytime","confined to bed")
#
#levels(prostate$hx) <- c("no","yes")
#
#levels(prostate$ekg) <- c("normal","benign", "rhythmic disturb & electrolyte ch",
#	"heart block or conduction def","heart strain",
#	"old MI","old MI")
#
#levels(prostate$bm) <- c("no","yes",NA)
#

prostate.file <- system.file("Examples","prostate.RData",package="moc",mustWork=TRUE)
load(prostate.file)
prostate <- transform(prostate,cons=1,log.ap=log(ap),sqrt.sz=sqrt(sz)) 

require(moc)    #load the essential

## We need the multinomial distribution for categorical variables

## The first version of the multinomial assumes that x is a set
## of indicator variables and p is the corresponding matrix
## of probabilities (one column for each category).
## It returns the joint density of the indicators.
dmnom <-                                                     
function(x,prob,...) {apply(prob*x+1-x,1,prod,na.rm=TRUE)}

mnom.mu <- function(p) {inv.glogit(p)}


## The second version of the multinomial assumes that x
## is the categorical variable coded from 1 to ncat.
## The extra parameter has ncat columns corresponding
## to each category.
dmnom.cat <- function(x,mu=1,shape=1,extra) {extra[,x]}


## We also need the univariate normal for continuous variables
dnormal <-
function(x,mu,sig,...) {dnorm(x,mu,sig)}

## and the multivariate normal too. Which is parametrized with the inverse of the Cholesky
## decomposition of the correlation matrix in the extra parameters.

mvnorm <- function(x,mu,sigma,extra)
{ 
  y <- (x-mu)/sigma
  ss <- diag(rep(1,dim(x)[[2]]))
  lind <- upper.tri(ss)
  for ( i in 1:dim(x)[1]) {
      ss[lind] <- extra[i,]
      na.ind <- is.na(y[i,])
      y[i,!na.ind] <- t(ss)[!na.ind,!na.ind]%*%y[i,!na.ind]
  }
  apply(dnorm(y)/sigma,1,prod,na.rm=TRUE)
}


mu.mvnorm <- function(pmu) rbind(pmu) #mvnorm mean base function

sig.mvnorm <- function(psh) rbind(exp(psh)) #mvnorm std.dev base function

extra.mvnorm <- function(pext) rbind(pext)  #mvnorm inverse of Cholesky (of correlation) base function

## The first model assumes local independence of all variables given the mixture group.

## We first construct the indicator variables for pf, hx, ekg and bm
prostate$pf.ind <- outer(unclass(prostate$pf),na.omit(unique(unclass(prostate$pf))),"==")*1
prostate$hx.ind <- outer(unclass(prostate$hx),na.omit(unique(unclass(prostate$hx))),"==")*1
prostate$ekg.ind <- outer(unclass(prostate$ekg),na.omit(unique(unclass(prostate$ekg))),"==")*1
prostate$bm.ind <- outer(unclass(prostate$bm),na.omit(unique(unclass(prostate$bm))),"==")*1

dmulti <-
function(x,mu,sig,...)  {apply(cbind(dmnom(x[,1:14],mu[,1:14]),
                                     dnorm(x[,-(1:14)],mu[,15:22],sig[,15:22])),1,prod,na.rm=TRUE)}

multimu <- list(
                G1 = function (p) 
            {
                rbind(c(inv.glogit(p[1:3]), inv.glogit(p[4]), 
                      inv.glogit(p[5:9]), inv.glogit(p[10]), 
                      p[11:18]))
            },
                G2 = function (p) 
            {
                rbind(c(inv.glogit(p[19:21]), inv.glogit(p[22]), 
                      inv.glogit(p[23:27]), inv.glogit(p[28]), 
                      p[29:36]))
            })

multisig <- list(
                 G1 = function (p) 
             {
                 rbind(c(rep(1, 14), exp(p[1:8])))
             },
                 G2 = function (p) 
             {
                 rbind(c(rep(1, 14), exp(p[9:16])))
             })

prost.mustart <-  c(-2.7,-4.49,-4.5,-0.14,-1.61, -1.8,0.14,-0.57,
                    -1.35,-5.66,71.79,100.23,14.42,8.22,137.22,9.05,3.13,1.63,
                    -2.267,-2.9,-4.51,-0.53,-1.98,-1.98,0.074,-0.84,-0.81,-0.49,
                    70.96,97.34,14.25,8.06,130.75,11.98,3.96,4.01)

prost.sigstart <- c( 1.872, 2.574, 0.950, 0.431, 2.916, 0.294, 0.322, -0.607 ,
                     2.021, 2.608, 0.801, 0.299, 2.985, 0.404, 0.523,  0.475 )


moc.prost <- moc(cbind(prostate$pf.ind, prostate$hx.ind,prostate$ekg.ind,prostate$bm.ind,
                 as.matrix(prostate[,c("age", "wti", "sbp", "dbp", "hg", "sg","sqrt.sz","log.ap")])),
                 density = dmulti, groups = 2, joint = TRUE, 
                 gmu = multimu, gshape = multisig, pgmu = prost.mustart, 
                 pgshape = prost.sigstart, pgmix = 0,gradtol=1e-3)

save(moc.prost,file="prost_moc1.RData",compress=TRUE)  #it takes a long time to estimate
                                                       #so you probably want to save it.

## As we can see only the three last continuous variables and the proportions of bm.ind
## substantially differs between those two groups.
plot(moc.prost,scale=TRUE)

## The most important remaining correlation is between SBP and DBP:
tmp <- na.omit(prostate[,c("age", "wti", "sbp", "dbp", "hg", "sg","sqrt.sz","log.ap")])
cov.wt(tmp,wt=post(moc.prost)[-na.action(tmp),1],cor=TRUE)$cor
cov.wt(tmp,wt=post(moc.prost)[-na.action(tmp),2],cor=TRUE)$cor

## We can thus incorporate this in a bivariate normal for those two variables.

dmulti.mv <-
function(x,mu,sig,extra)  {apply(cbind(dmnom(x[,1:14],mu[,1:14]),
                                 dnorm(x[,-c(1:14,17:18)],mu[,c(15:16,19:22)],sig[,c(15:16,19:22)]),
                                 mvnorm(x[,17:18],mu[,17:18],sig[,17:18],extra)),1,prod,na.rm=TRUE)}
                                                                        #correlation between sbp dbp

multiextra <-  list(G1 = function (p) { rbind(p[1]) },
                    G2 = function (p) { rbind(p[2]) })

moc.prost.2 <- moc(cbind(prostate$pf.ind, prostate$hx.ind,prostate$ekg.ind,prostate$bm.ind,
                 as.matrix(prostate[,c("age", "wti", "sbp", "dbp", "hg", "sg","sqrt.sz","log.ap")])),
                 density = dmulti.mv, groups = 2, joint = TRUE, 
                 gmu = multimu, gshape = multisig, gextra=multiextra, check.length=FALSE,
                 pgmu = moc.prost$coef[1:36], pgshape = moc.prost$coef[37:52],
                 pgextra = c(-0.2,-0.2) ,pgmix = moc.prost$coef[53],gradtol=1e-3)

save(moc.prost.2,file="prost_moc2.RData",compress=TRUE) 

## The corresponding correlations in each mixture are
cov2cor(solve(t(matrix(c(1,0,moc.prost.2$coef[53],1),2,2))%*% (matrix(c(1,0,moc.prost.2$coef[53],1),2,2))))
cov2cor(solve(t(matrix(c(1,0,moc.prost.2$coef[54],1),2,2))%*% (matrix(c(1,0,moc.prost.2$coef[54],1),2,2))))
## taking account of those correlations gives a better likelihood

AIC(moc.prost,moc.prost.2,k="BIC")

## but does not appear to have much effect (or benefits) on other aspects
## of the fit to the data. The resulting classification are also very
## similar.

source(system.file("Utils","combine.moc.R",package="moc"))  #load some utility functions
combine.prob(moc.prost,moc.prost.2)

## As a last example, we add a correlation between wti and hg and a dependence on bm.ind,
## since bm.ind are indicators of a categorical variable the dependence of wti and hg on bm
## is modeled by a location factor on the continuous variables (mean differences)

dmulti.mv.2 <-
function(x,mu,sig,extra)  {apply(cbind(dmnom(x[,1:14],mu[,1:14]),
                           dnorm(x[,-c(1:14,16:19)],mu[,c(15,20:22)],sig[,c(15,20:22)]),
                           mvnorm(x[,17:18],mu[,17:18],sig[,17:18],
                                  cbind(extra[,1])), #correlation between sbp dbp
                           mvnorm(x[,c(16,19)],apply(mu[,c(16,19)],2,"+",
                           cbind(extra[,3])*x[,14]), #corr. wti and hg and location with bm
                                  sig[,c(16,19)],cbind(extra[,2]))),
                           1,prod,na.rm=TRUE)}
multiextra.2 <-  list(G1 = function (p) { rbind(p[1:3]) },
                      G2 = function (p) { rbind(p[4:6]) })

## thus p[1] and p[4]  holds the parameters for the correlation between sbp dbp in group 1 and 2,
## p[2], p[5] are for the correlation between wti and hg in their respective group,
## while p[3] and p[4] stands for the mean differences of wti and hg between the two
## categories of bm in each mixture group.

## Be aware that this last model can take a while to estimate!

moc.prost.2.3 <- moc(cbind(prostate$pf.ind, prostate$hx.ind,prostate$ekg.ind,prostate$bm.ind,
                 as.matrix(prostate[,c("age", "wti", "sbp", "dbp", "hg", "sg","sqrt.sz","log.ap")])),
                 density = dmulti.mv.2, groups = 2, joint = TRUE, 
                 gmu = multimu, gshape = multisig, gextra=multiextra.2, check.length=FALSE,
                 pgmu = moc.prost.2$coef[1:36], pgshape = moc.prost.2$coef[37:52],
                 pgextra = c(moc.prost.2$coef[53],0,0,moc.prost.2$coef[54],0,0) ,pgmix = moc.prost.2$coef[55],gradtol=1e-3)

save(moc.prost.2.3,file="prost_moc23.RData",compress=TRUE) #keep the results

moc.prost.2.3
