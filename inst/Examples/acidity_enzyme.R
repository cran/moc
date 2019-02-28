## The following examples use the acidity and enzyme
## data available at
##  http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html
## They were used to illustrate bootstrapping of mixture models in
## McLachlan & Peel (2000) pp. 194-196.


## First we load the required libraries
## (we also use the boot library for bootstrapping)

require(moc)
require(boot)

## Then, we set the required density, gmu and gshape function 
## for a mixture of normals within MOC.

normal <-
function(x,mu,sig,...) {dnorm(x,mu,sig)}

norm.mu <- list(
G1 = function (p) {rbind(p[1])}, 
G2 = function (p) {rbind(p[2])},
G3 = function (p) {rbind(p[3])},
G4 = function (p) {rbind(p[4])}              
) #up to four groups.

norm.sig <- list(
G1 = function (p) {rbind(exp(p[1]))}, 
G2 = function (p) {rbind(exp(p[2]))},
G3 = function (p) {rbind(exp(p[3]))},
G4 = function (p) {rbind(exp(p[4]))}
)

## We load the first example data set:
## Acidity index in a sample of 155 lakes in north-central Wisconsin.
# (These data are available at  http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html)

	# If you have an internet connection
        # acidity <- scan(file="http://www.stats.bris.ac.uk/~peter/mixdata",skip=50,nmax=155)

acidity.file <- system.file("Examples","acidity.RData",package="moc",mustWork=TRUE)
load(acidity.file)


acd1 <-  moc(acidity,density=normal,groups=1,gmu=norm.mu[1],gshape=norm.sig[1],pgmu=c(5),pgshape=c(0.04))
acd2 <-  moc(acidity,density=normal,groups=2,gmu=norm.mu[1:2],gshape=norm.sig[1:2],pgmu=c(4,6),
             pgshape=c(-0.8,-0.8),pgmix=-0.5)
acd3 <-  moc(acidity,density=normal,groups=3,gmu=norm.mu[1:3],gshape=norm.sig[1:3],pgmu=c(4,6,8),
             pgshape=c(-1,-1,-1),pgmix=c(0.5,0.5))
acd4 <-  moc(acidity,density=normal,groups=4,gmu=norm.mu,gshape=norm.sig,pgmu=c(3,4.2,6,7),
             pgshape=rep(-1,4),pgmix=c(0,0,0))

AIC(acd1,acd2,acd3,acd4,k="BIC")
entropy(acd2,acd3,acd4)

## The ICL-BIC criterion clearly identifies the mixture of 2 normals as the best model.
## We can see that -2*logLik decreases drastically when going from a single groups to 2 groups.
## Adding more groups than only slightly decreases the likelihood at the price of a big
## increase of the entropy (lower separation between groups).
## As can also be seen from the following graphs, the separation between the groups
## is much better for the 2 groups model.

profilesplot(acd2,type="posterior")
profilesplot(acd3,type="posterior")

## Inspection of the residuals does not reveal any miss-fitting except for the lowest point
## which is a candidate outlier.

plot(residuals(acd2,within=TRUE))
plot(residuals(acd3,within=TRUE))

## The fitted density of the 2 normal mixture and the histogram closely agree.
## Both show the two principal modes in the data.

density.moc(acd2,var=1,plot="density",ylim=c(0,1))
hist(acidity,breaks="Sturges",prob=TRUE,add=TRUE)

## The following graph closely reproduce the density plot found in McLachlan & Peel p.195
## Strangely they present the density of the 3 groups model while arguing that the 2 groups
## model is preferable. This is probably due to the fact that they used a histogram with a
## lot of breakpoints that is more closely matched by the density estimate of the 3 groups
## mixture. However, you should be careful when using a histogram with too many breakpoints,
## since most cells then contains very few data points which can make the histogram highly
## unstable and unreliable.

density.moc(acd3,var=1,plot="density",ylim=c(0,1))
hist(acidity,breaks=30,prob=TRUE,add=TRUE)


## Though it must now be clear that a mixture of 2 normals best describe
## this data set, there are no formal test of this since the likelihood
## ratio test statistic (LRTS) of two models with different number of groups
## follows an unknown law. So we proceed by bootstrapping the LRTS as described
## in McLachlan & Peel pp. 192-201. However, there are no clear advantages in
## using the bootstrap to assess p-values for mixture models over good use
## of information criterion and other tools available in MOC to asses the fit of the model.
## In fact bootstrapping mixture models can be misleading since it completely ignores
## the dependence on starting values, the existence of local maxima's and other properties
## that would normally be assessed with a single data set.

## First we construct a function that will fit 2 MOC models: one with g groups, the other
## with g+1 groups and then return the LRTS.

boot.normloglike.moc <- function(data,gr,st1,st2) 
{ 
assign("boot.tmp",data,envir=.GlobalEnv)
 l1 <-  try(moc(boot.tmp,density=normal,groups=gr,gmu=norm.mu[1:gr],gshape=norm.sig[1:gr],
                pgmu=st1$mu,pgshape=st1$shape,pgmix=st1$mix,print.level=0)$loglike)
 l2 <- try(moc(boot.tmp,density=normal,groups=gr+1,gmu=norm.mu[1:(gr+1)],gshape=norm.sig[1:(gr+1)],
               pgmu=st2$mu,pgshape=st2$shape,pgmix=st2$mix,print.level=0)$loglike) 
rm("boot.tmp",envir=.GlobalEnv) 
diff(c(l1,l2)*2)}

## We also need a function that simulate a data set following a given mixture model.

gen.moc.norm <- function(d,mle)
{
    gk <- apply(rmultinom(mle$n,1,mle$pmix)==1,2,which)
    rnorm(mle$n,mle$mu[gk],mle$sig[gk])
}

if(interactive()) {
if(menu(c("YES","NO"),title="Proceed with time consuming bootstrap?")==1){

## We set the needed values to test H0: g=1 vs H1: g=2

mle1 <- list(n=155,mu=acd1$coef[1],sig=exp(acd1$coef[2]),pmix=1)
start1 <- list(mu=acd1$coef[1],shape=acd1$coef[2],mix=0)
start2 <- list(mu=acd2$coef[1:2],shape=acd2$coef[3:4],mix=acd2$coef[5])

acidity.boot1 <- boot(acidity, sim="parametric",boot.normloglike.moc,R=100,ran.gen=gen.moc.norm,mle=mle1,gr=1,st1=start1,st2=start2)

## We set the needed values to test H0: g=2 vs H1: g=3

mle2 <- list(n=155,mu=acd2$coef[1:2],sig=exp(acd2$coef[3:4]),pmix=inv.glogit(acd2$coef[5]))
start21 <- list(mu=acd2$coef[1:2],shape=acd2$coef[3:4],mix=acd2$coef[5])
start22 <- list(mu=acd3$coef[1:3],shape=acd3$coef[4:6],mix=acd3$coef[7:8])


acidity.boot2 <- boot(acidity, sim="parametric",boot.normloglike.moc,
                      R=100,ran.gen=gen.moc.norm,mle=mle2,gr=2,st1=start21,st2=start22)

mle3 <- list(n=155,mu=acd3$coef[1:3],sig=exp(acd3$coef[4:6]),pmix=inv.glogit(acd3$coef[7:8]))
start31 <- list(mu=acd3$coef[1:3],shape=acd3$coef[4:6],mix=acd3$coef[7:8])
start32 <- list(mu=acd4$coef[1:4],shape=acd4$coef[5:8],mix=acd4$coef[9:11])

## We set the needed values to test H0: g=3 vs H1: g=4

acidity.boot3 <- boot(acidity, sim="parametric",boot.normloglike.moc,
                      R=100,ran.gen=gen.moc.norm,mle=mle3,gr=3,st1=start31,st2=start32)
}}

## The next data set used by McLachlan & Peel to illustrate the bootstrap concerns
## Enzymatic activity (related to the metabolism of carcinogenic substance) 
## in the blood of 245 unrelated individuals.
## It is available at  http://www.maths.uq.edu.au/~gjm


	# If you have an internet connection
        # enzyme <- scan(file="http://www.stats.bris.ac.uk/~peter/mixdata",skip=18,nmax=245)

enzyme.file <- system.file("Examples","enzyme.RData",package="moc",mustWork=TRUE)
load(enzyme.file)

enzy1 <-  moc(enzyme,density=normal,groups=1,gmu=norm.mu[1],gshape=norm.sig[1],pgmu=c(5),pgshape=c(0.04))
enzy2 <-  moc(enzyme,density=normal,groups=2,gmu=norm.mu[1:2],gshape=norm.sig[1:2],
              pgmu=c(0.2,1.2),pgshape=c(0.04,0.04),pgmix=0)
enzy3 <-  moc(enzyme,density=normal,groups=3,gmu=norm.mu[1:3],gshape=norm.sig[1:3],
              pgmu=c(0.2,1.2,2),pgshape=c(0.04,0.04,0.04),pgmix=c(0,0))
enzy4 <-  moc(enzyme,density=normal,groups=4,gmu=norm.mu[1:4],gshape=norm.sig[1:4],
              pgmu=c(0.2,1.2,2,2.6),pgshape=c(0.04,0.04,0.04,0.04),pgmix=c(0,0,0))

## The negative log-likelihood decrease drastically when going from 1 group to 2 groups,
## while adding further groups shows less important decrease. The entropy of the 2 groups
## model is also smaller and its reduction in entropy greater than for the other models.
## The preferred model based on the ICL-BIC would be the 2 groups model.

AIC(enzy1,enzy2,enzy3,enzy4,k="BIC")
entropy(enzy2,enzy3,enzy4)

## The separation between group 2 and 3 in the 3 groups model is poor.

profilesplot(enzy2,type="posterior")
profilesplot(enzy3,type="posterior")

## Some residuals for the 2 groups model are still high particularly in the second group,
## but the situation is even worse for the 3 groups model.

plot(residuals(enzy2,within=TRUE))
plot(residuals(enzy3,within=TRUE))

## This is probably due to the fact that the normal distribution does not capture completely the
## asymmetry and longer tail distribution and would require more groups to achieve this.

density.moc(enzy2,var=1,plot="density")
hist(enzyme,breaks=25,prob=TRUE,add=TRUE)

density.moc(enzy2,var=1,plot="pq-plot")
density.moc(enzy3,var=1,plot="pq-plot")

## We would probably do better with a skewed and long tail distribution like the Gamma,
## which often do pretty well in modeling cells concentration.

GAMMA <- function(x,a,s,...) {dgamma(x,a,scale=s)}

gamma.shape <-  list(
G1 = function (p) {rbind(exp(p[1]))}, 
G2 = function (p) {rbind(exp(p[2]))},
G3 = function (p) {rbind(exp(p[3]))},
G4 = function (p) {rbind(exp(p[4]))}
)

gamma.scale <- list(
G1 = function (p) {rbind(exp(p[1]))}, 
G2 = function (p) {rbind(exp(p[2]))},
G3 = function (p) {rbind(exp(p[3]))},
G4 = function (p) {rbind(exp(p[4]))}
)

## The expected value of the Gamma distribution here is shape*scale and we have to compute
## this explicitly to compare directly the observed and predicted mean values.

expected.gamma1 <- list(G1=function(p) gamma.shape[[1]](p[1])*gamma.scale[[1]](p[2]))
enzy.gam1 <- moc(enzyme, density = GAMMA, groups = 1, gmu = gamma.shape[1], gshape = gamma.scale[1],
                 expected=expected.gamma1,pgmu = 0, pgshape = 0)

expected.gamma2 <- list(G1=function(p) gamma.shape[[1]](p[1])*gamma.scale[[1]](p[3]),
                        G2=function(p) gamma.shape[[1]](p[2])*gamma.scale[[1]](p[4]))
enzy.gam2 <- moc(enzyme, density = GAMMA, groups = 2, gmu = gamma.shape[1:2], gshape = gamma.scale[1:2],
                 expected=expected.gamma2,pgmu = log(c(0.2, 3)) , pgshape = c(0,0), pgmix = 0)

expected.gamma3 <- list(G1=function(p) gamma.shape[[1]](p[1])*gamma.scale[[1]](p[4]),
                        G2=function(p) gamma.shape[[1]](p[2])*gamma.scale[[1]](p[5]),
                        G3=function(p) gamma.shape[[1]](p[3])*gamma.scale[[1]](p[6]))
enzy.gam3 <- moc(enzyme, density = GAMMA, groups = 3, gmu = gamma.shape[1:3], gshape = gamma.scale[1:3],
                 expected=expected.gamma3,pgmu = log(c(0.2, 2, 4)) , pgshape = c(0,0,0), pgmix = c(0,0))

expected.gamma4 <- list(G1=function(p) gamma.shape[[1]](p[1])*gamma.scale[[1]](p[5]),
                        G2=function(p) gamma.shape[[1]](p[2])*gamma.scale[[1]](p[6]),
                        G3=function(p) gamma.shape[[1]](p[3])*gamma.scale[[1]](p[7]),
                        G4=function(p) gamma.shape[[1]](p[4])*gamma.scale[[1]](p[8]))
enzy.gam4 <- moc(enzyme, density = GAMMA, groups = 4, gmu = gamma.shape[1:4], gshape = gamma.scale[1:4],
                 expected=expected.gamma4,pgmu = log(c(0.2, 1, 2, 4)) , pgshape = c(0,0,0,0), pgmix = c(0,0,0))


## In this case the ICL-BIC  and entropy clearly favor the 2 groups model
## which reveals a pretty clear group separation and small residuals.
## Since the gmu function does not correspond to the expected values,
## the default definition of deviance residuals is not well defined in this case.
## We thus look at response and gradient residuals

AIC(enzy.gam1,enzy.gam2,enzy.gam3,enzy.gam4,k="BIC")
entropy(enzy.gam2,enzy.gam3,enzy.gam4)

profilesplot(enzy.gam2,type="posterior")
profilesplot(enzy.gam3,type="posterior")
plot(residuals(enzy.gam2,type="response",within=TRUE))
plot(residuals(enzy.gam2,type="gradient",within=TRUE))


## This model also do better at capturing the longer tail distribution for 
## the second group but does not fully account for the skewness.

density.moc(enzy.gam2,var=1,plot="density")
hist(enzyme,breaks=25,prob=TRUE,add=TRUE)
density.moc(enzy.gam2,var=1,plot="pq-plot")

## An alternative approach often used by most practitioner is to use transformations to
## "normalize" the distribution of the data (like log or sqrt).

logenzy1 <-  moc(log(enzyme),density=normal,groups=1,gmu=norm.mu[1],gshape=norm.sig[1],
                 pgmu=c(0),pgshape=c(0.04))
logenzy2 <-  moc(log(enzyme),density=normal,groups=2,gmu=norm.mu[1:2],gshape=norm.sig[1:2],
                 pgmu=c(-2,0.2),pgshape=c(-0.04,-0.04),pgmix=0)
logenzy3 <-  moc(log(enzyme),density=normal,groups=3,gmu=norm.mu[1:3],gshape=norm.sig[1:3],
                 pgmu=c(-3.2,-1.5,0.3),pgshape=c(-0.04,-0.04,-0.04),pgmix=c(0,0))
logenzy4 <-  moc(log(enzyme),density=normal,groups=4,gmu=norm.mu[1:4],gshape=norm.sig[1:4],
                 pgmu=c(-3.2,-2,-1,0.5),pgshape=rep(-0.08,4),pgmix=c(0,0,0))

sqrtenzy1 <-  moc(sqrt(enzyme),density=normal,groups=1,gmu=norm.mu[1],gshape=norm.sig[1],
                  pgmu=c(0.8),pgshape=c(0.04))
sqrtenzy2 <-  moc(sqrt(enzyme),density=normal,groups=2,gmu=norm.mu[1:2],gshape=norm.sig[1:2],
                  pgmu=c(0.4,1.1),pgshape=c(-0.04,-0.04),pgmix=0)
sqrtenzy3 <-  moc(sqrt(enzyme),density=normal,groups=3,gmu=norm.mu[1:3],gshape=norm.sig[1:3],
                  pgmu=c(0.4,1,1.5),pgshape=c(-0.04,-0.04,-0.04),pgmix=c(0,0))
sqrtenzy4 <-  moc(sqrt(enzyme),density=normal,groups=4,gmu=norm.mu[1:4],gshape=norm.sig[1:4],
                  pgmu=c(0.4,0.8,1.2,1.6),pgshape=rep(-0.08,4),pgmix=c(0,0,0))


## For the log-transformed data the ICL-BIC does not clearly favors the 3 groups
## model over the 2 groups model. Taking a closer look at the fitted densities
## and data reveals that in fact the log-transformation does a pretty good job
## at making two groups very close to normal distribution. But then there are
## three isolated observations (< -3) that becomes like outliers (or influential)
## which would require specific treatment or inclusion of a specific group for them.
## Note that this happens quite often with blood measurements when the values approach
## the detectable threshold (sensitivity). Further information about the data and
## measurement process would be needed here.

AIC(logenzy1,logenzy2,logenzy3,logenzy4,k="BIC")
density.moc(logenzy2,var=1,plot="density",ylim=c(0,0.8))
hist(log(enzyme),breaks=20,prob=TRUE,add=TRUE)
density.moc(logenzy3,var=1,plot="density",ylim=c(0,0.8))
hist(log(enzyme),breaks=20,prob=TRUE,add=TRUE)

## The square-root transform clearly prefers the 2 groups model over the 3 groups.
## However adding a fourth group does better than the 3 groups model which suggest
## that adding further groups could help us find a better model. But inspection
## of the fitted density indicates that this would only help capture the asymmetric and
## longer tail distribution of the second group that is not resolved by fitting
## normal distribution to the transformed data. So there are no real benefits in
## using transformed data over the use of a better fitting distribution 
## like the gamma on the original data.

AIC(sqrtenzy1,sqrtenzy2,sqrtenzy3,sqrtenzy4,k="BIC")

density.moc(sqrtenzy2,var=1,plot="density")
hist(sqrt(enzyme),breaks=20,prob=TRUE,add=TRUE)


## If you still think that bootstrapping the LRTS is useful, here it is.

if(interactive()) {
if(menu(c("YES","NO"),title="Proceed with time consuming bootstrap?")==1){

mle1 <- list(n=245,mu=enzy1$coef[1],sig=exp(enzy1$coef[2]),pmix=1)
start1 <- list(mu=enzy1$coef[1],shape=enzy1$coef[2],mix=0)
start2 <- list(mu=enzy2$coef[1:2],shape=enzy2$coef[3:4],mix=enzy2$coef[5])

enzyme.boot1 <- boot(enzyme, sim="parametric",boot.normloglike.moc,R=100,
                     ran.gen=gen.moc.norm,mle=mle1,gr=1,st1=start1,st2=start2)

mle2 <- list(n=245,mu=enzy2$coef[1:2],sig=exp(enzy2$coef[3:4]),pmix=inv.glogit(enzy2$coef[5]))
start21 <- list(mu=enzy2$coef[1:2],shape=enzy2$coef[3:4],mix=enzy2$coef[5])
start22 <- list(mu=enzy3$coef[1:3],shape=enzy3$coef[4:6],mix=enzy3$coef[7:8])


enzyme.boot2 <- boot(enzyme, sim="parametric",boot.normloglike.moc,R=100,
                     ran.gen=gen.moc.norm,mle=mle2,gr=2,st1=start21,st2=start22)

mle3 <- list(n=245,mu=enzy3$coef[1:3],sig=exp(enzy3$coef[4:6]),pmix=inv.glogit(enzy3$coef[7:8]))
start31 <- list(mu=enzy3$coef[1:3],shape=enzy3$coef[4:6],mix=enzy3$coef[7:8])
start32 <- list(mu=enzy4$coef[1:4],shape=enzy4$coef[5:8],mix=enzy4$coef[9:11])


enzyme.boot3 <- boot(enzyme, sim="parametric",boot.normloglike.moc,R=100,
                     ran.gen=gen.moc.norm,mle=mle3,gr=3,st1=start31,st2=start32)
}}
