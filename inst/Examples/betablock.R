## The following data concerns a 22-center clinical trial on beta-blockers
## effectiveness in reducing mortality after myocardial infarction.
## These data are available at  http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html
## and also in the examples of BUGS (http://www.mrc-bsu.cam.ac.uk/bugs/)
## and of Latent Gold (http://www.statisticalinnovations.com/products/latentgold_datasets.html).

## This data set is used to illustrate a kind of meta-analysis that is more
## convincing than meta-analysis based on literature which often just
## reflects publishing bias. In the present case, it is based on a proper trial
## with randomization and is a good example of the use of mixture as
## a non-parametric random effect logistic regression.

	#If you have an internet connection
        #bbloc <- read.table("http://www.maths.uq.edu.au/~gjm/DATA/Beta.dat",header=TRUE,skip=1)

require(moc)
beta.file <- system.file("Examples","bbloc.RData",package="moc",mustWork=TRUE)
load(beta.file)

## First, we simply reproduce the fixed treatment effect with a random intercept 
## as in McLachlan & Peel (2000) pp. 165-166

binom <- function(x,prob,sh=1,size){dbinom(x,size,prob)}

## Here p[1] stands for the fixed treatment and p[2:4] for the random effects that
## varies across group.

gmu3 <-list(
 Group1 = function (p) {cbind(1/(1 + exp(p[2])), 1/(1 + exp(p[2] + p[1])))},
 Group2 = function (p) {cbind(1/(1 + exp(p[3])), 1/(1 + exp(p[3] + p[1])))},
 Group3 = function (p) {cbind(1/(1 + exp(p[4])), 1/(1 + exp(p[4] + p[1])))})
attr(gmu3,"parameters")<-c("-log.OR.trt","-logodd.G1","-logodd.G2","-logodd.G3")

gext3 <- list(
 Group1 = function(p) cbind(bbloc[,3],bbloc[,5]),
 Group2 = function(p) cbind(bbloc[,3],bbloc[,5]),
 Group3 = function(p) cbind(bbloc[,3],bbloc[,5]))

bin.fitted <- list(
 Group1 = function(p) (gmu3[[1]](p)%x%rep(1,22))*gext3[[1]](p),
 Group2 = function(p) (gmu3[[2]](p)%x%rep(1,22))*gext3[[2]](p),
 Group3 = function(p) (gmu3[[3]](p)%x%rep(1,22))*gext3[[3]](p))
 

bbloc.moc3 <-  moc(bbloc[,c(2,4)],density=binom,group=3,gmu=gmu3,gextra=gext3,
                   expected=bin.fitted,pgmu=c(0.25,1.61,2.25,2.8),pgmix=c(0,0))

bbloc.moc2 <-  moc(bbloc[,c(2,4)],density=binom,group=2,gmu=gmu3[-3],gextra=gext3[-3],
                   expected=bin.fitted[-3],pgmu=c(0.25,1.61,2.25),pgmix=0)

## Since we have used gmu=gmu3[-3] the parameters attribute won't be retained
## in the print. Fortunately since the number of parameters is not the same.
## However we can add this attribute directly in the fitted moc model.
attr(bbloc.moc2$gmu,"parameters")<-c("-log.OR.trt","-logodd.G1","-logodd.G2")

bbloc.moc1 <-  moc(bbloc[,c(2,4)],density=binom,group=1,gmu=gmu3[-(2:3)],gextra=gext3[-(2:3)],
                   expected=bin.fitted[-(2:3)],pgmu=c(0.25,2.25))
attr(bbloc.moc1$gmu,"parameters")<-c("-log.OR.trt","-logodd.G1")

AIC(bbloc.moc1,bbloc.moc2,bbloc.moc3,k="BIC")

bbloc.moc3
plot(bbloc.moc3)

## Note that we report the log-likelihood instead of the deviance of the binomial.
## The non-parametric estimate of the random intercept is a three-point distribution
## located at:

bbloc.moc3$coef[2:4]

## with respective masses:

inv.glogit(bbloc.moc3$coef[5:6])

## This is a nearly symmetric distribution with mean and standard deviation:

m3 <- weighted.mean(bbloc.moc3$coef[2:4],w=inv.glogit(bbloc.moc3$coef[5:6]))
sd3 <- sqrt( weighted.mean((bbloc.moc3$coef[2:4]-m3)^2,w=inv.glogit(bbloc.moc3$coef[5:6])))

m3
sd3

## This could be compared to the usual normal random effect. To make such a model in MOC
## we need the quadrature points and weights for the normal which are available in the
## file system.file("Utils","rasch.R",package="moc"). A five point quadrature is usually enough
## to estimate the normal integral with sufficient precision.

source(system.file("Utils","rasch.R",package="moc"))
norm.quad <- gherm.quad(5)

## The first parameter is still for the treatment effect but the second and third are for the
## mean and log(standard.deviation) of the normal random effect.

gmu.nq5 <- lapply(norm.quad$nodes,function(nod) eval(parse(text=
           paste("function (p) {pq <- p[2]+",nod,
                 "*exp(p[3]); cbind(1/(1 + exp(pq[1])), 1/(1 + exp(pq[1] + p[1])))}"))))

## This  way of writing the list gmu.nq5 will generate the right list for any length of quadrature
## chosen above and is shorter then the explicit programming:
# gmu.nq5 <-list(
# Group1 = function (p) {pq <- p[2]+norm.quad$nodes*exp(p[3]);
#                        cbind(1/(1 + exp(pq[1])), 1/(1 + exp(pq[1] + p[1])))},
# Group2 = function (p) {pq <- p[2]+norm.quad$nodes*exp(p[3]);
#                        cbind(1/(1 + exp(pq[2])), 1/(1 + exp(pq[2] + p[1])))},
# Group3 = function (p) {pq <- p[2]+norm.quad$nodes*exp(p[3]);
#                        cbind(1/(1 + exp(pq[3])), 1/(1 + exp(pq[3] + p[1])))},
# Group4 = function (p) {pq <- p[2]+norm.quad$nodes*exp(p[3]);
#                        cbind(1/(1 + exp(pq[4])), 1/(1 + exp(pq[4] + p[1])))},
# Group5 = function (p) {pq <- p[2]+norm.quad$nodes*exp(p[3]);
#                        cbind(1/(1 + exp(pq[5])), 1/(1 + exp(pq[5] + p[1])))})

gext.nq5 <- list(
 Group1 = function(p) cbind(bbloc[,3],bbloc[,5]),
 Group2 = function(p) cbind(bbloc[,3],bbloc[,5]),
 Group3 = function(p) cbind(bbloc[,3],bbloc[,5]),
 Group4 = function(p) cbind(bbloc[,3],bbloc[,5]),
 Group5 = function(p) cbind(bbloc[,3],bbloc[,5]))

bin.fitted.nq5 <- list(
 Group1 = function(p) (gmu.nq5[[1]](p)%x%rep(1,22))*gext.nq5[[1]](p),
 Group2 = function(p) (gmu.nq5[[2]](p)%x%rep(1,22))*gext.nq5[[2]](p),
 Group3 = function(p) (gmu.nq5[[3]](p)%x%rep(1,22))*gext.nq5[[3]](p),
 Group4 = function(p) (gmu.nq5[[3]](p)%x%rep(1,22))*gext.nq5[[3]](p),
 Group5 = function(p) (gmu.nq5[[3]](p)%x%rep(1,22))*gext.nq5[[3]](p))

## The mixture function is then only the normal quadrature point mass.

norm.mix5 <- function(p) rbind(norm.quad$weights)

bbloc.nq5 <-  moc(bbloc[,c(2,4)],density=binom,group=5,gmu=gmu.nq5,gextra=gext.nq5,
gmixture=norm.mix5,expected=bin.fitted.nq5,pgmu=c(0.25,2.25,log(0.4)))

## As we can see: the likelihood and entropy of the 3 mixtures and of the normal random effect
## are quite close. The latter has a better ICL-BIC because it achieved this with less
## estimated parameters.

AIC(bbloc.moc3,bbloc.nq5,k="BIC")

## However the normal random effect tends to overestimate the number of deaths
## for the lowest groups.

plot(bbloc.nq5)

## We now proceed with the semi-parametric approach by adding a random slope
## as in McLachlan & Peel (2000) pp. 165-166 and Aitkin (1999)
## This means that we allow the treatment effect to vary between mixture groups,
## the treatment parameters are p[2], p[4] and p[6].

gmu3.rs <-list(
 Group1 = function (p) {cbind(1/(1 + exp(-p[1])), 1/(1 + exp(-p[1] - p[2])))},
 Group2 = function (p) {cbind(1/(1 + exp(-p[3])), 1/(1 + exp(-p[3] - p[4])))},
 Group3 = function (p) {cbind(1/(1 + exp(-p[5])), 1/(1 + exp(-p[5] - p[6])))})
attr(gmu3.rs,"parameters") <- c("log.odd.G1","log.OR.trt1","log.odd.G2",
                                "log.OR.trt2","log.odd.G3","log.OR.trt3")

bin.fitted.rs <- list(
 Group1 = function(p) (gmu3.rs[[1]](p)%x%rep(1,22))*gext3[[1]](p),
 Group2 = function(p) (gmu3.rs[[2]](p)%x%rep(1,22))*gext3[[2]](p),
 Group3 = function(p) (gmu3.rs[[3]](p)%x%rep(1,22))*gext3[[3]](p))

bbloc.moc3.rs <-  moc(bbloc[,c(2,4)],density=binom,group=3,gmu=gmu3.rs,gextra=gext3,
                      expected=bin.fitted.rs,pgmu=c(-1.58,-0.32,-2.25,-0.25,-2.9,-0.08),
                      pgmix=bbloc.moc3$coef[5:6])

bbloc.moc3.rs

## The gain in likelihood and in entropy is very small with respect
## to the number of parameters added.
## So the ICL-BIC still favors the simple random intercept model.

AIC(bbloc.moc3,bbloc.moc3.rs,k="BIC")

