## This example concerns 121 patients having treatment for breast cancer.
## The original data are from Boag (1949) but are directly available from
## Example 9.2 of McLachlan & Peel pp. 280-282
## (or http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html)
## This example illustrates the use of mixture to estimate cure rate
## (survival) in the presence of competing risk.

require(moc)

#bcs <- read.table("DataBank/breast_cancer.dat",header=TRUE)#replace to read data on your computer
bcs <-
structure(list(
delta = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
time = c(0.3, 12.2, 17.5, 28.2, 41, 78, 5, 12.3, 17.9, 29.1, 
42, 80, 5.6, 13.5, 19.8, 30, 44, 84, 6.2, 14.4, 20.4, 31, 46, 
87, 6.3, 14.4, 20.9, 31, 48, 89, 6.6, 14.8, 21, 32, 48, 90, 6.8, 
15.7, 21, 35, 51, 97, 7.5, 16.2, 21.1, 35, 51, 98, 8.4, 16.3, 
23, 38, 52, 100, 8.4, 16.5, 23.6, 39, 54, 114, 10.3, 16.8, 24, 
40, 56, 126, 11, 17.2, 24, 40, 60, 131, 11.8, 17.3, 27.9, 41, 
78, 174, 0.3, 110, 4, 111, 7.4, 112, 15.5, 112, 23.4, 162, 46, 
46, 51, 65, 68, 83, 88, 96, 111, 136, 112, 141, 113, 143, 114, 
167, 114, 177, 117, 179, 121, 189, 123, 201, 129, 203, 131, 203, 
133, 213, 134, 228, 134)), 
.Names = c("delta", "time"), class = "data.frame", row.names = c(NA, -121L))

## We first try to estimate the model of McLachlan & Peel using formula (10.9) p. 271. 
## In fact we use the following identity for the log-likelihood in formula (10.9)
#  I[1](d) log(pi_1*f1) + I[2](d) log(pi_2*f2) + I[0](d) log(S)
#  = log ( (pi_1*f1)^I[1](d) ) + log ( (pi_2*f2)^I[2](d) ) + log(S^I[0](d))
#  = log( (pi_1*f1)^I[1](d)*(pi_2*f2)^I[2](d)*S^I[0](d) )
#  = log( pi_1*f1*I[1](d) + pi_2*f1*I[2](d) + S*I[0](d) )
#  Where I[j](d) represents the indicator function such that
#  I[j](d) = 1 when d = j and is 0 otherwise. In this example
#  d corresponds to bcs[,1] which equals 0 for right-censored data (survivor),
#  2 for death due to breast cancer and 1 for death of other causes.
## This way, it can now be easily formulated as a standard mixture
## model in MOC.

## First we define the density, cumulative distribution and quantile
## functions for the Gompertz distribution.

pgompertz <- function(x,la,ka) {1-exp(-la*(exp(ka*x)-1)/ka)}
dgompertz <- function(x,la,ka) {(exp(-la * (exp(ka * x) - 1)/ka) * la * exp(ka * x) )}
qgompertz <-  function(prob,la,ka) {log(1-ka*log(1-prob)/la)/ka}

## There is no closed form for the expected value of the Gompertz distribution
## so we must use numerical integration over an infinite domain.
## If you use this, you should verify that there is no numerical instabilities or convergence
## problems afterward. It is probably safer to avoid this function.

gompertz.mean <- function(la,ka) integrate(function(x) x*dgompertz(x,la,ka),0,Inf)$value

## Note that both parameters of the Gompertz should be greater than 0 (la > 0 and ka > 0).
## McLachlan & Peel allowed negative values for ka but then it is not a density
## anymore. In the following example we allowed negative values for ka just to
## show how to obtain the same results as McLachlan & Peel,
## but negative values should not be allowed.

## Then we define the Gompertz survival function which is coded in its general
## form such that it would accept uncensored and censored (right or left) data.
## Here the second column of x is used to indicate the kind of censoring
## while the parameter cens will be used as an indicator for the cause
## of death.

Gompertz.Surv <- function(x,la,ka,cens)
  {
          (x[,2]==0)*(1-pgompertz(x[,1],la[,1],ka[,1]))*cens[,2]+
          (x[,2]<0)*pgompertz(x[,1],la[,1],ka[,1])*cens[,2]+
          (x[,2]>0)*dgompertz(x[,1],la[,1],ka[,1])*cens[,2]
  }


la.surv <- list(
                G1=function(p) {cbind(exp(p[1]),c(1,1,0)[bcs[,1]+1])},
                G2=function(p) {cbind(exp(p[2]),c(1,0,1)[bcs[,1]+1])}
                )
attr(la.surv,"parameters") <- c("log.lambda.1","log.lambda.2")

## Note also that we model the exponential distribution as the limit of a
## Gompertz as ka tends to 0 which we approximate by ka = .Machine$double.eps
## Thus the survival for the exponential distribution is
#  S2(t) = exp( -exp(la_2)*t ) instead of exp( -la_2*t) as in the
## usual parametrization. But this is a simple change of parameter that
## makes the parameters of the exponential and Gompertz mode comparable.

ka.surv <- list(
                G1=function(p) {cbind(p[1],0)},
                G2=function(p) {cbind(.Machine$double.eps,0)}
                )
attr(ka.surv,"parameters") <- c("kappa.1")

cens.surv <- list(
                  G1=function(p) {cbind(0,c(1,1,0)[bcs[,1]+1])},
                  G2=function(p) {cbind(0,c(1,0,1)[bcs[,1]+1])}
                  )


## The expected value for the Gompertz distribution is hard to compute since its mean value
## does not have a closed form. So we use its median instead which is easier to compute.
## Remember that the expected value does not affect the computation of the model itself,
## but only some specific residuals and the posterior expected mean which is hard to compare 
## to the mean observed value in the presence of censoring. (We could also use the 60th
## percentile which is close to mean value in the case of the Gompertz.)

gompertz.expected <- list(
                          G1=function(p) {cbind(qgompertz(0.5,exp(p[1]),p[3]),inv.glogit(p[4])[1])},
                          G2=function(p) {cbind(qgompertz(0.5,exp(p[2]),sqrt(.Machine$double.eps)),
                             inv.glogit(p[4])[2])}
                          )

## We could also use the 60th percentile which is close to mean value in the case of the Gompertz.

bcs.moc <- moc(bcs[,2:1],density=Gompertz.Surv,joint=TRUE,groups=2,
               gmu=la.surv,gshape=ka.surv,gextra=cens.surv,expected=gompertz.expected,
               pgmu=c(-5.4,-3),pgshape=c(-0.001),pgmix=0.7,gradtol=1e-6)


bcs.moc

## We use confint to obtain the parameters lamda_2 and the proportion pi_1
## instead of the log-odds, as reported by McLachlan & Peel. 
## Also note that McLachlan & Peel report the variances instead of the standard errors.

confint(bcs.moc,parm=list(~exp(p2),~1/(1+exp(p4)))) 

## We can obtain a graph of the survival function (in fact the cumulative death function)
## easily:

try(density.moc(bcs.moc,var=1:2,along=bcs[,1],plot="pq-plot"))

## The third graph correspond to 1-S_2(t) which is the conditional death function
## from breast cancer. Note that the first graph is for censored data so the empirical
## distribution stands over the estimated values since the observed times are lower
## bounds.

## As noted earlier the solution that is found has a ka < 0, then the corresponding
## Gompertz model is not a density and the survival function has
## a lower asymptote when t goes to infinity which is exp(la/ka).
## Thus the solution found is not a proper one, but ka is near zero
## and if we restrict ka to be > 0 we find that ka ~= 0 which give an exponential
## distribution.

## The preceding model assumes that the mixture corresponds exactly to the
## measured risk variable. That is, the form of the survival is exactly exponential
## when death occurs due to breast cancer, it is exactly Gompertz when the
## death is from another cause and it is a censored mixture of those when the death
## did not occur at the end of the study. This is a strong assumption that should
## not me made a priori, without further investigation, in most applications.
## In fact even in cases where we investigate a single cause of death, the
## mixture model still provide a powerful way of doing a semi-parametric
## survival analysis that can compete with the nonparametric estimator.
## Thus, we proceed by relaxing this assumption and show how
## to construct a model where death due to each cause follow a
## different mixture model. Note that in most applications of this
## kind of model the number of mixtures does not have to be the
## same as the number of considered causes of death since the
## mixtures are there to fit the survival functions and the
## different causes mix these components in different proportions.

## We now consider a mixture of two Gompertz for both causes (remember that
## when the death is not due to breast cancer the variable only indicate that
## the cause of death is somewhat else).

## Since we don't want to allow improper solutions we restrict ka > 0.

ka.surv2 <- list(
                 G1=function(p) {cbind(exp(p[1])+2*.Machine$double.eps,0)},
                 G2=function(p) {cbind(exp(p[2])+2*.Machine$double.eps,0)}
                 )
attr(ka.surv2,"parameters") <- c("log.kappa.1","log.kappa.2")

## Contrary to the preceding example we mix the Gompertz in different unknown
## proportions (to be estimated) in both groups.

cens.surv2 <- list(
                   G1=function(p) {p1 <- inv.glogit(p[1]); cbind(0,c(1,p1)[bcs[,1]+1])},
                   G2=function(p) {p2 <- inv.glogit(p[2]); cbind(0,c(1,p2)[bcs[,1]+1])}
                   )
attr(cens.surv2,"parameters") <- c("logit.1","logit.2")

gompertz.expected2 <- list(
                           G1=function(p) {cbind(qgompertz(0.5,exp(p[1]),exp(p[3])),
                             inv.glogit(p[5])%*%c(1,2))},
                           G2=function(p) {cbind(qgompertz(0.5,exp(p[2]),exp(p[4])),
                             inv.glogit(p[6])%*%c(1,2))}
                           )



bcs.free2 <- 
moc(bcs[, 2:1], density = Gompertz.Surv, joint = TRUE, groups = 2, 
    gmu = la.surv, gshape = ka.surv2, gextra = cens.surv2,expected=gompertz.expected2,
    pgmu = c(-5, -2), pgshape = c(-5, -3), pgextra=c(0.5,-0.5),pgmix = 0, 
    gradtol = 1e-06, iterlim = 500)


bcs.free2

## As can be seen from the shape parameter, the first group has a
## ka = exp(log.kappa.1) =  5.27e-08  which is in fact an exponential
## model. It is also seen from the extra parameters that the second
## Gompertz applies solely to death due to cancer since
## inv.glogit(logit.2 = 20.99) = c(0,1).
## Again we can obtain the graphs of 1-S1(t) and 1-S2(t) easily

try(density.moc(bcs.free2,var=1:2,along=bcs[,1],plot="pq-plot"))

## The Gompertz distribution is highly sensitive in its parameters, it
## change rapidly of shape. This make this distribution really hard to fit
## with many local maximums. Moreover the moments are hard to compute
## since there is no closed form expression for them. A more tractable
## distribution showing a broad range of shapes and forms widely used
## in survival analysis is the Weibull distribution (which also include
## the exponential as a special case). We now fit the same
## model as the previous one but with a Weibull distribution.

Weibull.Surv <- function(x,la,ka,cens)
  {
      (x[,2]==0)*(1-pweibull(x[,1],shape=la[,1],scale=ka[,1]))*cens[,2]+
          (x[,2]<0)*pweibull(x[,1],shape=la[,1],scale=ka[,1])*cens[,2]+
              (x[,2]>0)*dweibull(x[,1],shape=la[,1],scale=ka[,1])*cens[,2]
  }

weibull.mean <- function(a,b) {b*gamma(1+1/a)}

ka.wb <- list(G1=function(p) {cbind(exp(p[1]),0)},
              G2=function(p) {cbind(exp(p[2]),0)})
attr(ka.wb,"parameters") <- c("log.kappa.1","log.kappa.2")

expected.wb.bcs <- list(
                        G1=function(p) {cbind(weibull.mean(la.surv[[1]](p[1:2])[,1],ka.wb[[1]](p[3:4])[,1]),
                        inv.glogit(p[5])%*%c(1,2))},
                        G2=function(p) {cbind(weibull.mean(la.surv[[2]](p[1:2])[,1],ka.wb[[2]](p[3:4])[,1]),
                        inv.glogit(p[6])%*%c(1,2))}
                        )

bcs.moc.wb <- moc(bcs[,2:1],density=Weibull.Surv,joint=TRUE,groups=2,
                  gmu=la.surv, gshape=ka.wb, gextra=cens.surv2, expected=expected.wb.bcs,
                  pgmu=c(-1,1), pgshape=c(6,8), pgextra=c(0,0), pgmix=0,
                  gradtol=1e-6, iterlim=200)


bcs.moc.wb

try(density.moc(bcs.moc.wb,var=1:2,along=bcs[,1],plot="pq-plot"))


