## Survival mixture Stanford heart valve transplant example
## These data are available at  http://lib.stat.cmu.edu/datasets/
## and can also be found in the package survival as heart and stanford2
## in a diiferent form not suitable for the example here.
## We consider the data at http://lib.stat.cmu.edu/datasets/
## from Crowley and Hu, Covariance Analysis of Heart Transplant Survival
## Data, J. Amer. Stat. Assoc, vol 72 (1977), p 27-36.
## from which we have constructed the extended survival status
## using the variable dead and reject giving a status variable
## 0 for censored, 1 dead rejection, 2 dead other cause.
## The first analysis here mimics the results obtained in
## MacLachlan et al. http://www.jstatsoft.org/v02/i07/0010.pdf
## used to illustrate their program http://www.jstatsoft.org/v02/i07/code/mgomp.f


heart <- read.table("DataBank/heart_extended.dat",header=TRUE)  ## change to suit your data location and format
heart2 <- na.omit(heart)   ## not necessary in moc but that's what MacLachlan et al. did
heart2 <- transform(heart2,age.std=scale(age),mismat.std=scale(mismat))

library(moc)    ## load the essential

## See breast_cancer.R for more coding details.

pgompertz <- function(x,la,ka) {1-exp(-la*(exp(ka*x)-1)/ka)}
dgompertz <- function(x,la,ka) {(exp(-la * (exp(ka * x) - 1)/ka) * la * exp(ka * x) )}
qgompertz <-  function(prob,la,ka) {log(1-ka*log(1-prob)/la)/ka}

gompertz.mean <- function(la,ka) integrate(function(x) x*dgompertz(x,la,ka),0,Inf)$value

Gompertz.Surv <- function(x,la,ka,cens)
  {
          (x[,2]==0)*(1-pgompertz(x[,1],la[,1],ka[,1]))*cens[,2]+
          (x[,2]<0)*pgompertz(x[,1],la[,1],ka[,1])*cens[,2]+
          (x[,2]>0)*dgompertz(x[,1],la[,1],ka[,1])*cens[,2]
  }

la.heart <- list(G1=function(p) {cbind(exp(p[1]+p[2]*heart2$age.std),0)},
                 G2=function(p) {cbind(exp(p[3]+p[4]*heart2$age.std),0)})
attr(la.heart,"parameters") <- c("lambda_1","age_1","lambda_2","age_2")

ka.heart <- list(G1=function(p) {cbind(p[1],0)},
                 G2=function(p) {cbind(p[2],0)})

## Note that ka is not restricted here as in McLachlan et al. but it should be ka > 0
## (see the comments in breast_cancer.R for more information)

attr(ka.heart,"parameters") <- c("kappa_1","kappa_2")

cens.heart <- list(G1=function(p) {cbind(0,c(1,1,0)[heart2$status+1])},
                   G2=function(p) {cbind(0,c(1,0,1)[heart2$status+1])})

mix.heart <- function(p) {t(apply(cbind(p[1]+p[2]*heart2$mismat.std+p[3]*heart2$age.std),1,inv.glogit))}
attr(mix.heart,"parameters") <- c("Cons","mismatch","age")

heart.expected2 <- list(
                           G1=function(p) {cbind(qgompertz(0.5,la.heart[[1]](p[1:4])[,1],ka.heart[[1]](p[5:6])[,1]),
                           1)},
                           G2=function(p) {cbind(qgompertz(0.5,la.heart[[2]](p[1:4])[,1],ka.heart[[2]](p[5:6])[,1]),
                           2)}
                           )


heart.moc <- moc(heart2[,c("ftime","status")],density=Gompertz.Surv,joint=TRUE,groups=2,
                 gmu=la.heart,gshape=ka.heart,gextra=cens.heart,gmixture=mix.heart,expected=heart.expected2,
                 pgmu=c(-6,1,-4,0.3),pgshape=c(-0.0015,-0.0055),pgmix=c(1.4,0.4,0.2),
                 gradtol=1e-6,iterlim=200)

heart.moc
## You should have obtained exactly the same results as in
## http://www.jstatsoft.org/v02/i07/0010.pdf
## which I do.

## Then we can make a model with unrestricted mixtures,
## each death cause corresponding to a mixture.

cens.heart2 <- list(G1=function(p) {p1 <- inv.glogit(p[1]); cbind(0,c(1,p1)[heart2$status+1])},
                    G2=function(p) {p2 <- inv.glogit(p[2]); cbind(0,c(1,p2)[heart2$status+1])})
attr(cens.surv2,"parameters") <- c("OR_1","OR_2")

heart.moc2 <- moc(heart2[,c("ftime","status")],density=Gompertz.Surv,joint=TRUE,groups=2,
                  gmu=la.heart,gshape=ka.heart,gextra=cens.heart2,gmixture=mix.heart,
                  pgmu=heart.moc$coef[1:4],pgshape=heart.moc$coef[5:6],pgmix=heart.moc$coef[7:9],pgextra=c(0,0),
                  gradtol=1e-6,iterlim=200)

heart.moc2

## Weibull is easier to fit and less sensitive.

Weibull.Surv <- function(x,la,ka,cens)
  {
      (x[,2]==0)*(1-pweibull(x[,1],shape=la[,1],scale=ka[,1]))*cens[,2]+
          (x[,2]<0)*pweibull(x[,1],shape=la[,1],scale=ka[,1])*cens[,2]+
              (x[,2]>0)*dweibull(x[,1],shape=la[,1],scale=ka[,1])*cens[,2]
  }

weibull.mean <- function(a,b) {b*gamma(1+1/a)}

ka.wb.heart <- list(G1=function(p) {cbind(exp(p[1]+p[2]*heart2$age.std),0)},
                 G2=function(p) {cbind(exp(p[3]+p[4]*heart2$age.std),0)})
attr(ka.wb.heart,"parameters") <- c("log(beta_1)","age_1","log(beta_2)","age_2")

la.wb.heart <- list(G1=function(p) {cbind(exp(p[1]),0)},
                 G2=function(p) {cbind(exp(p[2]),0)})
attr(la.wb.heart,"parameters") <- c("log(alpha_1)","log(alpha_2)")

expected.wb.heart <- list(
                          G1=function(p) {cbind(weibull.mean(la.wb.heart[[1]](p[1:2])[,1],ka.wb.heart[[1]](p[3:6])[,1]),
                          heart2$status)},
                          G2=function(p) {cbind(weibull.mean(la.wb.heart[[2]](p[1:2])[,1],ka.wb.heart[[2]](p[3:6])[,1]),
                          heart2$status)}
                          )
heart2.tmp <- heart2
heart2.tmp$ftime[which(heart2.tmp$ftime==0)] <- 0.1  ## survival time of 0 are problematic ( strange data )

heart.moc.wb <- moc(heart2.tmp[,c("ftime","status")],density=Weibull.Surv,joint=TRUE,groups=2,
                  gmu=la.wb.heart,gshape=ka.wb.heart,gextra=cens.heart2,gmixture=mix.heart,
                  expected=expected.wb.heart,pgmu=c(-1,1),pgshape=c(6,0,8,0),pgmix=c(0,0,0),pgextra=c(0,0),
                  gradtol=1e-6,iterlim=200)

heart.moc.wb
