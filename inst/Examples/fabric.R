## Fitting the model in McLachlan & Peel pp.155-157 on fabric data
## available at http://www.maths.uq.edu.au/~gjm/DATA/mmdata.html

library(moc)

fabric <-  read.table("DataBank/Fabric.dat",header=TRUE)  #change with your own reading command

## It is used here to illustrate the mixture of Poisson regression
## with a random intercept that is
#  log mu_i_j = b0_i + b1 * log l_j
## or equivalently
#  mu_i_j = exp(  b0_i + b1 * log l_j )

poiss <-
function(x,la,...) {dpois(x,la)}

gmup <-list(
Group1=function(pmu) {exp(pmu[2]+pmu[1]*cbind(log(fabric[,1])))},
Group2=function(pmu) {exp(pmu[3]+pmu[1]*cbind(log(fabric[,1])))},
Group3=function(pmu) {exp(pmu[4]+pmu[1]*cbind(log(fabric[,1])))}
)
attr(gmup,"parameters") <- c("b1.log_length","b0_1","b0_2","b0_3")
     
fab1 <-
moc(fabric[, 2], density = poiss, groups = 1, gmu = gmup[-(2:3)], 
    pgmu = c(0.99, -4.17), gradtol = 1e-05)
attr(fab1$gmu,"parameters") <- c("b1.log_length","b0_1")

fab2 <-
moc(fabric[, 2], density = poiss, groups = 2, gmu = gmup[-3],
    pgmu =  c(0.8, -2.3, -2.7), pgmix=0, gradtol = 1e-05)
attr(fab2$gmu,"parameters") <- c("b1.log_length","b0_1","b0_2")

fab3 <-
moc(fabric[, 2], density = poiss, groups = 3, gmu = gmup,
    pgmu =  c(0.8, -2.3, -2.7, -3), pgmix = c(0,0), gradtol = 1e-05)

AIC(fab1,fab2,fab3,k="BIC")

## A 2 groups mixture is clearly needed in this case as supported by the ICL-BIC
## and inspection of the residuals. Adding another group however is unnecessary.

plot(residuals(fab1))
plot(residuals(fab2))

fab2

## We can obtain the parameters u1, u2 u3 of table 5.5 in McLachlan & Peel p. 156
## with the help of the function confint.moc. The standard deviation are different
## because they are computed differently.

confint(fab2,parm=list(~1/(1+exp(p4))*p2+exp(p4)/(1+exp(p4))*p3,
             ~(p2-p3)*exp(p4)/(1+exp(p4)),~(p3-p2)/(1+exp(p4))),
        profiling="simple",offscal=seq(-6,6,0.5))

## Note that we use the internal option offscal which determine the search range
## for the parameters.
