## The following example concerns the data of Thall and Vail (1990)
## on a randomized trial comparing drug to placebo in epileptic
## seizures counts during 4 successive two weeks periods.
## Those data are available at http://www.mrc-bsu.cam.ac.uk/bugs/
## in the examples of BUGS software and also in many R packages:
## faraway (epilepsy), geepack (seizure), GLMMGibbs (epil.bugs), MASS (epil).

library(moc)    ## load the essential
data(epil,package="MASS")

## The data are organized on a time basis (one row for each clinic visit),
## we first organized the data on a subject basis.

epil.y <- t(data.frame(split(as.integer(epil$y),epil$subject)))
epil.base <- t(data.frame(split(epil$base,epil$subject)))
epil.age <- t(data.frame(split(epil$age,epil$subject)))
epil.trt <- t(data.frame(split(unclass(factor(epil$trt)),epil$subject)))
epil.V4 <- t(data.frame(split(epil$V4,epil$subject)))

## We consider a mixture of Poisson distribution for the seizures count.

poiss <-
function(x,la,...) {dpois(x,la)}

## We code the same covariates as used by previous authors on these data,
## except for the special indicator for the 4th visit. Since we are
## using semi-parametric mixture model this effect will be incorparated
## later in time specific components.

epil.lbase.m <- scale(log(epil.base/4),scale=FALSE)
epil.lage.m <- scale(log(epil.age),scale=FALSE)
epil.trt.m <- scale(epil.trt-1,scale=FALSE)
epil.lbasetrt.m <- scale(log(epil.base/4)*(epil.trt==2),scale=FALSE)

poiss.reg <- function(p) exp(p[1]*epil.lbase.m+p[2]*epil.trt.m+
                             p[3]*epil.lbasetrt.m+p[4]*epil.lage.m+p[5])


mu.reg.4 <- list(G1=function(p) poiss.reg(p[1:5]),
                 G2=function(p) poiss.reg(p[c(1:4,6)]),
                 G3=function(p) poiss.reg(p[c(1:4,7)]),
                 G4=function(p) poiss.reg(p[c(1:4,8)]))


epil.1 <- moc(epil.y,density=poiss,gmu=mu.reg.4[1],pgmu=c(0,0,0,0,3))

epil.2 <- moc(epil.y,density=poiss,groups=2,gmu=mu.reg.4[1:2],pgmu=c(epil.1$coef,2),pgmix=0)

epil.3 <- moc(epil.y,density=poiss,groups=3,gmu=mu.reg.4[1:3],pgmu=c(epil.2$coef[1:6],0.5),pgmix=c(0,0))

epil.4 <- moc(epil.y,density=poiss,groups=4,gmu=mu.reg.4,pgmu=c(epil.2$coef[1:7],0),pgmix=c(0,0,0),iterlim=200)

AIC(epil.1,epil.2,epil.3,epil.4,k="BIC")

## A 3 group model is sufficient to capture the within subject dependance (between subject random effect).

## The previous models are the finite mixture counterparts of the subject random effects model.
## For these data, the models (as the standard random effect model) are highly sensitive to the starting values
## since the likelihood contains a lot of local maximums. Depending on the starting values we find
## similar parameters estimates as previously reported by other authors which doesn't appear
## to be maximum likelihood. Torturing the likelihood help us find a better solution.

epil.3.torture <- confint(epil.3,profiling="complete",iterlim=15) ## torturing the likelihood is computer intensive
                                                                  ## (about 10 min. on a fast computer)

## The preceeding torture found 19 better solutions (with respect to the likelihood) than the
## one previously found. The best solution being: (found in the torture of the trt*lbase interaction parameter) 

epil.3.fin <- update(epil.3,parm=c(0.9784716,-1.4367440,0.6278039,0.681641,1.620017,2.573327,0.7142217,-1.449593,-1.378657),evaluate=TRUE)

plot(epil.3.fin)

## As we can see the preceeding model does not capture all variations in time, the
## finite mixture equivalent of a full random effect model (subject random effect and
## time within subject random effect) is obtained by simply adding time specific parameters
## to each mixture group.

epil.time <- t(data.frame(split(epil$period,epil$subject)))

poiss.reg.time <- function(p) exp(p[1]*epil.lbase.m+p[2]*epil.trt.m+
                                p[3]*epil.lbasetrt.m+p[4]*epil.lage.m+
                             p[5]*(epil.time==1)+p[6]*(epil.time==2)+p[7]*(epil.time==3)+p[8]*(epil.time==4))

epil.3.times <-  moc(epil.y,density=poiss,groups=3,
                     gmu=list(G1=function(p) poiss.reg.time(p[1:8]),
                              G2=function(p) poiss.reg.time(p[c(1:4,9:12)]),
                              G3=function(p) poiss.reg.time(p[c(1:4,13:16)])),
                     pgmu=c(epil.3.fin$coef[1:4],rep(epil.3.fin$coef[5],4),rep(epil.3.fin$coef[6],4),rep(epil.3.fin$coef[7],4)),
                     pgmix=epil.3.fin$coef[8:9])

 AIC(epil.3.fin,epil.3.times,k="BIC")

## Including the time specific parameters moderately decreased the likelihood and
## the entropy. As we can see it helps essentially to capture some big increase
## at the 3rd visit ( which is due to some few influential observations).
## The 4th visit (as specially included by other authors) does not appear to differ substantially
## from the other points (except the 3rd which is the one that differs for a small group).

plot(epil.3.times)

## Note that the parameters estimates of the treatment effect and the interaction
## between treatment and lbase  are highly correlated

 cov2cor(epil.3.times$cov)[2,3]

## the inclusion of this interaction in conjonction with some influential data at the 3rd
## visit may explain the high sensitivity of the parameters estimates discussed above.

## It would be possible to use a normal quadrature as in the example betablock.R to compare the
## preceeding models to a standard normal random effect model. However in the case of the Poisson
## distribution, the random effect being multiplicative with respect to the mean parameter of
## the Poisson, the model then easily becomes unstable and require a large number of quadrature
## points. The semiparametric model used here represents a more flexible and stable approach.

## The treatment variable here is a true fixed variable that was fixed by design while
## the other variables used as covariates are in fact observational variables
## subject to sampling variations like the response variables. As such the
## treatment variable would be more accurately treated if included in the
## mixture probability function.

poiss.reg.time.mix.trt <- function(p) exp(p[1]*epil.lbase.m+p[2]*epil.lage.m+
                             p[3]*(epil.time==1)+p[4]*(epil.time==2)+p[5]*(epil.time==3)+p[6]*(epil.time==4))

mix.trt <- function(p) {t(apply(cbind(epil.trt[,1]),1,function(trt) inv.glogit(c(p[1]+p[2]*trt,p[3]+p[4]*trt))))}

epil.3.times.mix.trt <-  moc(epil.y,density=poiss,groups=3,
                     gmu=list(G1=function(p) poiss.reg.time.mix.trt(p[1:6]),
                              G2=function(p) poiss.reg.time.mix.trt(p[c(1:2,7:10)]),
                              G3=function(p) poiss.reg.time.mix.trt(p[c(1:2,11:14)])),
                     pgmu=epil.3.times$coef[c(1,4,5:16)],
                     gmixture=mix.trt,pgmix=c(epil.3.times$coef[17],0,epil.3.times$coef[18],0))

## As we can see now have 2 lower groups almost constant in time and a much higher group
## which shows more variations in time.

epil.3.times.mix.trt
plot(epil.3.times.mix.trt)

## The prior and posterior proportions falling in each group according to treatment
## can be obtained easily by the following function

epil.3.timetrt.obsfit <- obsfit.moc(epil.3.times.mix.trt,along=list(trt=epil$trt[epil$period==1]))
epil.3.timetrt.obsfit$"Mean Prior Probabilities"
epil.3.timetrt.obsfit$"Mean Posterior Probabilities"

## This clearly shows that the highest group (about 12%) is totally unaffected by treatment
## and that while 81% of the placebo subjects fall in the middle group this proportion
## is only 49% for progabide the other 37% of treated patients showing a reduced level (group 3)
## compared to 9% for the untreated.


## Finally we consider an alternative way of accounting for the dependance between the
## measurements on the same subject.
## This example is used solely to show how to implement a hidden markov chain
## in MOC, a model that incorporate a form of relationship on previous observations.
## That is we consider that between occasion the subjects have some probabilities of
## doing a transition from one mixture group to another. In the case of 4 occasions with
## 3 basic mixture groups, this makes 3*3*3*3=81 possible combinations of groups by occasion.
## However those 81 groups combinations have structured probabilities that depends on few
## parameters.

mu.hmc <- list()     ## this will hold the list of mu functions for the 81 combinations


## The poisson distribution of each observations needs to be computed only for 3 different states (group).
## The 81 combinations use those same values in different arrangements, so to speed up the computations
## we use a global variable that will hold the values for all groups and then compute the likelihood
## of each specific combination.

hmc.dens.tmp <- array(NA,c(dim(epil.y),3)) ## holds the density values for each subject, time and state.

## hmc.dens receive an index (as mu parameter) representing the time by state combination
## for example index=c(1,3,2,3) means that time 1 is in state 1, time 2 in state 3 ans so on.

hmc.dens <- function(x,index,...) {cbind(hmc.dens.tmp[,1,index[1,1]],hmc.dens.tmp[,2,index[1,2]],
                                   hmc.dens.tmp[,3,index[1,3]],hmc.dens.tmp[,4,index[1,4]])}

## mu.hmc simply returns the index combination mentioned above.

for(i in 1:3) for(j in 1:3) for(k in 1:3) for(l in 1:3) {
    mu.hmc[[paste("G",i,j,k,l,sep="")]] <-  evalq(substitute(function(p) {rbind(c(i,j,k,l))},list(i=i,j=j,k=k,l=l)))}

 mu.hmc <- lapply(mu.hmc,function(ll) eval(ll[-4]))  ## make the elements functions

## However we need to replace the first function (which is going to be called first) such that it will compute
## and assign the poisson density (hmc.dens.tmp) for all states at once in .GlobalEnv

mu.hmc[["G1111"]] <- function(p) {
    assign("p.tmp",p,envir=.GlobalEnv)  ## the parameters also need to be available in .GlobalEnv, then the density for the 3 states
    eval(expression(hmc.dens.tmp[,,1] <- poiss(epil.y,poiss.reg.time.mix.trt(c(p.tmp[1:2],p.tmp[rep(3,4)])))),envir=.GlobalEnv) 
    eval(expression(hmc.dens.tmp[,,2] <- poiss(epil.y,poiss.reg.time.mix.trt(c(p.tmp[1:2],p.tmp[rep(4,4)])))),envir=.GlobalEnv)
    eval(expression(hmc.dens.tmp[,,3] <- poiss(epil.y,poiss.reg.time.mix.trt(c(p.tmp[1:2],p.tmp[rep(5,4)])))),envir=.GlobalEnv)
    rbind(c(1,1,1,1))}         ## still returns the index




## The mixture probability function then construct the probablities of the combinations by
## first assigning  probabilities to the 3 states to the first observation and then
## multiplying by the appropriate transition probabilities for subsequent observations.
## The treatment effect could be included in the mixture, but we don't include it for simplicity.

mix.hmc <- function(p) { prob <- inv.glogit(p[1:2])                                                 ## 3 groups prob. for first observations
                         p.trans <- rbind(inv.glogit(p[3:4]),inv.glogit(p[5:6]),inv.glogit(p[7:8])) ## transition probabilities
                         p.joint <- numeric()
                         for(i in 1:3) for(j in 1:3) for(k in 1:3) for(l in 1:3) p.joint <-
                             c(p.joint,prob[,i]* p.trans[i,j]*p.trans[j,k]*p.trans[k,l])            ## joint probabilities of group 
                         rbind(p.joint)}                                                            ## by occasion combinations

epil.3.hmc <- moc(epil.y,density=hmc.dens,groups=81,gmu=mu.hmc,pgmu=c(epil.3.fin$coef[c(1,4,5:7)]),gmixture=mix.hmc,pgmix=c(epil.3.fin$coef[8:9],0,0,0,0,0,0),iterlim=200)    ## takes about 3 minutes on a fast computer

## The default plot method of MOC is not usefull for HMC models, 
## instead we make a graph of the poisson mean corrected levels with arrows
## representing the transitions ( the length of the array beeing proportional
## to the transition probabilities ).

 matplot(t(matrix(1:4,3,4)),t(matrix(exp(epil.3.hmc$coef[3:5]),3,4)),xlab="Time",ylab="Mean level",type="p")
for(i in (1:3)+0.01){
arrows(rep(i,3),rep(exp(epil.3.hmc$coef[3]),3),i+inv.glogit(epil.3.hmc$coef[8:9]),
       exp(epil.3.hmc$coef[3])+inv.glogit(epil.3.hmc$coef[8:9])*(exp(epil.3.hmc$coef[3:5])-exp(epil.3.hmc$coef[3])),
       length=.05,col=c("black","red","green"))

 arrows(rep(i,3),rep(exp(epil.3.hmc$coef[4]),3),i+inv.glogit(epil.3.hmc$coef[10:11]),
        exp(epil.3.hmc$coef[4])+inv.glogit(epil.3.hmc$coef[10:11])*(exp(epil.3.hmc$coef[3:5])-exp(epil.3.hmc$coef[4])),
        length=.05,col=c("black","red","green"))

 arrows(rep(i,3),rep(exp(epil.3.hmc$coef[5]),3),i+inv.glogit(epil.3.hmc$coef[12:13]),
        exp(epil.3.hmc$coef[5])+inv.glogit(epil.3.hmc$coef[12:13])*(exp(epil.3.hmc$coef[3:5])-exp(epil.3.hmc$coef[5])),
        length=.05,col=c("black","red","green"))}


