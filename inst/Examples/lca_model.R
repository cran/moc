## The following examples illustrate how we can
## perform latent class analysis (LCA) with MOC.
## For this purpose, we use some data consisting of 4 dichotomous
## items tabulated with their frequency in
## Hagenaars (1990) Categorical longitudinal data -- log-linear analysis
## of panel, trend and cohort data. Newbury Park : Sage.,

library(moc)    ## load the essential

hag90 <-
structure(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 
1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 
2, 2, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 
2, 1, 2, 1, 2, 59, 56, 14, 36, 7, 15, 4, 23, 75, 161, 22, 115, 
8, 68, 22, 123), .Dim = as.integer(c(16, 5)), .Dimnames = list(
    NULL, c("eq.right.gender", "good.edu", "good.med", "eq.right.work", 
    "freq")))

dummy.code <- function(vec,exclude=NA) {
    tmp <- sapply(levels(factor(vec,exclude=exclude)),"==",factor(vec,,exclude=exclude))*1
    dimnames(tmp) <- list(NULL,paste(make.names(deparse(substitute(vec))),levels(factor(vec,exclude=exclude)),sep="."))
    tmp}                               ## an helper function to dummy indicator coding



dmnom.joint <-                                               ## this version of the multinomial assumes that x is a set
function(x,prob,...) {apply(prob*x+1-x,1,prod,na.rm=TRUE)}   ## of indicator variables and p is the corresponding matrix
                                                             ## of probabilities ( one colum for each category).
                                                             ## It returns the joint density of the indicators.

dmnom <-                                                     ## the same that can be used without the joint option of moc
function(x,prob,...) {prob*x+1-x,1}

mnom.mu <- function(p) {inv.glogit(p)}

dmnom.cat <- function(x,mu=1,shape=1,extra) {extra[,x]}   ## this version of the multinomial assumes that x
                                                          ## is the categorical variable coded from 1 to ncat
                                                          ## The extra parameter has ncat columns corresponding
                                                          ## to the probability of the category

## First, we code the items with indicator variables for each category.

hag90.dummy <- cbind(dummy.code(hag90[,1]),dummy.code(hag90[,2]),dummy.code(hag90[,3]),dummy.code(hag90[,4]))

## A two classes (mixtures) saturated model is obtained by letting the category probabilities
## (item loadings) vary with the classes. Since we have tabulated data with frequencies
## corresponding to the response patterns, we use the frequency as a weight in MOC: wt=hag90[,5].

hag90.sat <- moc(hag90.dummy,density=dmnom,groups=2,wt=hag90[,5],
    gmu=list(G1=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[4])),
             G2=function(p) cbind(inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[7]),inv.glogit(p[8]))),
    pgmu=c(2,2,2,2,-2,-2,-2,-2),pgmix=c(0))

hag90.sat
AIC(hag90.sat,k="BIC")
entropy(hag90.sat)

## Constraining loading equalities for item 1 and 4 and item 2 and 3 is straightforward

hag90.eqcon <- moc(hag90.dummy,density=dmnom,groups=2,wt=hag90[,5],   ## don't forget the frequency weights
    gmu=list(G1=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[2]),inv.glogit(p[2]),inv.glogit(p[1])),
             G2=function(p) cbind(inv.glogit(p[3]),inv.glogit(p[4]),inv.glogit(p[4]),inv.glogit(p[3]))),
    pgmu=c(2,2,-2,-2),pgmix=c(0))

AIC(hag90.sat,hag90.eqcon,k="BIC")
entropy(hag90.sat,hag90.eqcon)

## We could obtain the same result by updating the saturated model

update(hag90.sat,parm=list(mu=hag90.sat$coef[c(1,2,2,1,5,6,6,5)],shape=NULL,extra=NULL,mix=hag90.sat$coef[9]),
                 what=list(mu=c(1,2,2,1,3,4,4,3),shape=NULL,extra=NULL,mix=0), evaluate=TRUE)

## Performing a 2 latent variables with 2 classes each is just the same as doing a 4 classes
## analysis. Applying the same parameters on the loadings makes the correspondance
## between items and latent variables (latent 1 for items 1,4 and latent 2 for items 2,3)

hag90.4g.rest <- moc(hag90.dummy,density=dmnom,groups=4,wt=hag90[,5],   ## don't forget the frequency weights
    gmu=list(G1=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[4])),
             G2=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[4])),
             G3=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[8])),
             G4=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[8]))),
    pgmu=c(2,2,2,2,-2,-2,-2,-2),pgmix=c(0,0,0))

t(matrix(inv.glogit(hag90.4g.rest$coef[9:11]),2,2)) ## show the 4 classes probabilities as a cross-tabulation
                                                    ## of two two classes variables


## We can impose some class probability to be equal to 0 simply by using the appropriate
## mixture function

mix.zero <- function(p) inv.glogit(c(p[1],-Inf,p[2]))

hag90.4g.rest.0 <- moc(hag90.dummy, density=dmnom, groups=4, wt=hag90[,5],  ## don't forget the frequency weights
    gmu=list(G1=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[4])),
             G2=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[4])),
             G3=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[8])),
             G4=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[8]))),
    pgmu=c(2,2,2,2,-2,-2,-2,-2), gmixture=mix.zero,pgmix=c(0,0))

AIC(hag90.4g.rest,hag90.4g.rest.0, k="BIC")

## Specific mixture probabilities corresponding to independance of the 2 latent variables is easy too

mix.ind <- function(p) inv.glogit(c(p[1],p[2],p[1]+p[2]))  ## the same as doing inv.glogit(p[1])%x%inv.glogit(p[2])

hag90.4g.rest.ind <-  moc(hag90.dummy,density=dmnom,groups=4,wt=hag90[,5],   ## don't forget the frequency weights
    gmu=list(G1=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[4])),
             G2=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[4])),
             G3=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[8])),
             G4=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[8]))),
    pgmu=hag90.4g.rest$coef[1:8],gmixture=mix.ind,pgmix=c(0,0))

matrix(mix.ind(hag90.4g.rest.ind$coef[9:10]),2,2)

hag90.4g.rest.ind

## Other mixture probabilities patterns are easily obtained, like symmetric probabilities
## between the latent variables

mix.symm <- function(p) {rbind(inv.glogit(p[1:2])[,c(1,2,2,3)]/c(1,2,2,1))} ## we need 3 base prob. the second is
                                                                            ## divided by 2 and repeated 2 times

hag90.4g.rest.sym <- moc(hag90.dummy,density=dmnom,groups=4,wt=hag90[,5],  ## don't forget the frequency weights
    gmu=list(G1=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[4])),
             G2=function(p) cbind(inv.glogit(p[1]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[4])),
             G3=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[2]),inv.glogit(p[3]),inv.glogit(p[8])),
             G4=function(p) cbind(inv.glogit(p[7]),inv.glogit(p[5]),inv.glogit(p[6]),inv.glogit(p[8]))),
    pgmu=c(2,2,2,2,-2,-2,-2,-2), gmixture=mix.symm, pgmix=c(0,0))

t(matrix(mix.symm(hag90.4g.rest.sym$coef[9:10]),2,2))

AIC(hag90.4g.rest.sym,k="BIC")
