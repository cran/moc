#
#  moc : Function to fit general mixture of curves models 
#         
#  Copyright (C) 2000 Bernard Boulerice
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public Licence as published by
#  the Free Software Foundation; either version 2 of the Licence, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public Licence for more details.
#
#  You should have received a copy of the GNU General Public Licence
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     moc(y,density=NULL,joint=FALSE, groups=1,
#       gmu=NULL, gshape=NULL, gextra=NULL, gmixture=glogit, expected = NULL,
#       pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, wt=1,
#       ndigit=10, gradtol=0.0001, 
#       steptol=gradtol,iterlim=100,...)
#
#  DESCRIPTION
#
# Function to fit general trajectory models with mixture
#
#
moc<- function(y,density=NULL,joint=FALSE, groups=1,
       gmu=NULL, gshape=NULL, gextra=NULL, gmixture=glogit, expected = NULL,
       pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, wt=1,
       ndigit=10, gradtol=0.0001, 
       steptol=gradtol,iterlim=100,...)
{
glogit<- function(gmix) {rbind(c(1,exp(gmix)))/(1+sum(exp(gmix)))}
call <- sys.call()
resp<-as.matrix(y)
ndim<-dim(resp)
n<-ndim[1]
nt<-ndim[2]
ng<-groups
#
# check density
#
if(!is.function(density)) stop("density must be a function")
if(length(formals(density))>4) stop("density must not use more than 4 arguments")
dens<-density
#
# count the number of parameters
#
npl <- length(pgmu)
nps <- length(pgshape)
npext <- length(pgextra)
npmix <- length(pgmix)
np <- npl+nps+npmix+npext
#
# create local versions of functions and
# check that the functions return correct values
#
repn<-deparse(substitute(rep(1,n),list(n=n)))
if(length(dim(gmu(pgmu)))!=2||dim(gmu(pgmu))[2]!=ng*nt) 
      stop("gmu must return a matrix with vectors of length times*groups=",ng*nt)
if(any(is.na(gmu(pgmu)))) stop("The gmu function returns NAs")
if(dim(gmu(pgmu))[1]==1){
   fcall<-paste(deparse(gmu)[1],collapse="")
   fbody<-paste(deparse(gmu)[-1],collapse="")
   ftext<-paste(fcall,"{",fbody,"%x%",repn,"}")
   gmu1<-eval(parse(text=ftext))} else
   if(dim(gmu(pgmu))[1]==n) gmu1<-gmu else
      stop("gmu should return a matrix of length 1 or",n)
#
if(is.null(gshape) && is.null(pgshape)) gshape1<- function(...) 1 else
  {
  if(length(dim(gshape(pgshape)))!=2||dim(gshape(pgshape))[2]!=ng*nt) 
     stop("gshape must return a matrix with vectors of length times*groups=",ng*nt)
  if(any(is.na(gshape(pgshape)))) stop("The shape function returns NAs")
  if(dim(gshape(pgshape))[1]==1) {
   fcall<-paste(deparse(gshape)[1],collapse="")
   fbody<-paste(deparse(gshape)[-1],collapse="")
   ftext<-paste(fcall,"{",fbody,"%x%",repn,"}")
   gshape1<-eval(parse(text=ftext))} else
   if(dim(gshape(pgshape))[1]==n ) gshape1<-gshape else
        stop("gshape should return a matrix of length 1 or",n)
   }
#
if(is.null(gextra) && is.null(pgextra)) gextra1<- function(...) 1 else
  {
  if(length(dim(gextra(pgextra)))!=2||dim(gextra(pgextra))[2]!=ng*nt) 
     stop("gextra must return a matrix with vectors of length times*groups=",ng*nt)
  if(any(is.na(gextra(pgextra)))) stop("The extra parameters function returns NAs")
  if(dim(gextra(pgextra))[1]==1) {
   fcall<-paste(deparse(gextra)[1],collapse="")
   fbody<-paste(deparse(gextra)[-1],collapse="")
   ftext<-paste(fcall,"{",fbody,"%x%",repn,"}")
   gextra1<-eval(parse(text=ftext))} else
   if(dim(gextra(pgextra))[1]==n ) gextra1<-gextra else
        stop("gextra should return a matrix of length 1 or",n)
   }
#
if(is.null(gmixture) && is.null(pgmix)) { gmixt<-function(...) 1 } else
   {
  if(dim(gmixture(pgmix))[2]!=ng) stop("gmixture must return a matrix with vectors of length groups=",ng)
  if(any(is.na(gmixture(pgmix)))) stop("The mixture function returns NAs")
  if(any(gmixture(pgmix)<0)||any(apply(gmixture(pgmix),1,sum))!=1) 
  stop("The mixture function components must be >=0 and sum to 1")
  if(dim(gmixture(pgmix))[1]==1 ) {
    fcall<-paste(deparse(gmixture)[1],collapse="")
    fbody<-paste(deparse(gmixture)[-1],collapse="")
    ftext<-paste(fcall,"{",fbody,"%x%",repn,"}")
    gmixt<-eval(parse(text=ftext))} else
    if(dim(gmixture(pgmix))[1]==n) gmixt<-gmixture else
      stop("gmixture should return a matrix of length 1 or",n)
    }

if (is.null(expected))
  {
    environment(gmu1)<-globalenv()
    expected<- eval(expression(function(p) {gmu(p[1:npl])}),list(gmu=gmu1,npl=npl))
  }
else {
  ptot <- c(pgmu,pgshape,pgextra) 
  if(length(dim(expected(ptot)))!=2||dim(expected(ptot))[2]!=ng*nt) 
   stop("expected must return a matrix with vectors of length times*groups=",ng*nt)
  if(any(is.na(expected(ptot)))) stop("The expected function returns NAs")
  if(dim(expected(ptot))[1]==1){
  fcall<-paste(deparse(expected)[1],collapse="")
  fbody<-paste(deparse(expected)[-1],collapse="")
  ftext<-paste(fcall,"{",fbody,"%x%",repn,"}")
  expected<-eval(parse(text=ftext))} else
  if(dim(expected(ptot))[1]!=n)  stop("gmu should return a matrix of length 1 or",n)
   environment(expected)<-globalenv()
     }
    
#
# check weights 
#
if(length(wt)==1) wt <- rep(wt,n)
if(length(wt)!=n) stop("wt must be the same length as the other variables")
if(min(wt)<0) stop("All weights must be non-negative")
#
#define the likelihood functions
#
tmp<-array(resp,dim=c(n,nt,ng))
likelihood<-if(joint){
  function(p){
   cm<-array(gmu1(p[1:npl]),dim=c(n,nt,ng))
   sc<-array(gshape1(p[(npl+1):(npl+nps)]),dim=c(n,nt,ng))
   extra<-array(gextra1(p[(npl+nps+1):(npl+nps+npext)]),dim=c(n,nt,ng))
   dens(tmp,cm,sc,extra)}} else
   {function(p){
     cm<-array(gmu1(p[1:npl]),dim=c(n,nt,ng))
     sc<-array(gshape1(p[(npl+1):(npl+nps)]),dim=c(n,nt,ng))
     extra<-array(gextra1(p[(npl+nps+1):(npl+nps+npext)]),dim=c(n,nt,ng))
     apply(dens(tmp,cm,sc,extra),c(1,3),prod,na.rm=TRUE)}}
#   
#   
loglike<-function(p)
{
 mix<-gmixt(p[(npl+nps+npext+1):(npl+nps+npext+npmix)])
 sum(-wt*log(apply(likelihood(p[1:(npl+nps+npext)])*mix,1,sum)))
}
#
# check that the likelihood returns an appropriate value and minimize
#
p<-c(pgmu,pgshape,pgextra,pgmix)
if(is.na(loglike(p)))
	stop("Likelihood returns NAs: probably invalid initial values")
if(np>0){
	dt<-system.time(z0 <- nlm(loglike,p=p,hessian=TRUE,
		ndigit=ndigit,gradtol=gradtol,steptol=steptol,
		iterlim=iterlim,...))[3]} else
         {dt<-system.time(z0 <- list(minimum=fscale,estimate=p,code=0,iterations=0))[3]}
cat("\n Estimation took ",dt," seconds.\n")
#
# compute cov and se's
#
if(np==0)cov <- NULL else 
if(np==1)cov <- 1/z0$hessian else 
{
	a <- if(any(is.na(z0$hessian))||any(abs(z0$hessian)==Inf)) 0 else 
               qr(z0$hessian)$rank
	if(a==np)cov <- solve(z0$hessian) else 
           cov <- matrix(NA,ncol=np,nrow=np)}
se <- sqrt(diag(cov))
#
# compute mixture and posterior probabilities
#
pmix<-gmixt(z0$estimate[(npl+nps+npext+1):(npl+nps+npext+npmix)])
post<-likelihood(z0$estimate)*pmix
post<-post/apply(post,1,sum)
#
# compute fitted and obserserved values 
#
mpost<-apply(post,2,mean)
fitted<-array(expected(z0$estimate),dim=c(n,nt,ng))
fitted.mean<-apply(aperm(fitted,c(1,3,2))*rep(post,nt),c(2,3),mean,na.rm=TRUE)/mpost
dimnames(fitted.mean)<-list(paste("Group",1:ng,sep=""), ifelse(is.na(nchar(temp<-colnames(resp))[1:nt]) ,paste("Time",1:nt,sep=""),temp))
observed.mean<-t(apply(array(apply(post,2,"*",resp),dim=c(n,nt,ng)),c(2,3),mean,na.rm=TRUE))/mpost
dimnames(observed.mean)<-list(paste("Group",1:ng,sep=""), ifelse(is.na(nchar(temp<-colnames(resp))[1:nt]) ,paste("Time",1:nt,sep=""),temp))
#
# clean functions environment and
# return a list of class moc
#
environment(dens)<-globalenv()
environment(gmu1)<-globalenv()
environment(gshape1)<-globalenv()
environment(gextra1)<-globalenv()
environment(gmixt)<-globalenv()

z1 <- list(
	call=call,
        resp=call[[2]],
	density=dens,
        joint=joint,
        nsubject=n,
        ntimes=nt,
        nobs=sum(!is.na(resp)),
        groups=ng,
        npar=c(npl,nps,npext,npmix),
	gmu=gmu1,
	gshape=gshape1,
	gextra=gextra1,
	gmixture=gmixt,
	expected = expected,
        prior.weights=wt,
	loglikelihood=-z0$minimum,
	df=sum((!is.na(resp))*wt)-np,
	AIC=2*z0$minimum+2*np,
        BIC=2*z0$minimum+np*log(sum((!is.na(resp))*wt)),
	coefficients=z0$estimate,
	cov=cov,
        fitted.mean=fitted.mean,
        observed.mean=observed.mean,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "moc"
return(z1) }

post<-function(object,...) UseMethod("post")

post.moc<-function(object,...)
{
 n<-object$nsubject
 dim1<-c(n,object$ntimes,object$groups)
 nm<-object$npar[1]
 ns<-object$npar[2]
 nextra<-object$npar[3]
 nmix<-object$npar[4]
 pm<-object$coefficients[1:nm]
 ps<-object$coefficients[(1+nm):(nm+ns)]
 pext<-object$coefficients[(1+nm+ns):(nm+ns+nextra)]
 pmix<-object$coefficients[(1+nm+ns+nextra):(nm+ns+nextra+nmix)]
 if(object$gmixture(pmix)==1) 1 else
   {
   y<-array(as.matrix(eval(object$resp)),dim=dim1)
   m<-array(object$gmu(pm),dim=dim1)
   s<-array(object$gshape(ps),dim=dim1)
   extra<-array(object$gextra(pext),dim=dim1)
   mix<-object$gmixture(pmix)
   if(object$joint) {post<-object$density(y,m,s,extra)*mix} else
   {post<-apply(object$density(y,m,s,extra),c(1,3),prod,na.rm=TRUE)*mix}
   dimnames(post)<-list(NULL,paste("Group",1:dim1[3],sep=""))
   post/apply(post,1,sum)
   }
}

fitted.moc<-function(object,...)
{
 n<-object$nsubject
 dim1<-c(n,object$ntimes,object$groups)
 fit<-array(object$expected(object$coefficients),dim=dim1)
 dimnames(fit)<-list(NULL, ifelse(is.na(nchar(temp<-colnames(eval(object$resp)))[1:dim1[2]]) ,paste("Time",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
 return(fit)
}

residuals.moc<-function(object,...,type="deviance",post.weight=TRUE,within=T)
{
 choices<-c("deviance","response")
 type<-match.arg(type,choices)
 n<-object$nsubject
 nt<-object$ntimes
 ng<-object$groups
 dim1<-c(n,nt,ng)
 nm<-object$npar[1]
 ns<-object$npar[2]
 nextra<-object$npar[3]
 pm<-object$coefficients[1:nm]
 ps<-object$coefficients[(1+nm):(nm+ns)]
 pext<-object$coefficients[(1+nm+ns):(nm+ns+nextra)]
 y<-array(as.matrix(eval(object$resp)),dim=dim1)
 m<-array(object$gmu(pm),dim=dim1)
 s<-array(object$gshape(ps),dim=dim1)
 extra<-array(object$gextra(pext),dim=dim1)
 wts<-object$prior.weights
 wpost<-post(object)
 if (within) wpost<-t(t(wpost)/apply(wpost,2,mean))
 response<-(y-fitted(object))
 switch(type,
        response= res<-response,
        deviance= res<- if(!object$joint) {sqrt(2*wts*(log(object$density(y,y,s,extra)/object$density(y,m,s,extra))))*sign(response)}
        else { tmp<-array(NA,dim1)
               for (i in 1:nt)
                 {mi<-y; mi[,i,]<-m[,i,]
                  tmp[,i,]<-sqrt(2*wts*(log(object$density(y,y,s,extra)/object$density(y,mi,s,extra))))*sign(response[,i,])
                } 
            tmp }
       )

 dimnames(res)<-list(NULL,ifelse(is.na(nchar(temp<-colnames(eval(object$resp)))[1:dim1[2]]) ,paste("Time",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
 if(post.weight) res<-res*array(wpost[,rep(1:ng,rep(nt,ng))],dim=dim1)
 return(res)
}


print.moc<-function(x,...)
{
 object<-x
 cat("\n\t\t",object$groups,"Mixture of Curves\n\n\n")
 if(object$joint) cat("joint ")
 cat("Density:\n",gsub("[{} ]","",paste(deparse(object$density)[-1],
    sep="",collapse=" ")),"\n") 
 cat("\nLocation function:\n",gsub("[{} ]","",paste(deparse(object$gmu)[-1],
    sep="",collapse=" ")),"\n") 
 cat("\nExpectation function:\n",gsub("[{} ]","",paste(deparse(object$expected)[-1],
    sep="",collapse=" ")),"\n") 
 cat("\nShape function:\n",gsub("[{} ]","",paste(deparse(object$gshape)[-1],
    sep="",collapse=" ")),"\n")
 cat("\nExtra parameter function:\n",gsub("[{} ]","",paste(deparse(object$gextra)[-1],
    sep="",collapse=" ")),"\n")
 cat("\nMixture function:\n",gsub("[{} ]","",paste(deparse(object$gmixture)[-1],
    sep="",collapse=" ")),"\n")
 cat("\n\n\t\t\tMaximum Likelihood Estimates\n\n")
 if(object$code>2) {cat("\n\nWARNING - CONVERGENCE NOT ACHIEVED - ERROR",object$code,"\n\n\n")}
 cat(formatC(c("",""," Standard","Wald Chi-Sq:"),width=13),"\n")
 cat(formatC(c("Parameter    ","Estimate"," Error"," Param = 0"," Prob > Wald"),width=13),"\n")
 object.se<-sqrt(diag(object$cov))
 cat("\nLocation:\n\n")
 nl<-object$npar[1];ns<-object$npar[2];nextra<-object$npar[3];nm<-object$npar[4]
 param<-attr(object$gmu,"parameters")
 if (length(param) != nl) param<-"       "
 coeftable<-formatC(cbind(est<-object$coef[1:nl],se<-object.se[1:nl],
                         w<-(est/se)^2,(1-pchisq(w,1))),digits=5,width=13)
 cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
  if(object$npar[2]>0)
 {
   cat("\nShape:\n\n")
 param<-attr(object$gshape,"parameters")
 if (length(param) != ns) param<-"       "
   coeftable<-formatC(cbind(est<-object$coef[(nl+1):(nl+ns)],se<-object.se[(nl+1):(nl+ns)],
                w<-(est/se)^2,(1-pchisq(w,1))),digits=5,width=13)
   cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
 }
 if(object$npar[3]>0)
 {
   cat("\nExtra Parameters:\n\n")
 param<-attr(object$gextra,"parameters")
 if (length(param) != nextra) param<-"       "
   coeftable<-formatC(cbind(est<-object$coef[(nl+ns+1):(nl+ns+nextra)],se<-object.se[(nl+ns+1):(nl+ns+nextra)],
                w<-(est/se)^2,(1-pchisq(w,1))),digits=5,width=13)
   cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
 }
 if(object$npar[4]>0)
 {
   cat("\nMixture:\n\n")
 param<-attr(object$gmixture,"parameters")
 if (length(param) != nm) param<-"       "
   coeftable<-formatC(cbind(est<-object$coef[(nl+ns+nextra+1):(nl+ns+nextra+nm)],se<-object.se[(nl+ns+nextra+1):(nl+ns+nextra+nm)],
                w<-(est/se)^2,(1-pchisq(w,1))),digits=5,width=13)
   cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
 
   cat("\nMean Mixture Probabilities:\n")
   mmix<-apply(object$gmixture(est),2,mean)
   names(mmix)<-paste("Group",1:object$groups,sep="")
   print(mmix,digits=5)
 }
 cat("\n",formatC(c("-2 Log Likelihood","df","AIC","BIC"),width=11),"\n")
 cat("      ",formatC(c(-2*object$loglikelihood,object$df,object$AIC,object$BIC),digits=8,width=11),"\n")
 cat("\nMean predicted values:\n")
 print(object$fitted.mean,digits=5)
 cat("\nMean observed values:\n")
 print(object$observed.mean,digits=5)
 invisible(x)
}


AIC.moc <- function(object,...,k=2) 
{
  mlist<-  list(object,...)
  ml<-length(mlist)
  nobslist<-sapply(mlist,function(tmp) tmp$nobs)
  nglist<-sapply(mlist,function(tmp) tmp$groups)
  if((k==0) && ((max(nglist)!=min(nglist)) || (max(nobslist)!=min(nobslist))))
    warning("log Likelihood should only be used to compare nested models on the same observations with the same number of mixture.\nTry using AIC or BIC instead.\n") else
  {if((k!="BIC") && (max(nobslist)!=min(nobslist)))
    warning("AIC like statistics should not be used to compare models with differing number of observations, use BIC instead.\n")}
  cnames<-as.character(match.call()[-1])[1:ml]
  if(k=="BIC")
  val<-as.data.frame(t(sapply(mlist,function(tmp) c(-2*tmp$loglikelihood+log(tmp$nobs)*sum(tmp$npar),tmp$nobs-sum(tmp$npar)))))
  else
  val<-as.data.frame(t(sapply(mlist,function(tmp) c(-2*tmp$loglikelihood+k*sum(tmp$npar),tmp$nobs-sum(tmp$npar)))))
  names(val)<-c(switch(as.character(k),"0"="-2*logLik","2"="AIC",BIC="BIC","generalized AIC"),"Df")
  row.names(val)<-cnames                     
  val
}

logLik.moc <- function(object,...)
  {
    val <- object$loglikelihood
    attr(val,"nobs") <- object$nobs
    attr(val,"df") <- object$df
    class(val) <- "logLik"
    val
  }
   
  
".First.lib"<-
function(lib, pkg)
{
        cat("\n\tMOC library version 0.5-9\n")
        cat("\tThis library is provided by Bernard Boulerice <Bernard.Boulerice@umontreal.ca>\n\n")
        if (!exists('AIC',envir=.GlobalEnv)){
    eval(expression(AIC<-function(object,...,k=2) UseMethod('AIC')),envir=.GlobalEnv)} 
}
