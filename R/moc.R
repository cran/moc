#
#  moc : Function to fit general mixture of curves models 
#         
#  Copyright (C) 2000-2002 Bernard Boulerice
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
#       gmu=NULL, gshape=NULL, gextra=NULL, gmixture=inv.glogit, expected = NULL,
#       pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, wt=1, scale.weight=FALSE, data=NULL,
#       ndigit=10, gradtol=0.0001, 
#       steptol=gradtol,iterlim=100,print.level=2,...)
#
#  DESCRIPTION
#
# Function to fit general nonlinear mixture models
#
#
moc<- function(y,density=NULL,joint=FALSE, groups=1,
       gmu=NULL, gshape=NULL, gextra=NULL, gmixture=inv.glogit, expected = NULL,
       pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, scale.weight=FALSE, wt=1, data=NULL,
       ndigit=10, gradtol=0.0001, 
       steptol=gradtol,iterlim=100,print.level=2,...)
{
  if(!is.null(data)) { attach(data); on.exit(detach(data))  }
  call <- sys.call()
  resp<-as.matrix(eval(y))
  ndim<-dim(resp)
  n<-ndim[1]
  nt<-ndim[2]
  ng<-groups
  inv.glogit<- if(ng==1) {function(gmix) {cbind(1)}} else {function(gmix) {rbind(c(1,exp(gmix)))/(1+sum(exp(gmix)))}}
  if(groups > 1) attr(inv.glogit,"parameters")<-paste("  G",2:groups," vs G1",sep="")
                                        #
                                        # check density
                                        #
  if(!is.function(density)) stop("density must be a function")
  if(length(formals(density))>4) stop("density must not use more than 4 arguments")
  .density<-density
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
                                        #repn<-deparse(substitute(rep(1,n),list(n=n)))
  if(length(dim(gmu(pgmu)))!=2||dim(gmu(pgmu))[2]!=ng*nt) 
    stop("gmu must return a matrix with vectors of length times*groups=",ng*nt)
  if(any(is.na(gmu(pgmu)))) stop("The gmu function returns NAs")
  if(dim(gmu(pgmu))[1]==1){
    fcall<-paste(deparse(gmu)[1],collapse="")
    fbody<-paste(deparse(gmu)[-1],collapse="\n")
    ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",ng*nt,",byrow=TRUE)}"),collapse="")
    .gmu<-eval(parse(text=ftext))
    environment(.gmu)<-environment(gmu)
  } else
  if(dim(gmu(pgmu))[1]==n) .gmu<-gmu else
  stop("gmu should return a matrix of length 1 or",n)
                                        #
  if(is.null(gshape) && is.null(pgshape)) .gshape<- function(...) 1 else
  {
    if(length(dim(gshape(pgshape)))!=2||dim(gshape(pgshape))[2]!=ng*nt) 
      stop("gshape must return a matrix with vectors of length times*groups=",ng*nt)
    if(any(is.na(gshape(pgshape)))) stop("The shape function returns NAs")
    if(dim(gshape(pgshape))[1]==1) {
      fcall<-paste(deparse(gshape)[1],collapse="")
      fbody<-paste(deparse(gshape)[-1],collapse="\n")
      ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",nt*ng,",byrow=TRUE)}"),collapse="")
      .gshape<-eval(parse(text=ftext))
      environment(.gshape)<-environment(gshape)
    } else
    if(dim(gshape(pgshape))[1]==n ) .gshape<-gshape else
    stop("gshape should return a matrix of length 1 or",n)
  }
                                        #
  if(is.null(gextra) && is.null(pgextra)) .gextra<- function(...) 1 else
  {
    if(length(dim(gextra(pgextra)))!=2||dim(gextra(pgextra))[2]!=ng*nt) 
      stop("gextra must return a matrix with vectors of length times*groups=",ng*nt)
    if(any(is.na(gextra(pgextra)))) stop("The extra parameters function returns NAs")
    if(dim(gextra(pgextra))[1]==1) {
      fcall<-paste(deparse(gextra)[1],collapse="")
      fbody<-paste(deparse(gextra)[-1],collapse="\n")
      ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",nt*ng,",byrow=TRUE)}"),collapse="")
      .gextra<-eval(parse(text=ftext))
      environment(.gextra)<-environment(gextra)
    } else
    if(dim(gextra(pgextra))[1]==n ) .gextra<-gextra else
    stop("gextra should return a matrix of length 1 or",n)
  }
                                        #
  if(is.null(gmixture) && is.null(pgmix)) { .gmixture<-function(...) 1 } else
  {
    if(dim(gmixture(pgmix))[2]!=ng) stop("gmixture must return a matrix with vectors of length groups=",ng)
    if(any(is.na(gmixture(pgmix)))) stop("The mixture function returns NAs")
    if(any(gmixture(pgmix)<0)||any(apply(gmixture(pgmix),1,sum))!=1) 
      stop("The mixture function components must be >=0 and sum to 1")
    if(dim(gmixture(pgmix))[1]==1 ) {
      fcall<-paste(deparse(gmixture)[1],collapse="")
      fbody<-paste(deparse(gmixture)[-1],collapse="\n")
      ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",ng,",byrow=TRUE)}"),collapse="")
      .gmixture<-eval(parse(text=ftext))
      environment(.gmixture)<-environment(gmixture)
    } else
    if(dim(gmixture(pgmix))[1]==n) .gmixture<-gmixture else
    stop("gmixture should return a matrix of length 1 or",n)
  }

  if (is.null(expected))
    {
      environment(.gmu)<-globalenv()
      .expected<- eval(expression(function(p) {gmu(p[1:npl])}),list(gmu=.gmu,npl=npl))
    } else
  {
    ptot <- c(pgmu,pgshape,pgextra) 
    if(length(dim(expected(ptot)))!=2||dim(expected(ptot))[2]!=ng*nt) 
      stop("expected must return a matrix with vectors of length times*groups=",ng*nt)
    if(any(is.na(expected(ptot)))) stop("The expected function returns NAs")
    if(dim(expected(ptot))[1]==1){
      fcall<-paste(deparse(expected)[1],collapse="")
      fbody<-paste(deparse(expected)[-1],collapse="\n")
      ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",ng*nt,",byrow=TRUE)}"),collapse="")
      .expected<-eval(parse(text=ftext))
      environment(.expected)<-environment(expected)
    } else
    if(dim(expected(ptot))[1]!=n)  stop("gmu should return a matrix of length 1 or",n)
    environment(.expected)<-globalenv()
  }
    
                                        #
                                        # check and scale weights if necessary  
                                        #
  wt<-as.vector(wt)
  if((length(wt)==1) && (wt==1)) wt <- rep(1,n)
  if(length(wt)!=n) stop("wt must be the same length as the other variables")
  if(min(wt)<0) stop("All weights must be non-negative")
  if((sum(wt)!=n) & scale.weight) {
    warning("\nThe weights have been rescaled to sum to the the sample size ",n,"\n")
  wt<-wt/mean(wt)
  }
                                        #
                                        #define the likelihood functions
                                        #
  tmp<-array(resp,dim=c(n,nt,ng))
  likelihood<-if(joint){
    function(p){
      cm<-array(.gmu(p[1:npl]),dim=c(n,nt,ng))
      sc<-array(.gshape(p[(npl+1):(npl+nps)]),dim=c(n,nt,ng))
      extra<-array(.gextra(p[(npl+nps+1):(npl+nps+npext)]),dim=c(n,nt,ng))
      .density(tmp,cm,sc,extra)}} else
  {function(p){
    cm<-array(.gmu(p[1:npl]),dim=c(n,nt,ng))
    sc<-array(.gshape(p[(npl+1):(npl+nps)]),dim=c(n,nt,ng))
    extra<-array(.gextra(p[(npl+nps+1):(npl+nps+npext)]),dim=c(n,nt,ng))
    apply(.density(tmp,cm,sc,extra),c(1,3),prod,na.rm=TRUE)}}
                                        #   
                                        #   
  loglike<-function(p)
    {
      mix<-.gmixture(p[(npl+nps+npext+1):(npl+nps+npext+npmix)])
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
                              iterlim=iterlim,print.level=print.level,...))[3]} else
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
  pmix<-.gmixture(z0$estimate[(npl+nps+npext+1):(npl+nps+npext+npmix)])
  post<-likelihood(z0$estimate)*pmix
  post<-post/apply(post,1,sum)
                                        #
                                        # compute fitted and obserserved values 
                                        #
  mpost<-apply(post,2,mean)
  fitted<-array(.expected(z0$estimate),dim=c(n,nt,ng))
  fitted.mean<-apply(aperm(fitted,c(1,3,2))*rep(post,nt),c(2,3),mean,na.rm=TRUE)/mpost
  dimnames(fitted.mean)<-list(paste("Group",1:ng,sep=""), ifelse(is.na(nchar(temp<-dimnames(resp)[[2]])[1:nt]) ,paste("Time",1:nt,sep=""),temp))
  observed.mean<-t(apply(array(apply(post,2,"*",resp),dim=c(n,nt,ng)),c(2,3),mean,na.rm=TRUE))/mpost
  dimnames(observed.mean)<-list(paste("Group",1:ng,sep=""), ifelse(is.na(nchar(temp<-dimnames(resp)[[2]])[1:nt]) ,paste("Time",1:nt,sep=""),temp))
                                        #
                                        # clean functions environment and
                                        # makes the weights callable
                                        #
  environment(.density)<-globalenv()
  environment(.gmu)<-globalenv()
  environment(.gshape)<-globalenv()
  environment(.gextra)<-globalenv()
  environment(.gmixture)<-globalenv()
  attributes(.gmu)<-attributes(gmu)
  attributes(.gshape)<-attributes(gshape)
  attributes(.gextra)<-attributes(gextra)
  attributes(.gmixture)<-attributes(gmixture)
  attributes(.density)<-attributes(density)
  wt<-match.call()$wt
  if(is.null(wt)) {wt<-call("rep",1,n)} else
  {
  if(scale.weight) {wt<-substitute(as.vector(a)/mean(as.vector(a)),list(a=wt)) } else
  wt<-substitute(as.vector(a),list(a=wt))
  }
                                        #
                                        # return a list of class moc
                                        #
  moc.out <- list(
                  call=call,
                  data=match.call()$data,
                  resp=match.call()$y,
                  density=.density,
                  joint=joint,
                  nsubject=n,
                  ntimes=nt,
                  nobs=sum(!is.na(resp)),
                  groups=ng,
                  npar=c(npl,nps,npext,npmix),
                  gmu=.gmu,
                  gshape=.gshape,
                  gextra=.gextra,
                  gmixture=.gmixture,
                  expected = .expected,
                  prior.weights=wt,
                  loglikelihood=-z0$minimum,
                  df=n*nt-np,
                  AIC=2*z0$minimum+2*np,
                  BIC=2*z0$minimum+np*log(sum(eval(wt))*nt),
                  coefficients=z0$estimate,
                  cov=cov,
                  fitted.mean=fitted.mean,
                  observed.mean=observed.mean,
                  iterations=z0$iterations,
                  code=z0$code)
  class(moc.out) <- "moc"
  return(moc.out) }

post<-function(object,...) UseMethod("post")

post.moc<-function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data)))
  {
    attach(eval(object$data),name=deparse(object$data))
    on.exit(eval(substitute(detach(a),list(a=deparse(object$data)))))
  }
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
    invisible(post/apply(post,1,sum))
  }
}

fitted.moc<-function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){ attach(eval(object$data),name=deparse(object$data));
                                   on.exit(eval(substitute(detach(a),list(a=deparse(object$data)))))}
  n<-object$nsubject
  dim1<-c(n,object$ntimes,object$groups)
  fit<-array(object$expected(object$coefficients),dim=dim1)
  dimnames(fit)<-list(NULL, ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:dim1[2]]) ,paste("Time",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
  invisible(fit)
}

residuals.moc<-function(object,...,type="deviance",post.weight=TRUE,within=FALSE)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){ attach(eval(object$data),name=deparse(object$data));
                                   on.exit(eval(substitute(detach(a),list(a=deparse(object$data)))))}
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
  wts<-eval(object$prior.weights)
  wpost<-post(object)
  if (within) wpost<-t(t(wpost)/apply(wpost,2,mean))
  response<-(y-fitted(object))
  switch(type,
         response= res<-response,
         deviance= res<- if(!object$joint) {sqrt(2*wts*(log(object$density(y,y,s,extra)/object$density(y,m,s,extra))))*sign(response)} else
         { tmp<-array(NA,dim1)
           for (i in 1:nt)
             {mi<-y; mi[,i,]<-m[,i,]
              tmp[,i,]<-sqrt(2*wts*(log(object$density(y,y,s,extra)/object$density(y,mi,s,extra))))*sign(response[,i,])
            } 
           tmp }
         )

  dimnames(res)<-list(NULL,ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:dim1[2]]) ,paste("Time",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
  if(post.weight) res<-res*array(wpost[,rep(1:ng,rep(nt,ng))],dim=dim1)
  class(res)<-c("residuals.moc","residuals")
  attr(res,"type")<-type
  attr(res,"post.weight")<-post.weight
  attr(res,"within")<-within
  invisible(res)
}


print.moc<-function(x,digits=5,...)
{
  object<-x
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){ attach(eval(object$data),name=deparse(object$data));
                                   on.exit(eval(substitute(detach(a),list(a=deparse(object$data)))))}

  cat("\n\t\t",object$groups,"Mixture of Curves\n\n\n")
  cat("Response is ",deparse(object$resp),"\n\n")
  if(object$joint) cat("joint ")
  cat("Density:\n");print(object$density)
  cat("\nLocation:\n");print(object$gmu )
  cat("\nExpectation:\n");print(object$expected )
  cat("\nShape:\n");print(object$gshape)
  cat("\nExtra parameter:\n");print(object$gextra)
  cat("\nMixture:\n");print(object$gmixture)
  cat("\n\n\t\t\tMaximum Likelihood Estimates\n\n")
  if(object$code>2) {cat("\n\nWARNING - CONVERGENCE NOT ACHIEVED - ERROR",object$code,"\n\n\n")}
  cat(formatC(c("",""," Standard","Wald Chi-Sq:"),width=13),"\n")
  cat(formatC(c("Parameter    ","Estimate"," Error"," Param = 0"," Prob > Wald"),width=13),"\n")
  object.se<-sqrt(diag(object$cov))
  if(object$npar[1]>0)
    {
      cat("\nLocation:\n\n")
      nl<-object$npar[1];ns<-object$npar[2];nextra<-object$npar[3];nm<-object$npar[4]
      param<-attr(object$gmu,"parameters")
      if (length(param) != nl) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[1:nl],se<-object.se[1:nl],w<-(est/se)^2,(1-pchisq(w,1))),nl,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),nl,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Location<-coeftable
    }
  if(object$npar[2]>0)
    {
      cat("\nShape:\n\n")
      param<-attr(object$gshape,"parameters")
      if (length(param) != ns) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[(nl+1):(nl+ns)],se<-object.se[(nl+1):(nl+ns)],w<-(est/se)^2,(1-pchisq(w,1))),ns,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),ns,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Shape<-coeftable
    }
  if(object$npar[3]>0)
    {
      cat("\nExtra Parameters:\n\n")
      param<-attr(object$gextra,"parameters")
      if (length(param) != nextra) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[(nl+ns+1):(nl+ns+nextra)],se<-object.se[(nl+ns+1):(nl+ns+nextra)],
                              w<-(est/se)^2,(1-pchisq(w,1))),nextra,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),nextra,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Extra<-coeftable
    }
  if(object$npar[4]>0)
    {
      cat("\nMixture:\n\n")
      param<-attr(object$gmixture,"parameters")
      if (length(param) != nm) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[(nl+ns+nextra+1):(nl+ns+nextra+nm)],se<-object.se[(nl+ns+nextra+1):(nl+ns+nextra+nm)],
                              w<-(est/se)^2,(1-pchisq(w,1))),nm,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),nm,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Mixture<-coeftable
      cat("\nMean Mixture Probabilities:\n")
      mmix<-apply(object$gmixture(est),2,mean)
      names(mmix)<-paste("Group",1:object$groups,sep="")
      print(mmix,digits=5)
      cat("\n")
      object$MixtProb<-coeftable
    }
  modelfit<-cbind("-2*logLikelihood"=-2*object$loglik, AIC=object$AIC,AIC(object,k="BIC"))
  dimnames(modelfit)[[1]]<-" "
  print(modelfit)
  object$ModelFit<-modelfit
  cat("\n\nMean predicted values:\n")
  print(object$fitted.mean,digits=digits)
  cat("\nMean observed values:\n")
  print(object$observed.mean,digits=digits)
  invisible(object)
}


AIC.moc <- function(object,...,k=2) 
{
  mlist<-  list(object,...)
  if(!all(sapply(mlist,inherits,"moc"))) stop("\n\tAll objects must be of class moc\n")
  ml<-length(mlist)
  cnames<-as.character(match.call()[-1])
  names(mlist)<- cnames[1:ml]
  nobslist<-sapply(mlist,function(tmp) tmp$nobs)
  nglist<-sapply(mlist,function(tmp) tmp$groups)
  if((k==0) && ((max(nglist)!=min(nglist)) || (max(nobslist)!=min(nobslist))))
  warning("log Likelihood should only be used to compare nested models on the same observations with the same number of mixture.\nTry using AIC or BIC or ICL-BIC instead.\n") else
  {if((k!="BIC") && (max(nobslist)!=min(nobslist)))
     warning("AIC like statistics should not be used to compare models with differing number of observations, use BIC or ICL-BIC instead.\n")}
  if(k=="BIC")
    val<-as.data.frame(t(sapply(mlist,function(tmp) {bic<-tmp$BIC
                                                     po<-post(tmp); entropy<--sum(ifelse(po==0,0,po*log(po)))
                                                     c(bic,entropy,bic+2*entropy,tmp$df)}))) else
  val<-as.data.frame(t(sapply(mlist,function(tmp)
                              c(-2*tmp$loglikelihood+k*sum(tmp$npar),tmp$df))))
  names(val)<-c(switch(as.character(k),"0"="-2*logLik","2"="AIC",BIC=c("BIC","Entropy","ICL-BIC"),"generalized AIC"),"Df")
  row.names(val)<-cnames[1:ml]                     
  val
}

logLik.moc <- function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  val <- object$loglikelihood
  attr(val,"df") <- object$df
  attr(val,"nobs")<-object$ntimes*object$nsubject
  class(val) <- "logLik"
  val
}


obsfit.moc<-function(object,along=NULL,FUN=function(x) x)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){ attach(eval(object$data),name=deparse(object$data));
                                   on.exit(eval(substitute(detach(a),list(a=deparse(object$data)))))}
  n<-object$nsubject
  nt<-object$ntimes
  ng<-object$groups
  post.obj<-post.moc(object)
  mpost<-by(post.obj,along,mean,na.rm=TRUE)
  fitted.mean<-by(FUN(object$expected(object$coef))*array(apply(post.obj,2,rep,nt),c(n,nt*ng)),
                 along,mean,na.rm=TRUE)
  observed.mean<-by(array(apply(post.obj,2,"*",FUN(as.matrix(eval(object$resp)))),c(n,nt*ng)),along,mean,na.rm=TRUE)
  nlist<-dim(mpost)
  for (i in 1:prod(nlist))
    {
      if(!is.null(fitted.mean[[i]]))
        {
          fitted.mean[[i]] <- t(array(fitted.mean[[i]],c(nt,ng)))
          fitted.mean[[i]]<-(fitted.mean[[i]]/mpost[[i]])
          dimnames(fitted.mean[[i]])<-list(paste("Group",1:ng,sep=""),
                                           ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:nt]) ,paste("Time",1:nt,sep=""),temp))
        }
      if(!is.null(observed.mean[[i]]))
        {
          observed.mean[[i]] <- t(array(observed.mean[[i]],c(nt,ng)))
          observed.mean[[i]]<-(observed.mean[[i]]/mpost[[i]])
          dimnames(observed.mean[[i]])<-list(paste("Group",1:ng,sep=""),
                                             ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:nt]) ,paste("Time",1:nt,sep=""),temp))
        }
    }

  val<-list("Mean Posterior Probabilities"=mpost,
            "Function"=substitute(FUN),
            "Mean function Expected Values"=fitted.mean,
            "Mean function Observed Values"=observed.mean)
  print(val)
  invisible(val)
}


plot.moc<-function(x,times=1:x$ntimes,main="",xlab="",ylab="",...)
{
  if(!inherits(x,"moc")) stop("Not an object of class moc")
  if(dim(as.matrix(times))[2]==1) {w<-times} else {w<-cbind(times,times)}
  matplot(w,t(rbind(x$observed.mean,x$fitted.mean)),type="n",main=main,xlab=xlab,ylab=ylab,...)
  matpoints(times,t(x$observed.mean),...)
  matlines(times,t(x$fitted.mean),...)
  invisible(list(times=times,fitted=x$fitted,observed=x$observed))
}

plot.residuals.moc<-function(x,against="Index",sunflower=FALSE,...)
  {
    if(sunflower) thisplot<-sunflowerplot else thisplot<-plot
    if(!inherits(x,"residuals.moc")) stop("Not an object of class residuals.moc")
    dim1<-dim(x)
    again<-against
    res.min<-min(x,na.rm=TRUE);res.max<-max(x,na.rm=TRUE)
    oldpar<-par(mfrow=c(dim1[3],1),oma=c(0,0,2,0));on.exit(par(oldpar))
    vname<-deparse(substitute(against))
    if(against=="Subject") against<-c(rep(1:dim1[1],dim1[2]))
    if(against=="Observation") against<-c(rep(1:dim1[2],rep(dim1[1],dim1[2])))
    if(against=="Index") against<-1:(dim1[1]*dim1[2])
    res.apply<-sapply(1:dim1[3],function(i){
      z<-na.omit(cbind(c(against),c(x[,,i])))
      dimnames(z)<-list(NULL,c(vname,attr(x,"type")))
      thisplot(z,main=paste("Group",i),...)
      if(again=="Index") {for(i in 1:(dim1[2]-1)) lines(rep(i,2)*dim1[1]+0.5,c(res.min,res.max))}
    })
     mtext(paste("MOC ",attr(x,"type")," residuals"),side=3,outer=TRUE,font=2)
     invisible(res.apply)
  }
    

  
".First.lib"<-
function(lib, pkg)
{
  cat("\n",readLines(system.file("DESCRIPTION",package="moc")),sep="\n","\n")
  if (!exists('AIC',envir=.GlobalEnv)){
  eval(expression(AIC<-function(object,...,k=2) UseMethod('AIC')),envir=.GlobalEnv)} 
}

















