#
#  moc : Function to fit general multivariate mixture models 
#         
#  Copyright (C) 2000-2003 Bernard Boulerice
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
#  moc(y, density=NULL, joint=FALSE, groups=1,
#               gmu=NULL, gshape=NULL, gextra=NULL, gmixture=inv.glogit, expected = NULL,
#               pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, check.length=TRUE,
#               scale.weight=FALSE, wt=1, data=NULL,
#               ndigit=10, gradtol=0.0001, steptol=gradtol, iterlim=100, print.level=2,...)
#
#  DESCRIPTION
#
# Function to fit general nonlinear multivariate mixture models
#
#
moc<- function(y, density=NULL, joint=FALSE, groups=1,
               gmu=NULL, gshape=NULL, gextra=NULL, gmixture=inv.glogit, expected = NULL,
               pgmu=NULL, pgshape=NULL, pgextra=NULL, pgmix=NULL, check.length=TRUE,
               scale.weight=FALSE, wt=1, data=NULL,
               ndigit=10, gradtol=0.0001, steptol=gradtol, iterlim=100, print.level=2,...)
{
  if (!is.null(data)) {attach(data,pos=2);on.exit(detach(data),add=TRUE)}
                                        #   thisEnv <- environment()
                                        #   for( i in names( data ) ) {
                                        #        assign( i, data[[i]], envir = thisEnv )
                                        #    }

  call <- sys.call()
  resp<-as.matrix(eval(substitute(y)))
  ndim<-dim(resp)
  n<-ndim[1]
  nt<-ndim[2]
  ng<-groups
  inv.glogit<- if(ng==1) {function(gmix) {cbind(1)}} else
  {function(gmix) {rbind(c(1,exp(gmix)))/(1+sum(exp(gmix)))}}
  if(groups > 1) attr(inv.glogit,"parameters")<-paste("  G",2:groups," vs G1",sep="")
                                        #
                                        # check density
                                        #
  if(!is.function(density)) stop("\ndensity must be a function\n")
  if(length(formals(density))>4) stop("\ndensity must not use more than 4 arguments\n")
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

  if(!is.list(gmu) & (length(gmu) != ng) & !is.null(gmu))
    stop(paste("\ngmu must be a list of functions of length ",ng))
  if(!is.list(gshape) & (length(gshape) != ng) & !is.null(gshape))
    stop(paste("\ngshape must be a list of functions of length ",ng))
  if(!is.list(gextra) & (length(gextra) != ng) & !is.null(gextra))
    stop(paste("\ngextra must be a list of functions of length ",ng))
  .gmu<-list()
  .gshape<-list()
  .gextra<-list()
  .expected<-list()  
  for( ig in 1:ng) {
    if(length(dim(gmu[[ig]](pgmu)))!=2 || dim(gmu[[ig]](pgmu))[2]!=nt) 
      stop(paste("\ngmu in group",ig,"must return a matrix with vectors of length nvar = ",nt))
    if(any(is.na(gmu[[ig]](pgmu)))) stop(paste("\nThe gmu function returns NAs in group ",ig))
    if(dim(gmu[[ig]](pgmu))[1]==1) {
      fcall<-paste(deparse(gmu[[ig]])[1],collapse="")
      fbody<-paste(deparse(gmu[[ig]])[-1],collapse="\n")
      ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",nt,",byrow=TRUE)}"),collapse="")
      .gmu[[ig]]<-eval(parse(text=ftext))
      environment(.gmu[[ig]])<-environment(gmu[[ig]])
    } else
    if(dim(gmu[[ig]](pgmu))[1]==n) .gmu[[ig]]<-gmu[[ig]] else
    stop(paste("\ngmu in group",ig,"should return a matrix of length 1 or ",n))
                                        #
    if(is.null(gshape) && is.null(pgshape)) .gshape[[ig]]<- function(...) 1 else
    {
      if(length(dim(gshape[[ig]](pgshape)))!=2||dim(gshape[[ig]](pgshape))[2]!=nt) 
        stop(paste("\ngshape in group",ig,"must return a matrix with vectors of length nvar = ",nt))
      if(any(is.na(gshape[[ig]](pgshape))))
        stop(paste("\nThe shape function returns NAs in group",ig))
      if(dim(gshape[[ig]](pgshape))[1]==1) {
        fcall<-paste(deparse(gshape[[ig]])[1],collapse="")
        fbody<-paste(deparse(gshape[[ig]])[-1],collapse="\n")
        ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",nt,",byrow=TRUE)}"),collapse="")
        .gshape[[ig]]<-eval(parse(text=ftext))
        environment(.gshape[[ig]])<-environment(gshape[[ig]])
      } else
      if(dim(gshape[[ig]](pgshape))[1]==n ) .gshape[[ig]]<-gshape[[ig]] else
      stop("\ngshape in group",ig,"should return a matrix of length 1 or",n)
    }
                                        #
    if(is.null(gextra) && is.null(pgextra)) .gextra[[ig]]<- function(...) 1 else
    {
      if(length(dim(gextra[[ig]](pgextra)))!=2||
         (((lgext<-dim(gextra[[ig]](pgextra))[2])!=nt)&check.length)) 
        stop(paste("\ngextra in group",ig,"must return a matrix with vectors of length nvar =",nt))
      if(any(is.na(gextra[[ig]](pgextra))))
        stop(paste("\nThe extra parameters function returns NAs in group",ig))
      if(dim(gextra[[ig]](pgextra))[1]==1) {
        fcall<-paste(deparse(gextra[[ig]])[1],collapse="")
        fbody<-paste(deparse(gextra[[ig]])[-1],collapse="\n")
        ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",ifelse(check.length,nt,lgext),
                       ",byrow=TRUE)}"),collapse="")
        .gextra[[ig]]<-eval(parse(text=ftext))
        environment(.gextra[[ig]])<-environment(gextra[[ig]])
      } else
      if(dim(gextra[[ig]](pgextra))[1]==n ) .gextra[[ig]]<-gextra[[ig]] else
      stop(paste("\ngextra in group",ig,"should return a matrix of length 1 or",n))
    }
    if (is.null(expected))
      {
        environment(.gmu[[ig]])<-globalenv()
        .expected[[ig]]<- eval(expression(function(p) {gmu[[k]](p[1:npl])}),
                               list(k=ig,gmu=.gmu,npl=npl))
      } else
    {
      ptot <- c(pgmu,pgshape,pgextra) 
      if(length(dim(expected[[ig]](ptot)))!=2 || dim(expected[[ig]](ptot))[2]!=nt) 
        stop(paste("\nexpected in group",ig,"must return a matrix with vectors of length nvar =",nt))
      if(any(is.na(expected[[ig]](ptot)))) stop(paste("\nThe expected function returns NAs",ig))
      if(dim(expected[[ig]](ptot))[1]==1){
        fcall<-paste(deparse(expected[[ig]])[1],collapse="")
        fbody<-paste(deparse(expected[[ig]])[-1],collapse="\n")
        ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",nt,",byrow=TRUE)}"),collapse="")
        .expected[[ig]]<-eval(parse(text=ftext))
        environment(.expected[[ig]])<-environment(expected[[ig]])
      } else
      if(dim(expected[[ig]](ptot))[1]==n)  .expected[[ig]]=expected[[ig]] else
      stop(paste("\nexpected for group",ig,"should return a matrix of length 1 or",n))
      environment(.expected[[ig]])<-globalenv()
    }
 
  }
                                        #
  if(is.null(gmixture) && is.null(pgmix)) { .gmixture<-function(...) 1 } else
  {
    if(dim(gmixture(pgmix))[2]!=ng)
      stop(paste("\ngmixture must return a matrix with vectors of length groups=",ng))
    if(any(is.na(gmixture(pgmix)))) stop("\nThe mixture function returns NAs")
    if(any(gmixture(pgmix)<0)||any(apply(gmixture(pgmix),1,sum))!=1) 
      stop("\nThe mixture function components must be >=0 and sum to 1")
    if(dim(gmixture(pgmix))[1]==1 ) {
      fcall<-paste(deparse(gmixture)[1],collapse="")
      fbody<-paste(deparse(gmixture)[-1],collapse="\n")
      ftext<-paste(c(fcall,"{matrix(",fbody,",",n,",",ng,",byrow=TRUE)}"),collapse="")
      .gmixture<-eval(parse(text=ftext))
      environment(.gmixture)<-environment(gmixture)
    } else
    if(dim(gmixture(pgmix))[1]==n) .gmixture<-gmixture else
    stop(paste("\ngmixture should return a matrix of length 1 or",n))
  }
                                        #
                                        # check and scale weights if necessary  
                                        #
  wt<-as.vector(wt)
  if((length(wt)==1) && (wt==1)) wt <- rep(1,n)
  if(length(wt)!=n) stop("\nwt must be the same length as the other variables")
  if(min(wt)<0) stop("\nAll weights must be non-negative")
  if(any(is.na(wt))) stop("\n weights contain NA\n")
  if((sum(wt)!=n) & scale.weight) {
    warning("\nThe weights have been rescaled to sum to the the sample size ",n,"\n")
    wt<-wt/mean(wt)
  }
                                        #
                                        #define the likelihood functions
                                        #

  loglike<-if(joint){
    function(p)
      {
        parm<-split(p,rep(c("mu","shape","extra","mix"),c(npl,nps,npext,npmix)))
        dens<- sapply(1:ng,
                      function(ind) .density(resp,.gmu[[ind]](parm$mu),
                                             .gshape[[ind]](parm$shape),.gextra[[ind]](parm$extra)))
        sum(-wt*log(apply(dens*.gmixture(parm$mix),1,sum)))
      }}else
  {
    function(p)
      {
        parm<-split(p,rep(c("mu","shape","extra","mix"),c(npl,nps,npext,npmix)))
        dens<- sapply(1:ng,function(ind)
                      apply(.density(resp,.gmu[[ind]](parm$mu),.gshape[[ind]](parm$shape),
                                     .gextra[[ind]](parm$extra)),1,prod,na.rm=TRUE))
        sum(-wt*log(apply(dens*.gmixture(parm$mix),1,sum)))
      }}
                                        #
                                        # check that the likelihood returns an appropriate value and minimize
                                        #
  p<-c(pgmu,pgshape,pgextra,pgmix)
  if(is.na(loglike(p)))
    stop("\nLikelihood returns NAs: probably invalid initial values")
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
    cov <- matrix(NA,ncol=np,nrow=np)
  }
  se <- sqrt(diag(cov))
                                        #
                                        # compute mixture and posterior probabilities
                                        #
  parm<-split(z0$estimate,rep(c("mu","shape","extra","mix"),c(npl,nps,npext,npmix)))
  pmix<-.gmixture(parm$mix)
  if(joint){
    post1<-sapply(1:ng,function(ind)
                 .density(resp,.gmu[[ind]](parm$mu),.gshape[[ind]](parm$shape),
                          .gextra[[ind]](parm$extra)))*pmix
  } else
  { post1<- sapply(1:ng,function(ind)
                  apply(.density(resp,.gmu[[ind]](parm$mu),.gshape[[ind]](parm$shape),
                                 .gextra[[ind]](parm$extra)),1,prod,na.rm=TRUE))*pmix }
  post1<-post1/apply(post1,1,sum)*wt
                                        #
                                        # compute fitted and obserserved values 
                                        #
  mpost<-apply(post1,2,mean)
  fitted.mean<-matrix(t(sapply(1:ng,function(ind)
                        apply(.expected[[ind]](z0$estimate)*post1[,ind],2,mean,na.rm=TRUE)))/mpost,ng,nt)
  mpost<-sapply(1:ng,function(ind) apply(post1[,ind]*ifelse(is.na(resp),NA,1),2,mean,na.rm=TRUE))
  observed.mean<-matrix(t(sapply(1:ng,function(ind) apply(post1[,ind]*resp,2,mean,na.rm=TRUE))/mpost),ng,nt)
  dimnames(fitted.mean)<-list(paste("Group",1:ng,sep=""),
                              ifelse(is.na(nchar(temp<-dimnames(resp)[[2]])[1:nt]) ,
                                     paste("Time",1:nt,sep=""),temp))
  dimnames(observed.mean)<-list(paste("Group",1:ng,sep=""),
                                ifelse(is.na(nchar(temp<-dimnames(resp)[[2]])[1:nt]) ,
                                       paste("Time",1:nt,sep=""),temp))
                                        #
                                        # clean functions environment and
                                        # makes the weights callable
                                        #
  environment(.density)<-globalenv()
  for(ig in 1:ng)
    {
      environment(.gmu[[ig]])<-globalenv()
      environment(.gshape[[ig]])<-globalenv()
      environment(.gextra[[ig]])<-globalenv()
    }
  environment(.gmixture)<-globalenv()
  attributes(.gmu)<-attributes(gmu)
  attributes(.gshape)<-attributes(gshape)
  attributes(.gextra)<-attributes(gextra)
  attributes(.gmixture)<-attributes(gmixture)
  attributes(.density)<-attributes(density)
  gname<-paste("Group",1:ng,sep="")
  names(.gmu)<-gname
  names(.gshape)<-gname
  names(.gextra)<-gname
  names(.expected)<-gname

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
                  resp=substitute(as.matrix(txt),list(txt=match.call()$y)),
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
                  hessian=z0$hessian,
                  fitted.mean=fitted.mean,
                  observed.mean=observed.mean,
                  iterations=z0$iterations,
                  code=z0$code,
                  execution.time=dt)
  class(moc.out) <- "moc"
  return(moc.out) }


post<-function(object,...) UseMethod("post")

post.moc<-function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data)))
    {
      attach(eval(object$data),name=deparse(object$data),pos=2)
      on.exit(detach(2),add=TRUE)
    }
  parm<-split(object$coefficients,rep(c("mu","shape","extra","mix"),object$npar))
  mix<-object$gmixture(parm$mix)
  if(mix==1) cbind(1) else
  {
    if(object$joint) {
      post1<-sapply(1:object$groups,function(ind)
                   object$density(as.matrix(eval(object$resp)),object$gmu[[ind]](parm$mu),
                                  object$gshape[[ind]](parm$shape),object$gextra[[ind]](parm$extra)))*mix
    } else
    {post1<-sapply(1:object$groups,function(ind)
                  apply(object$density(as.matrix(eval(object$resp)),object$gmu[[ind]](parm$mu),
                                       object$gshape[[ind]](parm$shape),object$gextra[[ind]](parm$extra)),
                        1,prod,na.rm=TRUE))*mix}
    dimnames(post1)<-list(NULL,paste("Group",1:object$groups,sep=""))
    invisible(post1/apply(post1,1,sum))
  }
}

fitted.moc<-function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){ attach(eval(object$data),name=deparse(object$data),pos=2);
                                   on.exit(detach(2),add=TRUE)}
  n<-object$nsubject
  dim1<-c(n,object$ntimes,object$groups)
  fit<-array(NA,dim=dim1)
  for(ig in 1:object$groups) fit[,,ig]<-object$expected[[ig]](object$coefficients)
  dimnames(fit)<-list(NULL,
                      ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:dim1[2]]),
                             paste("Time",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
  invisible(fit)
}

residuals.moc<-function(object,...,type="deviance",post.weight=TRUE,within=FALSE)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){
    attach(eval(object$data),name=deparse(object$data),pos=2)
    on.exit(detach(2),add=TRUE)}
  choices<-c("deviance","response")
  type<-match.arg(type,choices)
  n<-object$nsubject
  nt<-object$ntimes
  ng<-object$groups
  dim1<-c(n,nt,ng)
  parm<-split(object$coefficients,rep(c("mu","shape","extra","mix"),object$npar))
  wts<-eval(object$prior.weights)
  wpost<-post.moc(object)
  if (within) wpost<-t(t(wpost)/apply(wpost,2,mean))
  y<-as.matrix(eval(object$resp))
  res<-array(y,dim=dim1)-fitted(object)
  for(ig in 1:object$groups)
    {
      m<-object$gmu[[ig]](parm$mu)
      s<-object$gshape[[ig]](parm$shape)
      extra<-object$gextra[[ig]](parm$extra)
      switch(type,
             response= res[,,ig]<-res[,,ig],
             deviance= res[,,ig]<-
             sqrt(2*wts*(log(object$density(y,y,s,extra)/object$density(y,m,s,extra))))*
             sign(res[,,ig])
             )
    }

  dimnames(res)<-list(NULL,
                      ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:dim1[2]]),
                             paste("Time",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
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
  if(!is.null(eval(object$data))){
    attach(eval(object$data),name=deparse(object$data),pos=2)
    on.exit(detach(2),add=TRUE)}

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
  nl<-object$npar[1];ns<-object$npar[2];nextra<-object$npar[3];nm<-object$npar[4]
  if(object$npar[1]>0)
    {
      cat("\nLocation:\n\n")
      param<-attr(object$gmu,"parameters")
      if (length(param) != nl) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[1:nl],se<-object.se[1:nl],w<-(est/se)^2,
                              (1-pchisq(w,1))),nl,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),nl,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Location<-coeftable
    }
  if(object$npar[2]>0)
    {
      cat("\nShape:\n\n")
      param<-attr(object$gshape,"parameters")
      if (length(param) != ns) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[(nl+1):(nl+ns)],se<-object.se[(nl+1):(nl+ns)],
                              w<-(est/se)^2,(1-pchisq(w,1))),ns,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),ns,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Shape<-coeftable
    }
  if(object$npar[3]>0)
    {
      cat("\nExtra Parameters:\n\n")
      param<-attr(object$gextra,"parameters")
      if (length(param) != nextra) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[(nl+ns+1):(nl+ns+nextra)],
                              se<-object.se[(nl+ns+1):(nl+ns+nextra)],w<-(est/se)^2,
                              (1-pchisq(w,1))),nextra,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),nextra,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Extra<-coeftable
    }
  if(object$npar[4]>0)
    {
      cat("\nMixture:\n\n")
      param<-attr(object$gmixture,"parameters")
      if (length(param) != nm) param<-"             "
      coeftable<-matrix(cbind(est<-object$coef[(nl+ns+nextra+1):(nl+ns+nextra+nm)],
                              se<-object.se[(nl+ns+nextra+1):(nl+ns+nextra+nm)],
                              w<-(est/se)^2,(1-pchisq(w,1))),nm,4)
      coeftable<-matrix(apply(coeftable,2,formatC,digits=digits,width=13),nm,4)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Mixture<-coeftable
      cat("\nMean Mixture Probabilities:\n")
      mmix<-apply(object$gmixture(est),2,mean)
      names(mmix)<-paste("Group",1:object$groups,sep="")
      print(mmix,digits=5)
      object$MixtProb<-coeftable
    }
  modelfit<-cbind("-2*logLikelihood"=-2*object$loglik, AIC=object$AIC,AIC(object,k="BIC"))
  dimnames(modelfit)[[1]]<-" "
  cat("\n")
  print(modelfit)
  object$ModelFit<-modelfit
  coeftable<-apply(post.moc(object),2,mean,na.rm=TRUE)
  object$PostMixtProb<-coeftable
  if(object$groups>1) {cat("\n\nMean Posterior Mixture Probabilities:\n")
  print(coeftable,digits=5)}
  cat("\n")
  cat("\nPosterior mean predicted values:\n")
  print(object$fitted.mean,digits=digits)
  cat("\nPosterior mean observed values:\n")
  print(object$observed.mean,digits=digits)
  cat("\n")
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
    val<-as.data.frame(t(sapply(mlist,function(tmp){
      attach(eval(tmp$data),name=deparse(tmp$data),pos=2)
      on.exit(detach(2),add=TRUE)
      bic<-tmp$BIC
      po<-post.moc(tmp); entropy<--sum(ifelse(po==0,0,eval(tmp$prior.weight)*po*log(po)))
      c(bic,entropy,bic+2*entropy,tmp$df)}))) else
  val<-as.data.frame(t(sapply(mlist,function(tmp)
                              c(-2*tmp$loglikelihood+k*sum(tmp$npar),tmp$df))))
  names(val)<-c(switch(as.character(k),"0"="-2*logLik","2"="AIC",BIC=c("BIC","Entropy","ICL-BIC"),
                       "generalized AIC"),"Df")
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
  if(!is.null(eval(object$data))){
    attach(eval(object$data),name=deparse(object$data),pos=2)
    on.exit(detach(2),add=TRUE)}
  n<-object$nsubject
  if(is.null(along)) along<-rep(1,n)
  nt<-object$ntimes
  ng<-object$groups
  wts<-eval(object$prior.weight)
  post.obj<-post.moc(object)*wts
  mpost<-by(post.obj,along,mean,na.rm=TRUE)
  wts<-by(wts,along,mean)
  tmp<-matrix(NA,n,nt*ng)
  for (ig in 1:object$groups) tmp[,(1:nt)+(ig-1)*nt]<-FUN(object$expected[[ig]](object$coef))
  fitted.mean<-by(tmp*array(apply(post.obj,2,rep,nt),c(n,nt*ng)), along,mean,na.rm=TRUE)
  mpost.fitted<-by(ifelse(is.na(tmp),NA,1)*array(apply(post.obj,2,rep,nt),c(n,nt*ng)),
                   along,mean,na.rm=TRUE)
  tmp<-FUN(array(eval(object$resp),c(n,nt*ng)))
  observed.mean<-by(tmp*array(apply(post.obj,2,rep,nt),c(n,nt*ng)),along,mean,na.rm=TRUE)
  mpost.observed<-by(ifelse(is.na(tmp),NA,1)*array(apply(post.obj,2,rep,nt),c(n,nt*ng)),
                     along,mean,na.rm=TRUE)
  nlist<-dim(mpost)
  for (i in 1:prod(nlist))
    {
      if(!is.null(fitted.mean[[i]])){
        fitted.mean[[i]] <- t(array(fitted.mean[[i]]/mpost.fitted[[i]],c(nt,ng)))
        dimnames(fitted.mean[[i]])<-
          list(paste("Group",1:ng,sep=""),
               ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:nt]) ,
                      paste("Time",1:nt,sep=""),temp))
      }
      if(!is.null(observed.mean[[i]]))
        {
          observed.mean[[i]] <- t(array(observed.mean[[i]]/mpost.observed[[i]],c(nt,ng)))
          dimnames(observed.mean[[i]])<-
            list(paste("Group",1:ng,sep=""),
                 ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:nt]) ,
                        paste("Time",1:nt,sep=""),temp))
        }
    }

  val<-list("Mean Posterior Probabilities"=mpost,
            "Function"=substitute(FUN),
            "Mean function Expected Values"=fitted.mean,
            "Mean function Observed Values"=observed.mean)
  print(val)
  invisible(val)
}


plot.moc<-function(x,against=1:x$ntimes,main="",xlab="",ylab="",prob.legend=TRUE,scale=FALSE,group.colors=rainbow(x$groups),...)
{
  if(!inherits(x,"moc")) stop("Not an object of class moc")
  if(!is.null(eval(x$data))){
      attach(eval(x$data),name=deparse(x$data),pos=2)
      on.exit(detach(2),add=TRUE)}
  if(prob.legend) {
    oldpar<-par("oma"=c(0,0,0,8),"las"=1)
    legend.text<-paste("Group",1:x$groups,"=",
                       formatC(apply(x$gmixture(x$coef[(sum(x$npar[1:3])+1):sum(x$npar[1:4])]),
                                     2,mean),digit=4))
    on.exit(par(oldpar),add=TRUE)
  }
  if(dim(as.matrix(against))[2]==1) {w<-cbind(against)} else {w<-cbind(against,against)}
  if(scale) {
    center<-apply(eval(x$resp),2,mean,na.rm=TRUE)
    scale<-apply(eval(x$resp),2,sd,na.rm=TRUE)
  } else {center<-rep(0,x$ntimes);scale<-rep(1,x$ntimes)}
  matplot(w,cbind((t(x$observed.mean)-center)/scale,(t(x$fitted.mean)-center)/scale),
          type="n",main=main,xlab=xlab,ylab=ylab,col=group.colors,...)
  matpoints(against,(t(x$observed.mean)-center)/scale,col=group.colors,...)
  matlines(against,(t(x$fitted.mean)-center)/scale,type="o",pch=20,cex=.75,col=group.colors,...)
  if(prob.legend) { mtext(legend.text,side=4,outer=TRUE,at=((1:x$groups)+3)/(x$groups+6),col=group.colors) }
  invisible(list(against=against,fitted=x$fitted,observed=x$observed,center=center,scale=scale))
}

profilesplot<-function(x,...) UseMethod("profilesplot")

profilesplot.moc<-function(x,against=1:x$ntimes,main="",xlab="",ylab="",col.legend=TRUE,scale=FALSE,group.colors=rainbow(x$groups),...)
{
  if(!inherits(x,"moc")) stop("Not an object of class moc")
  if(!is.null(eval(x$data)))
    {
      attach(eval(x$data),name=deparse(x$data),pos=2)
      on.exit(detach(2),add=TRUE)
    }
  group.colors<-col2rgb(group.colors)/256
  if(col.legend) {
    oldpar<-par("oma"=c(0,0,0,6),"las"=1);on.exit(par(oldpar),add=TRUE)
    legend.text<-paste("Group",1:x$groups)
  }
  if(scale) {
    center<-apply(eval(x$resp),2,mean,na.rm=TRUE)
    scale<-apply(eval(x$resp),2,sd,na.rm=TRUE)
  } else {center<-rep(0,x$ntimes);scale<-rep(1,x$ntimes)}
  group.rgb<-t(apply(post.moc(x),1,function(y) (group.colors%*%y)))
  matplot((against),(t(eval(x$resp))-center)/scale,
          type="o",pch=20,cex=0.75,main=main,xlab=xlab,ylab=ylab,
          col=rgb(group.rgb[,1],group.rgb[,2],group.rgb[,3]))
  if(col.legend) { mtext(legend.text,side=4,outer=TRUE,at=((1:x$groups)+3)/(x$groups+6),
                         col=rgb(group.colors[1,],group.colors[2,],group.colors[3,])) }
}

plot.residuals.moc<-function(x,against="Index",groups=1:dim(x)[3],sunflower=FALSE,...)
{
  if(sunflower) thisplot<-sunflowerplot else thisplot<-plot
  if(!inherits(x,"residuals.moc"))
    stop("Not an object of class residuals.moc !")
  if (!all(groups %in% (1:dim(x)[3])))
    stop("You requested residuals for non-existing groups !")
  dim1<-dim(x)
  again<-against
  oldpar<-par(mfrow=c(length(groups),1),oma=c(0,0,2,0));on.exit(par(oldpar),add=TRUE)
  vname<-deparse(substitute(against))
  if(against=="Subject") against<-c(rep(1:dim1[1],dim1[2]))
  if(against=="Observation") against<-c(rep(1:dim1[2],rep(dim1[1],dim1[2])))
  if(against=="Index") against<-1:(dim1[1]*dim1[2])
  res.apply<-sapply(groups,function(i){
    z<-na.omit(cbind(c(against),c(x[,,i])))
    dimnames(z)<-list(NULL,c(vname,attr(x,"type")))
    thisplot(z,main=paste("Group",i),...)
    if(again=="Index")  abline(v=(1:dim1[2])*dim1[1]+0.5)
  })
  mtext(paste("MOC ",attr(x,"type")," residuals"),side=3,outer=TRUE,font=2)
  invisible(res.apply)
}


# Generalized logit and inverse logit with respect to a reference group

inv.glogit<-function(gmix,ref=1) {rbind(append(exp(gmix),1,ref-1))/(1+sum(exp(gmix)))}
glogit<-function(p,ref=1) {log(p[-ref]/p[ref])}

# Mix group colors according to posterior mixture probabilities

mix.colors.moc<-function(object,group.colors=rainbow(object$groups))
{
  if(!inherits(object,"moc")) stop("Not an object of class moc")
    if(length(group.colors)!=object$groups)
        stop("Oooups ! Wrong number of colors: should be the same as number of groups.")
	  group.rgb<-t(apply(post(object),1,function(y) ((col2rgb(group.colors)/256)%*%y)))
	    invisible(rgb(group.rgb[,1],group.rgb[,2],group.rgb[,3]))
	    }
	    

".First.lib"<-
function(lib, pkg)
{
  cat("\n",readLines(system.file("DESCRIPTION",package="moc")),sep="\n","\n")
  cat("\nSupplementary utility functions to help combine MOC models can be found in\n",
      system.file("Utils","combine.moc.R",package="moc"),
      ".\nYou will need to source that file to use them. See the file for simple documentation.\n\n")
}
    
