.packageName <- "moc"
#
#  moc : Library to fit general multivariate mixture models 
#         
#  Copyright (C) 2000-2004 Bernard Boulerice
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
  if (!is.null(data)) {attach(data,pos=2);on.exit(detach(2),add=TRUE)}
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
                                        # check gmu
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
                                        # check gshape
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
                                        # check gextra
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
                                        #
                                        # check expected
                                        #
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
                                        # check the returned values of the mixture function
                                        #
   if(is.null(gmixture) && is.null(pgmix)) { .gmixture<-function(...) 1 } else
  {
    if(dim(gmixture(pgmix))[2]!=ng)
      stop(paste("\ngmixture must return a matrix with vectors of length groups=",ng))
    if(any(is.na(gmixture(pgmix)))) stop("\nThe mixture function returns NAs")
    if(any(gmixture(pgmix)<0)||any(abs(apply(gmixture(pgmix),1,sum)-1)>.Machine$double.eps^0.5))
      warning("\nThe mixture function components probabilities must be >=0 and sum to 1")
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
                                        # define the likelihood functions
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
  post1<-post1/apply(post1,1,sum)
  dimnames(post1)<-list(NULL,paste("Group",1:ng,sep=""))
                                        #
                                        # compute posterior prob., fitted and obserserved values 
                                        #
  mpost<-apply(post1*wt,2,mean)
  fitted.mean<-matrix(t(sapply(1:ng,function(ind)
                        apply(.expected[[ind]](z0$estimate)*post1[,ind]*wt,2,mean,na.rm=TRUE)))/mpost,ng,nt)
  mpost<-sapply(1:ng,function(ind) apply(post1[,ind]*wt*ifelse(is.na(resp),NA,1),2,mean,na.rm=TRUE))
  observed.mean<-matrix(t(sapply(1:ng,function(ind) apply(post1[,ind]*wt*resp,2,mean,na.rm=TRUE))/mpost),ng,nt)
  dimnames(fitted.mean)<-list(paste("Group",1:ng,sep=""),
                              ifelse(is.na(nchar(temp<-dimnames(resp)[[2]])[1:nt]) ,
                                     paste("V",1:nt,sep=""),temp))
  dimnames(observed.mean)<-list(paste("Group",1:ng,sep=""),
                                ifelse(is.na(nchar(temp<-dimnames(resp)[[2]])[1:nt]) ,
                                       paste("V",1:nt,sep=""),temp))
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
  attributes(.expected)<-attributes(expected)
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
     if(any(.gmixture(parm$mix)<0)||any(abs(apply(rbind(.gmixture(parm$mix)),1,sum)-1)>.Machine$double.eps^0.5)) 
      warning("\nThe final mixture probablities are not all >=0 or don't sum to 1.\n")
 
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
                  nvar=nt,
                  nobs=sum(!is.na(resp)),
                  groups=ng,
                  npar=c(npl,nps,npext,npmix),
                  gmu=.gmu,
                  gshape=.gshape,
                  gextra=.gextra,
                  gmixture=.gmixture,
                  expected = .expected,
                  prior.weights=wt,
                  post.prob=post1,
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
  if(!is.null(object$post.prob)) {return(invisible(object$post.prob))} else
  {
  if(!is.null(eval(object$data)))
    {
      attach(eval(object$data),name=deparse(object$data),pos=2)
      on.exit(detach(2),add=TRUE)
    }
  parm<-split(object$coefficients,rep(c("mu","shape","extra","mix"),object$npar))
  mix<-object$gmixture(parm$mix)
  if(object$groups==1) cbind(1) else
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
    dimnames(post1)<-list(paste("Group",1:object$groups,sep=""))
    post1 <- post1/apply(post1,1,sum)
    obj.name <- paste(substitute(object))
    object$post.prob <- post1
    eval(substitute(assign(obj.name,object,pos=sys.frame())))
    attr(post1,"moc.name") <- deparse(substitute(object))
    invisible(post1)
  }
}
}

fitted.moc<-function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){ attach(eval(object$data),name=deparse(object$data),pos=2);
                                   on.exit(detach(2),add=TRUE)}
  n<-object$nsubject
  dim1<-c(n,object$nvar,object$groups)
  fit<-array(NA,dim=dim1)
  for(ig in 1:object$groups) fit[,,ig]<-object$expected[[ig]](object$coefficients)
  dimnames(fit)<-list(NULL,
                      ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:dim1[2]]),
                             paste("V",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
  attr(fit,"moc.name") <- deparse(substitute(object))
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
  nt<-object$nvar
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
                             paste("V",1:dim1[2],sep=""),temp),paste("Group",1:dim1[3],sep=""))
  if(post.weight) res<-res*array(wpost[,rep(1:ng,rep(nt,ng))],dim=dim1)
  class(res)<-c("residuals.moc","residuals")
  attr(res,"type")<-type
  attr(res,"post.weight")<-post.weight
  attr(res,"within")<-within
  attr(res,"moc.name") <- deparse(substitute(object))
  invisible(res)
}


print.moc<-function(x,digits=5,...)
{
  object<-x
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){
    attach(eval(object$data),name=deparse(object$data),pos=2)
    on.exit(detach(2),add=TRUE)}

  cat("\n\t\t",object$groups,"mixtures MOC model\n\n\n")
  cat("Response: ",deparse(object$resp),"\n\n")
  if(object$joint) cat("joint ")
  cat("Density: ");cat(deparse(object$call$density),"\n");print(object$density)
  if(!is.null(object$call$gmu)) {cat("\nLocation: ");cat(deparse(object$call$gmu),"\n");print.listof(object$gmu )}
  if(!is.null(object$call$expected)) {cat("\nExpectation: ");cat(deparse(object$call$expected),"\n");print.listof(object$expected )}
  if(!is.null(object$call$gshape)) {cat("\nShape: ");cat(deparse(object$call$gshape),"\n");print.listof(object$gshape)}
  if(!is.null(object$call$gextra)) {cat("\nExtra parameter: ");cat(deparse(object$call$gextra),"\n");print.listof(object$gextra)}
  cat("\nMixture: ")
  if(!is.null(object$call$gmixture)) {cat(deparse(object$call$gmixture),"\n")} else {cat("inv.glogit\n")}
  print(object$gmixture)  
  cat("\n\n\t\t\tMaximum Likelihood Estimates\n\n")
  cat(formatC(c("",""," Standard","Wald Chi-Sq:"),width=13),"\n")
  cat(formatC(c("Parameter    ","Estimate"," Error"," Param = 0"," Prob > Wald"),width=13),"\n")
  object.coef<-split(object$coef,rep(c("mu","shape","extra","mix"),object$npar))
  object.se<-split(sqrt(diag(object$cov)),rep(c("mu","shape","extra","mix"),object$npar))
  nl<-object$npar[1];ns<-object$npar[2];nextra<-object$npar[3];nm<-object$npar[4]
  if(nl>0)
    {
      cat("\nLocation:\n\n")
      param<-attr(object$gmu,"parameters")
      if (length(param) != nl) param<-"             "
      coeftable<-matrix(cbind(est<-object.coef$mu,se<-object.se$mu,w<-(est/se)^2,
                              (1-pchisq(w,1))),nl,4)
      coeftable<-formatC(matrix(apply(coeftable,2,format,digits=digits,width=13),nl,4),width=13)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Location<-coeftable
    }
  if(ns>0)
    {
      cat("\nShape:\n\n")
      param<-attr(object$gshape,"parameters")
      if (length(param) != ns) param<-"             "
      coeftable<-matrix(cbind(est<-object.coef$shape,se<-object.se$shape,
                              w<-(est/se)^2,(1-pchisq(w,1))),ns,4)
      coeftable<-formatC(matrix(apply(coeftable,2,format,digits=digits,width=13),ns,4),width=13)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Shape<-coeftable
    }
  if(nextra>0)
    {
      cat("\nExtra Parameters:\n\n")
      param<-attr(object$gextra,"parameters")
      if (length(param) != nextra) param<-"             "
      coeftable<-matrix(cbind(est<-object.coef$extra,
                              se<-object.se$extra,w<-(est/se)^2,
                              (1-pchisq(w,1))),nextra,4)
      coeftable<-formatC(matrix(apply(coeftable,2,format,digits=digits,width=13),nextra,4),width=13)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Extra<-coeftable
    }
  if(nm>0)
    {
      cat("\nMixture:\n\n")
      param<-attr(object$gmixture,"parameters")
      if (length(param) != nm) param<-"             "
      coeftable<-matrix(cbind(est<-object.coef$mix,
                              se<-object.se$mix,
                              w<-(est/se)^2,(1-pchisq(w,1))),nm,4)
      coeftable<-formatC(matrix(apply(coeftable,2,format,digits=digits,width=13),nm,4),width=13)
      cat(paste(formatC(param,width=-13),apply(coeftable,1,paste,collapse=" "),collapse="\n"),"\n")
      object$Mixture<-coeftable
      cat("\nMean Mixture Probabilities:\n")
      mmix<-apply(object$gmixture(est),2,weighted.mean,eval(object$prior.weights),na.rm=TRUE)
      names(mmix)<-paste("Group",1:object$groups,sep="")
      print(mmix,digits=5)
      object$MixtProb<-coeftable
    }
  modelfit<-cbind("-2*logLikelihood"=-2*object$loglik, AIC=object$AIC,AIC(object,k="BIC"))
  dimnames(modelfit)[[1]]<-" "
  cat("\n")
  print(modelfit)
  object$ModelFit<-modelfit
  coeftable<-apply(post.moc(object),2,weighted.mean,eval(object$prior.weights),na.rm=TRUE)
  object$PostMixtProb<-coeftable
  if(object$code>2) {cat("\n\nWARNING - CONVERGENCE NOT ACHIEVED - ERROR",object$code,"\n\n\n")}
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
      c(-2*tmp$loglikelihood,bic,entropy,bic+2*entropy,tmp$df)}))) else
  val<-as.data.frame(t(sapply(mlist,function(tmp)
                    c(-2*tmp$loglikelihood+k*sum(tmp$npar),tmp$df))))
  names(val)<-c(switch(as.character(k),"0"="-2*logLik","2"=c("AIC"),
                       BIC=c("-2*logLik","BIC","Entropy","ICL-BIC"),
                       "generalized AIC"),"Df")
  row.names(val)<-cnames[1:ml]                     
  val
}

entropy <- function(object,...) UseMethod("entropy")

entropy.default <- function(object,...)
  {
    obj2 <- as.matrix(object)
    if(!all((obj2 >= 0) & (obj2 <= 1)) || any(abs(apply(obj2,1,sum)-1) > .Machine$double.eps^0.5))
      stop("The probabilities must lies in [0,1] ans sum to 1 ")
    en <- apply(obj2,1,function(pro) -sum(ifelse(pro==0,0,pro*log(pro))))
    en <- cbind(en,en/log(dim(obj2)[2]))
    dimnames(en) <-  list(dimnames(obj2)[[1]],c("entropy","std.entropy"))
    return(en)
  }
                                          

entropy.moc <- function(object,...) 
{
  mlist<-  list(object,...)
  if(!all(sapply(mlist,inherits,"moc"))) stop("\n\tAll objects must be of class moc\n")
  ml<-length(mlist)
  cnames<-as.character(match.call()[-1])
  names(mlist)<- cnames[1:ml]
    val<-as.data.frame(t(sapply(mlist,function(tmp){
      attach(eval(tmp$data),name=deparse(tmp$data),pos=2)
      on.exit(detach(2),add=TRUE)
      parm<-split(tmp$coef,rep(c("mu","shape","extra","mix"),tmp$npar)) 
      pri <- tmp$gmixture(parm$mix)
      po<-post.moc(tmp)
      pri.entropy<--sum(ifelse(pri==0,0,eval(tmp$prior.weight)*pri*log(pri)))
      post.entropy<--sum(ifelse(po==0,0,eval(tmp$prior.weight)*po*log(po)))
      c(tmp$groups,pri.entropy,post.entropy,
        pri.entropy/log(tmp$groups)/sum(eval(tmp$prior.weight)),post.entropy/log(tmp$groups)/sum(eval(tmp$prior.weight)),
        1-post.entropy/pri.entropy)}))) 
  names(val)<-c("Groups","Total Prior Entropy","Total Posterior Entropy",
                "Mean Prior Standardized Entropy","Mean Posterior Standardized Entropy",
                "% Reduction")
  row.names(val)<-cnames[1:ml]                     
  val
}


logLik.moc <- function(object,...)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  val <- object$loglikelihood
  attr(val,"df") <- object$df
  attr(val,"nobs")<-object$nvar*object$nsubject
  class(val) <- "logLik"
  attr(val,"moc.name") <- deparse(substitute(object))
  val
}


obsfit.moc<-function(object,along=list(cons=rep(1,object$nsubject)),FUN=function(x) x)
{
  if(!inherits(object,"moc")) stop("\n\tObject must be of class moc\n")
  if(!is.null(eval(object$data))){
    attach(eval(object$data),name=deparse(object$data),pos=2)
    on.exit(detach(2),add=TRUE)}
  along.name <- substitute(along)
  if(!is.list(along)) eval(parse(text=paste("along <-list(",substitute(along),"=along)")))
  n<-object$nsubject
  nt<-object$nvar
  ng<-object$groups
  parm<-split(object$coef,rep(c("mu","shape","extra","mix"),object$npar))
  wts<-eval(object$prior.weight)
  post.obj<-post.moc(object)*wts
  mpost<-by(post.obj,along,mean,na.rm=TRUE)
  wts<-by(wts,along,mean)
  tmp <- object$gmixture(parm$mix)
  dimnames(tmp) <- list(NULL,paste("Group",1:object$groups,sep=""))
  gmix.mean <-by(tmp,along,mean)
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
                      paste("V",1:nt,sep=""),temp))
      }
      if(!is.null(observed.mean[[i]]))
        {
          observed.mean[[i]] <- t(array(observed.mean[[i]]/mpost.observed[[i]],c(nt,ng)))
          dimnames(observed.mean[[i]])<-
            list(paste("Group",1:ng,sep=""),
                 ifelse(is.na(nchar(temp<-dimnames(eval(object$resp))[[2]])[1:nt]) ,
                        paste("V",1:nt,sep=""),temp))
        }
      gmix.mean[[i]]<-gmix.mean[[i]]/wts[[i]]
    }

  val<-list("Mean Prior Probabilities"=gmix.mean,
            "Mean function Expected Values"=fitted.mean,
            "Mean function Observed Values"=observed.mean,
            "Mean Posterior Probabilities"=mpost)
  structure(val,moc.name=deparse(substitute(object)),FUN=substitute(FUN),along=deparse(along.name))
}


plot.moc<-function(x,against=1:x$nvar,main=paste(substitute(x)),xlab="",ylab="",prob.legend=TRUE,scale=FALSE,group.colors=rainbow(x$groups),...)
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
    center<-apply(eval(x$resp),2,function(v)
                  weighted.mean(v,eval(x$prior.weights),na.rm=TRUE))
    scale<-sqrt(apply(eval(x$resp),2,function(v)
                 mean(eval(x$prior.weights)*v[nav<-!is.na(v)]^2)/
                 mean(eval(x$prior.weights)[nav]))-center^2)
  } else {center<-rep(0,x$nvar);scale<-rep(1,x$nvar)}
  matplot(w,cbind((t(x$observed.mean)-center)/scale,(t(x$fitted.mean)-center)/scale),
          type="n",main=main,xlab=xlab,ylab=ylab,col=group.colors,...)
  matpoints(against,(t(x$observed.mean)-center)/scale,col=group.colors,...)
  matlines(against,(t(x$fitted.mean)-center)/scale,type="o",pch=20,cex=.75,col=group.colors,...)
  if(prob.legend) { mtext(legend.text,side=4,outer=TRUE,at=((1:x$groups)+3)/(x$groups+6),col=group.colors) }
  invisible(list(moc.name=deparse(substitute(x)),against=against,fitted=x$fitted,observed=x$observed,center=center,scale=scale))
}

profilesplot<-function(x,...) UseMethod("profilesplot")

profilesplot.moc<-function(x,against=1:x$nvar,main=NULL,xlab="",ylab="",col.legend=TRUE,scale=FALSE,group.colors=rainbow(x$groups),type="subject",...)
{
  if(!inherits(x,"moc")) stop("Not an object of class moc")
  if(!is.null(eval(x$data)))
    {
      attach(eval(x$data),name=deparse(x$data),pos=2)
      on.exit(detach(2),add=TRUE)
    }
  if (is.null(main)) main <- paste(paste(x$resp,collapse=" "),"with",substitute(x))
  type <- match.arg(type,c("subject","variable","posterior"))
  group.colors<-col2rgb(group.colors)/256
  if(col.legend) {
    oldpar<-par("oma"=c(0,0,0,6),"las"=1);on.exit(par(oldpar),add=TRUE)
    legend.text<-paste("Group",1:x$groups)
  }
  if(scale) {
    center<-apply(eval(x$resp),2,mean,na.rm=TRUE)
    scale<-apply(eval(x$resp),2,sd,na.rm=TRUE)
  } else {center<-rep(0,x$nvar);scale<-rep(1,x$nvar)}
  group.rgb<-t(apply(post.moc(x),1,function(y) (group.colors%*%y)))
  if(type=="subject") {
    matplot((against),(t(eval(x$resp))-center)/scale,type="o",pch=20,cex=0.75,
            main=main,xlab=xlab,ylab=ylab,col=rgb(group.rgb[,1],group.rgb[,2],group.rgb[,3]))
  } else
  { if(type=="variable") pairs(eval(x$resp),col=rgb(group.rgb[,1],group.rgb[,2],group.rgb[,3]),upper.panel=NULL,main=main)
  else if(type=="posterior")  pairs(post.moc(x),col=rgb(group.rgb[,1],group.rgb[,2],group.rgb[,3]),upper.panel=NULL,main=main) }
  if(col.legend) { mtext(legend.text,side=4,outer=TRUE,at=((1:x$groups)+3)/(x$groups+6),
                         col=rgb(group.colors[1,],group.colors[2,],group.colors[3,])) }
}

entropyplot <- function(x,...) UseMethod("entropyplot")

entropyplot.moc <- function(x,main=NULL,std=TRUE,lwd=1.5,col=c("red3","green3","gray95"),legend=TRUE,...)
  {
    if(!inherits(x,"moc")) stop("Not an object of class moc")
    if(!is.null(eval(x$data)))
      {
        attach(eval(x$data),name=deparse(x$data),pos=2)
        on.exit(detach(2),add=TRUE)
      }
    parm<-split(x$coef,rep(c("mu","shape","extra","mix"),x$npar))
    prior.entropy <- entropy(x$gmixture(parm$mix))
    posterior.entropy <- entropy(post.moc(x))
    max.ent <- ifelse(std,1,log(x$groups))
    if(std) {col.ind <- 2} else { col.ind <- 1}
    order.pripost <- order(prior.entropy[,col.ind],posterior.entropy[,col.ind])
    n <- length(order.pripost)
    if (is.null(main)) main <- paste("Prior and posterior",ifelse(std,"standardized",""),"entropy for",substitute(x))
    plot(c(0,1),c(0,max.ent),xaxt="n",xlim=c(0,1),ylim=c(0,max.ent),type="n",main=main,xlab="",ylab="")
    polygon(cbind((0:(n-1))/(n-1),((n-1):0)/(n-1)),
            cbind(prior.entropy[order.pripost,col.ind], posterior.entropy[order.pripost[n:1],col.ind]),
                  col=col[3],border=NA)
    lines((0:(n-1))/(n-1),prior.entropy[order.pripost,col.ind],col=col[1],lwd=lwd)
    lines(((n-1):0)/(n-1),posterior.entropy[order.pripost[n:1],col.ind],col=col[2],lwd=lwd)
    if(legend) legend(1,0,c("Prior","Posterior"),lwd=lwd,col=col[1:2],ncol=2,xjust=1,bty="n")
  }
    
    
plot.residuals.moc<-function(x,against="Index",groups=1:dim(x)[3],sunflower=FALSE,group.colors=NULL,...)
{
  if(sunflower) thisplot<-sunflowerplot else thisplot<-plot
  if(!inherits(x,"residuals.moc"))
    stop("Not an object of class residuals.moc !")
  if (!all(groups %in% (1:dim(x)[3])))
    stop("You requested residuals for non-existing groups !")
  if(is.null(group.colors)) group.colors=rainbow(eval(as.name(attr(x,"moc.name")))$groups)
  group.rgb <- mix.colors.moc(eval(as.name((attr(x,"moc.name")))),group.colors=group.colors)
  dim1<-dim(x)
  again<-against
  oldpar<-par(mfrow=c(length(groups),1),oma=c(0,0,2,0));on.exit(par(oldpar),add=TRUE)
  vname<-deparse(substitute(against))
  if(against=="Subject") against<-c(rep(1:dim1[1],dim1[2]))
  if(against=="Observation") against<-c(rep(1:dim1[2],rep(dim1[1],dim1[2])))
  if(against=="Index") against<-1:(dim1[1]*dim1[2])
  res.apply<-sapply(groups,function(i){
    z<-na.omit(cbind(group.rgb,c(against),c(x[,,i])))
    tmp.rgb <- z[,1]
    z <- z[,-1]
    dimnames(z)<-list(NULL,c(vname,attr(x,"type")))
    thisplot(z,main=paste("Group",i),col=tmp.rgb,...)
    if(again=="Index")  abline(v=(1:dim1[2])*dim1[1]+0.5)
  })
  mtext(paste("MOC ",attr(x,"type")," residuals of",attr(x,"moc.name")),side=3,outer=TRUE,font=2)
  attr(res.apply,"moc.name") <- deparse(substitute(x))
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
  cat("\n   Supplementary utility functions to help combine MOC models can be found in\n",
      system.file("Utils","combine.moc.R",package="moc"),"\n   See the file for simple documentation.",
      "\n   You must source that file to use those functions.\n\n")
}
    
