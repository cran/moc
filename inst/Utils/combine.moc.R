   # Copyright (C) 2002-2019 Bernard Boulerice

require(moc)

## This file contains the following functions:
##    find.unique.pattern
##    combine.density
##    combine.parfun
##    combine.prob

# The function combine.density combines two densities for two data sets
# to be used in a joint MOC analysis. The function returns a function that
# computes the joint density of the combined vector and must be assigned.
# The two functions name must be quoted.
# The length of the data vector (number of variables)
# to which they apply is specified through nvar1 and nvar2.
# If dname1 and/or dname2 already computes the
# joint density of the vector to which they apply this must be
# specified with is.joint1 and is.joint2. If you want to combine more
# than two data sets you can use the function recursively.

## For example if you have 2 fitted moc objects (moc1,moc2): you can use
## density.combined <- combine.density(moc1$call$density,moc2$call$density,
## moc1$nvar, moc2$nvar,is.joint1=moc1$joint,is.joint2=moc2$joint)
## BE CAREFULL HOWEVER that moc{1,2}$call$density are in fact names,
## you can have big probelms when this is not the case.
## Also when moc have been called with a data argument, the density names
## needs that data frame to be resolved.

"combine.density" <-
function(dname1,dname2,nvar1,nvar2,is.joint1=FALSE,is.joint2=FALSE)
{
if (!exists(dname1) ) stop(paste("The density ",dname1," must be defined."))
if (!exists(dname2) ) stop(paste("The density ",dname2," must be defined."))
if ((length(c(nvar1,nvar2))!=2) | (nvar1<1) | (nvar2<1) | any(as.integer(c(nvar1,nvar2))!=c(nvar1,nvar2)))
  stop("nvar1 and nvar2 must be 2 integers greater or equal to 1")
range1<-paste("[,1:",nvar1,"]")
range2<-paste("[,",nvar1+1,":",nvar1+nvar2,"]")
fbody<-paste("function(x,mu,shape,extra){\n")
if(is.joint1)
  {fbody<-paste(fbody,"d1<-",dname1,"(x",range1,",mu",range1,",shape",range1,",extra",range1,")\n")} else
{fbody<-paste(fbody,"d1<-apply(",dname1,"(x",range1,",mu",range1,",shape",range1,",extra",range1,")",",1,prod,na.rm=TRUE)\n")}
if(is.joint2)
  {fbody<-paste(fbody,"d2<-",dname2,"(x",range2,",mu",range2,",shape",range2,",extra",range2,")\n")} else
{fbody<-paste(fbody,"d2<-apply(",dname2,"(x",range2,",mu",range2,",shape",range2,",extra",range2,")",",1,prod,na.rm=TRUE)\n")}
fbody<-paste(fbody,"d1*d2\n}\n")
eval(parse(text=fbody),envir=.GlobalEnv)
}


# The function combine.parfun can be used to combine the parameter
# functions for MOC (gmu, gshape and gextra) for two data sets,
# it returns a function that must be assigned.
#
# The function names fname1 and fname2 must be supplied quoted,
# one of the function name can be the empty string "" but not both.
#
# You must supply the parameter indexes of each function in np1 and np2,
# the returned number of rows in n (both function should return the same number of rows)
# the vector length of each data set in nvar1 and nvar2, and the number of groups
# in ng1 and ng2.
#
# The set parameter must be a two columns matrix that specifies which groups of the
# first data set should be combined with which groups of the second data set,
# by default it assumes that a full cross-classification of the groups is requested.
# A subset can be specified with commands like merge(1:3,1:4)[c(2,3,6,8,11,12),]

## For example you can combine the gmu function of 2 fitted moc objects (moc1,moc1) with:
## gmu.combined <- combine.parfun(moc1$call$gmu,moc2$call$gmu,
##                 np1=1:moc1$npar[1],np2=1:moc2$npar[1],
##                 n=unique(c(moc1$nsubject,moc2$nsubject)),
##                 ng1=moc1$groups,ng2=moc2$groups,set=merge(1:ng1,1:ng2))

"combine.parfun" <-
function(fname1,fname2,np1,np2,nvar1,nvar2,n,ng1,ng2,set=merge(1:ng1,1:ng2))
{
if(dim(set)[2]!=2) stop("set must be a two column matrix")
cat("\nThe combined functions is for a ",dim(set)[1]," groups model.\n\n")
if ( any(np1<0) | any(np2<0) | any(as.integer(c(np1,np2))!=c(np1,np2)))
  stop("np1 and np2 must be sets of integers greater or equal to 0")
if ((length(c(nvar1,nvar2))!=2) | (nvar1<1) | (nvar2<1) | any(as.integer(c(nvar1,nvar2))!=c(nvar1,nvar2)))
  stop("nvar1 and nvar2 must be 2 integers greater or equal to 1")
if ((length(c(ng1,ng2))!=2) | (ng1<1) | (ng2<1) | any(as.integer(c(ng1,ng2))!=c(ng1,ng2)))
  stop("ng1 and ng2 must be 2 integers greater or equal to 1")
if ((fname1=="") & (fname2=="")) stop("At least one of the function must be non-empty")
x<-list();y<-list()
if (fname1=="")
  { np1<-0
    for(i in 1:ng1) x[[i]]<-paste("x<-array(1,c(",paste(n,nvar1,sep=","),"))\n")
    n1<-n
  } else
  {
    if (!exists(fname1)) stop(paste("The list of functions ",fname1," must exists"))
    if (eval(parse(text=paste("length(",fname1,")")))!=ng1) stop("Wrong number of groups ng1.\n")
    for(i in 1:ng1) {
      dim1<- dim(eval(parse(text=paste(fname1,"[[",i,"]](rep(0,",length(np1),"))",sep=""))))
      n1<- dim1[1]
      if((n1>1) & (n1!=n)) stop(paste(fname1,"[[",i,"]] must return 1 row or ",
                     n," rows"))
      if(dim1[2]!=(nvar1)) stop(paste("The number of columns returned by ",fname1,
                                      "[[",i,"]] is incompatible with nvar1"))
      if((n1==1) & (n>1)) temp<-paste("%x%rep(1,",n,")") else temp=""
    x[[i]]<- paste("x<-array(",fname1,"[[",i,"]](p[",paste(deparse(np1),collapse=""),"])",temp,",c(",
                   paste(n,nvar1,sep=","),"))\n")
    }
  }
if (fname2=="")
  { np2<-0
     for(i in 1:ng1) y[[i]]<-paste("y<-array(1,c(",paste(n,nvar2,sep=","),"))\n")
    n2<-n
  } else
  {
    if (!exists(fname2)) stop(paste("The function ",fname2," must exists"))
    if (eval(parse(text=paste("length(",fname2,")")))!=ng2) stop("Wrong number of groups ng2.\n")
     for(i in 1:ng2) {
    dim2<- dim(eval(parse(text=paste(fname2,"[[",i,"]](rep(0,",length(np2),"))",sep=""))))
    n2<-dim2[1]
    if((n2>1) & (n2!=n)) stop(paste(fname2,"[[",i,"]] must return 1 row or ",n," rows"))
    if(dim2[2]!=(nvar2)) stop(paste("The number of columns returned by ",fname2,
                                        "[[",i,"]]is incompatible with nvar2"))
    if((n2==1) & (n>1)) temp<-paste("%x%rep(1,",n,")") else temp=""
    y[[i]]<- paste("y<-array(",fname2,"[[",i,"]](p[",paste(deparse(np2),collapse=""),"])",temp,",c(",paste(n,nvar2,sep=","),"))\n")
  }
  }
f.list<-list()
for(i in 1:dim(set)[1]) {
  res<-paste("cbind(x,y)\n")
  fbody<-paste("function(p){\n",x[[set[i,1]]],y[[set[i,2]]],res,"}")
  f.list[[ paste("Group", paste(set[i, ],collapse=""),sep="")]]<-
    eval(parse(text=fbody),envir=.GlobalEnv)
}
f.list
}


## The function combine.prob combines the mixture probability 
## for two moc models. This function is useful to compute starting values for
## the full joint model. For example glogit(combine.prob(moc1,moc2))

"combine.prob" <-
  function(moc1,moc2)
  {
    if(!inherits(moc1, "moc") | !inherits(moc2, "moc")) stop("moc1 and moc2 must be moc objects")
    post1<-post(moc1)
    post2<-post(moc2)
    if(dim(post1)[1] != dim(post2)[1]) stop("moc1 and moc2 should apply to the same subjects")
    crossprob<-apply(post2,2,function(x) apply(post1*x,2,mean,na.rm=TRUE))
    mix.logit<-c(crossprob)
    mix.logit<-log(t(mix.logit)[-1]/mix.logit[1])
    list(CrossTab=crossprob,logit=mix.logit)
  }

## The following function find unique patterns in the data X
## and returns a matrix with those pattern and a column of total weight
## The resulting pattern can be used as the response for a moc model
## with the weights as the wt argument in moc.

find.unique.pattern<-function(X,w=rep(1,dim(X)[1]))
{
  if(dim(X)[1]!=dim(as.matrix(w))[1]) stop("X and w must have the same number of rows")
  sx<-apply(X, 1, paste, collapse = " ")
  sxu<-unique(sx)
  nrow<-length(sxu)
  ncol<-dim(X)[2]
  index<-sapply(sxu,function(x) which(x==sx))
  names(index)<-NULL
  new.w<-unlist(lapply(index,function(z) sum(w[z])))
  sxu<-t(matrix(as.numeric(unlist(lapply(sxu,strsplit,split=" "))),ncol,nrow))
  cbind(sxu,new.w)
}
