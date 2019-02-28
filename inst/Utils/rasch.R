## author: Bernard Boulerice
## functions to construct item response theory models
## with MOC version 0.8 and higher.
## Copyright (C) 2002-2019 Bernard Boulerice

require(moc)

# gherm.quad gives the nodes and weights of a quadrature to
# approximate a normal integral (Hermite quadrature).
# rasch returns a function that depends on a mean value and "thresholds" and compute
# the joint likelihood of a group of items. It suppose the existence of a
# variable of type "list" herm containing the points and weights of the quadrature.

gherm.quad <-
function (n) 
{
    i <- 1:n
    i1 <- 1:(n - 1)
    a <- rep(0, n)
    b <- sqrt(i1/2)
    A <- rep(0, n * n)
    A[(n + 1) * (i - 1) + 1] <- a
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))^2
    x <- rev(vd$values)*sqrt(2)
    list(nodes = x, weights = w)
}

# The functions "rasch.norm" and "rasch.logis" return the response probabilities 
# to each items over the different nodes, mean and thresholds combinations.
# It supposes the existence of the variables:
            # "nitems" containing the number of items
            # "nlevs"  a vector containing the number of levels for each item
            # "herm"  a list containing the Hermite quadrature

rasch.norm <-function(mu,p,loading)
{
thresh<-split(p,rep(1:nitems,(nlevs+1)))
norm.val<-lapply(1:nitems,function(it) sapply(herm$nodes,
          function(nd) apply(outer(thresh[[it]],mu,
                             function(trh,mu1) pnorm(trh-loading[it]*(nd+mu1))),2,diff)))
lapply(1:nitems,function(it) array(norm.val[[it]],c(nlevs[it],dim(norm.val[[it]])/c(nlevs[it],1))))
}

rasch.logis <-function(mu,p,loading)
{
  thresh<-split(p,rep(1:nitems,(nlevs+1)))
norm.val<-lapply(1:nitems,function(it) sapply(herm$nodes,
          function(nd) apply(outer(thresh[[it]],mu,
                             function(trh,mu1) plogis(trh-loading[it]*(nd+mu1))),2,diff)))
lapply(1:nitems,function(it) array(norm.val[[it]],c(nlevs[it],dim(norm.val[[it]])/c(nlevs[it],1))))
}


## The function rasch.density only returns mu, however when we call MOC
## resp (x) should be an arbitrary matrix of dimension n.subjects by n.variables.
## We can use the standardized scores of the scale corresponding to the items
## such that the observed.mean make sens.

## The function rasch.gmu uses the real matrix containing the items responses
## which should be called "resp.item".
## An example of a function  rasch.gmu is given which supposes the existence
## of the variables "ntimes" and "ngroups" and of a function "norm.gmu"
## which gives the normal curve for each subject at each time and group
## like for any call to MOC.
## We have to assign one of the function  "rasch.norm" or "rasch.logis"
## to "rasch", example: "rasch <- rasch.norm".

##The first parameters "p" stands for the thresholds followed by the
## parameters of "norm.gmu".


rasch.density <- function(x,mu,...) mu

rasch.gmu.group <- function(p,group)
{
  nquad<-length(herm$weights)
  nsj<-dim(resp.item)[1]
  th.parm<-split(p[1:sum(nlevs-1)],rep(1:nitems,(nlevs-1)))
  loading<-p[(sum(nlevs-1)+1):(sum(nlevs-1)+nitems)]
  mu.parm<-p[-(1:(sum(nlevs-1)+nitems))]
  thresh<-c(unlist(sapply(th.parm,function(itpar) qnorm(c(1E-10,cumsum(inv.glogit(itpar))-1E-10)))))
  mu<-norm.gmu[[group]](mu.parm)
  prob.val<-rasch(mu,thresh,c(loading))
  tmp<-array(NA,c(nsj,ntimes))
   for(t in 1:ntimes)
    {
      tmp2<-array(1,c(nsj,nquad))
      for (i in 1:nitems)
        {
        nana<-which(!is.na(resp.item[,(t-1)*nitems+i]) )
        tmp2[nana,]<-tmp2[nana,]*prob.val[[i]][resp.item[nana,(t-1)*nitems+i],t,]
        }
      tmp[,t]<-apply(tmp2,1,function(x) x%*%herm$weights)
    }
  return(matrix(tmp,nsj,ntimes))
}

rasch.expected.group <- function(p,group)
{
  nquad<-length(herm$weights)
  nsj<-dim(resp.item)[1]
  th.parm<-split(p[1:sum(nlevs-1)],rep(1:nitems,(nlevs-1)))
  loading<-p[(sum(nlevs-1)+1):(sum(nlevs-1)+nitems)]
  mu.parm<-p[-(1:(sum(nlevs-1)+nitems))]
  thresh<-c(unlist(sapply(th.parm,function(itpar) qnorm(c(1E-10,cumsum(inv.glogit(itpar))-1E-10)))))
  mu<-norm.gmu[[group]](mu.parm)
  if (dim(rbind(mu))[1]==1) mu<-matrix(mu,nsj,ntimes,byrow=TRUE)
  prob.val<-rasch(mu,thresh,c(loading))
  tmp<-array(NA,c(nsj,ntimes))
   for(t in 1:ntimes)
    {
      tmp2<-array(1,c(nsj,nquad))
      for (i in 1:nitems)
        {
        nana<-which(!is.na(resp.item[,(t-1)*nitems+i]) )
        tmp2[nana,]<-tmp2[nana,]*prob.val[[i]][resp.item[nana,(t-1)*nitems+i],t,]
        }
      tmp[,t]<-apply(tmp2,1,function(x) ((x*herm$nodes)%*%herm$weights)/x%*%herm$weights)
    }
  return(matrix(tmp,nsj,ntimes)+mu)
}


rasch.item.prob <- function(object)
{
  if (!inherits(object, "moc"))
        stop("\n        Object must be of class moc\n")
    if (!is.null(eval(object$data))) {
        attach(eval(object$data), pos=2,name="rasch.item.data")
        on.exit(detach(2))
    }
  p<-object$coef[1:sum(object$npar[-4])]
  nquad<-length(herm$weights)
  nsj<-dim(resp.item)[1]
  th.parm<-split(p[1:sum(nlevs-1)],rep(1:nitems,(nlevs-1)))
  loading<-p[(sum(nlevs-1)+1):(sum(nlevs-1)+nitems)]
  mu.parm<-p[-(1:(sum(nlevs-1)+nitems-1))]
  thresh<-c(unlist(sapply(th.parm,function(itpar) qnorm(c(1E-10,cumsum(inv.glogit(itpar))-1E-10)))))
  mu<-lapply(1:ngroups,function(gr) norm.gmu[[gr]](mu.parm))
  prob.val<-lapply(mu,function(x) {rasch(x,thresh,loading)})
  prob.val<-lapply(1:ngroups,function(gr) {lapply(1:nitems,function(it) array(apply(prob.val[[gr]][[it]],c(1,2),
         function(vec) vec%*%herm$weights ),c(nlevs[it],ntimes)))})
  names(prob.val)<-(paste("Group",1:ngroups))
  for (g in 1:ngroups) {names(prob.val[[g]])<-(paste("Item",1:nitems))
  for(it in 1:nitems) dimnames(prob.val[[g]][[it]])<-list(NULL,paste("Time",1:ntimes))}
  prob.val
}


moc.rasch<-function(rasch.data=list(resp.item=NULL,rasch=rasch.norm,herm=gherm.quad(5),
                    norm.gmu=NULL,ngroups=NULL,ntimes=NULL,nitems=NULL,nlevs=NULL),
                    gmixture = inv.glogit,                   
                    mu.start=NULL,thresh.start=NULL,loading.start=NULL,pgmix=NULL,
                    scale.weight = FALSE,wt=1,data=NULL,
                    gradtol = 1e-04, steptol = gradtol, iterlim = 100)
  {
    if(dim(rasch.data$resp.item)[2]!=rasch.data$nitems*rasch.data$ntimes) stop(paste("resp.item should have",nitems*ntimes,"columns."))
    if(length(rasch.data$nlevs)!=rasch.data$nitems) stop(paste("nlevs should have",nitems,"elements."))
    if (!all(sapply(1:dim((rasch.data$resp.item))[2],function(i)
                   all(rasch.data$resp.item[,i]%in%c(1:rep(rasch.data$nlev,rasch.data$ntimes)[i],NA))))){
    stop(paste("Column",i,"of resp.item should contain only the numbers 1 to",
               rep(rasch.data$nlev,rasch.data$ntimes)[i],"or NA."))}
    rasch.gmu<-lapply(1:rasch.data$ngroups,function(i)
                      eval(parse(text=paste("function(p) rasch.gmu.group(p,",i,")"))))
    names(rasch.gmu)<-paste("Group",1:rasch.data$ngroups,sep="")
    rasch.expected<-lapply(1:rasch.data$ngroups,function(i)
                 eval(parse(text=paste("function(p) rasch.expected.group(p,",i,")"))))
    names(rasch.expected)<-paste("Group",1:rasch.data$ngroups,sep="")
    y<-(apply(rasch.data$resp.item,1,function(x) apply(matrix(x,rasch.data$nitems,rasch.data$ntimes),2,mean,na.rm=TRUE)))
    wt2 <- wt
    if(length(wt2)==1) wt2 <- rep(wt2,dim(y)[2])
    y<-y-apply(y,1,weighted.mean,wt2,na.rm=TRUE)
    y<-t(y/apply(y^2,1,weighted.mean,wt2,na.rm=TRUE))
     eval(substitute(expression(dat$resp.dumm <- y),list(dat=(substitute(rasch.data)))),envir = parent.frame())
    val<- eval(substitute(moc(y,density=function(x,mu,...) {mu},gmu=rasch.gmu,
               groups=ngroups,expected=rasch.expected,
               pgmu=c(start1,start2,start3),
               pgmix=pgmix,gmixture=gmixture,scale.weight = FALSE,wt=wt,
               data=c(data,rasch.data),gradtol = gradtol, steptol = steptol, iterlim = iterlim)),
                          list(ngroups=rasch.data$ngroups,wt=substitute(wt),
                               rasch.data=substitute(rasch.data),gmixture=substitute(gmixture),
                               data=substitute(data),start1=thresh.start,
                               start2=loading.start,start3=mu.start,pgmix=pgmix,
			       gradtol = gradtol, steptol = steptol, iterlim = iterlim))
    val$item.prob.expected<-with(rasch.data,{rasch.item.prob(val)})
    post.prob<-post(val)
    val$item.prob.observed <-
      lapply(1:rasch.data$nitems,function(it) array(apply(post.prob*wt2,2,function(y)
           apply(outer(rasch.data$resp.item[,(0:(rasch.data$ntimes-1))*rasch.data$nitems+it],
                       1:rasch.data$nlevs[it],"==")*1,c(2,3),
                 function(x) weighted.mean(x,y,na.rm=T))),
                             c(rasch.data$ntimes,rasch.data$nlevs[it],rasch.data$ngroups)))
    names(val$item.prob.observed) <- paste("Item",1:rasch.data$nitems)
    tmp <- lapply(1:rasch.data$nitems,function(it) names(val$item.prob.observed[[it]])<-
           list(paste("Time",1:rasch.data$ntimes),paste("level",1:rasch.data$nlevs[it]),
                paste("Group",rasch.data$ngroups)))
    val$resp <- as.name("resp.dumm")
    val
  }
