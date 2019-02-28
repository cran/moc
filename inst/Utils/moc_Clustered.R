    #Copyright (C) Bernard Boulerice 2015-2019
## The utility functions included here are used to construct 
## moc models for clustered sampling as in family data or litters
## data in animals studies. Note however that these are not intended
## to be used for large clusters.

require(moc)

# used to make the product likelihooh within cluster
clus.prod <- function(V,clus,na.rm=TRUE) {
  nc <- dim(V)
  apply(V,2,function(cc){tapply(cc,clus,prod,na.rm=na.rm)})  
} 


density.clus <- function(x,mu,sig=1,extra,Density=dnorm) {
  nc <- dim(x)
  tmp <- array(1,nc)
  tmp2 <- clus.prod(Density(x,mu+extra[,seq(2,2*nc[2],2)],extra[,seq(3,2*nc[2]+1,2)]),
                    extra[,1],na.rm=TRUE)
  tmp[1:dim(tmp2)[1],] <- tmp2
  tmp
}

normal.clus <- function(x,mu,sig=1,extra) {
  nc <- dim(x)
  tmp <- array(1,nc)
  tmp2 <- clus.prod(dnorm(x,mu+extra[,seq(2,2*nc[2],2)],extra[,seq(3,2*nc[2]+1,2)]),
                    extra[,1],na.rm=TRUE)
  tmp[1:dim(tmp2)[1],] <- tmp2
  tmp
}


    #Returns the nodes and weights of an Hermitian Quadrature
    #corresponding to a standard normal. Used when inter-cluster density
    #is assumed normally distrited (Gaussian).

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

    #Construct the clustered parameter functions
make.clus.extra <- function(cluster=NULL,param,labels="",covariates="",quad=gherm.quad(7),data=.GlobalEnv)
  {
    if(!is.character(param)) stop("param must be character.\n")
    if(!is.character(covariates)) stop("covariates must be character vector.\n")
    if(!exists("cluster") && !("cluster" %in% names(data)))
      stop("The cluster variable is not defined.\n")
    if(!all(c( "nodes"  , "weights") %in% names(quad)))
      stop("quad should be a list containing the vectors nodes and weights.\n")
    ncov <- length(covariates)
    if(all(covariates!="")) covstring <- paste(covariates,collapse=",") else {ncov <- 1;covstring="Cons"}
    tmp <- eval(substitute( lapply(.nodes,function(nod) eval(parse(text=
           paste("function (p) {
                  Cons <- cbind(rep(1,length(",cluster,")))
                  pmat <- matrix(",param,"*rep(c(",nod,",1),c(",ncov,",",ncov,")),",ncov,",2)\n
                  cbind(",cluster,",cbind(",covstring,")%*%pmat)}")),
                   envir=.data)), list(.nodes=quad$nodes,.data=data)))
    if(paste(labels,collapse="")!="") attr(tmp,"parameters") <- paste(labels)
    attr(tmp,"quad") <- quad
    attr(tmp,"ncov") <- ncov
    tmp
  }


make.clus.mu <- function(cluster=NULL,param,labels="",covariates="",quad=gherm.quad(7),data=.GlobalEnv)
  {
    if(!is.character(param)||length(param)!=1) stop("param must be character.\n")
    if(!is.character(labels)) stop("labels must be character vector.\n")
    if(!is.character(covariates)) stop("covariates must be character vector.\n")
    if(!exists("cluster") && !("cluster" %in% names(data)))
      stop("The cluster variable is not defined.\n")
    if(!all(c( "nodes"  , "weights") %in% names(quad)))
      stop("quad should be a list containing the vectors nodes and weights.\n")
    ncov <- length(covariates)
    if(all(covariates!="")) covstring <- paste(covariates,collapse=",") else {ncov <- 1;covstring="Cons"}
    tmp <- eval(substitute( lapply(.nodes,function(nod) eval(parse(text=
           paste("function (p) {
                       Cons <- cbind(rep(1,length(",cluster,")))
                       pmat <- matrix(",param,",",ncov,",1)\n
                       cbind(",covstring,")%*%pmat}")
                  ),envir=.data)),
                  list(.nodes=quad$nodes,.data=data)))
    if(paste(labels,collapse="")!="") attr(tmp,"parameters") <- paste(labels)
    attr(tmp,"ncov") <- ncov
    tmp
  }



make.clus.expected <- function(cluster=NULL,mu.param,mu.covariates="",
                               extra.param,extra.covariates="",
                               quad=gherm.quad(9),data=.GlobalEnv)
{
  if(!is.character(mu.param)||length(mu.param)!=1) stop("mu.param must be character.\n")
  if(!is.character(extra.param)||length(extra.param)!=1)
    stop("extra.param must be character.\n")
    if(!all(c( "nodes"  , "weights") %in% names(quad)))
      stop("quad should be a list containing the vectors nodes et weights.\n")

  mu.ncov <- length(mu.covariates)
  
  if(all(mu.covariates!="")) mu.covstring <- paste(mu.covariates,collapse=",") else {
    mu.ncov <- 1
    mu.covstring="Cons"}

  extra.ncov <- length(extra.covariates)

  if(all(extra.covariates!="")) extra.covstring <- paste(extra.covariates,collapse=",") else {
    extra.ncov <- 1
    extra.covstring="Cons"}

  tmp <- eval(substitute( lapply(.nodes,function(nod) eval(parse(text=
           paste("function (p) {
                       Cons <- cbind(rep(1,length(",cluster,")))
                       mu.pmat <- matrix(",mu.param,",",mu.ncov,",1)\n
                       p <- p[-(1:(",mu.ncov,"))]
                       extra.pmat <- matrix(",extra.param,"*rep(c(",nod,",1),c(",
                       extra.ncov,",",extra.ncov,")),",extra.ncov,",2)[,1]\n
                       cbind(",mu.covstring,")%*%mu.pmat+
                       cbind(",extra.covstring,")%*%extra.pmat}")
                  ),envir=.data)),
                  list(.nodes=quad$nodes,.data=data)))
  tmp
}


make.clus <- function(cluster=NULL,
                      mu.param, mu.labels="", mu.covariates="",
                      extra.param, extra.labels="", extra.covariates="",
                      quad=gherm.quad(7),data=.GlobalEnv,expected=FALSE)
  {
    result <- list()
    result$groups <- length(quad$nodes) 
    result$density <- normal.clus
    result$gmixture <- eval(substitute(function(p) rbind(quad$weights)))
    result$gmu <- make.clus.mu(cluster,mu.param,mu.labels,mu.covariates,quad=quad,data=data)
    result$gextra <- make.clus.extra(cluster,extra.param,extra.labels,extra.covariates,quad=quad,data=data)

    if(expected) {
      result$expected <- make.clus.expected(cluster,mu.param, mu.covariates, extra.param,extra.covariates,
                                            quad=quad,data=data)
    }

    result$data <- as.name(paste(substitute(data)))
    attr(result,"cluster") <- as.name(cluster)
    result
  }

fit.mocclus <- function(y,moc.list,pgmu,pgextra,...)
  {
    clus.list <- c(moc.list,...)
    clus.list$y <- substitute(y)
    attach(as.environment(eval(moc.list$data)),pos=2,name="fitmocclus")
    on.exit(detach("fitmocclus"))
    if(missing(pgmu)) {
      clus.list$pgmu <- c(mean(y),rep(0,attr(moc.list$gmu,"ncov")-1))
    } else clus.list$pgmu <- pgmu
    if(missing(pgextra)) {
      tmp <- c(var(by(y,eval(attr(moc.list,"cluster")),mean)),mean(by(y,eval(attr(moc.list,"cluster")),var)))
      pgextra <- log(tmp)/2
      clus.list$pgextra <- rep(pgextra,attr(moc.list$gextra,"ncov"))
    } else clus.list$pgextra <- pgextra
#    detach(2)
    clus.list$joint <- TRUE
    clus.list$check.length <- FALSE
    tmp <- do.call("moc",clus.list)
    if(!is.null(moc.list$expected)) tmp$call$expected <- as.name(paste(substitute(moc.list),"$expected",sep=""))
    tmp$call$density <- as.name("normal.clus")
    tmp$call$gmu <- as.name(paste(substitute(moc.list),"$gmu",sep=""))
    tmp$call$gextra <-as.name(paste(substitute(moc.list),"$gextra",sep=""))
    tmp$call$gmixture <- as.name(paste(substitute(moc.list),"$gmixture",sep=""))
    tmp
  }


mocclus.qqnormplot <- function(object,Cluster=NULL,quad=attr(object$gextra,"quad"))
  { if(!("moc" %in% class(object))) stop("Object should be of class moc.\n")
    if(is.null(Cluster)) stop("Cluster should be supplied there is no default.\n")
    if(is.null(quad)) stop("quad does not exist !\n")
    nfam <- length(unique(Cluster))
    tmp <- post(object)[1:nfam,]
    expected.clus <- apply(tmp,1,function(row) mean(row*quad$node))
    qqnorm(expected.clus)
    qqline(expected.clus)
    invisible(expected.clus)
  }


dens.clus <- function (dens, clus) 
{
    clus.tmp <- unclass(factor(clus))
    tmp <- function(x, mu, sig, extra) {
        dens.joint <- apply(dens(x, mu, sig, extra), 1, prod, 
            na.rm = TRUE)
        dens.clus <- tapply(dens.joint, clus.tmp, sum, na.rm = TRUE)
        dens.joint * dens.clus[clus]
    }
    tmp
}

cat("\nmoc_Clustered loaded!\n")
