##Copyright (C) 2003-2008 Bernard Boulerice

moc.gui <- function() {
                                        # Load required packages
require(moc) || warning("You won't be able to run your model: install MOC")
require(tcltk) || stop("You should install TclTk to run this program")

                                        #  Set the tcl global variables
## Set the number of groups and hold the current group function setting
Groups <- tclVar("1")
mu.Groups <- tclVar("1")
shape.Groups <- tclVar("1")
extra.Groups <- tclVar("1")

## Set various MOC options
is.joint <- tclVar("0")
wtvar <- tclVar("Weight")
scale.weight <- tclVar("0")
chklen <- tclVar("1")
maxit.value <- tclVar("100")
gradtol.value <- tclVar("1e-5")
steptol.value <- tclVar("1e-5")
printlevel.value <- tclVar("1")

## Hold lists for menus
datalst <- tclVar()
resplst <- tclVar()
moclst <- tclVar()
funlst <- tclVar()
mixfunlst <- tclVar()
mufunlst <- tclVar()
shapefunlst <- tclVar()
extrafunlst <- tclVar()
expectedfunlst <- tclVar()

## Hold selections in menus
data.sel <- tclVar(0)
resp.sel <- tclVar()
mix.sel <- tclVar()
dist.sel <- tclVar()
mu.sel <- tclVar()
shape.sel <- tclVar()
extra.sel <- tclVar()
expected.sel <- tclVar()
dist.varsel <- tclVar("ALL")
mu.link <- tclVar("Identity")
shape.link <- tclVar("Identity")
extra.link <- tclVar("Identity")
mu.fun <- tclVar("None")
shape.fun <- tclVar("None")
extra.fun <- tclVar("None")
degree <- tclVar(3)  ## degree of the polynomial function
distfun <- tclVar("User")
mixfun <- tclVar("inv.glogit")
expectedfun <- tclVar("Mu")

## Set starting values constraints and optional labels
mu.start.val <- tclVar()
mu.label.val <- tclVar()
mu.const.val <- tclVar()
shape.start.val <- tclVar()
shape.label.val <- tclVar()
shape.const.val <- tclVar()
extra.start.val <- tclVar()
extra.label.val <- tclVar()
extra.const.val <- tclVar()
mix.start.val <- tclVar()
mix.label.val <- tclVar()
mix.const.val <- tclVar()

## Message box showing polynomial degree
mu.mess <- tclVar("  ")
shape.mess <- tclVar("  ")
extra.mess <- tclVar("  ")
dist.mess <- tclVar("  ")
run.mess <- tclVar("Ready")
R.val <- NULL

if("courier" %in% as.character(tkfont.names())) tkfont.delete("courier")
courier.font <- tkfont.create("courier",size="10",family="terminal") #used by the Tk pager

# R global variables

moc.gui.env <- environment() #hold the following lists
mocid <- 1                   # to give a number to fitted models
gmu.rlst <- list()
gshape.rlst <- list()
gextra.rlst <- list()
dist.rlst <- list()


base<-tktoplevel()
tkwm.title(base,"Mixtures with MOC")
moc.menu <- tkmenu(base)
tkconfigure(moc.menu,type="menubar")
tkconfigure(base,menu=moc.menu)

data.frm <- tkframe( base, borderwidth=2, relief="groove")
distmix.frm <- tkframe( base, borderwidth=2, relief="sunken")
top.frm <- tkframe( base, borderwidth=2, relief="groove")
bot.frm <- tkframe( base, borderwidth=2, relief="groove")
optim.frm <- tkwidget( base, "labelframe",text="Algorithm",borderwidth=2, relief="groove")

run.message <- tkmessage(optim.frm,textvariable=run.mess,width=510,justify="right",anchor="e")


submit <- function()
{
  refresh()
  tclObj(run.mess) <- "Running"
  on.exit({tclObj(run.mess) <- "Error";tkmessageBox(type="ok",message=geterrmessage())},add=TRUE)
  moc.call <- formals(moc)
  moc.call <- moc.call[-which(names(moc.call)=="...")]
  moc.call <- as.call(c(quote(moc),moc.call))
  data.list <- list()
  y <-  as.character(tclObj(data.sel))
  base.name <- paste(y,".moctk",mocid,sep="")
  if(exists(base.name,envir=.GlobalEnv)) base.name <- paste(base.name,".new",sep="")
  data.name <- paste(base.name,".data",sep="")
  resp.name <- paste(y,".resp",sep="")
  data.list[[resp.name]] <- as.data.frame(eval(as.name(y)))
  moc.call$y <- parse(text=paste(resp.name,"[,",deparse(as.character(tclObj(resp.sel))),"]"))[[1]]
  nvar <- dim(eval(moc.call$y,env=data.list))[2]
  ns <- dim(eval(moc.call$y,env=data.list))[1]
  moc.call$data <- as.name(data.name)
  moc.call$groups <- as.numeric(tclObj(Groups))
  moc.call$check.length <- as.logical(tclObj(chklen))
  moc.call$scale.weight <- as.logical(tclObj(scale.weight))
  if(!is.na(as.character(tclObj(wtvar))[2])) {
    moc.call$wt <-  as.character(tclObj(wtvar))[2]
  } else {moc.call$wt <- NULL}
  moc.call$iterlim <- as.numeric(tclObj(maxit.value)) 
  moc.call$gradtol <- as.numeric(tclObj(gradtol.value))
  moc.call$steptol <- as.numeric(tclObj(steptol.value))
  moc.call$print.level <- as.numeric(tclObj(printlevel.value))
  mixture <- makemixfun()
  if(as.numeric(tclObj(Groups)) > 1) {
    moc.call$gmixture <- as.name(paste(base.name,".mixt",sep=""))
  } else {moc.call$gmixture <- NULL}
  data.list[[paste(base.name,".mixt",sep="")]] <- mixture$FUN
  moc.call$pgmix <- mixture$start.val
  dens.tk <- makedensity()
  moc.call$joint <- dens.tk$joint
  if(!dens.tk$joint) {
    if(names(dens.tk[-1])=="User") {moc.call$density <- dens.tk$User}
    else { moc.call$density <- as.name(names(dens.tk[-1]))
           data.list[[names(dens.tk[-1])]] <- dens.tk[[-1]]}
  } else {
    moc.call$density <- as.name("joint.dens")
    for(nm in names(dens.tk[-1])) data.list[[nm]] <- dens.tk[-1][[nm]]
  }
  gmu <- makegfun(gmu.rlst,nvar,ns)
  if(!is.null(gmu)) {
    moc.call$pgmu <- gmu$start.val
    moc.call$gmu <- as.name(paste(base.name,".gmu",sep=""))
    data.list[[paste(base.name,".gmu",sep="")]] <- gmu$FUN
  } else { moc.call$pgmu <- NULL; moc.call$gmu <- NULL}
  gshape <- makegfun(gshape.rlst,nvar,ns)
  if(!is.null(gshape)) {
    moc.call$pgshape <- gshape$start.val
    moc.call$gshape <- as.name(paste(base.name,".gshape",sep=""))
    data.list[[paste(base.name,".gshape",sep="")]] <- gshape$FUN
  }  else { moc.call$pgshape <- NULL; moc.call$gshape <- NULL}
  gextra <- makegfun(gextra.rlst,nvar,ns)
  if(!is.null(gextra)) {
    moc.call$pgextra <- gextra$start.val
    moc.call$gextra <- as.name(paste(base.name,".gextra",sep=""))
    data.list[[paste(base.name,".gextra",sep="")]] <- gextra$FUN
  }  else { moc.call$pgextra <- NULL; moc.call$gextra <- NULL}
  moc.call$expected <- makexpected()$FUN
  assign(data.name,data.list,envir=.GlobalEnv)
  assign(base.name,moc.call,envir=.GlobalEnv) # only assign a call, so if eval doesn't work this will be available
  constraints <- list("mu"=eval(parse(text=as.character(tclObj(mu.const.val)))),
                      "shape"=eval(parse(text=as.character(tclObj(shape.const.val)))),
                      "extra"=eval(parse(text=as.character(tclObj(extra.const.val)))),
                      "mix"=eval(parse(text=as.character(tclObj(mix.const.val)))))
  if(!is.null(unlist(constraints))) {
    tmp <- list()
    tmp$call <- eval(as.name(base.name),envir=.GlobalEnv)
    tmp$groups <- tmp$call$groups
    tmp$coefficients <- c(tmp$call$pgmu,tmp$call$pgshape,tmp$call$pgextra,tmp$call$pgmix)
    tmp$npar <- c(length(tmp$call$pgmu),length(tmp$call$pgshape),length(tmp$call$pgextra),
                  length(tmp$call$pgmix))
    class(tmp) <- "moc"
    for (j in 1:4) if((tmp$npar[j] > 0) & is.null(constraints[[j]])) constraints[[j]] <- rep(0,tmp$npar[j])
    moc.call <- with(eval(as.name(data.name),envir=.GlobalEnv),update.moc(tmp,what=constraints,evaluate=FALSE))
    assign(base.name,moc.call,envir=.GlobalEnv)
   rm(tmp)
  }
  on.exit()
  isok <- as.character(tclvalue(tkmessageBox(type="okcancel",message=paste("Run",base.name))))
  if(isok=="ok") tryCatch(assign(base.name,eval(eval(as.name(base.name))),envir=.GlobalEnv),
           error=function(e) tkmessageBox(type="ok",message=e[[1]]))
  tclObj(run.mess) <-"Ready"
  tmpfile <- tempfile("moc")
  sink(tmpfile)
  eval(substitute(print(tmp),list(tmp=as.name(base.name))),envir=.GlobalEnv)
  sink()
  tkpager(tmpfile,title="MOC",header=as.name(base.name),delete.file=TRUE)
  assign("mocid", mocid+1,envir=moc.gui.env)
  refresh()
}

refresh <- function()
{
    tclObj(data.sel) <- as.character(tclObj(datalst))[as.numeric(tkcurselection(data.listbox))+1]
    full.list <- makethelist()
    sel <- intersect(full.list$data,as.character(tclObj(data.sel)))
    tclObj(datalst) <- full.list$data
    tclObj(moclst) <- full.list$moc
#    tkdelete(resp.listbox,0.0,"end")
#    for( i in 1:ncol( get(sel))) 
#	tkinsert( resp.listbox, "end", paste(colnames( get(sel))[i]))
        tkselection.clear(data.listbox,0,"end")
    if(length(sel)>0) {
        tkselection.set(data.listbox,which(as.character(tclObj(datalst))%in%sel)-1)
        tclObj(data.sel) <-  as.character(tclObj(datalst))[as.numeric(tkcurselection(data.listbox))+1]
        tclObj(resp.sel) <-  as.character(tclObj(resplst))[as.numeric(tkcurselection(resp.listbox))+1]
        tclObj(resplst) <- colnames( eval(substitute(as.data.frame(vnam),list(vnam=as.name(sel)))))
        tkselection.clear(resp.listbox,0,"end")
        for(i in which(as.character(tclObj(resplst))%in%as.character(tclObj(resp.sel))))
            tkselection.set(resp.listbox,i-1)
        tclObj(resp.sel) <- as.character(tclObj(resplst))[as.numeric(tkcurselection(resp.listbox))+1]
    } else {tclObj(resplst)="";tclObj(resp.sel)=""}
    if(tclvalue(distfun)=="User") {
        tclObj(dist.sel) <- as.character(tclObj(funlst))[as.numeric(tkcurselection(dist.fun))+1]
        tclObj(funlst) <- full.list$fun
        tkselection.clear(dist.fun,0,"end")
        for(i in which(as.character(tclObj(funlst))%in%as.character(tclObj(dist.sel))))
            tkselection.set(dist.fun,i-1)
        tclObj(dist.sel) <- as.character(tclObj(funlst))[as.numeric(tkcurselection(dist.fun))+1]
    } else {
        tclObj(funlst) <- ""
        tclObj(dist.sel) <- ""
    }
        if(as.character(tclvalue(is.joint))=="0") {
        tkentryconfigure(dist.submenu,"MultiNormal",state="disabled")
        tkentryconfigure(dist.submenu,"Joint_Multinomial_Ind",state="disabled")
        if(as.character(tclvalue(distfun))=="MultiNormal") tclObj(distfun) <- ""
        tkconfigure(dist.spinbox,values="ALL")
    tclObj(dist.varsel) <- "ALL"}
       else {
           tkentryconfigure(dist.submenu,"MultiNormal",state="normal")
           tkentryconfigure(dist.submenu,"Joint_Multinomial_Ind",state="normal")
           tkconfigure(dist.spinbox,values=as.character(tclObj(resp.sel)))}
    tclObj(mix.sel) <- as.character(tclObj(mixfunlst))[as.numeric(tkcurselection(mix.fun))+1]
    if(tclvalue(mixfun)=="User") {
        tclObj(mixfunlst) <- full.list$fun
        tkconfigure(mix.fun,selectmode="single")
    }
    else {
        tclObj(mixfunlst) <- full.list$cov
        tkconfigure(mix.fun,selectmode="multiple")}
    tkselection.clear(mix.fun,0,"end")
    for(i in which(as.character(tclObj(mixfunlst))%in%as.character(tclObj(mix.sel))))
        tkselection.set(mix.fun,i-1)
    tclObj(mix.sel) <- as.character(tclObj(mixfunlst))[as.numeric(tkcurselection(mix.fun))+1]
        
    resetdist()
    resetgmu()
    resetgshape()
    resetgextra()
    if(tclvalue(expectedfun)=="User") {
        tclObj(expected.sel) <- as.character(tclObj(expectedfunlst))[as.numeric(tkcurselection(expected.fun))+1]
        tclObj(expectedfunlst) <- full.list$funlst
     tkselection.clear(expected.fun,0,"end")
        for(i in which(as.character(tclObj(expectedfunlst))%in%as.character(tclObj(expected.sel))))
            tkselection.set(expected.fun,i-1)
        tclObj(expected.sel) <- as.character(tclObj(expectedfunlst))[as.numeric(tkcurselection(expected.fun))+1]
    } else {
        tclObj(expectedfunlst) <- ""
        tclObj(expected.sel) <- ""
    }
}

endgui<-function() {
	tkdestroy(base)
	}

var.summ <- function()
{
    tksumm <- tktoplevel()
    dataname <- as.character(tclObj(datalst))[as.numeric(tkcurselection(data.listbox))+1]
    tkwm.title(tksumm,dataname)
    summtxt <- tktext(tksumm,wrap="word")
    ok.but <- tkbutton(tksumm,text="OK",command=function() tkdestroy(tksumm))
    if(all(is.na(dataname))) tmp <- cbind("Error"="You have to select some data first !")
    else if(all(is.na(respsel <- as.numeric(tkcurselection(resp.listbox))+1)))
        tmp <- summary(eval(as.name(dataname)))
    else    tmp <- summary(eval(as.name(dataname))[respsel])
    valtxt <- paste(apply(format(rbind(colnames(tmp),unclass(tmp))),
                           1,paste,collapse="  "),collapse="\n")
    tkinsert(summtxt,"end",valtxt)
    tkconfigure(summtxt,state="disabled")
    tkpack(summtxt,side="top",fill="both",expand=TRUE)
    tkpack(ok.but,side="top")
}

makegfun <- function(rlst,nvar,ns)
{
  if(length(rlst)==0) return(NULL)
  ngr <- as.numeric(tclObj(Groups))
  par.name <- strsplit(paste(substitute(rlst)),"\\.")[[1]][1]
  linktext <- list()
  linkarg <- list()
  functext <- list()
  curpar <- 1
  for (i in 1:ngr) {
    linktext[[i]] <- switch(rlst[[paste("Group",i,sep="")]]$link,
                       "Identity"="I",
                       "log"="exp",
                       "log-log"=c("exp","-exp"),
                       "logit"="local(function(lg) { exp(lg)/(1+exp(lg))})",
                       "exp"="log")
    type <- rlst[[paste("Group",i,sep="")]]$type
    ncov <- length(covlst <- rlst[[paste("Group",i,sep="")]]$cov)
    if(type=="User") {linkarg[[i]] <- rlst[[paste("Group",i,sep="")]]$cov
                      if(is.list(eval(parse(text=linkarg[[i]])))) linkarg[[i]] <- paste(linkarg[[i]],"[[",i,"]]",sep="")
                      linkarg[[i]] <- paste(linkarg[[i]],"(",par.name,")",sep="")}
    if(type=="Free") {if(ncov==0) {
      linkarg[[i]] <- paste("rbind(",par.name,"[",deparse(curpar:(curpar+nvar-1)),"]",")",sep="")
      curpar <- curpar+nvar }
    else {  linkarg[[i]] <- paste(par.name,"[",
                            split(curpar:(curpar+nvar*(ncov+1)-1),rep(1:(ncov+1),rep(nvar,ncov+1))),
                            "]", c("",paste("*t(matrix(",covlst,",",ns,",",nvar,"))")),sep="")
            linkarg[[i]] <- paste("t(",paste(linkarg[[i]],collapse=" + "),")")
             curpar <- curpar+nvar*(ncov+1)}
                    }
    if(type=="Constant") {if(ncov==0) {
      linkarg[[i]] <- paste("rbind(rep(",par.name,"[",curpar,"],",nvar,"))",sep="")
      curpar <- curpar+1 }
    else {  linkarg[[i]] <- paste(par.name,"[",
                                  split(curpar:(curpar+ncov+nvar-1),c(rep(1,nvar),2:(ncov+1))),
                 "]",c("",paste("*t(matrix(",covlst,",",ns,",",nvar,"))")),sep="")
            linkarg[[i]] <- paste("t(",paste(linkarg[[i]],collapse=" + "),")")
             curpar <- curpar+ncov+nvar}
                    }
    if(type=="Linear") {if(ncov==0) {
      linkarg[[i]] <- paste("rbind(rep(",par.name,"[",curpar,"],",nvar,"))",sep="")
      curpar <- curpar+1 }
    else {  linkarg[[i]] <- paste(par.name,"[",curpar:(curpar+ncov),
                 "]",c("",paste("*",covlst)),sep="")
            linkarg[[i]] <- paste(linkarg[[i]],collapse=" + ")
             curpar <- curpar+ncov+1}
                    }
    if(type=="Quadratic") {if(ncov==0) {
      linkarg[[i]] <- paste("rbind(rep(",par.name,"[",curpar,"],",nvar,"))",sep="")
      curpar <- curpar+1 }
    else {  linkarg[[i]] <- paste(par.name,"[",curpar:(curpar+2*ncov),
                 "]","*",c("1",covlst,paste(covlst,"^2",sep="")),sep="")
            linkarg[[i]] <- paste(linkarg[[i]],collapse=" + ")
             curpar <- curpar+2*ncov+1}
                    }
    if(type=="Polynomial") {if(ncov==0) {
      linkarg[[i]] <- paste("rbind(rep(",par.name,"[",curpar,"],",nvar,"))",sep="")
      curpar <- curpar+1 }
    else { polydeg <- as.numeric(rlst[[paste("Group",i,sep="")]]$degree)
      linkarg[[i]] <- paste(par.name,"[",curpar:(curpar+polydeg*ncov),
                            "]","*", c("1",sapply(covlst,paste,paste("^",1:polydeg,sep=""),sep="")),sep="")
      linkarg[[i]] <- paste(linkarg[[i]],collapse=" + ")
      curpar <- curpar+polydeg*ncov+1}
                    }
    functext[[paste("Group",i,sep="")]] <-
      eval(parse(text=paste("function(",par.name,"){",paste(linktext[[i]],"(",collapse=""),linkarg[[i]],
                   paste(rep(")",length(linktext[[i]])),collapse=""),"}")),env=.GlobalEnv)
  }
  label <- eval(substitute(as.character(tclObj(label.val)),
                             list(label.val=as.name(paste(substr(par.name,2,20),".label.val",sep="")))))
  if(length(label) != 0) attr(functext,"parameters") <- eval(parse(text=label))
  start.val <- eval(substitute(as.character(tclObj(st.val)),
                    list(st.val=as.name(paste(substr(par.name,2,20),".start.val",sep="")))))
  start.val <- eval(parse(text=start.val))
  return(list(FUN=functext,start.val=start.val))
}

makemixfun <- function()
  {ngr <- as.numeric(tclObj(Groups))
   if(ngr==1) return(list(FUN=NULL,start.val=NULL))
    type <- as.character(tclObj(mixfun))
   if(type=="User") FUN <- eval(parse(text=paste("function(pmix){ ",as.character(tclObj(mix.sel)),"(pmix) }",
                                        sep="")),env=.GlobalEnv)
   if(type=="inv.glogit") {
     covlst <-  as.character(tclObj(mix.sel))
     ncov <- length(covlst)
     largs <- sapply(1:(ngr-1),function(i)
                     paste(paste("pmix[",(1:(ncov+1))+(i-1)*(ncov+1),"]",sep=""),"*",
                           c("1",covlst),sep="",collapse="+"))
     largs <- paste("function(pmix) {t(apply(cbind(",paste(largs,collapse=","),"),1,inv.glogit))}")
     FUN <- eval(parse(text=largs),env=.GlobalEnv)
  }
   label <- as.character(tclObj(mix.label.val))
   if(length(label) != 0) attr(FUN,"parameters") <- eval(parse(text=label))
   start.val <- as.character(tclObj(mix.start.val))
   start.val <- eval(parse(text=start.val))
   return(list(FUN=FUN,start.val=start.val))                         
 }

makedensity <- function()
  { val.list <- list()
    fun.list <- character()
    density.list <- list("Normal"=function(x,mu,sd,...) {dnorm(x,mu,sd)},
                         "LogNormal"=function(x,mulog,sdlog,...) {dlnorm(x,mulog,sdlog)},
                         "CensoredNormal"=function(x,mu,sd,...) {
                           mi<-(x == min)*1
                           ma<-(x == max)*1
                           mi*pnorm((min-mu)/sd)+ma*(1-pnorm((max-mu)/sd))+
                             (1-mi-ma)*dnorm((x-mu)/sd)/sd},
                         "Student-t"=function(x,ncp,df,...) { dt(x,df,ncp) },
                         "Logistic"=function(x,mu,sc,...) {dlogis(x,mu,sc)},
                         "Gamma"=function(x,mu=1,sc,sh) {dgamma(x,sh,scal=sh)},
                         "Beta"=function(x,mu=1,sh1,sh2) {dbeta(x,sh1,sh2)},
                         "Exponential"=function(x,mu,...) {dexp(x,1/mu)},
                         "Weibull"=function(x,mu=1,sc,sh) {dweibull(x,sc,sh)},
                         "Binomial"=function(x, prob, sh = 1, size) { dbinom(x, size, prob) },
                         "Geometric"=function(x,prob,...) { dgeom(x,prob)},
                         "Hypergeometric"=function(x,m,n,k) {dhyper(x,m,n,k)},
                         "Poisson"=function(x,lda,...) {dpois(x,lda)},
                         "ZeroInflatedPoisson"=function(x,la,sh=1,zx)
                         { mix<- exp(zx)/(1+exp(zx))
                           mix*(x == 0)+(1-mix)*dpois(x,la) },
                         "NegBinomial"=function(x,mu,sh=1,size) {dnbinom(x,size,mu=mu)},
                         "Multinomial_Ind"=function(x,prob,...) {prob*x+1-x},
                         "Multinomial_Cat"=function(x,mu=1,shape=1,extra) {extra[,x]},
                         "Joint_Multinomial_Ind"=function(x,prob,...) {apply(prob*x+1-x,1,prod,na.rm=TRUE)},
                         "MultiNormal"=function(x,mu,sigma,extra)
                         { y <- (x-mu)/sigma
                           ss <- diag(rep(1,dim(x)[[2]]))
                           lind <- upper.tri(ss)
                           for ( i in 1:dim(x)[1]) {
                             ss[lind] <- extra[i,]
                             na.ind  <- is.na(y[i,])
                             y[i,!na.ind] <- t(ss)[!na.ind,!na.ind]%*%y[i,!na.ind] }
                           apply(dnorm(y)/sigma,1,prod,na.rm=TRUE) })
    if(!as.logical(tclObj(is.joint))) {
      type <- dist.rlst[["ALL"]][1]
      val.list[["joint"]] <- FALSE
      if(type=="User") { val.list[["User"]] <- as.name(dist.rlst[["ALL"]][2])}
      else if(type=="CensoredNormal") {
        minmax <-eval(parse(text=dist.rlst[["ALL"]][2]))
        val.list[[paste(type,dist.rlst[["ALL"]][2],sep=".")]] <- eval(substitute(function(x,mu,sd,...) {
                           mi<-(x == min)*1
                           ma<-(x == max)*1
                           mi*pnorm((min-mu)/sd)+ma*(1-pnorm((max-mu)/sd))+
                             (1-mi-ma)*dnorm((x-mu)/sd)/sd},
                                             list(min=minmax[1],max=minmax[2]))[-4])
      } else val.list[[type]] <- density.list[[type]]
    } else {
      val.list[["joint"]] <- TRUE
      type <- dist.rlst[as.character(tclObj(resp.sel))]
      for(i in 1:length(type)) {
        if(type[[i]][1]=="User") { val.list[["User"]] <- as.name(type[[i]][2]) }
        else if(type[[i]][1]=="CensoredNormal") {
        minmax <-eval(parse(text=type[[i]][2]))
        val.list[[paste(type[[i]],collapse=".")]] <- eval(substitute(function(x,mu,sd,...) {
                           mi<-(x == min)*1
                           ma<-(x == max)*1
                           mi*pnorm((min-mu)/sd)+ma*(1-pnorm((max-mu)/sd))+
                             (1-mi-ma)*dnorm((x-mu)/sd)/sd},
                                                     list(min=minmax[1],max=minmax[2]))[-4])
      } else {
        val.list[[type[[i]][1]]] <- density.list[[type[[i]][1]]] 
        fun.list <- c(fun.list,type[[i]][1])}
      }
      fun.index <- apply(outer(fun.list,unique(fun.list),"=="),2,function(cl) deparse(which(cl)))
      fun.list <- unique(fun.list)
      xtra.index <- fun.index
      tmp <- (fun.list%in% c("MultiNormal","Multinomial_Cat"))
      if(any(tmp)) {xtra.index[tmp] <- paste("-",deparse(unlist(lapply(xtra.index[!tmp],
                                       function(toc) eval(parse(text=toc))))),sep="")
                     tclvalue(chklen) <- "0"}
      tmp <- fun.list%in%"MultiNormal"
      if(any(tmp)) { nx <- choose(length(eval(parse(text=fun.index[tmp]))),2)
        xtra.index[tmp] <- paste(xtra.index[tmp],"][,",
                                 deparse(1:nx),sep="")}
      tmp2 <- fun.list%in%"Multinomial_Cat"
      if(any(tmp2)) {
        if(any(tmp))  xtra.index[tmp2] <- paste(xtra.index[tmp2],"][,-",
                                       deparse(1:nx),sep="")}
       fun.list <- paste(fun.list,"(","x[,",fun.index,"],mu[,",fun.index,"],sh[,",
                        fun.index,"],xtra[,",xtra.index,"])",collapse=",")
      fun.list <- paste("function(x,mu,sh,xtra) {\n apply(cbind(",fun.list,"),\n1,prod,na.rm=TRUE) }")
      val.list[["joint.dens"]] <- eval(parse(text=fun.list))
  }
    return(val.list)
  }

makexpected <- function()
  {
    type <- tclvalue(expectedfun)
    if(type=="User") {return(list("FUN"=as.name(as.character(tclObj(expected.sel)))))}
    else {return(list("FUN"=NULL))}
  }

setgmu <- function()
{
    gmu.rtmp <- gmu.rlst
    funcov <- as.character(tclObj(mu.sel))
    switch(tmp <- tclvalue(mu.fun),
           "None"={gmu.rtmp <- list();tclvalue(mu.mess) <- "  "},
           "User"={for(i in 1:as.numeric(tclObj(Groups))){
                   gmu.rtmp[[paste("Group",i,sep="")]] <- list(type="User",cov=paste(funcov),degree=NULL,link="Identity")
                 }
               tclvalue(mu.mess) <- "  "},
           {gmu.rtmp[[paste("Group",as.numeric(tclvalue(mu.Groups)),sep="")]] <- 
               list(type=tmp,cov=funcov,degree=tclvalue(degree),link=as.character(tclObj(mu.link)))
            if(tmp!="Polynomial") {tclvalue(mu.mess) <- "  "
            gmu.rtmp[[paste("Group",as.numeric(tclvalue(mu.Groups)),sep="")]]$degree <- "  "}})
    assign("gmu.rlst",gmu.rtmp,envir=moc.gui.env)
}

resetgmu <- function()
{
    gr <- paste("Group",as.numeric(tclvalue(mu.Groups)),sep="")
    if(length(gmu.rlst)==0)  {tclObj(mu.fun) <- "None"; tclObj(mufunlst) <- ""}
     else{   tclObj(mu.fun) <- paste("",gmu.rlst[[gr]][["type"]],sep="")
             tclvalue(mu.link) <- gmu.rlst[[gr]]$link
             if(paste("",gmu.rlst[[gr]][["type"]],sep="")=="User") {
               tclObj(mufunlst) <- makethelist()$funlst
               tkconfigure(mu.covlist,selectmode="single")}
             else {tclObj(mufunlst) <- makethelist()$cov
                   tkconfigure(mu.covlist,selectmode="multiple") }
             tkselection.clear(mu.covlist,0,"end")
             for(i in which(as.character(tclObj(mufunlst))%in%gmu.rlst[[gr]]$cov))
               tkselection.set(mu.covlist,paste("",i-1,sep=""))
         }
    tclvalue(mu.mess) <- gmu.rlst[[paste("Group",as.numeric(tclvalue(mu.Groups)),sep="")]]$degree
}

setgshape <- function()
{
    gshape.rtmp <- gshape.rlst
    funcov <- as.character(tclObj(shape.sel))
    switch(tmp <- tclvalue(shape.fun),
           "None"={gshape.rtmp <- list();tclvalue(shape.mess) <- "  "},
           "User"={for(i in 1:as.numeric(tclvalue(Groups))){
                   gshape.rtmp[[paste("Group",i,sep="")]] <- list(type="User",cov=paste(funcov),link="Identity")
               gshape.rtmp[[paste("Group",i,sep="")]]$degree <- NULL}
               tclvalue(shape.mess) <- "  "},
           {gshape.rtmp[[paste("Group",as.numeric(tclvalue(shape.Groups)),sep="")]] <- 
               list(type=tmp,cov=funcov,degree=tclvalue(degree),link=as.character(tclObj(shape.link)))
            if(tmp!="Polynomial") {tclvalue(shape.mess) <- "  "
                               gshape.rtmp[[paste("Group",as.numeric(tclvalue(shape.Groups)),sep="")]]$degree <- "  "}})
    assign("gshape.rlst",gshape.rtmp,envir=moc.gui.env)
}

resetgshape <- function()
{
    gr <- paste("Group",as.numeric(tclvalue(shape.Groups)),sep="")
    if(length(gshape.rlst)==0)  {tclObj(shape.fun) <- "None"; tclObj(shapefunlst) <- ""}
     else{   tclObj(shape.fun) <- paste("",gshape.rlst[[gr]][["type"]],sep="")
             tclvalue(shape.link) <- gshape.rlst[[gr]]$link
        if(paste("",gshape.rlst[[gr]][["type"]],sep="")=="User") {
            tclObj(shapefunlst) <- makethelist()$funlst
            tkconfigure(shape.covlist,selectmode="single")}
        else {tclObj(shapefunlst) <- makethelist()$cov
              tkconfigure(shape.covlist,selectmode="multiple") }
             tkselection.clear(shape.covlist,0,"end")
             for(i in which(as.character(tclObj(shapefunlst))%in%gshape.rlst[[gr]]$cov))
                 tkselection.set(shape.covlist,paste("",i-1,sep=""))
         }
             tclvalue(shape.mess) <- gshape.rlst[[paste("Group",as.numeric(tclvalue(shape.Groups)),sep="")]]$degree
}

setgextra <- function()
{
    gextra.rtmp <- gextra.rlst
    funcov <- as.character(tclObj(extra.sel))
    switch(tmp <- tclvalue(extra.fun),
           "None"={gextra.rtmp <- list();tclvalue(extra.mess) <- "  "},
           "User"={for(i in 1:as.numeric(tclvalue(Groups))){
                   gextra.rtmp[[paste("Group",i,sep="")]] <- list(type="User",cov=paste(funcov),link="Identity")
               gextra.rtmp[[paste("Group",i,sep="")]]$degree <- NULL}
               tclvalue(extra.mess) <- "  "},
           {gextra.rtmp[[paste("Group",as.numeric(tclvalue(extra.Groups)),sep="")]] <- 
               list(type=tmp,cov=funcov,degree=tclvalue(degree),link=as.character(tclObj(extra.link)))
            if(tmp!="Polynomial") {tclvalue(extra.mess) <- "  "
                               gextra.rtmp[[paste("Group",as.numeric(tclvalue(extra.Groups)),sep="")]]$degree <- "  "}})
    assign("gextra.rlst",gextra.rtmp,envir=moc.gui.env)
}

resetgextra <- function()
{
    gr <- paste("Group",as.numeric(tclvalue(extra.Groups)),sep="")
    if(length(gextra.rlst)==0)  {tclObj(extra.fun) <- "None"; tclObj(extrafunlst) <- ""}
     else{   tclObj(extra.fun) <- paste("",gextra.rlst[[gr]][["type"]],sep="")
             tclvalue(extra.link) <- gextra.rlst[[gr]]$link
        if(paste("",gextra.rlst[[gr]][["type"]],sep="")=="User") {
            tclObj(extrafunlst) <- makethelist()$funlst
            tkconfigure(extra.covlist,selectmode="single")}
        else {tclObj(extrafunlst) <- makethelist()$cov
              tkconfigure(extra.covlist,selectmode="multiple") }
             tkselection.clear(extra.covlist,0,"end")
             for(i in which(as.character(tclObj(extrafunlst))%in%gextra.rlst[[gr]]$cov))
                 tkselection.set(extra.covlist,paste("",i-1,sep=""))
         }
             tclvalue(extra.mess) <- gextra.rlst[[paste("Group",as.numeric(tclvalue(extra.Groups)),sep="")]]$degree
}


getpolydegree <- function()
{
    poly.win <- tktoplevel(parent=base)
    tkconfigure(poly.win,bd=0)
    tkwm.geometry(poly.win,paste("",tclvalue(tkwinfo("rootx",shape.fun.menu)),tclvalue(tkwinfo("rooty",shape.fun.menu)),sep="+"))
    tkfocus(poly.win)
    tkwm.title(poly.win,"Polynomial degree")
    degree.win <- tkwidget(poly.win,"spinbox",from=0,to=20,width=3,format="%3.0f",state="readonly",textvariable=degree,
                           readonlybackground="#FFFFFF")
    ok.but <- tkbutton(poly.win,text="OK",command=function() {
        tkgrab.release(poly.win)
        tkdestroy(poly.win)
        tkfocus(base)
        })
    tkpack(degree.win,ok.but,side="left",fill="x")
    tkgrab.set(poly.win)
    tkwait.window(poly.win)
}
    

tkadd(moc.menu,"command",label="Run",underline="0",command=submit)
tkadd(moc.menu,"command",label="Refresh",underline="0",command=refresh)
tkadd(moc.menu,"command",label="Summary",underline="0",command=var.summ)
tkadd(moc.menu,"command",label="Quit",underline="0",command=endgui)
help.menu <- tkmenu(moc.menu,tearoff="0")
tkadd(moc.menu,"cascade",label="Help",underline="0",menu=help.menu)
tkadd(help.menu,"command",label="MOC",underline="0",command=function() tkpager( system.file("help","moc",package="moc"),title="Help",header="moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Print Methods",underline="0",command=function() tkpager( system.file("help","print.moc",package="moc"),title="Help",header="print.moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Plot Methods",underline="0",command=function() tkpager( system.file("help","plot.moc",package="moc"),title="Help",header="plot.moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Information Criterions",underline="0",command=function() tkpager( system.file("help","AIC.moc",package="moc"),title="Help",header="AIC.moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Residuals & Diagnostics",underline="0",command=function() tkpager( system.file("help","residuals.moc",package="moc"),title="Help",header="residuals.moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Profiling & Density",underline="0",command=function() tkpager( system.file("help","confint.moc",package="moc"),title="Help",header="confint.moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Utilities",underline="0",command=function() tkpager( system.file("help","utils.moc",package="moc"),title="Help",header="utils.moc",delete.file=FALSE))
tkadd(help.menu,"command",label="Readme",underline="0",command=function() tkpager( system.file("GUI","moc-gui.Readme",package="moc"),title="Help",header="moc-gui.Readme",delete.file=FALSE))
tkinsert(help.menu,7,"separator")

## Choose which data file to load
load.file <- function() {
                        fname <- (tclvalue(tkgetOpenFile(filetypes="{{RData} {.RData .Rda .rdata .rda}} {{All files} {*}}")))
                        if(fname!="") {if(regexpr(".rda",tolower(fname))>0) try(load(fname,env=.GlobalEnv))
                        else assign(fname,import.data(fname),env=.GlobalEnv)}
                        oldfname <- tclvalue(tkget(data.file.name))
                        tkdelete(data.file.name,0.0,"end")
                        tkinsert(data.file.name,0.0,paste(fname))
                        refresh()
                        }

data.file.frm <- tkframe(data.frm, borderwidth=2)
data.file.name <- tkentry(data.file.frm,background="#FFFFFF")
tkbind(data.file.name,"<Return>",function() {
                        fname <- as.character(tclvalue(tkget(data.file.name)))
                        if(fname!="") {if(regexpr(".rda",tolower(fname))>0) try(load(fname,env=.GlobalEnv))
                        else assign(fname,import.data(fname),env=.GlobalEnv)}
                        oldfname <- tclvalue(tkget(data.file.name))
                        tkdelete(data.file.name,0.0,"end")
                        tkinsert(data.file.name,0.0,paste(fname))
                        refresh()
                        })
data.file <- tkbutton(data.file.frm,
                    text="Load data",
                    relief="groove", borderwidth="1",
                    command= load.file)

choose.wtvar <- function() {
    tmp <- tktoplevel(parent=base)
    tkwm.geometry(tmp,paste("",tclvalue(tkwinfo("rootx",weight.but)),tclvalue(tkwinfo("rooty",weight.but)),sep="+"))
    tkwm.title(tmp,"Weight")
    tkfocus(tmp)
    wtlst <- tclVar()
    wt.listbox <- tklistbox(tmp,
                            listvariable=wtlst,
                            yscrollcommand=function(...) tkset(wt.scroll,...),
                            selectmode="single",
                            ##width=20,
                            height=10,
                            exportselection=0,background="#FFFFFF")

    wt.scroll <- tkscrollbar( tmp, orient="vert",
                             command=function(...)tkyview(wt.listbox,...))
    tclObj(wtlst) <- makethelist()$cov
    tkpack(wt.listbox,wt.scroll,side="left",fill="both",expand="1")
       ok.but <- tkbutton(tmp,text="OK",command=function() {
           cursel <- as.character(tclObj(wtlst))[as.numeric(tclvalue(tkcurselection(wt.listbox)))+1]
           if(is.na(cursel)) cursel=""
           tclvalue(wtvar) <- paste("Weight",cursel)
           tkgrab.release(tmp)
           tkdestroy(tmp)
           tkfocus(base)
       })
        tkpack(ok.but,side="bottom")
     tkgrab.set(tmp)
 }

weight.but <- tkbutton(data.file.frm,textvariable=wtvar,command=choose.wtvar)
scale.weight.but <- tkcheckbutton(data.file.frm,text="scaled",variable=scale.weight)

data.list.frm <- tkwidget(data.frm, "labelframe",text="Data Object:", borderwidth=2)
data.listbox <- tklistbox(data.list.frm,
                          listvariable=datalst,
			yscrollcommand=function(...) tkset(data.scroll,...),
			selectmode="single",
#			width=20,
			height=4,
			exportselection=0,background="#FFFFFF")

data.scroll <- tkscrollbar( data.list.frm, orient="vert",
			command=function(...)tkyview(data.listbox,...))

moc.list.frm <- tkwidget(data.frm, "labelframe",#text="Moc Object:",
                         borderwidth=2)
mocplot.menubut <- tkmenubutton(moc.list.frm,text="Moc Object:",relief="raised",borderwidth=4)
mocplot.menu <- tkmenu(mocplot.menubut,tearoff="0")
tkadd(mocplot.menu,"command",label="Print",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       tmpfile <- tempfile("moc")
       sink(tmpfile)
       if (!is.na(mocsel)) eval(substitute(print(tmp),list(tmp=as.name(mocsel))))
       sink()
       tkpager(tmpfile,title="MOC",header=as.name(mocsel),delete.file=TRUE)
     })
mocplotchx.menu <- tkmenu(mocplot.menu,tearoff="0")
tkadd(mocplotchx.menu,"command",label="Plot",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(plot(tmp),list(tmp=as.name(mocsel))))})
tkadd(mocplotchx.menu,"command",label="Plot (scaled)",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(plot(tmp,scale=TRUE),list(tmp=as.name(mocsel))))})
tkadd(mocplot.menu,"cascade",label="Plot",menu=mocplotchx.menu)
mocplotres.menu <- tkmenu(mocplot.menu,tearoff="0")
tkadd(mocplotres.menu,"command",label="Deviance",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(plot(residuals(tmp,type="deviance")),list(tmp=as.name(mocsel)))) })
tkadd(mocplotres.menu,"command",label="Gradient",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel))  eval(substitute(plot(residuals(tmp,type="gradient")),list(tmp=as.name(mocsel))))})
tkadd(mocplotres.menu,"command",label="Mixture",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(plot(residuals(tmp,type="mixture")),list(tmp=as.name(mocsel)))) })
tkadd(mocplotres.menu,"command",label="Response",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(plot(residuals(tmp,type="response")),list(tmp=as.name(mocsel)))) })
tkadd(mocplot.menu,"cascade",label="Plot residuals",menu=mocplotres.menu)
mocplotprof.menu <- tkmenu(mocplot.menu,tearoff="0")
tkadd(mocplotprof.menu,"command",label="Subject",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) profilesplot(get(mocsel),type="subject")})
tkadd(mocplotprof.menu,"command",label="Variable",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) profilesplot(get(mocsel),type="variable")})
tkadd(mocplotprof.menu,"command",label="Posterior",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) profilesplot(get(mocsel),type="posterior")})

mocplotdens.menu <- tkmenu(mocplot.menu,tearoff="0")
tkadd(mocplotdens.menu,"command",label="Density",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(density.moc(tmp,var=1:tmp$nvar,plot="density",type="s"),list(tmp=as.name(mocsel)))) })
tkadd(mocplotdens.menu,"command",label="PP-plot",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(density.moc(tmp,var=1:tmp$nvar,plot="pp-plot",type="l"),list(tmp=as.name(mocsel)))) })
tkadd(mocplotdens.menu,"command",label="PQ-plot",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(density.moc(tmp,var=1:tmp$nvar,plot="pq-plot",type="l"),list(tmp=as.name(mocsel)))) })
tkadd(mocplot.menu,"cascade",label="Plot Profiles",menu=mocplotprof.menu)
tkadd(mocplot.menu,"cascade",label="Density Plot",menu=mocplotdens.menu)
tkadd(mocplot.menu,"command",label="Plot entropy",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(entropyplot(tmp),list(tmp=as.name(mocsel))))})
tkadd(mocplot.menu,"command",label="Delete",command=function()
      {mocsel <- as.character(tclObj(moclst))[as.numeric(tclvalue(tkcurselection(moc.listbox)))+1]
       if (!is.na(mocsel)) eval(substitute(remove(tmp,envir=.GlobalEnv),list(tmp=mocsel)))
       tclObj(moclst) <- makethelist()$moc})

tkinsert(mocplot.menu,1,"separator")
tkinsert(mocplot.menu,7,"separator")

tkconfigure(mocplot.menubut,menu=mocplot.menu)
tkconfigure(moc.list.frm,labelwidget=mocplot.menubut)

moc.listbox <- tklistbox(moc.list.frm,
                          listvariable=moclst,
			yscrollcommand=function(...) tkset(moc.scroll,...),
			selectmode="single",
#			width=20,
			height=4,
			exportselection=0,background="#FFFFFF")

moc.scroll <- tkscrollbar( moc.list.frm, orient="vert",
			command=function(...)tkyview(moc.listbox,...))



# initialize variables for data list.
# 'temp' is list of everything in global environment.
# 'full.list' will be list of all vector, matrix or data.frame objects in '.GlobalEnv'.
makethelist <- function()
{
    temp <- c(ls( name=".GlobalEnv"))
    full.list <- list()
    for( i in 1:length( temp)) {
        if((isvec <- is.vector(get(temp[i]),mode="numeric")) || is.matrix(get(temp[i])) || is.data.frame(get(temp[i])))
            {full.list$data <- c( full.list$data, temp[i])
            full.list$cov <- c( full.list$cov, temp[i],if(!isvec) paste(temp[i],
                    eval(substitute(colnames(as.data.frame(vnam)),list(vnam=as.name(temp[i])))),sep="$"))}
          if(is.function(get(temp[i])))  full.list$fun <- c(full.list$fun,temp[i])
        if(is.list(get(temp[i]))) if(all(sapply(get(temp[i]),is.function))) full.list$funlst <- c(full.list$funlst,temp[i])
        if(data.class(get(temp[i]))=="moc") full.list$moc <- c(full.list$moc,temp[i])
} # end of for i loop
    full.list
}

tclObj(datalst) <- makethelist()$data
tclObj(moclst) <- makethelist()$moc

resp.list.frm <- tkwidget(data.frm,"labelframe",text="Response:", borderwidth=2)
resp.listbox <-
	tklistbox(resp.list.frm,yscrollcommand=function(...)tkset(resp.scroll,...),
                  listvariable=resplst,
			selectmode="multiple",
#                  width=15,
                  height=4,exportselection=0,background="#FFFFFF")
resp.scroll <- tkscrollbar(resp.list.frm,orient="vert",
			command=function(...)tkyview(resp.listbox,...))

mocimage <- tkimage.create("photo","moc-logo",
                           file=system.file("GUI","moc-gui.ppm",package="moc"))
tkpack(tkbutton(data.frm,image="moc-logo",command=submit),side="left",expand="1")
tkpack(data.file,data.file.name,side="top",expand="1")
tkpack(scale.weight.but,weight.but,side="left",anchor="s",pady=15)
tkpack(data.file.frm,side="left",expand="1")
tkpack( data.listbox, side="left",fill="y",expand="1")
tkpack( data.scroll, side="right", fill="y")
tkpack( data.list.frm, after=data.file.frm,side="right",padx=4,fill="y",expand="1")
tkpack(resp.listbox,side="left",fill="y",expand="1")
tkpack(resp.scroll,side="right",fill="y")
tkpack( resp.list.frm, before=data.list.frm,side="right",fill="y",expand="1")

tkpack(moc.listbox,side="left",fill="y",expand="1")
tkpack(moc.scroll,side="right",fill="y")
tkpack( moc.list.frm,before= resp.list.frm,side="right",fill="y",expand="1")
tkpack( data.frm,expand="1",fill="x")


# place binding on data.listbox to reflect the chosen data from the list.
tkbind( data.listbox, "<Button-1>", "")
tkbind( data.listbox, "<ButtonRelease-1>", refresh)


# frame for distribution and mixture

getRvalue <- function(title="?",parent)
{
    Rval.win <- tktoplevel(parent=parent)
    tkconfigure(Rval.win,bd=0)
    tkwm.geometry(Rval.win,paste("",tclvalue(tkwinfo("rootx",distmix.entry)),tclvalue(tkwinfo("rooty",distmix.entry)),sep="+"))
    tkfocus(Rval.win)
    tkwm.title(Rval.win,title)
    val.win <- tktext(Rval.win,width=35,height=1,background="#FFFFFF")
    ok.but <- tkbutton(Rval.win,text="OK",command=function() {
        assign("R.val", eval(parse(text=tclvalue(tkget(val.win,"1.0","end")))),envir=moc.gui.env)
        tkgrab.release(Rval.win)
        tkdestroy(Rval.win)
        tkfocus(base)
        })
    tkpack(val.win,ok.but,side="left",fill="x")
    tkgrab.set(Rval.win)
    tkfocus(val.win)
    tkwait.window(Rval.win)
}
    
setdist <- function()
{ tmp.rlst <- dist.rlst
  tmp <- as.character(tclObj(distfun))
  varsel <- as.character(tclObj(dist.varsel))
  if(varsel%in%c("ALL",as.character(tclObj(resplst)))){
  if(tmp=="User") {tmp.rlst <- list()
                   tclvalue(dist.mess) <- "  "
                   for (i in c("ALL",as.character(tclObj(resplst)))) tmp.rlst[[i]] <-
      c("User",as.character(tclObj(dist.sel)))}
  else if(tmp=="CensoredNormal"){getRvalue("c(min,max) for censored normal",dist.menu)
                                 tmp <- deparse(R.val)
                                 tmp.rlst[[varsel]] <- c("CensoredNormal",tmp)
                                  tclvalue(dist.mess) <- tmp
                                 tclObj(funlst) <- ""}
  else {tmp.rlst[[varsel]] <- tmp;tclObj(funlst) <- ""; tclvalue(dist.mess) <- "  "}
} else tmp.rlst[[varsel]] <- NULL
  assign("dist.rlst",tmp.rlst,envir=moc.gui.env)
}

resetdist <- function()
{
   varsel <- as.character(tclObj(dist.varsel))
   tclObj(distfun) <- paste("",dist.rlst[[varsel]][1],sep="")
   if(paste("",dist.rlst[[varsel]][1],sep="")=="User") {
       tkselection.clear(dist.fun,0,"end")
       if(length(i <- which(as.character(tclObj(funlst))%in% dist.rlst[[varsel]][2]))>0) tkselection.set(dist.fun,i-1)
   } else tclObj(funlst) <- ""
   tmp <- dist.rlst[[varsel]][2]
   tclvalue(dist.mess) <- tmp
}
    
dist.frm <- tkwidget(distmix.frm,"labelframe",text="Density:",borderwidth=2,relief="raised")
dist.joint <- tkcheckbutton(dist.frm,text="joint",variable=is.joint,command=refresh)
dist.spinbox <- tkwidget(dist.frm,"spinbox",values="ALL",textvariable=dist.varsel,state="readonly",command=resetdist,
                         readonlybackground="#FFFFFF")
dist.menu <- tkwidget( dist.frm, "tk_optionMenu",distfun,
                      "User","Normal","LogNormal","CensoredNormal","Student-t","Logistic","Gamma","Beta",
                      "Exponential","Weibull",
                      "Binomial","Geometric","Hypergeometric",
                      "Poisson","ZeroInflatedPoisson","NegBinomial",
                      "Multinomial_Ind","Multinomial_Cat",
                      "MultiNormal","Joint_Multinomial_Ind")
tkconfigure(dist.menu,width=18)

dist.submenu <- .Tk.newwin(tkcget(dist.menu,"-menu"))
tkinsert(dist.submenu,1,"separator")
tkinsert(dist.submenu,11,"separator")
tkinsert(dist.submenu,20,"separator")
tkentryconfigure(dist.submenu,"User",command=function() {tclObj(funlst) <- makethelist()$fun;setdist()})
for(i in c("Normal","LogNormal","CensoredNormal","Student-t","Logistic","Gamma","Beta",
           "Exponential","Weibull",
           "Binomial","Geometric","Hypergeometric",
           "Poisson","ZeroInflatedPoisson","NegBinomial",
           "Multinomial_Ind","Multinomial_Cat",
           "MultiNormal","Joint_Multinomial_Ind"))
tkentryconfigure(dist.submenu,paste(i),command=setdist)
tkentryconfigure(dist.submenu,"MultiNormal",state="disabled")
tkentryconfigure(dist.submenu,"Joint_Multinomial_Ind",state="disabled")
dist.fun <- tklistbox( dist.frm,
                      listvariable=funlst,
		yscrollcommand=function(...) tkset( dist.scroll, ...),
		selectmode="single",
#		width=15,
		height=4,
		exportselection=0,background="#FFFFFF")

dist.scroll <- tkscrollbar( dist.frm, orient="vert",
		command=function(...) tkyview( dist.fun, ...))
dist.message <- tkmessage(parent=dist.frm,textvariable=dist.mess,width=210,justify="center",anchor="w")


mix.frm <- tkwidget(distmix.frm,"labelframe",text="Mixture:",borderwidth=2,relief="raised")
mix.selgr <- tkframe(mix.frm,relief="flat")
mix.group <- tkwidget(mix.selgr,"spinbox",from=1,to=20,width=3,format="%3.0f",state="readonly",
                            readonlybackground="#FFFFFF",textvariable=Groups,
                      command=function(){
                          gr <- as.numeric(tclvalue(Groups))
                          tkconfigure(mu.group,to=gr)
                          tkconfigure(shape.group,to=gr)
                          tkconfigure(extra.group,to=gr)
                          tclObj(mu.Groups) <- gr;tclObj(shape.Groups) <- gr;tclObj(extra.Groups) <- gr                    
                      })
mix.menu <- tkwidget( mix.selgr,"tk_optionMenu",mixfun,"User","inv.glogit")
tkconfigure(mix.menu,width=10)
mix.submenu <- .Tk.newwin(tkcget(mix.menu,"-menu"))
tkinsert(mix.submenu,1,"separator")
tkentryconfigure(mix.submenu,"User",command=function() {
    tclObj(mixfunlst) <-  makethelist()$fun
    tkconfigure(mix.fun,selectmode="single")})
tkentryconfigure(mix.submenu,"inv.glogit",command=function() {
    tclObj(mixfunlst) <-  makethelist()$cov
    tkconfigure(mix.fun,selectmode="multiple")})

mix.funscroll <- tkframe(mix.frm)
mix.fun <- tklistbox( mix.funscroll,
		yscrollcommand=function(...) tkset( mix.scroll, ...),
                     listvariable=mixfunlst,
		selectmode="single",
#		width=15,
		height=4,
		exportselection=0,background="#FFFFFF")

mix.scroll <- tkscrollbar( mix.funscroll, orient="vert",
		command=function(...) tkyview( mix.fun, ...))
tkpack(mix.fun,side="left",fill="both",expand="1")
tkpack(mix.scroll,side="left",fill="y")

tkbind(mix.fun,"<ButtonRelease-1>",function() tclObj(mix.sel) <-
       as.character(tclObj(mixfunlst))[as.numeric(tkcurselection(mix.fun))+1])

mix.start <- tkentry(mix.frm,textvariable=mix.start.val,background="#FFFFFF")
mix.const <- tkentry(mix.frm,textvariable=mix.const.val,background="#FFFFFF")
mix.label <- tkentry(mix.frm,textvariable=mix.label.val,background="#FFFFFF")

distmix.label <- tkwidget(distmix.frm,"labelframe",
                          text="R commands: C-RET eval window, MOUSE eval selection",
                          borderwidth=2, relief="ridge")
distmix.entry <- tktext(distmix.label,width=50,height=12,background="#FFFFFF")
tkpack(distmix.entry,side="top",fill="both",expand="1")

tkbind(distmix.entry,"<Button-1>",function() tkfocus(distmix.entry))
tkbind(distmix.entry,"<ButtonRelease-1>",
       function() {selection <- strsplit(tclvalue(tktag.ranges(distmix.entry,"sel"))," ")[[1]]
                   if(!is.na(selection[1])) eval(parse(text=tclvalue(tkget(distmix.entry,selection[1],selection[2]))),env=.GlobalEnv)
                                                     refresh()})
tkbind(distmix.entry,"<Control-Return>",function() { eval(parse(text=tclvalue(tkget(distmix.entry,"1.0","end"))),env=.GlobalEnv)
                                                     refresh()})

tkpack( tklabel(mix.selgr,text="Groups"),mix.group,side="left")
tkpack(mix.menu,side="left",after=mix.group)
tkpack(mix.selgr,side="top")
tkpack( mix.funscroll,side="top",fill="both",expand="1")
tkpack(mix.label,tklabel(mix.frm,text="Labels:"),fill="x",side="bottom")
tkpack(mix.const,tklabel(mix.frm,text="Constraints:"),fill="x",side="bottom")
tkpack(mix.start,tklabel(mix.frm,text="Starting values:"),fill="x",side="bottom")

tkpack(dist.joint,dist.spinbox)
tkpack(dist.menu,dist.message,anchor="s",after=dist.spinbox)
tkpack( dist.fun, side="left",fill="both",expand="1")
tkpack( dist.scroll, side="left",fill="y",after=dist.fun)

tkbind(dist.fun,"<ButtonRelease-1>",function() {tclObj(dist.sel) <-
       as.character(tclObj(funlst))[as.numeric(tkcurselection(dist.fun))+1];setdist()})
tkpack(mix.frm,distmix.label,dist.frm,fill="both",side="left",expand="1")
tkpack( distmix.frm, side="top",fill="both",expand="1")


# top frame for response variable

top.r <- tkframe(top.frm,borderwidth=2)
param.frm <- tkframe( top.r, borderwidth=2, relief="groove")
mu.frm <- tkwidget( param.frm, "labelframe",text="Mu parameter:",borderwidth=2, relief="groove")
shape.frm <- tkwidget( param.frm, "labelframe",text="Shape parameter:",borderwidth=2, relief="groove")
extra.frm <- tkwidget( param.frm, "labelframe",text="Extra parameter:", borderwidth=2, relief="groove")


tkpack(param.frm,side="left",fill="both")
tkpack(top.r,side="left",fill="both",expand="1")

# place binding on resp.listbox to eliminate the response from the 
# lists of covs.
tkbind(resp.listbox,"<ButtonRelease-1>",function() {
    tclObj(resp.sel) <-
       as.character(tclObj(resplst))[as.numeric(tkcurselection(resp.listbox))+1]
    if(as.character(tclvalue(is.joint))=="0") {
        tkentryconfigure(dist.submenu,"MultiNormal",state="disabled")
        tkentryconfigure(dist.submenu,"Joint_Multinomial_Ind",state="disabled")
        tkconfigure(dist.spinbox,values="ALL")
        tclObj(dist.varsel) <- "ALL"}
    else {
           tkentryconfigure(dist.submenu,"MultiNormal",state="normal")
           tkentryconfigure(dist.submenu,"Joint_Multinomial_Ind",state="normal")
          tkconfigure(dist.spinbox,values=as.character(tclObj(resp.sel)))}})

# mu frame

mu.l <- tkframe( mu.frm, borderwidth=2, relief="flat")
mu.r <- tkframe( mu.frm, borderwidth=2, relief="flat")
mu.covlist <- tklistbox( mu.l,listvariable=mufunlst,
		yscrollcommand=function(...) tkset( mu.scroll, ...),
		selectmode="multiple",
#		width=15,
		height=4,
		exportselection=0,background="#FFFFFF")

mu.scroll <- tkscrollbar( mu.l, orient="vert",
		command=function(...) tkyview( mu.covlist, ...))

tkbind(mu.covlist,"<ButtonRelease-1>",function() {
    tclObj(mu.sel) <- as.character(tclObj(mufunlst))[as.numeric(tkcurselection(mu.covlist))+1]
    setgmu()})

mu.link.menu <- tkwidget( mu.r, "tk_optionMenu",mu.link, "Identity", "log" ,"log-log", "logit","exp")
tkconfigure(mu.link.menu,width=8)
tkbind(mu.link.menu,"<ButtonRelease-1>",function() {setgmu()})
mu.fun.menu <- tkwidget( mu.r, "tk_optionMenu",mu.fun,"None","Free","Constant","Linear","Quadratic","Polynomial","User")
tkconfigure(mu.fun.menu,width=10)
mu.fun.submenu <- .Tk.newwin(tkcget(mu.fun.menu,"-menu"))
tkinsert(mu.fun.submenu,1,"separator")
tkinsert(mu.fun.submenu,7,"separator")

tkentryconfigure(mu.fun.submenu,"User",command=function() {
    tclObj(mufunlst) <- makethelist()$funlst
    tkconfigure(mu.covlist,selectmode="single")
    setgmu()})
tkentryconfigure(mu.fun.submenu,"None",command=function() {tclObj(mufunlst) <- "";setgmu()})

polydegree.mu <- function() {
    tmp.rlst <- gmu.rlst
    getpolydegree()
    tmp.rlst[[paste("Group",as.numeric(tclvalue(mu.Groups)),sep="")]]$degree <- tclvalue(degree)
    assign("gmu.rlst",tmp.rlst,envir=moc.gui.env)
    tclObj(mufunlst) <- makethelist()$cov
    tkconfigure(mu.covlist,selectmode="multiple")
    tclvalue(mu.mess) <- tclvalue(degree)
    setgmu()
}

tkentryconfigure(mu.fun.submenu,"Polynomial",command=polydegree.mu)
for(i in c("Free","Constant","Linear","Quadratic"))
    tkentryconfigure(mu.fun.submenu,i,command=function() {
        tclObj(mufunlst) <- makethelist()$cov
        tkconfigure(mu.covlist,selectmode="multiple")
        setgmu()})


mu.group <- tkwidget(mu.r,"spinbox",from=1,to=as.numeric(tclvalue(Groups)),width=3,
                     format="%3.0f",state="readonly",textvariable=mu.Groups, readonlybackground="#FFFFFF",
                     command=resetgmu)
mu.message <- tkmessage(parent=mu.r,textvariable=mu.mess,width=30,justify="left",anchor="w")

mu.start <- tkentry(mu.frm,textvariable=mu.start.val,background="#FFFFFF")
mu.label <- tkentry(mu.frm,textvariable=mu.label.val,background="#FFFFFF")
mu.const <- tkentry(mu.frm,textvariable=mu.const.val,background="#FFFFFF")

tkpack( mu.covlist, side="left",fill="both",expand="1")
tkpack( mu.scroll, side="left", fill="y")

tkpack( tklabel( mu.r, text="Link:"), side="left")
tkpack(mu.group,side="left")
tkpack(mu.link.menu,anchor="nw",padx=30)
tkpack(mu.fun.menu, anchor="s",side="left",fill="none")
tkpack(mu.message,anchor="w",after=mu.fun.menu,fill="none")

tkpack( mu.r, side="top",fill="x",expand="1")
tkpack( mu.l, side="top", after=mu.r,fill="both",expand="1")
tkpack(mu.label,tklabel(mu.frm,text="Labels:"),fill="x",side="bottom",expand="1")
tkpack(mu.const,tklabel(mu.frm,text="Constraints:"),fill="x",side="bottom",expand="1")
tkpack(mu.start,tklabel(mu.frm,text="Starting values:"),fill="x",side="bottom",expand="1")

# shape frame
 
shape.l <- tkframe( shape.frm,borderwidth=2)
shape.r <- tkframe( shape.frm,borderwidth=2)
shape.covlist <-
	tklistbox( shape.l,listvariable=shapefunlst,
                  yscrollcommand=function(...)tkset(shape.scroll,...),
		selectmode="multiple",
#                  width=15,
                  height=4,exportselection=0,background="#FFFFFF")
shape.scroll <- tkscrollbar( shape.l,orient="vert",
		command=function(...)tkyview(shape.covlist,...))


tkbind(shape.covlist,"<ButtonRelease-1>",function() {
    tclObj(shape.sel) <- as.character(tclObj(shapefunlst))[as.numeric(tkcurselection(shape.covlist))+1]
    setgshape()})

shape.link.menu <- tkwidget( shape.r, "tk_optionMenu",shape.link, "Identity", "log" ,"log-log", "logit","exp")
tkconfigure(shape.link.menu,width=8)
tkbind(shape.link.menu,"<ButtonRelease-1>",function() {setgshape()})
shape.fun.menu <- tkwidget( shape.r, "tk_optionMenu",shape.fun,"None","Free","Constant","Linear","Quadratic","Polynomial","User")
tkconfigure(shape.fun.menu,width=10)
shape.fun.submenu <- .Tk.newwin(tkcget(shape.fun.menu,"-menu"))
tkinsert(shape.fun.submenu,1,"separator")
tkinsert(shape.fun.submenu,7,"separator")
tkentryconfigure(shape.fun.submenu,"User",command=function() {
    tclObj(shapefunlst) <- makethelist()$funlst
    tkconfigure(shape.covlist,selectmode="single")
    setgshape()})
tkentryconfigure(shape.fun.submenu,"None",command=function() {tclObj(shapefunlst) <- "";setgshape()})

polydegree.shape <- function() {
    tmp.rlst <- gshape.rlst
    getpolydegree()
    tmp.rlst[[paste("Group",as.numeric(tclvalue(shape.Groups)),sep="")]]$degree <- tclvalue(degree)
    assign("gshape.rlst",tmp.rlst,envir=moc.gui.env)
    tclObj(shapefunlst) <- makethelist()$cov
    tkconfigure(shape.covlist,selectmode="multiple")
    tclvalue(shape.mess) <- tclvalue(degree)
    setgshape()
}

tkentryconfigure(shape.fun.submenu,"Polynomial",command=polydegree.shape)
for(i in c("Free","Constant","Linear","Quadratic"))
    tkentryconfigure(shape.fun.submenu,i,command=function() {
        tclObj(shapefunlst) <- makethelist()$cov
        tkconfigure(shape.covlist,selectmode="multiple")
        setgshape()})

shape.group <- tkwidget(shape.r,"spinbox",from=1,to=as.numeric(tclvalue(Groups)),width=3,
                         readonlybackground="#FFFFFF",format="%3.0f",state="readonly",
                        textvariable=shape.Groups,command=resetgshape)
shape.message <- tkmessage(parent=shape.r,textvariable=shape.mess,width=30,justify="left",anchor="w")

shape.start <- tkentry(shape.frm,textvariable=shape.start.val,background="#FFFFFF")
shape.const <- tkentry(shape.frm,textvariable=shape.const.val,background="#FFFFFF")
shape.label <- tkentry(shape.frm,textvariable=shape.label.val,background="#FFFFFF")

tkpack( shape.covlist, side="left",fill="both",expand="1")
tkpack( shape.scroll, side="right", fill="y")

tkpack( tklabel( shape.r, text="Link:"), side="left")
tkpack(shape.group,side="left")
tkpack(shape.link.menu,anchor="nw",padx=30)
tkpack(shape.fun.menu, anchor="s",side="left",fill="none")
tkpack(shape.message,anchor="w",after=shape.fun.menu,fill="none")

tkpack( shape.r, side="top",fill="x",expand="1")
tkpack( shape.l, side="top", after=shape.r,fill="x",expand="1")
tkpack(shape.label,tklabel(shape.frm,text="Labels:"),fill="x",side="bottom")
tkpack(shape.const,tklabel(shape.frm,text="Constraints:"),fill="x",side="bottom")
tkpack(shape.start,tklabel(shape.frm,text="Starting values:"),fill="x",side="bottom")



# extra frame
 
extra.l <- tkframe(extra.frm,borderwidth=2)
extra.r <- tkframe(extra.frm,borderwidth=2)
extra.covlist <-
	tklistbox(extra.l,listvariable=extrafunlst,
                  yscrollcommand=function(...)tkset(extra.scroll,...),
		selectmode="multiple",
                  #width=15,
                  height=4,exportselection=0,background="#FFFFFF")
extra.scroll <- tkscrollbar(extra.l,orient="vert",
		command=function(...)tkyview(extra.covlist,...))


tkbind(extra.covlist,"<ButtonRelease-1>",function() {
    tclObj(extra.sel) <- as.character(tclObj(extrafunlst))[as.numeric(tkcurselection(extra.covlist))+1]
    setgextra()})

extra.link.menu <- tkwidget( extra.r, "tk_optionMenu",extra.link, "Identity", "log" ,"log-log", "logit","exp")
tkconfigure(extra.link.menu,width=8)
tkbind(extra.link.menu,"<ButtonRelease-1>",function() {setgextra()})
extra.fun.menu <- tkwidget( extra.r, "tk_optionMenu",extra.fun,"None","Free","Constant","Linear","Quadratic","Polynomial","User")
tkconfigure(extra.fun.menu,width=10)
extra.fun.submenu <- .Tk.newwin(tkcget(extra.fun.menu,"-menu"))
tkinsert(extra.fun.submenu,1,"separator")
tkinsert(extra.fun.submenu,7,"separator")
tkentryconfigure(extra.fun.submenu,"User",command=function() {
    tclObj(extrafunlst) <- makethelist()$funlst
    tkconfigure(extra.covlist,selectmode="single")
    setgextra()})
tkentryconfigure(extra.fun.submenu,"None",command=function() {tclObj(extrafunlst) <- "";setgextra()})

polydegree.extra <- function() {
    tmp.rlst <- gextra.rlst
    getpolydegree()
    tmp.rlst[[paste("Group",as.numeric(tclvalue(extra.Groups)),sep="")]]$degree <- tclvalue(degree)
    assign("gextra.rlst",tmp.rlst,envir=moc.gui.env)
    tclObj(extrafunlst) <- makethelist()$cov
    tkconfigure(extra.covlist,selectmode="multiple")
    tclvalue(extra.mess) <- tclvalue(degree)
    setgextra()
}

tkentryconfigure(extra.fun.submenu,"Polynomial",command=polydegree.extra)
for(i in c("Free","Constant","Linear","Quadratic"))
    tkentryconfigure(extra.fun.submenu,i,command=function() {
        tclObj(extrafunlst) <- makethelist()$cov
        tkconfigure(extra.covlist,selectmode="multiple")
        setgextra()})
        
extra.group <- tkwidget(extra.r,"spinbox",from=1,to=as.numeric(tclvalue(Groups)),width=3,
                         readonlybackground="#FFFFFF",format="%3.0f",state="readonly",
                        textvariable=extra.Groups,command=resetgextra)
extra.message <- tkmessage(parent=extra.r,textvariable=extra.mess,width=30,justify="left",anchor="w")

extra.start <- tkentry(extra.frm,textvariable=extra.start.val,background="#FFFFFF")
extra.const <- tkentry(extra.frm,textvariable=extra.const.val,background="#FFFFFF")
extra.label <- tkentry(extra.frm,textvariable=extra.label.val,background="#FFFFFF")

tkpack( extra.covlist, side="left",fill="both",expand="1")
tkpack( extra.scroll, side="right", fill="y")

tkpack( tklabel( extra.r, text="Link:"), side="left")
tkpack(extra.group,side="left")
tkpack(extra.link.menu,anchor="nw",padx=30)
tkpack(extra.fun.menu, anchor="s",side="left",fill="none")
tkpack(extra.message,anchor="w",after=extra.fun.menu,fill="none")

tkpack( extra.r, side="top",fill="x",expand="1")
tkpack( extra.l, side="top", after=extra.r,fill="x",expand="1")
tkpack(extra.label,tklabel(extra.frm,text="Labels:"),fill="x",side="bottom")
tkpack(extra.const,tklabel(extra.frm,text="Constraints:"),fill="x",side="bottom")
tkpack(extra.start,tklabel(extra.frm,text="Starting values:"),fill="x",side="bottom")



# expected frame


expected.frm <- tkwidget(param.frm,"labelframe",text="Expected function:",borderwidth=2,relief="groove")
expected.menu <- tkwidget( expected.frm, "tk_optionMenu",expectedfun,"Mu","User")
tkconfigure(expected.menu,width=4)
expected.submenu <- .Tk.newwin(tkcget(expected.menu,"-menu"))
tkentryconfigure(expected.submenu,"User",command=function() tclObj(expectedfunlst) <- makethelist()$funlst)
tkentryconfigure(expected.submenu,"Mu",command=function() tclObj(expectedfunlst) <- "" )    
expected.funscr <- tkframe(expected.frm)
expected.fun <- tklistbox( expected.funscr,
                      listvariable=expectedfunlst,
		yscrollcommand=function(...) tkset( expected.scroll, ...),
		selectmode="single",
#		width=15,
		height=4,
		exportselection=0,background="#FFFFFF")

expected.scroll <- tkscrollbar( expected.funscr, orient="vert",
		command=function(...) tkyview( expected.fun, ...))
tkpack(expected.menu,side="top")
tkpack(expected.fun,side="left",fill="both",expand="1")
tkpack(expected.scroll,fill="y",expand="1")
tkpack(expected.funscr,fill="both",expand="1")


chklen.but <- tkcheckbutton(expected.frm,text="check length",variable=chklen)
tkpack(chklen.but)

tkbind(expected.fun,"<ButtonRelease-1>",function() tclObj(expected.sel) <-
       as.character(tclObj(expectedfunlst))[as.numeric(tkcurselection(expected.fun))+1])

# maxit frame

maxit.entry <- tkentry( optim.frm, textvariable=maxit.value, width=8,background="#FFFFFF")
tkpack( tklabel( optim.frm, text="Max iter.", padx=4), maxit.entry,
		side="left")
gradtol.entry <- tkentry( optim.frm, textvariable=gradtol.value, width=8,background="#FFFFFF")
tkpack( tklabel( optim.frm, text="Grad tol.", padx=4), gradtol.entry,
		side="left")
steptol.entry <- tkentry( optim.frm, textvariable=steptol.value, width=8,background="#FFFFFF")
tkpack( tklabel( optim.frm, text="Step tol.", padx=4), steptol.entry,
		side="left")
printlevel.entry <- tkentry( optim.frm, textvariable=printlevel.value, width=2,background="#FFFFFF")
tkpack( tklabel( optim.frm, text="Print level", padx=4), printlevel.entry,
		side="left")
tkpack(run.message,side="left",fill="x",expand=TRUE)

tkpack( top.frm, side="top", fill="x")
tkpack( optim.frm, side="top", fill="both")
tkpack( mu.frm, shape.frm, extra.frm,side="left",fill="y")
tkpack( expected.frm,fill="both",expand="1")
tkpack( bot.frm, side="top", fill="x")
refresh()
} # end of moc.gui func


import.data <- function(file)
  {
    # library foreign should be loaded
    require(foreign) || stop("Install the foreign package to have access to this functionality")
    owarn<-getOption("warn")
    on.exit(options("warn"=owarn),add=TRUE)
    options("warn"=2)
    if (!file.exists(file)) stop(paste("The file",file,"is unreachable"))
    essai<-try(read.spss(file,to.data.frame=TRUE),silent=TRUE)
    if (!inherits(essai,"try-error")) return(essai)
    essai<-try(read.xport(file),silent=TRUE)
    if (!inherits(essai,"try-error")) return(essai)
    header.csv <- (unlist(read.csv(file,nrow=1,header=FALSE,fill=FALSE,
                      as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE)))
    header.csv2 <- (unlist(read.csv2(file,nrow=1,header=FALSE,fill=FALSE,
                      as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE)))
    header.delim <- (unlist(read.delim(file,nrow=1,header=FALSE,fill=FALSE,
                      as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE)))
    header.delim2 <- (unlist(read.delim2(file,nrow=1,header=FALSE,fill=FALSE,
                      as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE)))
    nfields.tab <- count.fields(file,sep="\t",blank.lines.skip = FALSE,skip=1)
    nfields.tab2 <- (is.character(header.delim2) & all(nfields.tab==length(header.delim2)))
    nfields.tab <- (is.character(header.delim) & all(nfields.tab==length(header.delim)))
    nfields.comma <- count.fields(file,sep=",",blank.lines.skip = FALSE,skip=1)
    nfields.comma <- (is.character(header.csv) & all(nfields.comma==length(header.csv)))
    nfields.column <- count.fields(file,sep=";",blank.lines.skip = FALSE,skip=1)
    nfields.column <- (is.character(header.csv2) & all(nfields.column==length(header.csv2)))
    essai<-try(read.csv(file,fill=FALSE,as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE),silent=TRUE)
    if (!inherits(essai,"try-error"))
      if( !(dim(essai)[1]==0) &
         (length(grep("[\t;\\\ [:cntrl:]]",paste(unlist(essai[,]))))==0) &
         nfields.comma) return(essai)
    essai<-try(read.delim(file,fill=FALSE,as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE),silent=TRUE)
    if (!inherits(essai,"try-error"))
         if (!(dim(essai)[1]==0) &
         length(grep("[,;\ \\[:cntrl:]]",paste(unlist(essai[,]))))==0 &
         nfields.tab) return(essai)
    essai<-try(read.delim2(file,fill=FALSE,as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE),silent=TRUE)
    if (!inherits(essai,"try-error"))
        if (!(dim(essai)[1]==0) &
        (length(grep("[;,\\\ [:cntrl:]]",paste(unlist(essai[,]))))==0) &
            nfields.tab) return(essai)
    essai<-try(read.csv2(file,fill=FALSE,as.is=TRUE,strip.white=FALSE,blank.lines.skip = FALSE),silent=TRUE)
    if (!inherits(essai,"try-error"))
        if (!(dim(essai)[1]==0) &
        (length(grep("[\t,\\\ [:cntrl:]]",paste(unlist(essai[,]))))==0) &
        nfields.tab) return(essai)
    return(-1)
  }

##open a Tk window describing the GUI
if(require(tcltk)) tkpager(system.file("GUI","moc-gui.Readme",package="moc"),
        title="MOC-GUI",header="Readme",delete.file=FALSE)
