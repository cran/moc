#  Copyright (C) 2002-2003 Bernard Boulerice
#  Sweave utilities for Mixture of Curves
#  This file is not loaded with the library moc, it has to be sourced separately.
#  For examples of use of these function see the file moc.Rnw


#  This function produces the entries of the maximum likelihood estimates
#  to be included in a LaTeX table.
#  digits: number of significant digits to be printed
#  headers: logical value specifying the use of the default column labels (Parameters & Estimates ...)
#  spacing: LaTeX command to include extra spacing before the coefficient labels
#           extracted from the functions parameters attributes


TeXCoefTable.moc<-function(object,digits=NULL,headers=TRUE,spacing="\\hspace{24pt}")
  {
    mle <- cbind(object$coef, se<-sqrt(diag(object$cov)) , w <- (object$coef/se)^2, (1 - pchisq(w, 1)))
    nl<-object$npar[1];ns<-object$npar[2];nx<-object$npar[3];nm<-object$npar[4]

    if(headers) cat("Parameters & Estimates & Standard Error & Wald & $P[\\chi_1^2 > \\textrm{Wald}]$ \\\\ \n \\hline \n")

   
    if(nl >0){
      cat("Location  \\\\ \n ")
    lblname<-paste(spacing,attr(object$gmu,"parameters"))
    coeftable<-apply(matrix(mle[1:nl,],nl,4),2,format,digits=digits)
    cat(paste(apply(cbind(lblname,matrix(coeftable,nl,4)),
                    1,paste,collapse=" & "),collape=" \\\\ \n "),"\n")
    cat("\\hline\n")}

    if(ns >0){
      cat("Scale \\\\ \n ")
    lblname<-paste(spacing,attr(object$gshape,"parameters"))
    coeftable<-apply(matrix(mle[(1+nl):(nl+ns),],ns,4),2,format,digits=digits)
    cat(paste(apply(cbind(lblname,matrix(coeftable,ns,4)), 1,paste,collapse=" & "),
              collape=" \\\\ \n "),"\n")
    cat("\\hline\n")}
    if(nx >0){
      cat("Extra \\\\ \n ")
    lblname<-paste(spacing,attr(object$gextra,"parameters"))
    coeftable<-apply(matrix(mle[(1+nl+ns):(nl+ns+nx),],nx,4),2,format,digits=digits)
   cat(paste(apply(cbind(lblname,matrix(coeftable,nx,4)),1,paste,collapse=" & "),
             collape=" \\\\ \n"),"\n")
    cat("\\hline\n")}

    if(nm >0){
      cat("Mixture  \\\\ \n ")
    lblname<-paste(spacing,attr(object$gmixture,"parameters"))
    coeftable<-apply(matrix(mle[(1+nl+ns+nx):(nl+ns+nx+nm),],nm,4),2,format,digits=digits)
    cat(paste(apply(cbind(lblname,matrix(coeftable,nm,4)),1,paste,collapse=" & "),
              collape=" \\\\ \n"),"\n")
    cat("\\hline\n")}

    invisible(mle)

  }

#  The following function produces the entries of the fit (BIC) of one or more moc objects
#  to be included in a LaTeX table.
#  ...: moc objects.
#  lbl: character vector of length equals to the  number of objects used as labels
#       defaults to the sequence number
#  digits: number of significant digits to be printed
#  headers: logical value specifying the use of the default headers

TeXFitTable.moc<-function(...,lbl=NULL,digits=NULL,headers=TRUE)
  {
    fit<-as.matrix(cbind(AIC(...,k=0)[,1],AIC(...,k="BIC")))
    nmoc<-dim(fit)[1]
    if(is.null(lbl)) lbl<-paste(1:nmoc)
    if(headers) cat("Model & $-2 \\log(\\textrm{Likelihood})$ & BIC & Entropy & ICL-BIC & Df  \\\\ \n ")
    cat("\\hline\n")
    cat(paste(apply(cbind(lbl,formatC(fit,digits=digits)),1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n")
    invisible(fit)
  }

#  TeXMixPTable.moc generates the entries for the mean mixture probabilities of a
#  fitted moc object
#  object: a fitted moc object
#  lbl: character vector used as label
#  digits: number of digits to be printed
#  headers: logical value, when TRUE use the default headers

TeXMixPTable.moc<-function(object,lbl=NULL,digits=NULL,headers=TRUE)
  {
    np<-cumsum(object$npar)
    if(!is.null(lbl) & headers) cat(" & ")
    if(headers) cat(paste("Group",1:object$groups,collapse=" & "),"\\\\ \n ")
    prob<-apply(object$gmixture(object$coef[(np[3]+1):np[4]]), 2, mean)
    if(!is.null(lbl)) cat(lbl," & ")
    cat(paste(format(prob,digits=digits),collapse=" & "),"\\\\ \n")
    invisible(prob)
  }

#  The function TeXMeanTable.moc generates the mean fitted and observed values
#  to be used within a LaTeX table
#  object: a fitted moc object
#  digits: number of digits to be printed
#  headers: logical value specifying the use of the default headers
#  spacing: LaTeX command to include extra spacings before the labels

TeXMeanTable.moc<-function(object,digits=NULL,headers=TRUE,spacing="\\hspace{24pt}")
  {
    if(headers) cat("Fitted \\\\  \n")
    cat(" & ",paste(dimnames(object$fitted.mean)[[2]],collapse=" & "),"\\\\ \n")
    cat(paste(spacing,apply(cbind(dimnames(object$fitted.mean)[[1]],formatC(object$fitted.mean,digits=digits))
                    ,1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n")
    cat("\\hline \n")
    if(headers) cat("Observed \\\\  \n")
    cat(" & ",paste(dimnames(object$observed.mean)[[2]],collapse=" & "),"\\\\ \n")
    cat(paste(spacing,apply(cbind(dimnames(object$observed.mean)[[1]]
                    ,formatC(object$observed.mean,digits=digits)),1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n")
    cat("\\hline \n")
    invisible(rbind(object$fitted.mean,object$observed.mean))
  }

# The following function creates  LaTeX table entries from a general list or
# any object that can be coerced to a matrix. The function is called recursively
# to process all elements.
# x: a list or an object that can be coerced to a matrix
# digits: number of digits to be printed
# headers: if TRUE use the list component names as labels
# spacing:  LaTeX command to include extra spacing before the coefficient labels
            
TeXListTable.moc<-function(x,digits=NULL,headers=TRUE,spacing="\\hspace{24pt}")
  {
    if(is.list(x)) {
      lbl<-names(x)
      sapply(1:length(x),function(ind)
             { if(headers) cat(lbl[ind]," \\\\    \n")
               TeXListTable.moc(x[[ind]],digits=digits,headers=headers)
             })} else
    {
      matx<-as.matrix(format(x,digits=digits))
      lbl<-dimnames(matx)
      if(!is.null(lbl[[2]])) cat(" & ",paste(lbl[[2]],collapse=" & ")," \\\\ \n")
      cat(paste(spacing,apply(cbind(lbl[[1]],matx),1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n")
      cat("\\hline \n")
    }
    invisible()
  }
