#  Copyright (C) 2002-2005 Bernard Boulerice
#  Sweave utilities for Mixture of Curves

#  This file is not loaded with the library moc, it has to be sourced separately.
#  For examples of use of these function see the file moc.Rnw

#  The following function produces the entries of the maximum likelihood estimates
#  to be included in a LaTeX table.
#  digits: number of significant digits to be printed
#  headers: character vector specifying the column labels (Parameters & Estimates ...)
#  spacing: LaTeX command to include extra spacing before the coefficient labels
#           extracted from the functions parameters attributes
#  align: alignment string fot tabular environment
#  ... : arguments to be passed to cat like file and append


TeXCoefTable.moc<-function(object,digits=NULL,headers=
                  c("Parameters","Estimates","Standard Error","Wald","$P[\\chi_1^2 > \\textrm{Wald}]$"),
                           spacing="\\hspace{6pt}",align="*{5}{l}",...)
  {
    mle <- cbind(object$coef, se<-sqrt(diag(object$cov)) , w <- (object$coef/se)^2, (1 - pchisq(w, 1)))
    nl<-object$npar[1];ns<-object$npar[2];nx<-object$npar[3];nm<-object$npar[4]
    cat("\n \\begin{tabular}{",align, "}\n\\hline\\hline \n ",...)
    if(length(headers)>0) cat(paste(headers,collapse=" & ")," \\\\ \n \\hline \n",...)
    if(nl >0){
      cat("Location  \\\\ \n ",...)
    lblname<-paste(spacing,attr(object$gmu,"parameters"))
    coeftable<-rbind(apply(matrix(mle[1:nl,],nl,4),2,format,digits=digits))
    cat(paste(apply(cbind(lblname,matrix(coeftable,nl,4)),
                    1,paste,collapse=" & "),collape=" \\\\ \n "),"\n",...)
    cat("\\hline\n",...)}

    if(ns >0){
      cat("Scale \\\\ \n ",...)
    lblname<-paste(spacing,attr(object$gshape,"parameters"))
    coeftable<-rbind(apply(matrix(mle[(1+nl):(nl+ns),],ns,4),2,format,digits=digits))
    cat(paste(apply(cbind(lblname,matrix(coeftable,ns,4)), 1,paste,collapse=" & "),
              collape=" \\\\ \n "),"\n",...)
    cat("\\hline\n",...)}
    if(nx >0){
      cat("Extra \\\\ \n ",...)
    lblname<-paste(spacing,attr(object$gextra,"parameters"))
    coeftable<-rbind(apply(matrix(mle[(1+nl+ns):(nl+ns+nx),],nx,4),2,format,digits=digits))
   cat(paste(apply(cbind(lblname,matrix(coeftable,nx,4)),1,paste,collapse=" & "),
             collape=" \\\\ \n"),"\n",...)
    cat("\\hline\n",...)}

    if(nm >0){
      cat("Mixture  \\\\ \n ",...)
    lblname<-paste(spacing,attr(object$gmixture,"parameters"))
    coeftable<-rbind(apply(matrix(mle[(1+nl+ns+nx):(nl+ns+nx+nm),],nm,4),2,format,digits=digits))
    cat(paste(apply(cbind(lblname,matrix(coeftable,nm,4)),1,paste,collapse=" & "),
              collape=" \\\\ \n"),"\n",...)
    cat("\\hline\n",...)}
    cat("\\hline \n \\end{tabular} \n \n ",...)
    invisible(mle)
  }

#  The following function produces the entries of the fit (BIC) of one or more moc objects
#  to be included in a LaTeX table.
#  moc.list: character vector of moc objects names.
#  lbl: character vector of length equals to the  number of objects used as labels
#       defaults to the moc.list names
#  digits: number of significant digits to be printed
#  headers: character vector specifying the column labels
#  align: alignment string fot tabular environment
#  ... : arguments to be passed to cat like file, append

TeXFitTable.moc<-function(moc.list,lbl=moc.list,digits=NULL,
                          headers=c("Model", "$-2 \\log(\\textrm{Likelihood})$","BIC","Entropy","ICL-BIC","Df","Groups"),
                          align="*{7}l", ...)
  {
    if(!is.character(moc.list)) stop("\n\tmoc.list should be a character vector of moc object names.\n")
    ng<-sapply(eval(parse(text=paste("list(",paste(moc.list,collapse=","),")"))),function(x) x$groups)
    fit<-as.matrix(cbind(eval(parse(text=paste("AIC(",paste(moc.list,collapse=","),",k=\"BIC\")"))),Groups=ng))
    nmoc<-dim(fit)[1]
    if(is.null(lbl)) lbl<-paste(1:nmoc)
    cat("\n\\begin{tabular}{",align," } \n\\hline\\hline \n",...)
    if(length(headers)>0) cat(paste(headers,collapse=" & ")," \\\\ \n ",...)
    cat("\\hline\n",...)
    cat(paste(apply(cbind(lbl,rbind(apply(fit,2,format,digits=digits))),1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n",...)
    cat("\\hline\\hline \n \\end{tabular} \n  ",...)
    invisible(fit)
  }

#  TeXMixPTable.moc generates the entries for the mean mixture probabilities of a
#  fitted moc object
#  object: a fitted moc object
#  lbl: character vector used as label
#  digits: number of digits to be printed
#  headers: character vector specifying the column labels
#  align: alignment string fot tabular environment
#  ... : arguments to be passed to cat like file and append

TeXMixPTable.moc<-function(object,lbl=NULL,digits=NULL,headers=paste("Group",1:object$groups),
                           align=paste("l*{",object$groups,"}{l}"),...)
  {
    np<-cumsum(object$npar)
    cat(paste("\n\\begin{tabular}{",align,"}\n\\hline\\hline \n "),...)
    if(!is.null(lbl) & length(headers)> 0) cat(" & ",...)
    if(length(headers)>0) cat(paste(headers,collapse=" & "),"\\\\ \n ",...)
    prob<-apply(object$gmixture(object$coef[(np[3]+1):np[4]]), 2, weighted.mean,eval(object$prior.weights))
    if(!is.null(lbl)) cat(lbl," & ",...)
    cat(paste(rbind(format(prob,digits=digits)),collapse=" & "),"\\\\ \n",...)
    cat("\n \\hline\\hline \n\\end{tabular} \n \n ",...)
    invisible(prob)
  }

#  The function TeXMeanTable.moc generates the mean fitted and observed values
#  to be used within a LaTeX table
#  object: a fitted moc object
#  digits: number of digits to be printed
#  headers: character vector specifying sub-table labels
#  spacing: LaTeX command to include extra spacings before the labels
#  transpose: logical value indicating whether to transpose or not the matrices
#  align: alignment string fot tabular environment
#  ... : arguments to be passed to cat like file and append

TeXMeanTable.moc<-function(object,digits=NULL,headers=c("Fitted","Observed"),spacing="\\hspace{6pt}",transpose=FALSE,align=NULL,...)
  {
    if(transpose){
      tran <- function(x) t(x)
      if(is.null(align)) align <- paste("l*{",dim(object$fitted.mean)[1],"}{l}")
    } else {
      tran <- function(x) x
      if(is.null(align)) align <- paste("l*{",dim(object$fitted.mean)[2],"}{l}")
    }      
    cat(paste("\n\\begin{tabular}{ ",align,"}\n\\hline\\hline \n "),...)
    if(length(headers)>0) cat(headers[1]," \\\\  \n",...)
    cat(" & ",paste(dimnames(tran(object$fitted.mean))[[2]],collapse=" & "),"\\\\ \n",...)
    cat(paste(spacing,apply(cbind(dimnames(tran(object$fitted.mean))[[1]],rbind(apply(tran(object$fitted.mean),2,format,digits=digits)))
                    ,1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n",...)
    cat("\\hline \n",...)
    if(length(headers)>0) cat(headers[2]," \\\\  \n",...)
    cat(" & ",paste(dimnames(tran(object$observed.mean))[[2]],collapse=" & "),"\\\\ \n",...)
    cat(paste(spacing,apply(cbind(dimnames(tran(object$observed.mean))[[1]]
                    ,rbind(apply(tran(object$observed.mean),2,format,digits=digits))),1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n",...)
    cat("\n \\hline\\hline \n \\end{tabular} \n \n ",...)
    invisible(rbind(object$fitted.mean,object$observed.mean))
  }


#  TeXObsfitMixPTable.moc generates the entries for the mean mixture probabilities of a
#  fitted moc object
#  object: a fitted moc object
#  lbl: character vector used as label
#  digits: number of digits to be printed
#  headers: character vector specifying the prior and posterior probability labels
#  align: alignment string fot tabular environment
#  ... : arguments to be passed to cat like file and append

TeXObsfitMixPTable.moc<-function(object,lbl="Probabilities",digits=NULL,headers=c("Prior","Posterior"),
                           fbox="\\fbox",align=NULL,...)
  {
    if(!all(sapply(object,inherits,"by"))) stop("\n\tProbably not an obsfit.moc returned values !\n")
    if(is.null(align)) align <- paste("l*{", max(sapply(object$"Mean Prior Probabilities",length)),"}{l}")
    tex.by <- function (x, dig=digits) 
      {
        d <- dim(x)
        dn <- dimnames(x)
        dnn <- names(dn)
        lapply(seq(along = x), function(i, x, dig) {
          ii <- i - 1
          if(!is.null(x[[i]]) && length(x[[i]])!=0 )  {
            if(length(x)!=1){
              cat("\\\\",fbox,"{\\vbox{",...)
            for (j in seq(along = dn)) {
            iii <- ii%%d[j] + 1
            ii <- ii%/%d[j]
           cat("\\hbox{",dnn[j], ": ", dn[[j]][iii], "} ",sep = "",...)
          }
            cat("}}",...)}
            cat(" & ",paste(names(x[[i]]),collapse=" & "),"\\\\ \n",...)
            cat(paste(paste(lbl," & "),paste(rbind(format(x[[i]],digits=dig)),collapse=" & "))," \\\\ \n",...)}
        }, x, dig=digits)
        invisible()
      }
    cat(paste("\n\\begin{tabular}{",align,"}\n\\hline\\hline \n "),...)
    if(length(headers)>0) cat(headers[1]," \\\\ \n",...)
    tmp <- object$"Mean Prior Probabilities"
    tex.by(tmp,spacing)
    cat("\\hline \n",...)
    tmp <- object$"Mean Posterior Probabilities"
    if(length(headers)>0) cat(headers[2]," \\\\ \n",...)
    tex.by(tmp,spacing)
    cat("\n \\hline\\hline \n \\end{tabular} \n \n ",...)
    invisible(list("Mean Prior Probabilities"=object$"Mean Prior Probabilities",
                  "Mean Posterior Probabilities"= object$"Mean Posterior Probabilities"))
 }

#  The function TeXObsfitMeanTable.moc generates the table of mean observed and fitted
#  values from an obsfit.moc output
#  object: the output of obsfit.moc of a moc object
#  digits: number of digits to be printed
#  headers: character vector specifying sub-table labels first for fitted second for observed
#  spacing: LaTeX command to include extra spacings before the labels
#  transpose: logical value indicating whether to transpose or not the matrices
#  align: alignment string fot tabular environment
#  ... : arguments to be passed to cat like file and append

TeXObsfitMeanTable.moc<-function(object,digits=NULL,headers=c("Fitted","Observed"),spacing="\\hspace{6pt}",
                                 fbox="\\fbox",transpose=FALSE,align=NULL,...)
  {
    if(!all(sapply(object,inherits,"by"))) stop("\n\tProbably not an obsfit.moc returned values !\n")
    tex.by <- function (x, dig=digits) 
      {
        d <- dim(x)
        dn <- dimnames(x)
        dnn <- names(dn)
        lapply(seq(along = x), function(i, x, dig) {
          ii <- i - 1
          if(!is.null(x[[i]]) )  {
            if(length(x)!=1){
              cat("\\\\",fbox,"{\\vbox{",...)
            for (j in seq(along = dn)) {
            iii <- ii%%d[j] + 1
            ii <- ii%/%d[j]
           cat("\\hbox{",dnn[j], ": ", dn[[j]][iii], "} ",sep = "",...)
          }
            cat("}}",...)}
            cat(" & ",paste(dimnames(tran(x[[i]]))[[2]],collapse=" & "),"\\\\ \n",...)
            cat(paste(spacing,apply(cbind(dimnames(tran(x[[i]]))[[1]],
                                    rbind(apply(tran(x[[i]]),2,format,digits=dig))),
                            1,paste,collapse=" & "),collapse=" \\\\ \n")," \\\\ \n",...)}
        }, x, dig=digits)
        invisible()
      }
    rowcol <- apply(sapply(object$"Mean function Observed Values",
                           function(xl) if(is.null(xl)) c(0,0) else dim(xl)),1,max)
    if(transpose){
      tran <- function(x) t(x)
      if(is.null(align)) align <- paste("l*{",rowcol[1],"}{l}")
    } else {
      tran <- function(x) x
      if(is.null(align)) align <- paste("l*{",rowcol[2],"}{l}")
    }      
    cat(paste("\n\\begin{tabular}{",align,"}\n\\hline\\hline \n "),...)
    if(length(headers)>0) cat(headers[1]," \\\\ \n",...)
    tmp <- object$"Mean function Expected Values"
    tex.by(tmp,spacing)
    cat("\\hline \n",...)
    tmp <- object$"Mean function Observed Values"
    if(length(headers)>0) cat(headers[2]," \\\\ \n",...)
    tex.by(tmp,spacing)
    cat("\n \\hline\\hline \n \\end{tabular} \n \n ",...)
    invisible(list("Mean function Expected Values"=object$"Mean function Expected Values",
                  "Mean function Observed Values"= object$"Mean function Observed Values"))
  }



# The following function creates  LaTeX table entries from a general list or
# any object that can be coerced to a matrix. The function is called recursively
# to process all elements.
# x: a list or an object that can be coerced to a matrix
# digits: number of digits to be printed
# headers: if TRUE use the list component names as labels
# spacing:  LaTeX command to include extra spacing before the coefficient labels
# align1: single character for alignment of headers
# align2: single character for alignment of matrix elements
#  ... : arguments to be passed to cat like file and append

TeXListTable.moc<-function(x,digits=NULL,headers=TRUE,spacing="\\hspace{6pt}",align1="l",align2="l",...)
  {
    if(is.list(x)) {
      lbl<-names(x)
     cat("\n  \\begin{tabular}{*{",1,"}{",align1,"}}   \n",...)
       sapply(1:length(x),function(ind)
             { if(headers) cat("\\\\ \\multicolumn{",1,"}{",align1,"}{",lbl[ind],"} \\\\ \n",...)
               TeXListTable.moc(x[[ind]],digits=digits,headers=headers,spacing=spacing)
               if(ind==length(x)) cat("\n \\end{tabular} \\\\ \n ",...)
             })} else
    {
      lbl<-dimnames(x)
      matx<-rbind(apply(as.matrix(x),2,format,digits=digits))
      cat("\\begin{tabular}{*{",dim(matx)[2]+1,"}{",align2,"}} \\hline\\hline \n",...)
      if(!is.null(lbl[[2]])) cat(" & ",paste(lbl[[2]],collapse=" & ")," \\\\ \n",...)
      cat(paste(spacing,apply(cbind(lbl[[1]],matx),1,paste,collapse="& "),collapse=" \\\\ \n"),
" \\\\ \n",...)
      cat("\\hline\\hline \n \\end{tabular} \\\\   \n",...)
    }
    invisible()
  }


