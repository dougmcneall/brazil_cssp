#### make canconical correlation emulator based on cancor

## have to be careful to allow for missing x components where ncol(X)
## < ncol(Y)

shapiro.pvalue <- function(X) {
  if (sd(X) == 0.0) {
    return(0.0)
  } else {  
    return(shapiro.test(X)$p.value)
  }  
}


ksnorm <- function(x) {
  sdx <- sd(x)
  if (sdx == 0.0) {
    return(0.0)
  } else {  
    xx <- (x - mean(x)) / sdx
    ks <- ks.test(xx, 'pnorm', 0, 1)
    return(ks$p.value)
  }
}

extractAICc <- function(fit, scale, k=2, ...) {
    bAIC <- extractAIC(fit, scale, k=k, ...) 
    edf <- bAIC[1]
    n <- length(residuals(fit))
    bAIC <- bAIC[2] + (2 * edf + (edf * 1)) / (n - edf - 1) 
    return(bAIC)
}    


library(DiceKriging)


makeGPEmStep <- function(x, xNames, y, e, nugget=NULL, nuggetEstim=FALSE, meanTerms=NULL,
                     noiseVar=NULL, seed=NULL,
                     trace=FALSE, maxit=100, REPORT=10, factr=1e7, pgtol=0.0, 
                     parinit=NULL, popsize=100) {

    response <- as.vector(y)
    design <- as.matrix(x)
    if (is.null(nugget)) {nugget = sd(e)}
    
    if (ncol(design) > length(response)) stop('Too many terms in design matrix. Must be less than number of experiments')

    ## We need to standardise the data to the unit hypercube to build the emulator - otherwise the parameter units dominate. There might be a more efficient way of doing this but this is quick and easy to interpret so I have stuck with it. 
    nexp = nrow(design)
    npar = ncol(design)
    nval = length(response)
    stopifnot(nval == nexp)

    cmax = rep(0.5, npar)
    cmin = rep(-0.5, npar)
#    stopifnot(cmax <=-0.5 & cmin >= 0.5)
    
    subinputs = design #- 0.5
    xvars = labels(terms(xNames[[2]]))      # xNames is a formula
    control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, factr=factr, pgtol=pgtol, pop.size=popsize)

    # Build the emulator. 
      data = data.frame(response=response, x=subinputs)
      colnames(data) <- c("response", xvars)
      
      m0 = km(y ~ 1, design=data.frame(subinputs), response=response, nugget=nugget,
             nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
             
      coefs0 = m0@covariance@range.val 
      print('coefs0')
      print(coefs0)
      
      start.form = as.formula(paste("y ~ ", paste(xvars, collapse= "+")))
      
# use BIC so that model is parsimonious and allow GP to pick up any other behaviour not
# explained by the key linear terms      
      startlm = lm(start.form, data=data)
#      print('before step')
#      print(startlm)
      steplm = step(startlm, direction="both", k=log(nval), trace=FALSE)
      print('after step')
      print(steplm)
      form = as.formula(steplm)
      print('Formula')
      print(form)
      data$response = NULL
      labels = labels(terms(steplm))
      labels = labels[!(labels %in% c('response'))]
      if (length(labels) > 0) {
        start.form = as.formula(paste("~ ", paste(labels, collapse= "+")))
      } else {
        start.form = as.formula("~ 1")
      }    
      print("Step has found formula:")
      print(start.form)
      if (!is.null(seed)) {set.seed(seed)}
      m = km(start.form, design=data, response=response, nugget=nugget, parinit=coefs0,
             nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
    
    return(list(x=design, y=response, e=e, nugget=nugget, nugget.estim=nuggetEstim, noise.var=noiseVar, emulator=m, cmax=cmax, cmin=cmin, xNames=xNames, seed=seed, coefs=m@covariance@range.val, trends=m@trend.coef, meanTerms=all.vars(start.form), 
coefs0=coefs0))
}





makeGPEm <- function(x, xNames, y, e, nugget=NULL, nuggetEstim=FALSE, meanTerms=NULL,
                     noiseVar=NULL, fitMeanResponse=TRUE, useStep=FALSE, seed=NULL,
                     trace=FALSE, maxit=100, REPORT=10, factr=1e7, pgtol=0.0, 
                     parinit=NULL, popsize=100, reltol=sqrt(.Machine$double.eps)) {

    response <- as.vector(y)
    design <- as.matrix(x)
    if (is.null(nugget)) {nugget = sd(e)}
    
    if (ncol(design) > length(response)) stop('Too many terms in design matrix. Must be less than number of experiments')

    ## We need to standardise the data to the unit hypercube to build the emulator - otherwise the parameter units dominate. There might be a more efficient way of doing this but this is quick and easy to interpret so I have stuck with it. 
    nexp = nrow(design)
    npar = ncol(design)
    nval = length(response)
    stopifnot(nval == nexp)

    cmax = rep(0.5, npar)
    cmin = rep(-0.5, npar)
#    stopifnot(cmax <=-0.5 & cmin >= 0.5)
    
    subinputs = design #- 0.5
    xvars = labels(terms(xNames[[2]]))      # xNames is a formula
    control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, reltol=reltol, factr=factr, pgtol=pgtol, pop.size=popsize)

    # Build the emulator. 
#    print(paste("sgadhagdhsdh", nugget, sep=" "))
#    print(nugget.estim)
#    print(noise.var)
    if (useStep) {
#      print(xNames)
      data = data.frame(response=response, x=subinputs)
      colnames(data) <- c("response", xvars)
      
      if (!is.null(meanTerms)) {
        start.form = as.formula(paste("y ~ ", paste(meanTerms, collapse= "+")))
      } else {
        start.form = as.formula("y ~ 1")
      }
      print("start formula")
      print(start.form)  
# use BIC so that model is parsimonious and allow GP to pick up any other behaviour not
# explained by the key linear terms      
      startlm = lm(start.form, data=data)
      steplm = step(startlm, direction="both", k=log(nval), trace=FALSE)
      form = as.formula(steplm)
      data$response = NULL
      labels = labels(terms(steplm))
      labels = labels[!(labels %in% c('response'))]
      if (length(labels) > 0) {
        start.form = as.formula(paste("~ ", paste(labels, collapse= "+")))
      } else {
        start.form = as.formula("~ 1")
      }    
      print("Step has found formula:")
      print(start.form)
      if (!is.null(seed)) {set.seed(seed)}
      m = km(start.form, design=data, response=response, nugget=nugget, parinit=parinit,
             nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
    } else {
      if (!is.null(meanTerms)) {
        start.form = as.formula(paste("y ~ ", paste(meanTerms, collapse= "+")))
      } else if (fitMeanResponse) {
        start.form = as.formula("~ .")
      } else {
        start.form = as.formula("~ 1")
      }
      print("start formula")
      print(start.form)  
      colnames(subinputs) <- xvars

#      if (fitMeanResponse) {
       if (!is.null(seed)) {set.seed(seed)}
        m = km(start.form, design=data.frame(subinputs), response=response, nugget=nugget, 
               nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list,             
               parinit=parinit)
#      } else {
#        m = km(~1, design=data.frame(subinputs), response=response, nugget=nugget, 
#               nugget.estim=nuggetEstim, noise.var=noiseVar, control=list(trace=FALSE))
#      }                  
    }
    
#    print(m@covariance)
#    print(m@trend.coef)
#    print(m@covariance@range.val)

    
    return(list(x=design, y=response, e=e, nugget=nugget, nugget.estim=nuggetEstim, noise.var=noiseVar, emulator=m, cmax=cmax, cmin=cmin, xNames=xNames, fitMeanResponse=fitMeanResponse, seed=seed, coefs=m@covariance@range.val, trends=m@trend.coef, meanTerms=all.vars(start.form)))
}


predGPEm <- function(gp, Xnew, hIndex=1:cc$ny, quiet=FALSE, sample=FALSE, leaveOutNoiseInSample=FALSE, return_var=NULL) {

km.type = 'UK'     # universal kriging

if (is.null(return_var)) {
   need.sd <- FALSE
} else {
  need.sd <- return_var
}


#print(min(Xnew))
#print(max(Xnew))
delta = 1e-3
stopifnot(min(Xnew) >= -0.5-delta & max(Xnew) <= 0.5+delta)

design = as.matrix(Xnew)
if (ncol(design) == 1 & nrow(design) == ncol(gp$x)) {design = matrix(Xnew, nrow=1)}

#print(ncol(gp$x))
#print(dim(design))
result = predict.km(gp$emulator, design, km.type, se.compute=need.sd, checkNames=FALSE)
#print(attributes(result))
samples = NA
if (need.sd) {
  rvar = result$sd**2
} else {
  rvar = NULL
}
result = list(mean=result$mean, var=rvar, samples=samples)
return(result)


}


jackKnifeGPEm <- function(GPobj) {
  km.type = 'UK'     # universal kriging
  result <- leaveOneOut.km(GPobj$emulator, km.type, trend.reestim=TRUE)
  errors <- as.matrix(GPobj$y - result$mean)
  var <- result$sd**2
#  print(result)
  return(list(errors=errors, var=var))
}



makeccEm <- function(x, xNames, y, e, xtol=0.99, ytol=0.99, etol=1.0, quiet=FALSE, 
                     use_svd_on_y=TRUE) {
#                     ,power=0.5alter.corr=FALSE, rescale.var=FALSE, ) {

#  names(x) <- xNames
  xsafe <- x
  ysafe <- y
  esafe <- e
  
  x <- as.matrix(x)
#  colnames(x) <- xNames
  
#  x <- x[, !(names(x) == "Intercept" | names(x) == "1")]# remove intercept and must be before as.matrix

  notIntercept <- (xNames == "Intercept" | xNames == "1")
  x <- x[, !notIntercept]            # remove intercept and must be before as.matrix
  
#      return((xNames == "Intercept" | xNames == "1"))

  
  y <- as.matrix(y)
  e <- as.matrix(e)
  
  if ((nr <- nrow(x)) != nrow(y)) 
      stop("unequal number of rows in 'cancor'")

  xcenter <- colMeans(x,)
  x <- x - rep(xcenter, rep.int(nrow(x), ncol(x)))

  ycenter <- colMeans(y,)
  y <- y - rep(ycenter, rep.int(nrow(y), ncol(y)))
    
  ysd <- apply(y, 2, sd)
  ysd[ysd == 0.0] <- 1.0
  y <-  sweep(y, 2, ysd, "/")
    
  ecenter <- colMeans(e,)
  e <- e - rep(ecenter, rep.int(nrow(e), ncol(e)))
  e <-  sweep(e, 2, ysd, "/")
  

    
# if TRUE rescale variances so that coherent variables do not dominate  
#  if (rescale.var) {
#    it <- seq(1,144)
#    ip <- seq(145,288)
#    cory <- cor(y)
#    diag(cory)<-0.0
#    mnCorrT <- mean(abs(cor(y[,it])))
#    mnCorrP <- mean(abs(cor(y[,ip])))
#    rescaleT <- (mnCorrT/((mnCorrT+mnCorrP)/2))	#^0.5 
#    rescaleP <- (mnCorrP/((mnCorrT+mnCorrP)/2))	#^0.5
#    rescale <- c(rep(rescaleT,144), rep(rescaleP,144))  
#    rescale <- colMeans(abs(cory))	# lower failed tests
#    rescale <- rescale / mean(rescale)
#    y <- sweep(y, 2, rescale^power, '/')
#    ysd <- ysd * rescale^power
#}
#  
  
#  cat("x",dim(x),"\n")
  stopifnot(is.matrix(x),
            is.matrix(y),
            nrow(x) == nrow(y))				#,            nrow(x) > ncol(x)) 

  xsvd <- svd(x)
  xpc <- xsvd$u %*% diag(xsvd$d) 
  xevals <- xsvd$d^2
  dx <- sum(cumsum(xevals) < xtol*sum(xevals))
#  xdash <- xsvd$u[,1:dx] %*% diag(xsvd$d[1:dx]) %*% t(xsvd$v[,1:dx])
#  xepl <- apply(xdash,2,var)/apply(x,2,var)
  
  ysvd <- svd(y)
  ypc <- ysvd$u %*% diag(ysvd$d)
  ypcvar <- apply(ypc, 2, var)
  yevals <- ysvd$d^2

# project the eigenvalues of Y onto the control run to see which ones are noise and store the result in YOK  
  epc <- e %*% ysvd$v
  epcvar <- apply(epc, 2, var)
  yok <- ypcvar*etol > epcvar
  temp <- yevals
  temp[!yok] <- 0.0
  ycumsum <- cumsum(temp)/sum(yevals)
#  print(ypcvar)
#  print(etol)
#  print(epcvar)
#  print(ycumsum)
#  print(ytol)
  yok <- ypcvar*etol >= epcvar & ycumsum < ytol
  dy <- sum(yok)
  if (!use_svd_on_y) dy = ncol(y)
#  ydash <- ysvd$u[,yok] %*% diag(ysvd$d[yok]) %*% t(ysvd$v[,yok])
  ydash <- ysvd$u[,yok] %*% tcrossprod(diag(ysvd$d[yok]), ysvd$v[,yok])
  yexpl <- apply(ydash,2,var)/apply(y,2,var)

  if (!quiet) cat("xtol=",xtol," ytol=",ytol,"dx=",dx," dy=",dy,"\n")

# do the canonical correlation. Centered once already  
  cc<-cancor(xpc[,1:dx], ypc[,yok], xcenter=FALSE, ycenter=FALSE)
  r <- length(cc$cor)
  iB <- solve(cc$ycoef)

#  if (alter.corr) {
#  ecoef <- epc[,yok] %*% cc$ycoef
#  evar <- apply(ecoef, 2, var)
#  yvar <- apply(ypc[,yok] %*% cc$ycoef, 2, var)	# seems to equal 1.0/(nrow(y)-1.0)
#  var.diff <- yvar-evar
#  var.diff[var.diff<0] <- yvar[var.diff<0]	# remove -ve diffs as probably due to sampling
#  rescale <- (yvar/var.diff)^0.5
#  cc$cor <- cc$cor * rescale[1:r]	#1.0 - (1.0-cc$cor)*1.0
#  cc$cor[cc$cor > 1] <- 1.0
#}
  
  if (r == nrow(cc$ycoef)) { # this one is easy!
 
    C <- cc$xcoef[,1:r,drop=FALSE] %*% (diag(cc$cor) %*% iB) # first r can cors
    Sig <- diag(1 - cc$cor^2, r)
	Sig2 <- (Sig / (nrow(x)-1.0))^0.5
    Sig <- crossprod(iB, Sig) %*% iB / (nrow(x)-1.0)

  } else if (r == nrow(cc$xcoef)) { # this one needs padding

    q <- nrow(cc$ycoef)
    C <- cc$xcoef %*% (diag(cc$cor) %*% iB[1:r,, drop=FALSE])
    Sig <- diag(1, q); diag(Sig)[1:r] <- 1 - cc$cor^2
	Sig2 <- (Sig / (nrow(x)-1.0))^0.5
    Sig <- crossprod(iB, Sig) %*% iB / (nrow(x)-1.0)	# is the denominator correct??
    
  } else stop('Never get here!')
  
#  chol <- Sig2 %*% iB %*% t(ysvd$v[,yok]) 
  chol <- Sig2 %*% tcrossprod(iB, ysvd$v[,yok]) 
#  var <- ysvd$v[,yok] %*% Sig %*% t(ysvd$v[,yok]) 
  var <- ysvd$v[,yok] %*% tcrossprod(Sig, ysvd$v[,yok]) 
#  print(paste("dims ", dim(chol), dim(var) ,"\n"))
#  print(paste("check ", range(crossprod(chol) - var) ,"\n"))
  
  var <- sweep(var, 1, ysd, '*')
  var <- sweep(var, 2, ysd, '*') 
  vardash <- var
  yresid <- y - ydash
  yresid <-  sweep(yresid, 2, ysd, "*")
  var <- var + cov(yresid)
#  print(dim(var))
#  print(diag(var)**0.5)
#  print(apply(yresid,2,sd))

  vc <- xsvd$v[,1:dx] %*% C 
#  vcvt <- vc %*% t(ysvd$v[,yok])
  vcvt <- tcrossprod(vc, ysvd$v[,yok])

  return(list(x=xsafe,y=ysafe,e=esafe,ny=ncol(y),xcenter=xcenter,ycenter=ycenter,ysd=ysd,xsvd=xsvd,ysvd=ysvd,C=C,Sig=Sig,rank=r,cc=cc,xtol=xtol,ytol=ytol,etol=etol,dx=dx,dy=dy,yexpl=yexpl,yok=yok,ydash=ydash,var=var,xNames=xNames,notIntercept=notIntercept,vcvt=vcvt,yresid=yresid,chol=chol))
}	
	
predccEm <- function(cc, Xnew, hIndex=1:cc$ny, quiet=FALSE, sample=FALSE, 
                     leaveOutNoiseInSample=FALSE, return_var=NULL) {

#	cat("dfg",hIndex ,"\n")

    if (is.null(return_var)) {
      need.sd <- FALSE
    } else {
      need.sd <- return_var
    }

	if (is.vector(Xnew)) {Xnew = matrix(Xnew, nrow=1) }
	
#	cat(dim(Xnew), "\n")
    stopifnot(is.matrix(Xnew),
              ncol(Xnew) == nrow(cc$xcoef))
			  
	if (all(hIndex == 1:cc$ny)) {		  
#	Xnew <- as.matrix(Xnew, ncol=length(cc$xNames))
#	cat(length(colnames(Xnew)), length(cc$xNames) ,"\n")
#    colnames(Xnew) <- cc$xNames  
#    notIntercept <- (cc$xNames == "Intercept" | cc$xNames == "1")
    Xnew <- Xnew[, !cc$notIntercept, drop=FALSE] # remove intercept but do not let it reduce to a vector with drop=FALSE 	
#	cat("dfs",dim(Xnew), length(cc$xcenter) ,"asf\n")
    Xcentered <- sweep(Xnew, 2, cc$xcenter, '-')
    Xpc <- Xcentered %*% cc$xsvd$v[,1:cc$dx]
    if (need.sd) {
      pred.scale <- 1.0 + apply(sweep(Xpc, 2, cc$xsvd$d[1:cc$dx], '/'), 1, sumsq)
    } else {
      pred.scale <- 1.0
    }
    Ypc <- Xpc %*% cc$C 
#    Yp <- Ypc %*% t(cc$ysvd$v[,cc$yok])
    Yp <- tcrossprod(Ypc, cc$ysvd$v[,cc$yok])
    EYp <- sweep(Yp, 2, cc$ysd, '*')
    EYp <- sweep(EYp, 2, cc$ycenter, '+')
	} else { 
#	cat(dim(Xnew), dim(cc$xsvd$v[,1:cc$dx]), dim(cc$C), dim(t(cc$ysvd$v[,cc$yok])), "\n")
    Xnew <- Xnew[, !cc$notIntercept, drop=FALSE] # remove intercept but do not let it reduce to a vector with drop=FALSE 	
#	cat("dfs",dim(Xnew), length(cc$xcenter) ,"asf\n")
    Xcentered <- sweep(Xnew, 2, cc$xcenter, '-')
##    Xpc <- Xcentered %*% cc$xsvd$v[,1:cc$dx]
##    Ypc <- Xpc %*% cc$C 
#	v <- cc$ysvd$v[,cc$yok]
#    Yp <- Xcentered %*% cc$vc %*% t(v[hIndex,])	
	Yp <- Xcentered %*% cc$vcvt[,hIndex]
    EYp <- sweep(Yp, 2, cc$ysd[hIndex], '*')
    EYp <- sweep(EYp, 2, cc$ycenter[hIndex], '+')
	}
	
	if (sample) {
	  samples = EYp + sampleEm(cc, nrow(Xnew), leaveOutNoise=leaveOutNoiseInSample)
	} else { samples=NA }  
	
	ans <- list(mean=EYp, samples=samples, scale=pred.scale)
	gc()
	
	return(ans)	#,covariance=cc$var))
}

sampleEm <- function(cc, n, leaveOutNoise=FALSE) {
  nr <- nrow(cc$chol)
  rand <- matrix(rnorm(nr*n), ncol=nr)
  samples <- rand %*% cc$chol
  samples <- sweep(samples, 2, cc$ysd, '*')
  if (!leaveOutNoise) {  
    library(mvtnorm)
    samples = samples + rmvnorm(n, sigma = cov(cc$yresid),
        method="svd")
  }
  return(samples)
}

jackKnifeEm <- function(cc) {
  errors <- matrix(0.0, nrow(cc$y), cc$ny)
  for (i in 1:nrow(cc$y)) {
    ccDrop1 <- makeccEm(cc$x[-i,], cc$xNames, cc$y[-i,], cc$e, xtol=cc$xtol, ytol=cc$ytol, etol=cc$etol, quiet=TRUE)      #, power=0.5, alter.corr=FALSE, rescale.var=FALSE)
	predDrop1 <- predccEm(ccDrop1, cc$x[i,], return_var=FALSE)
    errors[i,] = cc$y[i,] - predDrop1$mean
  }
  return(errors)
}

diagnoseEm <- function(cc, file1, names=colnames(cc$y)) {
  colnames(cc$y) <- names
  full.pred <- predccEm(cc,cc$x)
  yp <- full.pred$mean
  full.err.var1 <- var(cc$y-full.pred$mean)
  ydash <- sweep(cc$ydash, 2, cc$ysd, '*')
  ydash <- sweep(ydash, 2, cc$ycenter, '+')
  ySvdErr <- sweep(cc$y - cc$ydash, 2, cc$ysd, '*')


  ksnormality <- apply(ydash-yp, 2, ksnorm)
  normality <- apply(ydash-yp, 2, shapiro.pvalue)
  cat("%failed Shapiro tests=",100.0*sum(normality < 0.05)/ncol(cc$y) ,"\n")
  cat("%failed Kolsmir tests=",100.0*sum(ksnormality < 0.05)/ncol(cc$y) ,"\n\n")
  explained <- apply(yp, 2, var) / apply(cc$y, 2, var)
#  cat("T expl",mean(explained[it]), "P expl",mean(explained[ip]) ,"\n")
  explained2 <- apply(yp, 2, var) / apply(ydash, 2, var)
#  cat("T expl",mean(explained2[it]), "P expl",mean(explained2[ip]) ,"\n")
  
  pdf(file=file1)
  plot(normality, ylim=c(0,1))
  points(ksnormality, col="red")
  
  par(mfrow=c(4,2))
  for (i in 1:ncol(cc$y)) {
  plot(cc$y[,i],yp[,i],main=colnames(cc$y)[i], xlab="Actual", ylab="Prediction")	#,ylim=c(0,8),xlim=c(0,8))
#  lines(c(-1000,1000), c(-1000,1000))
  abline(a=0, b=1)
#  plot(y[,i],y[,i]-yp[,i])	#,ylim=c(0,8),xlim=c(0,8))
  residuals <- ydash[,i]-yp[,i]
  plot(yp[,i], residuals, xlab="Prediction", ylab="Residual", main=paste("sd = ", sd(residuals), sep=""))
  abline(h=0)
#  plot(y[,i], ySvdErr[,i])
  }
  if (sum(normality < 0.05)/ncol(cc$y) > 0.15) {
    for (i in 1:ncol(cc$y)) {
	  if (normality[i] < 0.05) {
	    qqnorm(ydash[,i]-yp[,i], main=colnames(cc$y)[i])
	    qqline(ydash[,i]-yp[,i])
	  }	
	}
  }

dev.off()
return(list(normality=normality,yp=yp,ydash=ydash,diff=ydash-yp))
}
