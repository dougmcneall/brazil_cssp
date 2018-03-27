# ensembleFunctions.R
# useful functions for working with ensemble data



# HadCM3 lats and longs
HadCM3.lats <- seq(from = -90, to = 90, by = 2.5)
HadCM3.lons <- seq(from = 0, by = 3.75, length = 96)


processTS <- function(file, rlfile, yearvec,exclude = FALSE){
  # read in and process an aggregated
  # timeseries matrix for plotting
  
  dat <- read.table(file, na.strings = c('0.0000000', 'NA'))
  runlist <- scan(file = rlfile, what = c('character'))
  year <- yearvec

  colnames(dat) <- year
  row.names(dat) <- runlist
  dat

  if(exclude){
    dat[-exclude, ]
  }
  
}


anomalizeTS <- function(dat, i){
  ## anomalize and ensemble relative to a part of the ensemble
  subdat <- dat[ , i]

  sweepstats <- apply(subdat,1, FUN = mean)

  anom <- sweep(dat, 1,sweepstats, FUN = '-')
  anom
}




percChange <- function(a,b){
  # percentage change
  
 out <- ((b-a)/a)*100

 out
}




percChangeTS <- function(dat, i){
  ## normalize wrt a period of a timeseries
  ## returns percentage change from the
  ## reference period (the vector i)

  subdat <- dat[ , i]
  
  sweepstats <- apply(subdat,1, FUN = mean)

  anom <- sweep(dat, 1,sweepstats, FUN = '-')

  percanom <- 100 * sweep(anom, 1,sweepstats, FUN = '/')

  percanom
  
}
  




ensCompTS <- function(dat1,dat2,grid = TRUE,colvec,rugvec,legvec,mainvec,...){
  # Plot comparison ensemble time series

  
  source('/home/h01/hadda/code/R/useful/dougplotpars.R')  
  par(dougpar_web)

  matplot(colnames(dat1), t(dat1), type = 'l',
          lty = 1,
          col = colvec[1],
          ...
          )
  
  points(rep(rugvec[1],nrow(dat1)), dat1[,ncol(dat1)], pch = '-', col = colvec[1])

  matlines(colnames(dat2), t(dat2), type = 'l',
           lty = 1,
           col = colvec[2]
           )
  
  points(rep(rugvec[2],nrow(dat2)), dat2[,ncol(dat2)], pch = '-', col = colvec[2])
  
  legend('topleft', legvec,
         col = colvec,
         text.col = colvec,
         lwd = 1.5
         )

  
  mtext(side = 3, line = 1, adj = 0, mainvec, cex = 2, col = 'black')
  if(grid) {grid(lty ='dashed',col = 'grey')}
  

}



ensCompTShist <- function(dat1,dat2,grid = TRUE,colvec,legvec,mainvec,...){
  # Plot comparison ensemble time series
  # add a histogram on the end

  
  source('/home/h01/hadda/code/R/useful/dougplotpars.R')  
  par(dougpar_web)
  par(mar = c(5,5,4,0), mgp = c(3.5,1,0))

  
  nf <- layout(matrix(c(1,2,3),1,3,byrow=TRUE),widths = c(10,1,1), TRUE)
  #layout.show(nf)
  matplot(colnames(dat1), t(dat1), type = 'l',
          lty = 1,
          col = colvec[1],
          ...
          )

  matlines(colnames(dat2), t(dat2), type = 'l',
           lty = 1,
           col = colvec[2]
           )


  if(grid) {grid(lty ='dashed',col = 'grey')}

  mtext(side = 3, line = 1, adj = 0, mainvec, cex = 2, col = 'black')
  
  legend('topleft', legvec,
         fill = colvec,
         text.col = colvec,
         bg = 'white',
         border = par()$fg,
         cex = 1.5
         )

  
  # Add the histograms
  datran <- range(dat1,dat2, na.rm = TRUE)

  breaks <- seq(from = datran[1], to = datran[2], length = 15)
  
  dat1Hist <- hist( dat1[,ncol(dat1)],breaks = breaks, plot = FALSE)
  dat2Hist <- hist( dat2[,ncol(dat2)],breaks = breaks,  plot = FALSE)
  
  xlim = c(0, max(dat1Hist$counts, dat2Hist$counts))
  par(mar = c(5,0,4,1), fg = 'white')
  barplot(dat1Hist$counts, horiz = TRUE, col = colvec[1], space = 0, axes = FALSE, xlim = xlim)
  barplot(dat2Hist$counts, horiz = TRUE, col = colvec[2], space = 0, axes = FALSE, xlim = xlim)
    
}



loessSmooth <- function(dat,f = 0.3){
  # basic loess smoother for time series (vector):  
  
  datc <- c(dat, recursive = TRUE)
  names(datc) <- NULL
  
  x <- 1:length(datc)
  l <- loess(datc~x)
  
  dat.df <- data.frame(x,datc)

  l <- predict(loess(datc~x, data = dat.df, span = f, na.action = na.exclude))
  
  l 
}


loessSmoothdf <- function(dat,f = 0.3){
  # apply the loess smoother to the rows of a data frame

  lsmooth <- data.frame(t(apply(dat,1,loessSmooth,f = f)))
  colnames(lsmooth) <- colnames(dat)
  lsmooth
}



findTurnover <- function(smoothedDat){
  # find the year in which a smoothed data set
  # 'turns over' to a negative gradient
  # Use on a smoothed timeseries
  
  datDiff <- t(apply(smoothedDat,1,diff))
  datSign <- sign(datDiff)

  nens <- nrow(smoothedDat)

  turnover.ix <- rep(NA,nens)
  
  for(i in 1:nens){
    
    vec <- datSign[i,]
    
    turnover.ix[i] <- match(-1, vec)
    
  }
  
  turnover.year <- as.numeric(colnames(datSign)[turnover.ix])
  turnover.year
  
}


vec2map <- function(vec, nr = 73, nc = 96){
  # outputs matrix of data that can be fed to image.plot,
  # or image for plotting using
  # Default HadCM3 rows (73) and columns (96)

  mapraw <- matrix(vec, nrow = nr, ncol = nc, byrow = TRUE)
  map <- t(mapraw[nrow(mapraw):1, ])
  map
  
  
}


matAggDec <- function(mat, ndeltat, FUN, cn = NULL, ...){

  # aggregate an ensemble matrix of timeseries
  # to non-overlapping blocks
  # not tested yet.
  
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  newnc <- (nc/ndeltat)
  
  longmat <- matrix(c(t(mat), recursive = TRUE), ncol = ndeltat, nrow = newnc * nr, byrow = TRUE)

  decblocks <- apply(longmat,1,FUN,...)
  
  aggout <- matrix(decblocks, ncol = newnc, nrow = nr, byrow = TRUE)
  colnames(aggout) <- cn

  aggout
}


swapRun <- function(tabfile = '/home/h01/hadda/umui_jobs/run_management/kitchenSinkRuntab.txt',
                    na.strings = '99',
                    RUNIDS,inscen, outscen)
       {
   # Translate RUNIDs from one scenario to 
   # another in the Earth system ensemble
   #
   # INPUTS:
   # tabfile ... Experiment translation table
   # RUNIDS  ... Vector if RUNIDS
   # inscen  ... Scenario of the input RUNIDS
   # outscen ... Scenario of the output RUNIDS


   # translation table
   tab <- read.table(file = tabfile, na.strings = na.strings, header = TRUE)

   ix <- match(RUNIDS, tab[, inscen])

   out <- tab[ix,outscen]

   out
  
}



yearSelect <- function(dat,yearvec){
  # select columns from an enemble data matrix
  # based on a vector of years
  ix <- colnames(dat) %in% yearvec
  dat.trunc <- dat[ , ix]
  dat.trunc

}
