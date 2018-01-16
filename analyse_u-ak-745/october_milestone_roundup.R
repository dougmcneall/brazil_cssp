# october_milestone_roundup.R
# A few bits of code to fill in the gaps for the 
# october 2017 milestone report (presentation)

setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
## Try and build an emulator
library(devtools)

#install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)
library(coefplot)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(9, "RdYlBu")
byr = rev(ryb)


# Here is the leaf area index data
setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
source("/Users/dougmcneall/Documents/work/R/brazil_cssp/per_pft.R")


filelist.global.vegfrac <- paste0('frac_area_mean/global_area_mean_PFT',0:16, '.txt')
filelist.wus.vegfrac <- paste0('frac_area_mean/WUS_area_mean_PFT',0:16, '.txt')
filelist.sam.vegfrac <- paste0('frac_area_mean/SAM_area_mean_PFT',0:16, '.txt')

filelist.global.lai <- paste0('lai_area_mean/global_area_mean_lai_PFT',0:12, '.txt')
filelist.wus.lai <- paste0('lai_area_mean/WUS_area_mean_lai_PFT',0:12, '.txt')
filelist.sam.lai <- paste0('lai_area_mean/SAM_area_mean_lai_PFT',0:12, '.txt')

c.df <- function(fn){
  out = c(read.table(fn, header = FALSE, skip = 1), recursive = TRUE)
  out
}

# creates a list of (PFTs) long, each element of which has (ensemble members) data points
global_area_means_vegfrac.list <- lapply(filelist.global.vegfrac, c.df)
global_area_means_vegfrac_standard <- as.numeric(
  readLines('frac_area_mean/global_area_mean_PFTs_standard.txt'))

wus_area_means_vegfrac.list <- lapply(filelist.wus.vegfrac, c.df)
wus_area_means_vegfrac_standard <- as.numeric(
  readLines('frac_area_mean/WUS_area_mean_PFTs_standard.txt'))

sam_area_means_vegfrac.list <- lapply(filelist.sam.vegfrac, c.df)
sam_area_means_vegfrac_standard <- as.numeric(
  readLines('frac_area_mean/SAM_area_mean_PFTs_standard.txt'))


# creates a list of (PFTs) long, each element of which has (ensemble members) data points
global_area_means_lai.list <- lapply(filelist.global.lai, c.df)
#global_area_means_standard <- as.numeric(
#  readLines('lai_area_mean/global_area_mean_lai_PFTs_standard.txt'))

wus_area_means_lai.list <- lapply(filelist.wus.lai, c.df)
#wus_area_means_standard <- as.numeric(
#  readLines('lai_area_mean/WUS_area_mean_lai_PFTs_standard.txt'))

sam_area_means_lai.list <- lapply(filelist.sam.lai, c.df)
#sam_area_means_standard <- as.numeric(
#  readLines('lai_area_mean/SAM_area_mean_lai_PFTs_standard.txt'))



# How are LAI and vegetation fraction related?
# Turns out, not very strongly.
plot(sam_area_means_lai.list[[3]], sam_area_means_vegfrac.list[[3]])

sam_frac = c(sam_area_means_vegfrac.list[[2]], recursive = TRUE)

pdf(file = 'sam_param_vs_tffraction.pdf', width = 10, height = 6)
par(mfrow = c(6,13), mar = c(2,0.5,0.5,0.5))

for(i in 1:d){
  plot(lhs[, i], sam_frac, axes=FALSE, xlab = '', ylab = '', pch = 20, cex = 0.8)
  mtext(side = 1, text = cn[i], line = 0.3, cex = 0.7)
}
dev.off()

# Build a stepwise model of SAM BLET and find where the "good"
# (or at least functional) forests are.

lhs <- read.table('lhs_u-ak745.txt', header = TRUE)
d <- ncol(lhs)
cn <- colnames(lhs)
X.norm <- normalize(lhs)

y = sam_frac
dat = data.frame(y=y, x=X.norm)
colnames(dat) = c('y', colnames(X.norm))

initfit = lm(y ~ ., data = dat)
stepfit = step(initfit, direction="both", k=log(length(y)), trace=TRUE)

# Sample from the emulator
nsamp = 100000
X.unif = data.frame(samp.unif(nsamp, mins = rep(0,d), maxes = rep(1,d)))
colnames(X.unif) = colnames(X.norm)

pred.unif = predict(stepfit, newdata = X.unif)

thres = 0.1

ix.kept = which(pred.unif > thres)

hist(pred.unif[ix.kept])
ix.rejected = which(pred.unif < thres)

stepvars = names(stepfit$coefficients)[-1]

stepvars.ix = which(colnames(X.norm)%in% stepvars)

# how much space is rejected when we remove those samples without a functioning
# forest?
prop.rejected = length(ix.rejected) / nsamp

kept = X.unif[ix.kept, stepvars.ix]


dfunc.up <- function(x,y,...){
  require(MASS)
  require(RColorBrewer)
  
  rb = brewer.pal(9, "RdBu")
  br  = rev(rb)
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde = kde2d(x,y)
  image(kde, col = rb, add = TRUE)
}

rb = brewer.pal(9, "RdBu")
br  = rev(rb)

pdf(file = 'sam_blet_stepwise_density.pdf', width = 10, height = 10)
#dev.new()
par(oma = c(0,0,0,3))
pairs(kept,
      gap = 0, lower.panel = NULL, xlim = c(0,1), ylim = c(0,1),
      cex.labels = 0.5,
      panel = dfunc.up)

image.plot(legend.only = TRUE,
           legend.shrink = 0.7,
           zlim = c(0,1),
           col = rb,
           legend.args = list(text = 'Density of model runs with a forest', side = 3, line = 1),
           horizontal = TRUE
)
dev.off()




