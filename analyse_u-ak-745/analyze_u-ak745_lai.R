# analyze_u-ak745_lai.R
setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
source("/Users/dougmcneall/Documents/work/R/brazil_cssp/per_pft.R")


filelist.global.lai <- paste0('lai_area_mean/global_area_mean_lai_PFT',0:12, '.txt')
filelist.wus.lai <- paste0('lai_area_mean/WUS_area_mean_lai_PFT',0:12, '.txt')
filelist.sam.lai <- paste0('lai_area_mean/SAM_area_mean_lai_PFT',0:12, '.txt')

c.df <- function(fn){
  out = c(read.table(fn, header = FALSE, skip = 1), recursive = TRUE)
  out
}

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


# ------------------------------------------------------
# Build emulators
# ------------------------------------------------------
lhs <- read.table('lhs_u-ak745.txt', header = TRUE)
d <- ncol(lhs)
cn <- colnames(lhs)
lhs.norm <- normalize(lhs)

global.sensmat.lai <- globalsens(lhs, global_area_means_lai.list)
global.sensmat.lai.co <- cutoff(global.sensmat, c(-0.11, 0.11))

wus.sensmat.lai <- globalsens(lhs, wus_area_means_lai.list)
wus.sensmat.lai.co <- cutoff(wus.sensmat.lai, c(-0.11, 0.11))

sam.sensmat.lai <- globalsens(lhs, sam_area_means_lai.list)
sam.sensmat.lai.co <- cutoff(sam.sensmat.lai, c(-0.11, 0.11))


# sensmat.summary <- function(sensmat, cutoff){
#   
#  summary.vec <- apply(abs(sensmat),1, sum, na.rm = TRUE)
#  summary.vec.ix <- order(summary.vec, decreasing = TRUE)
#  summary.ordered <- summary.vec[summary.vec.ix]
#   
# }

wus.sens.lai.summary <- apply(abs(wus.sensmat.lai),1, sum, na.rm = TRUE)
wus.sens.lai.ix <- order(wus.sens.lai.summary, decreasing = TRUE)

sam.sens.lai.summary <- apply(abs(sam.sensmat.lai),1, sum, na.rm = TRUE)
sam.sens.lai.ix <- order(sam.sens.lai.summary, decreasing = TRUE)

global.sens.lai.summary <- apply(abs(global.sensmat.lai),1, sum, na.rm = TRUE)
global.sens.lai.ix <- order(global.sens.lai.summary, decreasing = TRUE)

pdf(file = 'response_summary_step_global_lai_both_measures.pdf', width = 12, height = 8)
par(mar = c(4,7,2,7), mfrow = c(2,1), oma = c(3,0,0,0))
image(global.sensmat.co, col = br, axes = FALSE, main = 'Global')
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/12) , labels = pfts, las = 1)

image.plot(global.sensmat.lai.co, col = br, legend.only = TRUE)

plot(global.sens.lai.summary[global.sens.lai.ix], axes = FALSE, 
     xlab = '', ylab = 'Sensitivity summary', ylim = c(0,1.5))
axis(1, labels=(colnames(lhs))[global.sens.lai.ix], at = 1:length(global.sens.lai.summary), las = 3, cex.axis = 0.8 )
axis(2, las = 1)

dev.off()

pdf(file = 'response_summary_step_wus_lai_both_measures.pdf', width = 12, height = 8)
par(mar = c(4,7,2,7), mfrow = c(2,1), oma = c(3,0,0,0))
image(wus.sensmat.co, col = br, axes = FALSE, main = 'WUS')
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/12) , labels = pfts, las = 1)

image.plot(wus.sensmat.lai.co, col = br, legend.only = TRUE)

plot(wus.sens.lai.summary[wus.sens.lai.ix], axes = FALSE, 
     xlab = '', ylab = 'Sensitivity summary', ylim = c(0,1.5))
axis(1, labels=(colnames(lhs))[wus.sens.lai.ix], at = 1:length(wus.sens.lai.summary), las = 3, cex.axis = 0.8 )
axis(2, las = 1)

dev.off()

pdf(file = 'response_summary_step_sam_lai_both_measures.pdf', width = 12, height = 8)
par(mar = c(4,7,2,7), mfrow = c(2,1), oma = c(3,0,0,0))
image(sam.sensmat.co, col = br, axes = FALSE, main = 'SAM')
axis(1, at = seq(from = 0, to = 1, by = 1/(d-1)), labels = colnames(lhs), las = 3, cex.axis = 0.8)
axis(2, at = seq(from =0, to = 1, by = 1/12) , labels = pfts, las = 1)

image.plot(sam.sensmat.lai.co, col = br, legend.only = TRUE)

plot(sam.sens.lai.summary[sam.sens.lai.ix], axes = FALSE, 
     xlab = '', ylab = 'Sensitivity summary', ylim = c(0,1.5))
axis(1, labels=(colnames(lhs))[sam.sens.lai.ix], at = 1:length(sam.sens.lai.summary), las = 3, cex.axis = 0.8 )
axis(2, las = 1)

dev.off()


