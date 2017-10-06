# create_design_textfile.R

setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")

lhs <- read.table('lhs_u-ak745.txt', header = TRUE)

rownames(lhs) <- paste0("P", sprintf("%04d", 0:99))
