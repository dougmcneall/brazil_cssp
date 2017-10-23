# create_design_textfile.R

setwd("/Users/dougmcneall/Documents/work/R/brazil_cssp/analyse_u-ak-745")
library(MASS)

lhs <- read.table('lhs_u-ak745.txt', header = TRUE)

rownames(lhs) <- paste0("P", sprintf("%04d", 0:99))

write.table(lhs, file = 'u-ak745_design.txt')
write.csv(lhs, file = 'u-ak745_design.csv')


