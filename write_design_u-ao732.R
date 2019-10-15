
source('write_jules_design2.R')
source('jules_params_u-ao732.R')


# Easiest way to generate a design of the right size is to have a "fac" which takes
# the names from the parameter list, and then multiplies everything by 0.5 or 2

tf = 'l_vg_soil'
# we don't want anything that is TRUE/FALSE to be in fac
fac.init = names(paramlist)
not_tf.ix = which(names(paramlist)!=tf)
paramlist.trunc = paramlist[not_tf.ix]

fac = names(paramlist.trunc)

maxfac = lapply(paramlist.trunc,function(x) x$max[which.max(x$max)] / x$standard[which.max(x$max)])
minfac = lapply(paramlist.trunc,function(x) x$min[which.max(x$max)] / x$standard[which.max(x$max)])

write_jules_design2(paramlist, n = 300,
                    fac = fac, minfac = minfac, maxfac = maxfac,
                    tf = tf,
                    fnprefix = 'conf_files_u-ao732/param-perturb-P',
                    lhsfn = 'conf_files_u-ao732/lhs_u-ao732.txt',
                    stanfn = 'conf_files_u-ao732/stanparms_u-ao732.txt',
                    allstanfn = 'conf_files_u-ao732/allstanparms_u-ao732.txt',
                    rn = 12)

print(cbind(minfac, maxfac))


# Checking section
lapply(paramlist, function(x) length(x$standard))
lapply(paramlist, function(x) length(x$standard)==length(x$min) & length(x$standard)==length(x$max))



#lhs = read.table('conf_dummy/lhs.txt', head = TRUE)

#stan = read.table('conf_dummy/stanparms.txt', head = TRUE)
#allstan = read.table('conf_dummy/allstanparms.txt', head = TRUE)

# Check that what we've created is consistent with the design that we created before.

#lhs_u_ak745 <- read.table('analyse_u-ak-745/lhs_u-ak745.txt', head = TRUE)


#dev.new(width = 10, height = 10)
#par(mfrow = c(10, 8), mar = c(1,1,1,1))


#for(i in 1 : ncol(lhs_u_ak745)){

#  dat = lhs_u_ak745[, i]
#  plot(range(dat))
#  points(stan[i], col = 'red')
 # hist(lhs_u_ak745[, i], breaks = 10, axes = FALSE, xlab = '', ylab = '')
   
#}
