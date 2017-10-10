

source('write_jules_design2.R')
#source('../jules_params_u-ao732.R')


# display standard settings for all pfts
#lapply(paramlist, function(x) x$standard)

#fac = c('g_root_io', 'retran_l_io')
#minfac = c(0.5, 0.5)
#maxfac = c(2,2)
#tf = 'l_vg_soil'

#write_jules_design2(paramlist, n = 100, fac = fac, minfac = minfac, maxfac = maxfac, tf = tf, fnprefix = 'param-perturb-dummy')

source('jules_params_awiltshire.R')
fac = c('g_root_io', 'retran_l_io')
minfac = c(0.5, 0.5)
maxfac = c(2,2)
tf = 'l_vg_soil'

write_jules_design2(paramlist, n = 100,
                    fac = fac, minfac = minfac, maxfac = maxfac,
                    tf = tf,
                    fnprefix = 'conf_dummy/param-perturb-test',
                    lhsfn = 'conf_dummy/lhs.txt',
                    stanfn = 'conf_dummy/stanparms.txt',
                    allstanfn = 'conf_dummy/allstanparms.txt')

lhs = read.table('conf_dummy/lhs.txt', head = TRUE)

stan = read.table('conf_dummy/stanparms.txt', head = TRUE)
allstan = read.table('conf_dummy/allstanparms.txt', head = TRUE)

# Check that what we've created is consistent with the design that we created before.

lhs_u_ak745 <- read.table('analyse_u-ak-745/lhs_u-ak745.txt', head = TRUE)


dev.new(width = 10, height = 10)
par(mfrow = c(10, 8), mar = c(1,1,1,1))


for(i in 1 : ncol(lhs_u_ak745)){

  dat = lhs_u_ak745[, i]
  plot(range(dat))
  points(stan[i], col = 'red')
 # hist(lhs_u_ak745[, i], breaks = 10, axes = FALSE, xlab = '', ylab = '')
  
  
}
