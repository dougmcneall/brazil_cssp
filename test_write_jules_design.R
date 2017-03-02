# test_write_jules_design.R
# Tests and examples for write_jules_design.R
# Doug McNeall dougmcneall@gmail.com

fac = c('g_root_io', 'retran_l_io')
minfac = c(0.5, 0.5)
maxfac = c(2,2)
tf = 'l_veg_soil'

write_jules_design(paramlist, n = 10, fac = fac, minfac = minfac, maxfac = maxfac, tf = tf)

