# test_write_jules_design.R
# Tests and examples for write_jules_design.R
# Doug McNeall dougmcneall@gmail.com


# Source design writing function
#source('https://raw.githubusercontent.com/dougmcneall/brazil_cssp/master/write_jules_design.R')
# It's probably easier and safer for testing to source the local file.
source('write_jules_design.R')

# Create some "dummy" parameters that are easy to test.
param_a = list(
  'standard' = rep(1, 13),
  'min' = rep(0, 13),
  'max' = rep(10,13),
  'namelist' = 'jules_pftparm'
)

param_b = list(
  'standard' = c(0.001, 0.001, 0.001, 0.1, 0.1, 0.1, 1, 1, 1,10, 10, 10, 100),
  'min' = rep(0, 13),
  'max' = 10 * c(0.001, 0.001, 0.001, 0.1, 0.1, 0.1, 1, 1, 1,10, 10, 10, 100),
  'namelist' = 'jules_triffid'
)

param_c = list(
  'standard' = 1,
  'min' = 0,
  'max' = 10,
  'namelist' = 'jules_surface'
)

param_d = list(
  'standard' = '.true.',
  'min' = '.false.',
  'max' = '.true.',
  'namelist' = 'jules_soil'
)


paramlist = list(
  'param_a' = param_a,
  'param_b' = param_b,
  'param_c' = param_c,
  'param_d' = param_d
)

fac = c('param_b')
minfac = c(0.5)
maxfac = c(2)
tf = 'param_d'

write_jules_design(paramlist, n = 10, fac = fac, minfac = minfac, maxfac = maxfac, tf = tf, rn = 5)



# Load Andy Wiltshire's paramter choices
source('https://raw.githubusercontent.com/dougmcneall/brazil_cssp/master/jules_params_awiltshire.R')

fac = c('g_root_io', 'retran_l_io')
minfac = c(0.5, 0.5)
maxfac = c(2,2)
tf = 'l_veg_soil'

write_jules_design(paramlist, n = 100, fac = fac, minfac = minfac, maxfac = maxfac, tf = tf, fnprefix = 'param-perturb-test')

