# jules_params_ga7.R
# Jules parameters used in GA7, from D Sexton.

# Note, it appears that GA7 uses a different JULES setup and default parameters. These parameters
# are found by multiplying a standard by the factor used to multiply the GA7 parameters.

tupp_io_frac = c(25,41) / 36
f0_io_frac = c(0.65, 0.972) / 0.875
dz0v_dh_io_frac = c(0, 0.16) / 0.05
nl0_io_frac = c(0, 0.12) / 0.04
rootd_ft_io_frac = c(0, 8) /3

tupp_io = list(
  'standard' = c(43, 43, 43, 26, 32, 32, 32, 32, 45, 45, 45, 40, 36),
  'min' = tupp_io_frac[1] * c(43, 43, 43, 26, 32, 32, 32, 32, 45, 45, 45, 40, 36),
  'max' = tupp_io_frac[2] * c(43, 43, 43, 26, 32, 32, 32, 32, 45, 45, 45, 40, 36)
  'namelist' = 'jules_pftparm'
)

f0_io = list(
  'standard' = c(0.875, 0.875, 0.892, 0.875, 0.875, 0.931, 0.931, 0.931, 0.8, 0.8, 0.8, 0.875, 0.875),
  'min' =  f0_io_frac[1] * c(0.875, 0.875, 0.892, 0.875, 0.875, 0.931, 0.931, 0.931, 0.8, 0.8, 0.8, 0.875, 0.875),
  'max' = f0_io_frac[2] * c(0.875, 0.875, 0.892, 0.875, 0.875, 0.931, 0.931, 0.931, 0.8, 0.8, 0.8, 0.875, 0.875),
  'namelist' = 'jules_pftparm'
)

dz0v_dh_io = list(
  'standard' = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.1 , 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  'min' = rep(0, 13),
  'max' =  dz0v_dh_io_frac[2] * c(0.05, 0.05, 0.05, 0.05, 0.05, 0.1 , 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
  'namelist' = 'jules_pftparm'
)

nl0_io = list(
  'standard' = c(0.046, 0.046, 0.046, 0.033, 0.033, 0.073, 0.073, 0.073, 0.06, 0.06, 0.06, 0.06, 0.06),
  'min' = rep(0,13),
  'max' = nl0_io_frac[2] * c(0.046, 0.046, 0.046, 0.033, 0.033, 0.073, 0.073, 0.073, 0.06, 0.06, 0.06, 0.06, 0.06),
  'namelist' = 'jules_pftparm'
)

rootd_ft_io = list(
  'standard' = c(2, 3, 2, 2, 1.8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  'min' = rep(0, 13),
  'max' = rootd_ft_io_frac[2] * c(2, 3, 2, 2, 1.8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  'namelist' = 'jules_pftparm'
)


