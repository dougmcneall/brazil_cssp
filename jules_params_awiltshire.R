# Andy Wiltshire's JULES Parameter choices
# Doug McNeall dougmcneall@gmail.com


hw_sw_io = list(
  'standard' = rep(0.5, 13),
  'min' = rep(0, 13),
  'max' = rep(1,13),
  'namelist' = 'jules_pftparm'
)

knl_io = list(
  'standard' = rep(0.2, 13),
  'min' = rep(0.1, 13),
  'max' = rep(0.5,13),
  'namelist' = 'jules_pftparm'
)

a_wl_io = list(
  'standard' = c(0.78, 0.845, 0.78, 0.8, 0.65, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.2, 0.13),
  'min' = 0.5 * c(0.78, 0.845, 0.78, 0.8, 0.65, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.2, 0.13),
  'max' = 2* c(0.78, 0.845, 0.78, 0.8, 0.65, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.2, 0.13),
  'namelist' = 'jules_pftparm'
)

g_root_io = list(
  'standard' = c(0.15, 0.25, 0.25, 0.15, 0.15, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.15, 0.15),
  'min' = rep(0, 13),
  'max' = rep(0.5, 13),
  'namelist' = 'jules_triffid'
)

g_wood_io = list(
  'standard' = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05, 0.05),
  'min' = rep(0,13),
  'max'      = c(0.05, 0.05, 0.05, 0.05, 0.05, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1),
  'namelist' = 'jules_triffid' 
)

retran_l_io = list(
  'standard' = c(0.5, 0.5, 0.5, 0.77, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  'min' = rep(0,13),
  'max' = rep(1,13),
  'namelist' = 'jules_triffid'
)

retran_r_io = list(
  'standard' = rep(0.2,13),
  'min' = rep(0,13),
  'max' = rep(1,13),
  'namelist' = 'jules_triffid'
)

kaps_roth = list(
  'standard' = c(1.61E-07, 4.83E-09, 1.06E-08, 3.22E-10),
  'min' = 0.5 * c(1.61E-07, 4.83E-09, 1.06E-08, 3.22E-10),
  'max' = 2 * c(1.61E-07, 4.83E-09, 1.06E-08, 3.22E-10),
  'namelist' = 'jules_soil_biogeochem'
)

n_inorg_turnover = list(
  'standard' = 1,
  'min' = 0,
  'max' = 10,
  'namelist' = 'jules_soil_biogeochem'
)

sorp = list(
  'standard' = 10,
  'min' = 0,
  'max' = 20,
  'namelist' = 'jules_soil_biogeochem'
)

# test
l_vg_soil = list(
  'standard' = '.true.',
  'min' = '.false.',
  'max' = '.true.',
  'namelist' = 'jules_soil'
)

paramlist = list('g_root_io' = g_root_io,
                 'g_wood_io' = g_wood_io,
                 'retran_l_io' = retran_l_io,
                 'retran_r_io' = retran_r_io,
                 'hw_sw_io' = hw_sw_io, 
                 'knl_io' = knl_io, 
                 'a_wl_io' = a_wl_io,
                 'kaps_roth' = kaps_roth,
                 'n_inorg_turnover' = n_inorg_turnover,
                 'sorp' = sorp,
                 'l_vg_soil' = l_vg_soil
)
