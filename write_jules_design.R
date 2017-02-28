# write_jules_design.R
# Write a design matrix and push the parameters to 
# configuration files.

library(lhs)
library(MASS)
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")

# ----------------------------------------------------
# Create a list of parameters
# ----------------------------------------------------

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
  'namelist' = 'jules_surface'
)

n_inorg_turnover = list(
  'standard' = 1,
  'min' = 0,
  'max' = 10,
  'namelist' = 'jules_surface'
)

sorp = list(
  'standard' = 10,
  'min' = 0,
  'max' = 20,
  'namelist' = 'jules_surface'
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
                 'sorp' = sorp
)


write_jules_design = function(paramlist, n, fac, minfac, maxfac, fnprefix = 'test',
                              lhsfn = 'lhs.txt', rn = 3){
  # This code writes a design taking either a 'factor', min and max by which
  # to multiply all pfts, or perturbing each pft individually according to
  # their maximum and minimum in the parameter list.
  # fac is a character vector of names of variables that you would like to alter
  # by a factor. Everything else gets variaed by PFT
  # minfac and maxfac must correspond to fac - i.e. one value per parameter, in the
  # correct order.
  
  paramvec = names(paramlist)
  nmlvec = unlist(lapply(paramlist, FUN = function(x) x$namelist))
  
  # which parameters do we want as a parameter list?
  fac.ix = which(names(paramlist) %in% fac)
  paramfac = paramlist[fac.ix]
  parampft = paramlist[-fac.ix]
  pftvec = names(parampft)
  
  parampft_nml = unlist(lapply(parampft, FUN = function(x) x$namelist))
  paramfac_nml = unlist(lapply(paramfac, FUN = function(x) x$namelist))
  
  parampft_standard = unlist(lapply(parampft, FUN = function(x) x$standard))
  parampft_mins = unlist(lapply(parampft, FUN = function(x) x$min))
  parampft_maxes = unlist(lapply(parampft, FUN = function(x) x$max))
  
  all_names = c(names(parampft_standard), fac)
  k = length(all_names)
  
  # Do the pfts and then the factors
  lhs = unnormalize(
    maximinLHS(n = n, k = k, dup = 1),
    un.mins = c(parampft_mins , minfac),
    un.maxes = c(parampft_maxes, maxfac)
  )
  colnames(lhs) = all_names
  
  for(i in 1:nrow(lhs)){
    fn = paste0(fnprefix,i,'.txt')
    
    for(el in unique(nmlvec)){
      write(paste0('[namelist:',el,']'), file = fn, append = TRUE)
      
      # grab the parts of the list that match
      pft_elms = parampft[parampft_nml==el] # & statement
      pft_elms_vec = names(pft_elms)
      fac_elms = paramfac[paramfac_nml==el]
      fac_elms_vec = names(fac_elms)
      
      if(length(pft_elms_vec) > 0){
        for(j in 1:length(pft_elms_vec)){
          param = pft_elms_vec[j]
          colix = grep(param, colnames(lhs))
          values.out = lhs[i, colix]
          write(paste0(param,'=',(paste0(round(values.out,rn), collapse = ',')), collapse = ''),
                file = fn, append = TRUE)
        }
      }
      
      if(length(fac_elms_vec) > 0){
        for(k in 1: length(fac_elms_vec)){
          param = fac_elms_vec[k]
          colix = grep(param, colnames(lhs))
          lhs.factor = lhs[i, colix]
          values.out = lhs.factor * get(param, paramlist)$standard
          write(paste0(param,'=',(paste0(round(values.out,rn), collapse = ',')), collapse = ''),
                file = fn, append = TRUE)
        }
      }
      write('/', file = fn, append = TRUE)
    }
  }
  write.matrix(lhs, file = lhsfn)
}

# vary these parameters only by a factor
fac = c('g_root_io', 'retran_l_io')
minfac = c(0.5, 0.5)
maxfac = c(2,2)

write_jules_design(paramlist, n = 10, fac = fac, minfac = minfac, maxfac = maxfac)


# -----------------------------------------------------------
# Development code from here on
# -----------------------------------------------------------

write_jules_design_bypft = function(paramlist, n, fnprefix = 'test',
                                    lhsfn = 'lhs.txt', rn = 3){
  # A function for generating a latin hypercube design for individual PFTs in a number of parameters
  # paramlist .......... A list of named parameter lists, each sublist with standard, min and max,
  #           .......... and namelist
  
  paramvec = names(paramlist)
  all_standard = unlist(lapply(paramlist, FUN = function(x) x$standard))
  all_mins = unlist(lapply(paramlist, FUN = function(x) x$min))
  all_maxes = unlist(lapply(paramlist, FUN = function(x) x$max))
  nmlvec = unlist(lapply(paramlist, FUN = function(x) x$namelist))
  
  k = length(all_standard)
  
  lhs = unnormalize(
    maximinLHS(n = n, k = k, dup = 1),
    un.mins = all_mins,
    un.maxes = all_maxes
  )
  colnames(lhs) = names(all_standard)
  
  for(i in 1:nrow(lhs)){
    
    fn = paste0(fnprefix,i,'.txt')
    
    for(el in unique(nmlvec)){
      write(paste0('[namelist:',el,']'), file = fn, append = TRUE)
      # grab the parts of the list that match
      param_elms = paramlist[nmlvec==el]
      param_elms_vec = names(param_elms)
      
      for(j in 1:length(param_elms)){
        
        param = param_elms_vec[j]
        colix = grep(param, colnames(lhs))
        values.out = lhs[i, colix]
        write(paste0(param,'=',(paste0(round(values.out,rn), collapse = ',')), collapse = ''),
              file = fn, append = TRUE)
        
      }
      write('/', file = fn, append = TRUE)
    }
  }
  write.matrix(lhs, file = lhsfn)
}


write_jules_design_byparam = function(paramlist, n, minfac, maxfac, fnprefix = 'test',
                                      lhsfn = 'lhs.txt', rn = 3){
  
  paramvec = names(paramlist)
  nmlvec = unlist(lapply(paramlist, FUN = function(x) x$namelist))
  k = length(paramvec)
  
  lhs = unnormalize(
    maximinLHS(n = n, k = k, dup = 1),
    un.mins = minfac,
    un.maxes = maxfac
  )
  colnames(lhs) = paramvec
  
  for(i in 1:nrow(lhs)){
    
    fn = paste0(fnprefix,i,'.txt')
    
    for(el in unique(nmlvec)){
      write(paste0('[namelist:',el,']'), file = fn, append = TRUE)
      
      # grab the parts of the list that match
      param_elms = paramlist[nmlvec==el]
      param_elms_vec = names(param_elms)
      
      for(j in 1:length(param_elms)){
        
        param = param_elms_vec[j]
        colix = grep(param, colnames(lhs))
        
        lhs.factor = lhs[i, colix]
        values.out = lhs.factor * get(param, paramlist)$standard
        write(paste0(param,'=',(paste0(round(values.out,rn), collapse = ',')), collapse = ''),
              file = fn, append = TRUE)
      }
      write('/', file = fn, append = TRUE)
    }
  }
  write.matrix(lhs, file = lhsfn)
}

minfac = rep(0.5, length(paramlist))
maxfac = rep(2, length(paramlist))
write_jules_design_byparam(paramlist, minfac, maxfac, n = 10, rn = 3)

write_jules_design_bypft(paramlist, n = 100, rn = 10)

# -----------------------------------------------------------------
# This piece of code works if we are to alter all PFTs in a 
# parameter at the same time, by a single factor.
# 
# -----------------------------------------------------------------
# might have to shift everything up and down by a factor to keep
# it from getting crazy.

# Standard Parameter values
hw_sw_io = rep(0.5, 13)
knl_io = rep(0.2, 13)
a_wl_io = c(0.78, 0.845, 0.78, 0.8, 0.65, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.2, 0.13)

paramlist = list('hw_sw_io', 'knl_io', 'a_wl_io')

# build a latin hypercube of "factors", which you multiply the
# standard values by to get the ensemble member values

n = 9  # number of ensemble members 
# use unnormalize to translate between the lhs 0-1 range,
# and the appropriate range
lhs = unnormalize(
  maximinLHS(n = n, k = length(paramlist), dup = 1),
  un.mins = c(0.1,1,10),
  un.maxes = c(0.2,2,20)
)

colnames(lhs) = unlist(paramlist)

rn = 3  # rounding number
 
for(i in 1:nrow(lhs)){
  
  # Section header
  nml = "[namelist:jules_pftparm]"
  
  # output file name
  fn = paste0('test',i,'.txt')
  write(nml, file = fn )
  
  # loop over the parameters in paramlist
  for(j in 1:length(paramlist)){
    out = paste0(paramlist[[j]], '=',
                 paste0(round( (lhs[i, paramlist[[j]] ] * get(paramlist[[j]])), rn), collapse = ','),
                 collapse = '')
    
    write(out, file = fn, append = TRUE)
  }
  write('/', file = fn, append = TRUE)
}


# -------------------------------------------------------------
# Build a latin hypercube if we are changing individual
# parameters at a time
# 
# -------------------------------------------------------------

# How long would it take to build a latin hypercube?
ptm <- proc.time()
k = 110
lhs = maximinLHS(n = 10*k, k = k, dup = 1)
proc.time() - ptm

ens = c(60,70,80,90,100,110)
tim = c(5,10,16,26,40,58)
plot(ens,sqrt(tim))

slope = coef(lm(sqrt(tim)~ens))
((0.107*200) - 4.3)^2

# 200 parameters would take 300 seconds
# 300 around 12 minutes

getvarname <- function(x){
  # turns variable name into a string
  deparse(substitute(x))
}


# here are the standard parameters, along with some
# simple values for mins/maxes.
hw_sw_io = rep(0.5, 13)
hw_sw_io.min = hw_sw_io/2
hw_sw_io.max = hw_sw_io*2

knl_io = rep(0.2, 13)
knl_io.min = knl_io/2
knl_io.max = knl_io*2

a_wl_io = c(0.78, 0.845, 0.78, 0.8, 0.65, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.2, 0.13)
a_wl_io.min = a_wl_io/2
a_wl_io.max = a_wl_io*2

# This has the advantage of producing headers automatically!
params_standard = unlist(mget(unlist(paramlist)))

k = length(params_standard)

lhs = unnormalize(
  maximinLHS(n = 10*k, k = k, dup = 1),
  un.mins = c(hw_sw_io.min, knl_io.min, a_wl_io.min),
  un.maxes = c(hw_sw_io.max, knl_io.max, a_wl_io.max)
)

colnames(lhs) = names(params_standard)

# Now need to split the lhs back into its constituent
# parts, and write them individually to file.

paramvec = unlist(paramlist)

for(i in 1:nrow(lhs)){
  
  fn = paste0('test',i,'.txt')
  write(nml, file = fn )
  
  for(j in 1:length(paramvec)){
    
    param = paramvec[j]
    colix = grep(param, colnames(lhs))
    values.out = lhs[i, colix]
    write(paste0(param,'=',(paste0(round(values.out,rn), collapse = ',')), collapse = ''),
          file = fn, append = TRUE)
  }
  
  write('/', file = fn, append = TRUE)
  
}


write_jules_design_factor = function(paramlist, n, fnprefix = 'test', nml = "[namelist:jules_pftparm]",
                                     lhsfn = 'lhs.txt', rn = 3){
  # A function for generating a latin hypercube design for individual PFTs in a number of parameters
  # paramlist .......... A list of named parameter lists, each sublist with standard, min and max
  
  paramvec = names(paramlist)
  all_standard = unlist(lapply(paramlist, FUN = function(x) x$standard))
  all_mins = unlist(lapply(paramlist, FUN = function(x) x$min))
  all_maxes = unlist(lapply(paramlist, FUN = function(x) x$max))
  
  k = length(all_standard)
  
  lhs = unnormalize(
    maximinLHS(n = n, k = k, dup = 1),
    un.mins = all_mins,
    un.maxes = all_maxes
  )
  
  colnames(lhs) = names(all_standard)
  
  for(i in 1:nrow(lhs)){
    fn = paste0(fnprefix,i,'.txt')
    write(nml, file = fn )
    
    for(j in 1:length(paramvec)){
      param = paramvec[j]
      colix = grep(param, colnames(lhs))
      values.out = lhs[i, colix]
      write(paste0(param,'=',(paste0(round(values.out,rn), collapse = ',')), collapse = ''),
            file = fn, append = TRUE)
    }
    write('/', file = fn, append = TRUE)
  }
  write.matrix(lhs, file = lhsfn)
}

write_jules_design_factor(paramlist, n = 10)

