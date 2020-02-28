# es_ppe_ii_history_match_waves.R
# A more formal history match for the ESPPEii ensemble.
# Use A History matching to create candidates for a new wave
# of model runs.

source('es_ppe_ii_history_match_common.R')

# -----------------------------------------------------------------------------------------------
# Use History matching to create candidates for a new design
# -----------------------------------------------------------------------------------------------

# Which emulator to build with? All, level0 constraint or level 1 constraint?

# What algorithm to use for keeping new candidates?
# 1) Keep the n candidates with the smallest implausibility
# 2) Generate candidate points near the NROY runs, and keep a proportion of them (expand
# towards the edges).
# 3) Just go ahead and generate random points and keep n NROY ones.


# Initially, copy functions from waves/hmwaves.R

create.kmfit.list = function(X, Y){
  # create a list of km fits for multivariate output
  Y.list = as.list(as.data.frame(Y))
  
  fit.list = NULL
  for(i in 1:length(Y.list)){
    
    fit = km(formula = ~.,  design = X, response = Y.list[[i]])
    fit.list = c(fit.list, fit)
  }
  fit.list
}

# Are all the elements of a matrix row below a threshold?
all.bt = function(x, thres) all(x < thres)

ChooseMaximinNroy = function(n.app, waveobj, nreps){
  # Choose a set of NROY points with the largest minimum
  # distance
  ix.list = vector(mode='list', length = nreps)
  mindist.vec = rep(NA, nreps)
  
  for(i in 1:nreps){
    ix = sample(1:nrow(waveobj$X.nroy), n.app)
    X.cand = waveobj$X.nroy[ix, ]
    ix.list[[i]] = ix
    mindist = min(dist( X.cand))
    mindist.vec[i] = mindist
  }
  ix.chosen = ix.list[[which.max(mindist.vec)]]
  
  return(waveobj$X.nroy[ix.chosen, ])
}

# Improvements to the above function from Andrianakis et al. (2015)
# 1) Reduce the range the emulator is fit over, as the waves continue.
# 2) Sample candidate design points from NEAR the existing NROY design points.

add.nroy.design.points = function(X, Y, Y.target, n.aug,
                                  mins.aug, maxes.aug,
                                  thres = 3, disc.list,
                                  disc.sd.list, obs.sd.list){
  
  # Add NROY design points to a design, using uniform sampling from the
  # entire input space.
  # Inputs
  # X            ...       design matrix (output from maximinLHS)    
  # Y            ...       model output matrix
  # Y.target     ...       Target output, or "observation"
  # n.aug        ...       number of candidate points to augment the lhs
  # mins.aug
  # maxes.aug    ...       Limits on the candidate points
  
  # list of fitted km objects, one list element for each output
  fit.list = create.kmfit.list(X=X, Y=Y)
  
  # How good is the fit?
  loo.list = lapply(fit.list, FUN = leaveOneOut.km, type = 'UK', trend.reestim = TRUE)
  
  loo.mse.vec = rep(NA, length(loo.list))
  for(i in 1:length(loo.list)){
    
    loo.mse = mean((loo.list[[i]]$mean - Y.target[,i])^2)
    loo.mse.vec[i] = loo.mse
  }
  
  # create a new set of candidate points
  #X.aug = augmentLHS(X, n.aug) # this will add to the current design.
  X.aug = samp.unif(n.aug, mins = mins.aug, maxes = maxes.aug)
  colnames(X.aug) = colnames(X)
  
  # predicting the output at each design point
  pred.list = lapply(fit.list, FUN = 'predict', newdata = X.aug, type = 'UK')
  # calculate the implausibility at each design point
  
  impl.mat = NULL
  for(i in 1:length(pred.list)){
    
    pred.impl = impl(em = pred.list[[i]]$mean,
                     em.sd = pred.list[[i]]$sd,
                     disc = disc.list[[i]],
                     disc.sd = disc.sd.list[[i]],
                     obs = Y.target[[i]],
                     obs.sd = obs.sd.list[[i]])
    
    impl.mat = cbind(impl.mat, pred.impl)
  }
  
  # Which of the candidte design points are NROY?
  
  # create a matrix of the implausibility measures
  # find the indices of the matrix where all are below the threshold.
  nroy.tf = apply(impl.mat, 1, FUN = all.bt, thres = thres)
  nroy.ix = which(nroy.tf==TRUE)
  X.nroy = X.aug[nroy.ix, ]
  
  X.nroy.max = apply(X.nroy, 2, max)
  X.nroy.min = apply(X.nroy, 2, min)
  
  return(list(X.nroy = X.nroy,
              X.aug = X.aug, 
              impl.mat = impl.mat, 
              loo.mse.vec = loo.mse.vec,
              fit.list = fit.list,
              X.nroy.max = X.nroy.max,
              X.nroy.min = X.nroy.min)
  )
}

WithinRange = function(x, maxes, mins){
  # which elements of a vector are between
  # elements of the min and max vectors?
  
  all(x < maxes && x > mins)
}





allix = 1:(nrow(dat.norm))

# --------------------------------------------------------------------------------
# Level 0 constraint
# Model runs and has some kind of carbon cycle
# --------------------------------------------------------------------------------

# Constrain with global runoff and nbp, so that the emulator
# is not negatively affected by really bad points
level0.ix = which(dat.norm[,'runoff'] >0.5 & dat.norm[,'nbp'] > -10)

dat.level0  = dat.norm[level0.ix, ]     # (scaled to sensible) data that passes level 0
X.level0 = X[level0.ix, ]               # Normalised inputs that pass level 0
lhs.level0 = lhs[level0.ix, ]           # Non-normalised inputs that pass level 0

nlevel0.ix = setdiff(allix, level0.ix)  # These inputs DO NOT PASS
X.nlevel0 = X[nlevel0.ix, ]
lhs.nlevel0 = lhs[nlevel0.ix, ]

# These are the ensemble members which do not run (produce NA)
na.ix = which(is.na(dat.norm[,'runoff']))
X.na = X[na.ix, ]
lhs.na = lhs[na.ix, ]

# --------------------------------------------------------------------------------
# Level 1 constraint
# Minimally functioning carbon cycle
# --------------------------------------------------------------------------------

# Apply level 1 constraints
level1.ix = which(dat.norm[,'cs_gb'] > 750 & dat.norm[,'cs_gb'] < 3000 &
  dat.norm[,'cv'] > 300 & dat.norm[,'cv'] < 800 & 
  dat.norm[,'npp_n_gb'] > 35 &
  dat.norm[,'npp_n_gb'] < 80 &
  dat.norm[,'runoff'] >0.5 &
  dat.norm[,'nbp'] > -10
  )

X.level1 = X[level1.ix, ]

# Ranges of the inputs after a level 1 constraint
rn.l1 = apply(X.level1,2, range)

nlevel1.ix = setdiff(allix, level1.ix)
X.nlevel1 = X[nlevel1.ix, ]

dat.level1 = dat.norm[level1.ix, ]

# In order to use the history matching code, we'll need to specify targets -
# specifically "observations" with uncertainty thatapproximately match our
# "hard boundaries" from the existing constraints.

# Initial guess just choose the centre of the (implied) uniform distribution.
#cs_gb.target = (3000 - 750) / 2
#cv.target = (800 - 300) / 2
#npp_n_gb.target = (80 - 35) / 2
#runoff.target = 1
#nbp.target = 0
#gpp.target = 75

#Y.target = c(cs_gb.target, cv.target, npp_n_gb.target, runoff.target, nbp.target)

# I'm going to set it so that +3sd aligns approximately with the original limits
# given by Andy Wiltshire.
#cs_gb       cv    gpp_gb        nbp npp_n_gb    runoff
Y.lower = c(750, 300, 50, -10, 35, 0.5)
Y.upper = c(3000, 800,200, 10, 80, 1.5)
Y.target = (Y.upper - abs(Y.lower)) / 2 # abs() to fix the problem with negative numbers


# standard deviation is derived from the limits and the central target
# (this distance is assumed to be 3 standard deviations.
Y.sd = (Y.upper - Y.target) / 3
names(Y.sd) = colnames(dat.norm)

obs.sd.list = as.list(rep(0.01,p))
disc.list =  as.list(rep(0,p)) 
disc.sd.list =  as.list(Y.sd/2)
thres = 3

mins.aug = apply(X.level1, 2, FUN = min)
maxes.aug =apply(X.level1, 2, FUN = max)

# convert Y.target for ingestion into function
Y.target = matrix(Y.target, nrow = 1)


# Final output needs to be expressed in terms of original LHS, then put back out to conf files.

# This function adds n.aug potential design points, and finds their implausibility
# score in X.nroy
wave1 = add.nroy.design.points(X = X.level1, 
                               Y = dat.level1, 
                               Y.target = Y.target,
                               n.aug = 10000, 
                               mins.aug = mins.aug,
                               maxes.aug = maxes.aug,
                               thres = 3,
                               disc.list=disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

# Tighter restrictions on model outputs will lead to a smaller
# X.nroy object, as more potential inputs are ruled out.
dim(wave1$X.nroy)

# Exclude any design points outside of the ranges returned as NROY
keep1.ix = which(apply(X.level1, FUN = WithinRange,1,
                      maxes = wave1$X.nroy.max,
                      mins = wave1$X.nroy.min))

X.keep1 = X.level1[keep1.ix, ]

# choose the number of NROY design points that we require for the next wave
# of model runs.
n.app = 100

# The ChooseMaximinNroy algorithm picks a matrix which has space filling properties
X.aug1 = ChooseMaximinNroy(n.app = n.app, waveobj=wave1, nreps = 10000)


X2 = rbind(X.keep1, X.aug1)

# Express this matrix in terms of the original latin hypercube

# Visualise the various design matrices.
dev.new(width = 10, height = 10)
pairs(X.keep1, gap = 0, xlim = c(0,1), ylim = c(0,1)
)

# Here is the new part of the design.
# This has pleasing gaps with no design points, meaning we've ruled out parts of space.
dev.new(width = 10, height = 10)
pairs(X.aug1, gap = 0, xlim = c(0,1), ylim = c(0,1)
)


# end of analysis
stop()

# Here is code from write_jules_design_3.R
library(lhs)
library(MASS)
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")

# ----------------------------------------------------
# Create a list of parameters
# ----------------------------------------------------
# NOTE potential problem with input lhs

write_jules_design = function(X.aug, paramlist, fac, minfac, maxfac, tf,
  fnprefix = 'param-perturb-test',
                              lhsfn = 'lhs.txt',stanfn = 'stanparms.txt', allstanfn = 'allstanparms.txt', rn = 5, startnum=0){
  # This code writes a design taking either a 'factor', min and max by which
  # to multiply all pfts, or perturbing each pft individually according to
  # their maximum and minimum in the parameter list.
  # fac is a character vector of names of variables that you would like to alter
  # by a factor. Everything else gets variaed by PFT
  # minfac and maxfac must correspond to fac - i.e. one value per parameter, in the
  # correct order.
  # tf is a character vector containing the logical parameters
  
  stopifnot(
    all(
      length(fac) == length(minfac),
      length(fac) == length(maxfac)
    )
  )
  
  paramvec = names(paramlist)
  nmlvec = unlist(lapply(paramlist, FUN = function(x) x$namelist))
  
  # which parameters do we want as a parameter list?
  fac.ix = which(names(paramlist) %in% fac)
  tf.ix = which(names(paramlist) %in% tf)
  
  paramfac = paramlist[fac.ix]
  paramtf = paramlist[tf.ix]
  parampft = paramlist[-c(fac.ix, tf.ix)]
  pftvec = names(parampft)
  
  parampft_nml = unlist(lapply(parampft, FUN = function(x) x$namelist))
  paramfac_nml = unlist(lapply(paramfac, FUN = function(x) x$namelist))
  paramtf_nml = unlist(lapply(paramtf, FUN = function(x) x$namelist))
  
  parampft_standard = unlist(lapply(parampft, FUN = function(x) x$standard))
  parampft_mins = unlist(lapply(parampft, FUN = function(x) x$min))
  parampft_maxes = unlist(lapply(parampft, FUN = function(x) x$max))

  paramfac_standard = unlist(lapply(paramfac, FUN = function(x) x$standard))
  paramtf_standard = unlist(lapply(paramtf, FUN = function(x) x$standard))
  
  all_names = c(names(parampft_standard), fac, tf)
  k = length(all_names)

  # This writes out the 'standard' parameters, but sets the "factors"
  # used to multiply them to "1". It should produce the same
  # number of columns as the output design.
  #
  standard_matrix = matrix( c(parampft_standard, rep(1, length(fac)), paramtf_standard), nrow = 1)
  colnames(standard_matrix) = all_names
  write.matrix(standard_matrix, file = stanfn)

  standard_matrix_all = matrix( c(parampft_standard, paramfac_standard, paramtf_standard), nrow = 1)

  # Writes standard values for all pfts, even if they aren't used.
  colnames(standard_matrix_all) = c(names(parampft_standard),
            names(paramfac_standard),
            names(paramtf_standard))
  
  write.matrix(standard_matrix_all, file = allstanfn)



  # ------------------------------------------------------------------------
  # I think I just need to replace this section, putting in the
  # new design from the HM wave rather than generating a new maximin lhs
  
  # Do the pfts and then the factors and then the logical
#  lhs = unnormalize(
#    maximinLHS(n = n, k = k, dup = 1),
#    un.mins = c(parampft_mins , minfac, rep(0, length(tf)), recursive = TRUE),
#    un.maxes = c(parampft_maxes, maxfac, rep(1, length(tf)), recursive = TRUE)
#  )
#  colnames(lhs) = all_names

  #
  lhs = unnormalize(
    X.aug,
    un.mins = c(parampft_mins , minfac, rep(0, length(tf)), recursive = TRUE),
    un.maxes = c(parampft_maxes, maxfac, rep(1, length(tf)), recursive = TRUE)
    )
  colnames(lhs) = all_names
    
)

  # -------------------------------------------------------------------------
  
  for(i in 1:nrow(lhs)){
    fn = paste0(fnprefix,sprintf("%04d",(startnum+i)-1),'.conf')
    
    for(el in unique(nmlvec)){
      write(paste0('[namelist:',el,']'), file = fn, append = TRUE)
      
      # grab the parts of the list that match
      pft_elms = parampft[parampft_nml==el] # & statement
      pft_elms_vec = names(pft_elms)
      fac_elms = paramfac[paramfac_nml==el]
      fac_elms_vec = names(fac_elms)
      tf_elms = paramtf[paramtf_nml==el]
      tf_elms_vec = names(tf_elms)
      
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
      
      if(length(tf_elms_vec) > 0){
        for(l in 1: length(tf_elms_vec)){
          param = tf_elms_vec[l]
          colix = grep(param, colnames(lhs))
          lhs.factor = lhs[i, colix]
          logical.out = NA
          if (lhs.factor>=0.5){logical.out = '.true.'}
          else {logical.out = '.false.'}
          write(paste0(param,'=',logical.out, collapse = ''), file = fn, append = TRUE)
        }
      }
      write('\n', file = fn, append = TRUE)
    }
  }
  
  write.matrix(lhs, file = lhsfn)
}


# Here is the instance of the previous being called.

#source('write_jules_design3.R')
source('~/brazilCSSP/code/brazil_cssp/jules_params_u-ao732.R')


# Easiest way to generate a design of the right size is to have a "fac" which takes
# the names from the parameter list, and then multiplies everything by 0.5 or 2

tf = 'l_vg_soil'
# we don't want anything that is TRUE/FALSE to be in fac
fac.init = names(paramlist)
not_tf.ix = which(names(paramlist)!=tf)
paramlist.trunc = paramlist[not_tf.ix]

fac = names(paramlist.trunc)

# these just find max and min for the design and transleate it into multiplication factors.
maxfac = lapply(paramlist.trunc,function(x) x$max[which.max(x$max)] / x$standard[which.max(x$max)])
minfac = lapply(paramlist.trunc,function(x) x$min[which.max(x$max)] / x$standard[which.max(x$max)])

write_jules_design(X.aug = X2, paramlist,
                    fac = fac, minfac = minfac, maxfac = maxfac,
                    tf = tf,
                    fnprefix = '~/brazilCSSP/code/brazil_cssp/conf_files_u-ao732b/param-perturb-P',
                    lhsfn = '~/brazilCSSP/code/brazil_cssp/conf_files_u-ao732b/lhs_u-ao732a.txt',
                    stanfn = '~/brazilCSSP/code/brazil_cssp/conf_files_u-ao732b/stanparms_u-ao732a.txt',
                    allstanfn = '~/brazilCSSP/code/brazil_cssp/conf_files_u-ao732b/allstanparms_u-ao732a.txt',
                    rn = 12)

# everything is done
print(cbind(minfac, maxfac))


# Checking section
lapply(paramlist, function(x) length(x$standard))
lapply(paramlist, function(x) length(x$standard)==length(x$min) & length(x$standard)==length(x$max))












X2 = rbind(X.level1[keep1.ix, ] , ChooseMaximinNroy(n.app = n.app, waveobj=wave1, nreps = 10000))


Y2.raw = run.model(X2)
Y2 = normalize(Y2.raw, wrt = Y.raw)

wave2 = add.nroy.design.points(X = X2,
                               Y = Y2, 
                               Y.target = Y.target,
                               n.aug = 30000,
                               mins.aug = wave1$X.nroy.min,
                               maxes.aug = wave1$X.nroy.max,
                               thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

keep2.ix = which(apply(X2, FUN = WithinRange,1,
                       maxes = wave2$X.nroy.max,
                       mins = wave2$X.nroy.min))

X3 = rbind(X2[keep2.ix, ], ChooseMaximinNroy(n.app = n.app, waveobj=wave2, nreps = 10000))
