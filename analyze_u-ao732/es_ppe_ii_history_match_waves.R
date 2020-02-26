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




wave1 = add.nroy.design.points(X = X, 
                               Y = Y, 
                               Y.target = Y.target,
                               n.aug = 50000, 
                               mins.aug = mins.aug,
                               maxes.aug = maxes.aug,
                               thres = 3,
                               disc.list=disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

# Exclude any design points outside of the ranges returned as NROY



WithinRange = function(x, maxes, mins){
  # which elements of a vector are between
  # elements of the min and max vectors?
  
  all(x < maxes && x > mins)
}

keep1.ix = which(apply(X, FUN = WithinRange,1,
                      maxes = wave1$X.nroy.max,
                      mins = wave1$X.nroy.min))

X2 = rbind(X[keep1.ix, ] , ChooseMaximinNroy(n.app = n.app, waveobj=wave1, nreps = 10000))
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
