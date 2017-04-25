# summerNWtemperature.R
# Doug McNeall

# -------------------------------------------------------------------
# 0. Administration
# -------------------------------------------------------------------
library(devtools)
library(DiceKriging)
install_github(repo = "dougmcneall/hde")
library(hde)
library(fields)
library(sensitivity)

# some useful emulation and visualisation tools
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

# -------------------------------------------------------------------
# 1. Data loading
# -------------------------------------------------------------------

X = read.table(file = 'parameter_sets_input_170.txt', header = TRUE)
y = read.table(file = 'output_modeled_JJA_NW_temp_170.txt')
colnames(y) = 'SNWT'

X.norm = normalize(X)

# A pairs plot of the input space quickly shows that the entcoef
# parameter has a major impact on the probability of return
pdf(file = 'Xpairs.pdf', width = 7, height = 7)
pairs(X, pch = 19, cex = 0.3, gap = 0, lower.panel = NULL)
dev.off()

# Adding the output onto the end of the pairs plot indicates that
# entcoef, vf1 and v_crif_alpha will probably be the most important
# inputs on temperature
Xy = cbind(X,y)

pdf(file = 'Xypairs.pdf', width = 7, height = 7)
pairs(Xy, pch = 19, cex = 0.3, gap = 0, lower.panel = NULL)
dev.off()

# -------------------------------------------------------------------
# 2. Fit emulator
# -------------------------------------------------------------------

# linear response prior 
fit = km(formula=~., design=X.norm, response=y)

# flat response prior
#fit = km(formula=~1, design=X.norm, response=y)

# Emulator seems to work well for leave-one-out
pdf(file = 'fit.pdf', width = 3, height = 6)
plot(fit)
dev.off()

# -------------------------------------------------------------------
# 3. One at a time sensitivity analysis
# -------------------------------------------------------------------
n = 21
X.oaat= oaat.design(X.norm, n)
y.oaat = predict(fit, newdata = X.oaat, type = 'UK')

pdf(file = 'oaat.pdf', width = 12, height = 6)
par(mfrow = c(3,6), mar = c(4,3,1,1))
ylim = range(y.oaat$mean)

for(i in 1:17){
  ix = seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oaat[ix,i], y.oaat$mean[ix],
       type = 'o',
       xlab = colnames(X.norm)[i], ylab= '', ylim = ylim)
}
dev.off()

# -------------------------------------------------------------------
# 3. FAST sensitivity analysis
# -------------------------------------------------------------------

# generate the design to run the emulator at, using fast99
x = fast99(model = NULL, factors = colnames(X.norm), n = 1000,
            q = "qunif", q.arg = list(min = 0, max = 1))

fast.pred = predict(fit, newdata = x$X, type = 'UK')
fast <- tell(x, fast.pred$mean)
pdf(file = 'fast.pdf', width = 12, height = 6)
par(las = 2, mar = c(10,4,2,1))
plot(fast)
dev.off()



