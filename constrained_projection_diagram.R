# constrained_projection_diaram.R
# An R example of constrained projection


n = 25
intercept = rnorm(n, mean = 50, sd = 40)
beta = rnorm(n, mean = 1.2, sd = 0.8)

x = matrix(rep(1:150, n), nrow = n, byrow = TRUE)

y = intercept + (beta*x)

# add in some random noise
epsmat = x
for (i in 1:n){
  epsmat[i, ] = rnorm(ncol(y), mean = 0, sd = 8)
}

z = y + epsmat

yrs = 1901:2050

pdf(file = 'constraint_plot_1.pdf', width = 7, height = 5)
par(las = 1)
matplot(yrs, t(y), type = 'n', axes = FALSE, bty = 'n',
        ylab = 'impact', xlab = 'year')
matlines(yrs, t(z), col = 'black', lty = 'solid', lwd = 1.2)
axis(1)
axis(2)
dev.off()


#points(c(1900, 1900, 1900),obs_1900, col = 'red')
#points(c(1950, 1950, 1950),obs_1950, col = 'red')
#points(c(2000, 2000, 2000),obs_2000, col = 'red')

obsmean = 50 + (1.2*1:150)
#lines(yrs,obsmean, col = 'red')

obs_1900 = c(obsmean[1] - 40, obsmean[1],  obsmean[1] + 40)
obs_1950 = c(obsmean[50] - 50, obsmean[50],  obsmean[50] + 50)
obs_2000 = c(obsmean[100] - 50, obsmean[100],  obsmean[100] + 50)

#segments(x0 = c(1900, 1950, 2000), y0 = c(obs_1900[1],  obs_1950[1], obs_2000[1]), 
#         x1 = c(1900, 1950, 2000), y1 = c(obs_1900[3],  obs_1950[3], obs_2000[3]),
#col = 'darkorange', lwd = 3)

#points(c(1900, 1950, 2000),c(obsmean[1], obsmean[50], obsmean[100]),
#       col = 'darkorange',pch = 19)


ix.1900 = min(obs_1900) < z[,1] &  z[,1] < max(obs_1900)

lcol = rep('grey', n)
lcol[ix.1900] = 'black'

pdf(file = 'constraint_plot_2.pdf', width = 7, height = 5)
par(las = 1)
matplot(yrs, t(y), type = 'n', axes = FALSE, bty = 'n',
        ylab = 'impact', xlab = 'year')
matlines(yrs, t(z), lty = 'solid', lwd = 1.2, col = lcol)
axis(1)
axis(2)

segments(x0 = c(1900), y0 = c(obs_1900[1]), 
         x1 = c(1900), y1 = c(obs_1900[3]),
         col = 'darkorange', lwd = 3)

points(1900,c(obsmean[1]),
       col = 'darkorange',pch = 19)
legend('topleft', legend = 'observation', pch = 19, col = 'darkorange', bty = 'n')
dev.off()


ix.1950 = min(obs_1900) < z[,1] &  z[,1] < max(obs_1900) & min(obs_1950) < z[,50] &  z[,50] < max(obs_1950)
lcol = rep('grey', n)
lcol[ix.1950] = 'black'

pdf(file = 'constraint_plot_3.pdf', width = 7, height = 5)
par(las = 1)
matplot(yrs, t(y), type = 'n', axes = FALSE, bty = 'n',
        ylab = 'impact', xlab = 'year')
matlines(yrs, t(z), lty = 'solid', lwd = 1.2, col = lcol)
axis(1)
axis(2)

segments(x0 = c(1900, 1950), y0 = c(obs_1900[1],  obs_1950[1]), 
         x1 = c(1900, 1950), y1 = c(obs_1900[3],  obs_1950[3]),
         col = 'darkorange', lwd = 3)

points(c(1900, 1950), c(obsmean[1], obsmean[50]),
       col = 'darkorange',pch = 19)
legend('topleft', legend = 'observation', pch = 19, col = 'darkorange', bty = 'n')
dev.off()


ix.2000 = min(obs_1900) < z[,1] &  z[,1] < max(obs_1900) & min(obs_1950) < z[,50] &  z[,50] < max(obs_2000)& min(obs_2000) < z[,100] &  z[,100] < max(obs_2000)
lcol = rep('grey', n)
lcol[ix.2000] = 'black'

pdf(file = 'constraint_plot_4.pdf', width = 7, height = 5)
par(las = 1)
matplot(yrs, t(y), type = 'n', axes = FALSE, bty = 'n',
        ylab = 'impact', xlab = 'year')
matlines(yrs, t(z), lty = 'solid', lwd = 1.2, col = lcol)
axis(1)
axis(2)

segments(x0 = c(1900, 1950, 2000), y0 = c(obs_1900[1],  obs_1950[1], obs_2000[1]), 
         x1 = c(1900, 1950, 2000), y1 = c(obs_1900[3],  obs_1950[3], obs_2000[3]),
         col = 'darkorange', lwd = 3)

points(c(1900, 1950, 2000), c(obsmean[1], obsmean[50], obsmean[100]),
       col = 'darkorange',pch = 19)
legend('topleft', legend = 'observation', pch = 19, col = 'darkorange', bty = 'n')
dev.off()









  
  