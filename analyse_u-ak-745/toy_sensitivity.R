# toy_sensitivity.R


# Goal is to simultaneously estimate true climate sensitivity
# Se, and biases or offsets given data on sensitivities of Paleo, Historic and GCM.
# Assume each of those three has a bias

# Sp = Se + mu_b_p + epsilon_p
# Sh = Se + mu_b_h + epsilon_h
# Sg = Se + mu_b_g + epsilon_g
  
Se = 3

mu_b_p = -1
mu_b_h = -0.5
mu_b_g = 1

epsilon_p = 0.1
epsilon_h = 0.5
epsilon_g = 1

bp = rnorm(mu_b_p, epsilon_p)
bh = rnorm(mu_b_h,epsilon_h)
bg = rnorm(mu_b_g,epsilon_g)

dat = 

