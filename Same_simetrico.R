####### Simulación del modelo de Samejima (GRM) ################
####### Y ecuperación de parámetros por Bayes ##################

##########################################################################################
### Modelo modificado de dos parametros


same_repa_acum <- function(theta, alpha, delta, tau, cate){
  # El bloque anterior permite que la función se pueda aplicar a un vector de thetas.
    acum <- array(data = NA, dim = c(length(theta), length(alpha), cate))
    
    for (i in 1:length(theta)) {
      for (j in 1:length(delta)) {
        acum[i, j, 1] <- 1
        acum[i, j, 2] <- exp(alpha[j] * (theta[i] - (delta[j] - tau[i]))) / (1 + exp(alpha[j] * (theta[i] - (delta[j] - tau[i]))))
        acum[i, j, 3] <- exp(alpha[j] * (theta[i] - delta[j])) / (1 + exp(alpha[j] * (theta[i] - delta[j])))
        acum[i, j, 4] <- exp(alpha[j] * (theta[i] - (delta[j] + tau[i]))) / (1 + exp(alpha[j] * (theta[i] - (delta[j] + tau[i]))))
      }
    }
    
    prob <- list(cate1 = acum[,,1] - acum[,,2], 
                 cate2 = acum[,,2] - acum[,,3],
                 cate3 = acum[,,3] - acum[,,4],
                 cate4 = acum[,,4])
    
    return(prob)
}

#######################################################################################################
# Parámetros a recobrar.

# Deltas's item
set.seed(1)
delta <- rnorm(15, 0, 1)

print(delta)

set.seed(14)
thetas <- rnorm(1000, 0, 1)
hist(thetas)
summary(thetas)

set.seed(13)
alphas <- rgamma(15, 2, 2)
print(alphas)

set.seed(12)
taus <- rlnorm(1000, 0, 1)
print(taus)

cate = 4
#######################################################################################################
# Simulación de datos

acum = same_repa_acum(thetas, alphas, delta, taus, cate)


simulated_data_GRM <- function(thetas, alphas, delta, tau, num_cate){
  
  prob_cate <- same_repa_acum(theta = thetas, alpha = alphas, delta = delta, tau = tau, cate = num_cate)
  data = matrix(data = NA, nrow = 1000, ncol = 15)
  set.seed(10)
  
    for (item in 1:15) {
    for (person in 1:1000) {
      resp_item = sample(c(0, 1, 2, 3), size = 1, replace = FALSE, prob = c(prob_cate$cate1[person, item], 
                                                                            prob_cate$cate2[person, item], 
                                                                            prob_cate$cate3[person, item], 
                                                                            prob_cate$cate4[person, item]))   
      data[person, item] = resp_item
    }   
    
  }
  return(data)
}

generated_data <- simulated_data_GRM(thetas, alphas, delta, taus, cate)

###################################################################################################################
###################################################################################################################
# Bayesian Modelling.
n_per <- dim(generated_data)[1]
n_items <- dim(generated_data)[2]
n_cate <- 4 # Número de categorías es igual m-1.

generated_data = generated_data + 1
View(generated_data)

data_model = list(
  J = n_cate,
  n_per = n_per,
  n_items = n_items,
  DATA = generated_data
)

write("data{
  int<lower=1, upper=4> J; //número de categorías m-1
  int<lower = 0> n_per; //número de personas
  int<lower = 0> n_items; //número de ítems
  int<lower = 1, upper = J> DATA[n_per, n_items];
}
parameters{
    vector[n_per] theta;
    real<lower=0> alpha [n_items];
    real delta[n_items]; //item position
    //ordered[J-1] tau[n_items]; //psichological distance between categories
    vector<lower = 0>[n_per] tau;
}
model{
    alpha ~ lognormal(0, 2);
    theta ~ normal(0,1);
    delta~normal(0,5);
    tau ~ lognormal(0, 2);


    for (i in 1:n_per){
      for (j in 1:n_items){

        vector[3] beta;
        beta[1] = delta[j] - tau[i];
        beta[2] = delta[j];
        beta[3] = delta[j] + tau[i];

        DATA[i,j] ~ ordered_logistic(theta[i]*alpha[j], alpha[j]*beta[]);
    }}
}",

"model_samejima_simetrico.stan")

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit <- stan(file = "model_samejima_simetrico.stan", data = data_model, model_name = "Samejima_Model_Simetrico_Fit", 
            iter = 15000, warmup = 2000, chains = 4, thin = 2, verbose = T, cores = 6)

suma <- summary(fit) # Descripción de los datos 
desc <- suma$summary # Visualizar estadísticos importantes.

summary(desc[, 10]) # Rhat
summary(desc[,9]) # n.efe.

theta_mean = desc[1:1000, 1]
alpha_mean = desc[1001:1015, 1]
delta_mean = desc[1016:1030, 1]
tau_mean = desc[1031:2030, 1]
tau_median = desc[1031:2030, 6]


#############################################################################################################################
# Visualización
############################


person_plot = c()

set.seed(1408)
person = sample(x = 1:1000, size = 20, replace = T)
person = sort(person)

person_tau = sample(x = 1:1000, size = 20, replace = T)
person_tau = sort(person_tau)

for (i in 1:20) {
  stri = paste("theta[", round(person[i]), "]", sep = "")
  person_plot[i] = stri
}

for (i in 1:20) {
  stri = paste("tau[", round(person[i]), "]", sep = "")
  person_tau[i] = stri
}

plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = person_plot, point_est = "mean") # Densities thetas
plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = "alpha", point_est = "mean") # Densities alphas
plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = "delta", point_est = "mean") # Densities betas
plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = person_tau, point_est = "median") # Densities taus

########################################################################################################################
# Comparación de Resultados

###################################
# Theta de las personas
###########

col_theta = "#D78D50"
plot(NULL, ann = F, axes = F, xlim = c(-3, 3), ylim = c(-3, 3))
axis(1, at = c(-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4))
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(thetas, theta_mean, pch = 21, 
       bg = paste(col_theta,'aa',sep=''), cex=1.5, scatter = T)
mtext('\u03B8 Real', 1,cex=1.5,line=3)
mtext('\u03B8 Estimada',2,cex=1.5,line=2.5)
mtext('\u03B8 Real vs \u03B8 Estimada',3,cex=2, line = 1.15)
text(x = 1.5, y = -0.5, paste("\u03C1: ", round(cor(thetas, theta_mean), 3), sep = ""), cex = 2)


###################################
# Alphas de los ítems
###########


col_alpha = "#BD2629"
plot(NULL, ann = F, axes = F, xlim = c(0, 3), ylim = c(0, 3))
axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5))
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "#85CF00", lw = 2)
points(alphas, alpha_mean, pch = 21, 
       bg = paste(col_alpha,'aa',sep=''), cex=1.5, scatter = T)
mtext('\u03B1 Real', 1,cex=1.5,line=3)
mtext('\u03B1 Estimada',2,cex=1.5,line=2.5)
mtext('\u03B1 Real vs \u03B1 Estimada',3,cex=2, line = 1.15)
text(x = 1.5, y = 1, paste("\u03C1: ", round(cor(alphas, alpha_mean), 3), sep = ""), cex = 2)

###################################
# Deltas de los ítems
###########


col_delta = "#AC6A9F"
plot(NULL, ann = F, axes = F, xlim = c(-3, 2), ylim = c(-3, 2))
axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5))
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "#85CF00", lw = 2)
points(delta, delta_mean, pch = 21, 
       bg = paste(col_delta,'aa',sep=''), cex=1.5, scatter = T)
mtext('\u03B4 Real', 1,cex=1.5,line=3)
mtext('\u03B4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03B4 Real vs \u03B4 Estimada',3,cex=2, line = 1.15)
text(x = 1.5, y = 1, paste("\u03C1: ", round(cor(delta, delta_mean), 3), sep = ""), cex = 2)

###################################
# Taus de las personas
###########

col_tau = "#4F98C4"
plot(NULL, ann = F, axes = F, xlim = c(0, 10), ylim = c(0, 10))
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(taus, tau_mean, pch = 21, 
       bg = paste(col_tau,'aa',sep=''), cex=1.25, scatter = T)
mtext('\u03C4 Real', 1,cex=1.5,line=3)
mtext('\u03C4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03C4 Real vs \u03C4 Estimada (Valor Esperado)',3,cex=2, line = 1.15)

plot(NULL, ann = F, axes = F, xlim = c(10, 40), ylim = c(10, 40))
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(taus, tau_mean, pch = 21, 
       bg = paste(col_tau,'aa',sep=''), cex=1.25, scatter = T)
mtext('\u03C4 Real', 1,cex=1.5,line=3)
mtext('\u03C4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03C4 Real vs \u03C4 Estimada (Valor Esperado)',3,cex=2, line = 1.15)

############## Figura Completa

plot(NULL, ann = F, axes = F, xlim = c(0, 40), ylim = c(0, 40))
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(taus, tau_mean, pch = 21, 
       bg = paste(col_tau,'aa',sep=''), cex=1.25, scatter = T)
mtext('\u03C4 Real', 1,cex=1.5,line=3)
mtext('\u03C4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03C4 Real vs \u03C4 Estimada (Valor Esperado C.)',3,cex=2, line = 1.15)

#################################################### Point Estimate: Median

plot(NULL, ann = F, axes = F, xlim = c(0, 2), ylim = c(0, 2))
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(taus, tau_median, pch = 21, 
       bg = paste(col_tau,'aa',sep=''), cex=1.25, scatter = T)
mtext('\u03C4 Real', 1,cex=1.5,line=3)
mtext('\u03C4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03C4 Real vs \u03C4 Estimada (Mediana)',3,cex=2, line = 1.15)


text(x = 1.5, y = -0.5, paste("\u03C1: ", round(cor(thetas, theta_mean), 3), sep = ""), cex = 2)

col_tau = "#4F98C4"
plot(NULL, ann = F, axes = F, xlim = c(2, 25), ylim = c(2, 25))
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(taus, tau_median, pch = 21, 
       bg = paste(col_tau,'aa',sep=''), cex=1.25, scatter = T)
mtext('\u03C4 Real', 1,cex=1.5,line=3)
mtext('\u03C4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03C4 Real vs \u03C4 Estimada (Mediana)',3,cex=2, line = 1.15)

###### Fig Completa

plot(NULL, ann = F, axes = F, xlim = c(0, 25), ylim = c(0, 25))
axis(1)
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(taus, tau_median, pch = 21, 
       bg = paste(col_tau,'aa',sep=''), cex=1.25, scatter = T)
mtext('\u03C4 Real', 1,cex=1.5,line=3)
mtext('\u03C4 Estimada',2,cex=1.5,line=2.5)
mtext('\u03C4 Real vs \u03C4 Estimada (Mediana C.)',3,cex=2, line = 1.15)


summary(taus)
