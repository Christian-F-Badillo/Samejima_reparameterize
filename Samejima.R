####### Simulación del modelo de Samejima (GRM) ################
####### Y ecuperación de parámetros por Bayes ##################

##########################################################################################
### Modelo de dos parámetros
Prob_2PLM <- function(theta, alpha, beta){
    # Revisamos que los parámetros tengan valores válidos
    
    if (length(theta) > 1){
        result <- sapply(theta, FUN = Prob_2PLM, alpha = alpha, beta = beta)
        return(result)
    }
    # El bloque anterior permite que la función se pueda aplicar a un vector de thetas.
    
    tmp <- exp(alpha * (theta - beta))
    return(tmp / (1 + tmp))    
}

############################################################################################
# Modelo de Samejima (GRM)

ProbCum_GRM <- function(theta, alpha, beta, cat = NA){
    if (any(sort(beta) != beta))
        stop("Los elementos en el vector beta deben estar ordenados")
    # Revisamos que los parámetros tengan valores válidos
    
    if (length(cat) == 1)
        if (is.na(cat))
            cat <- seq(0, length(beta))
    # Si cat = NA, entonces cat será un vector (0, 1, ..., m - 1); es decir, la función entregará la probabilidad de todas
    #      las categorías.
    
    if (length(theta) > 1){
        result <- sapply(theta, FUN = ProbCum_GRM, alpha = alpha, beta = beta, cat = cat)
        if (length(cat) > 1)
            return(t(result))
        return(result)
    }
    # El bloque anterior permite que la función se pueda aplicar a un vector de thetas.
    
    result <- c()
    for (j in 1:length(cat)){
        if (cat[j] == 0) {
            result <- c(result, 1)
        }
        else {
            result <- c(result, Prob_2PLM(theta = theta, alpha = alpha, beta = beta[cat[j]]))
        }
    }
    return(result)  
}

#######################################################################################################
# Calcular las probabilidades de responder en cada categoría

Prob_GRM <- function(theta, alpha, beta, cat = NA){
    if (any(sort(beta) != beta))
        stop("Los elementos en el vector beta deben estar ordenados")
    # Revisamos que los parámetros tengan valores válidos
    
    if (length(cat) == 1)
        if (is.na(cat))
            cat <- seq(0, length(beta))
    # Si cat = NA, entonces cat será un vector (0, 1, ..., m - 1); es decir, la función entregará la probabilidad de todas
    #      las categorías.
    
    if (length(theta) > 1){
        result <- sapply(theta, FUN = Prob_GRM, alpha = alpha, beta = beta, cat = cat)
        if (length(cat) > 1)
            return(t(result))
        return(result)
    }
    ProbCums <- ProbCum_GRM(theta = theta, alpha = alpha, beta = beta)
    # Calculamos primero las probabilidades acumuladas
    result <- c()
    for (j in 1:length(beta)){
        result <- c(result, ProbCums[j] - ProbCums[j + 1])
        # Calculamos para cada categoría la diferencia entre la probabilidad acumulada de esta categoría y de la
        #      siguiente categoría.
    }
    result <- c(result, ProbCums[length(ProbCums)])
    # La probabilidad de responder en la categoría más alta es igual a la probabilidad acumulada de esta categoría.
    return(result[cat + 1])
}

#######################################################################################################
# Parámetros a recobrar.
item_betas <- matrix(data = NA, nrow = 15, ncol = 3)

for (i in 1:15) {
    set.seed(i)
    betas <- rnorm(3, 0, 2)
    betas <- sort(betas)
    
    for (j in 1:3) {
        item_betas[i, j] = betas[j]
    }
}

print(item_betas)

set.seed(1)
thetas <- rnorm(1000, 0, 1)
hist(thetas)
summary(thetas)

set.seed(1)
alphas <- rgamma(15, 2, 2)
print(alphas)
#######################################################################################################
# Simulación de datos

simulated_data_GRM <- function(thetas, alphas, item_betas){
    data = matrix(data = NA, nrow = 1000, ncol = 15)
    set.seed(10)
    for (item in 1:15) {
        prob_cate <- Prob_GRM(theta = thetas, alpha = alphas[item], beta = as.vector(item_betas[item, ]), cat = NA)
        for (person in 1:1000) {
            prob = as.vector(prob_cate[person, ])
            resp_item = sample(c(0, 1, 2, 3), size = 1, replace = FALSE, prob = prob)   
            data[person, item] = resp_item
        }   
        
    }
    return(data)
}

generated_data <- simulated_data_GRM(thetas, alphas, item_betas)

write.csv(item_betas, file = "betas.csv")
write.csv(alphas, file = "alphas.csv")
write.csv(thetas, file = "thetas.csv")
write.csv(generated_data, file = "sim_data.csv")


#########################################################################################################
# Bayes Modelling.

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

write("
data{
  int<lower=1, upper=5> J; //número de categorías m-1
  int<lower = 0> n_per; //número de personas
  int<lower = 0> n_items; //número de ítems
  int<lower = 1, upper = J> DATA[n_per, n_items];
}
parameters{
    vector[n_per] theta;
    real<lower=0> alpha [n_items];
    ordered[J-1] beta[n_items]; //category difficulty
}
model{
    alpha ~ gamma(1.5,2);
    theta ~ normal(0,1);
    for (i in 1:n_items){
      for (j in 1:(J-1)){
        beta[i,j] ~ normal(0,5);
      }}
    for (i in 1:n_per){
      for (j in 1:n_items){
        DATA[i,j] ~ ordered_logistic(theta[i]*alpha[j],beta[j]);
      }}
}" 
,"model_samejima.stan")

library(rstan)

fit <- stan(file = "model_samejima.stan", data = data_model, model_name = "Samejima_Model_Fit", 
     iter = 8000, warmup = 1000, chains = 4, thin = 2, verbose = T, cores = 6)

suma <- summary(fit) # Descripción de los datos 
desc <- suma$summary # Visualizar estadísticos importantes.

summary(desc[1:1060, 10]) # Rhat
summary(desc[1:1060, 9]) # n_effect

theta_means = desc[1:1000, 1] # Expected value Betas.
alpha_means = desc[1001:1015, 1] # Expected value Alphas.
beta_means = matrix(data = desc[1016:1060], nrow = 15, ncol = 3, byrow = T) # Expected value Betas.


#############################################################################################################################
# Visualización
############################


person_plot = c()

set.seed(25)
person = sample(x = 1:1000, size = 20, replace = T)
person = sort(person)

for (i in 1:20) {
  stri = paste("theta[", round(person[i]), "]", sep = "")
  person_plot[i] = stri
}

plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = "alpha", point_est = "mean") # Densities alphas
plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = "beta", point_est = "mean") # Densities betas
plot(fit, show_density = T, ci_level = 0.95, fill_color = "#AFDEA2", pars = person_plot, point_est = "mean") # Densities thetas

########################################################################################################################
# Comparación de Resultados
###################################
# Theta de las personas
###########

col_theta = "#D78D50"
plot(NULL, ann = F, axes = F, xlim = c(-3.5, 4), ylim = c(-4, 4))
axis(1, at = c(-3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4))
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "red", lw = 2)
points(thetas, theta_means, pch = 21, 
       bg = paste(col_theta,'aa',sep=''), cex=1.5, scatter = T)
mtext('\u03B8 Real', 1,cex=1.5,line=3)
mtext('\u03B8 Estimada',2,cex=1.5,line=2.5)
mtext('\u03B8 Real vs \u03B8 Estimada',3,cex=2, line = 1.15)
text(x = 1.5, y = -0.5, paste("\u03C1: ", round(cor(thetas, theta_means), 3), sep = ""), cex = 2)


###################################
# Alphas de los ítems
###########


col_alpha = "#BD2629"
plot(NULL, ann = F, axes = F, xlim = c(0, 2.5), ylim = c(0, 2.5))
axis(1, at = c(0, 0.5, 1, 1.5, 2, 2.5))
axis(2, las = 1)
abline(a = 0, b = 1, lty = 'dashed', col = "#85CF00", lw = 2)
points(alphas, alpha_means, pch = 21, 
       bg = paste(col_alpha,'aa',sep=''), cex=1.5, scatter = T)
mtext('\u03B1 Real', 1,cex=1.5,line=3)
mtext('\u03B1 Estimada',2,cex=1.5,line=2.5)
mtext('\u03B1 Real vs \u03B1 Estimada',3,cex=2, line = 1.15)
text(x = 1.5, y = 1, paste("\u03C1: ", round(cor(alphas, alpha_means), 3), sep = ""), cex = 2)


###################################
# Betas de las Categorías
###########


col_beta = "#0078BD"

par(mfrow = c(3, 5))

for (i in 1:15) {
  plot(NULL, ann = F, axes = F, xlim = c(-4, 5), ylim = c(-5, 7), new = T)
  axis(1)
  axis(2, las = 1)
  abline(a = 0, b = 1, lty = 'dashed', col = "#E78140", lw = 2)
  points(item_betas[i, ], beta_means[i, ], pch = 21, 
         bg = paste(col_beta,'aa',sep=''), cex=3, scatter = T)
  mtext('\u03B2 Real', 1,cex=1.5,line=3)
  mtext('\u03B2 Estimada',2,cex=1.5,line=2.5)
  mtext(paste('Betas Item: ', i, sep = ""),3,cex=2, line = 1.15)
}

