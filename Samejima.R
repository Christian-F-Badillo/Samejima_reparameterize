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

write.csv(generated_data, file = "datos.csv")


#########################################################################################################
# Bayes Modelling.

n_per <- dim(generated_data)[1]
n_items <- dim(generated_data)[2]
n_cate <- 4 # Número de categorías es igual m-1.


obs_var <- list("n_per", "n_items", "n_cate", "generated_data")

unobs_var <- c("theta_prior_p", "alpha_prior_i",
               "theta_pos_p", "alpha_pos_i", 
               "beta_prior_ij", "beta_pos_ij")

write(
    "
    model{
    
    # Priors sobre personas.
    
        for (p in 1:n_per){
            theta_prior_p[p] ~ dnorm(0, 5)
            theta_pos_p[p] ~ dnorm(0, 5)
        }
    
    # Priors sobre alphas items
        for (a in 1:n_items){
            alpha_prior_i[a] ~ dgamma(1.5, 1)  T(0,)
            alpha_pos_i[a] ~ dgamma(1.5, 1) T(0,)
        }
    
    # Priors sobre betas de categorías e ítems
        for(i in 1:n_items){
            for(j in 1:n_cate){
    
                beta_prior_ij[i, j] ~ dnorm(0, 5)
                beta_pos_ij[i, j] ~ dnorm(0, 5)
            }
        }
        
       # for (i in n_items){
        #    betas[1:n_cate] <- beta_pos_ij[i, ]
         #   betas[1:n_cate] <- sort(betas)
    #
     #       for (j in 1:n_cate){
      #          betas_pos[i, j] <- betas[j]
       # } }
    
    # Probabilidades acumuladas
    
    for (persona in 1:n_per) {
            for (item in 1:n_items) {
                
                generated_data[persona,item] ~ dordered.logit(alpha_pos_i[item] * theta_pos_p[persona], beta_pos_ij[item,])
            
            }
        }


    
    }    
", "model_samejima.bug")

library(R2jags)

initsvalues = list(list(theta_prior_p = runif(1, -1,1), alpha_prior_0= runif(15, 0.5,2), theta_pos_p=runif(1, -1,1), beta_prior_ij=runif(1, -2,2), beta_pos_ij=matrix(0.01, ncol = 4, nrow = 15)), 
                   list(theta_prior_p = runif(1, -1,1), alpha_prior_i=runif(15, 0.5,2), theta_pos_p=runif(1, -1,1), beta_prior_ij=runif(1, -2,2), beta_pos_ij=matrix(0.1, ncol = 4, nrow = 15)), 
                   list(theta_prior_p = runif(1, -1,1), alpha_prior_i=runif(15, 0.5,2), theta_pos_p=runif(1, -1,1), beta_prior_ij=runif(1, -2,2), beta_pos_ij=matrix(0.1, ncol = 4, nrow = 15)), 
                   list(theta_prior_p = runif(1, -1,1), alpha_prior_i=runif(15, 0.5,2), theta_pos_p=runif(1, -1,1), beta_prior_ij=runif(1, -2,2), beta_pos_ij=matrix(0.01, ncol = 4, nrow = 15)))

set.seed(10156542)
bayes_samejima <- jags(
    data = obs_var, inits = initsvalues,
    parameters.to.save = unobs_var,
    model.file = "model_samejima.bug",
    n.chains = 4,
    n.iter = 25000,
    n.burnin = 2000,
    n.thin = 4 ) 

