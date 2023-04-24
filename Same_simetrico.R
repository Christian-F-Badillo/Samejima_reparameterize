####### Simulación del modelo de Samejima (GRM) ################
####### Y ecuperación de parámetros por Bayes ##################

##########################################################################################
### Modelo modificado de dos parametros


same_repa_acum <- function(theta, alpha, delta, tau, num_cate){
  
  if (length(theta) > 1){
    result <- sapply(theta, FUN = same_repa_acum, alpha = alpha, delta = delta, tau = tau, num_cate = num_cate)
    return(t(result))
  }
  # El bloque anterior permite que la función se pueda aplicar a un vector de thetas.
      
    acum <- c()  
  
    acum[2] <- exp(alpha * (theta - delta)) / (1 + exp(alpha * (theta - delta)))
    acum[1] <- exp(alpha * (theta - (delta - tau))) / (1 + exp(alpha * (theta - (delta - tau))))
    acum[3] <- exp(alpha * (theta - (delta + tau))) / (1 + exp(alpha * (theta - (delta + tau))))
    
    return(acum)
  
}

#######################################################################################################
# Calcular las probabilidades de responder en cada categoría

Prob_SR <- function(theta, alpha, delta, tau, num_cate){
  #if (any(sort(beta) != beta))
   # stop("Los elementos en el vector beta deben estar ordenados")
  # Revisamos que los parámetros tengan valores válidos
  
  if (length(theta) > 1){
    result <- sapply(theta, FUN = Prob_SR, alpha = alpha, delta = delta, tau = tau, num_cate = num_cate)
    return(result)
  }
  
  cat <- seq(0, length(num_cate) - 1)
  
  ProbCums <- same_repa_acum(theta = theta, alpha = alpha, delta = delta, tau = tau, num_cate = num_cate)
  # Calculamos primero las probabilidades acumuladas
  result <- c()
  for (j in 1:(num_cate-1)){
    result <- c(result, ProbCums[j] - ProbCums[j + 1])
    # Calculamos para cada categoría la diferencia entre la probabilidad acumulada de esta categoría y de la
    #      siguiente categoría.
  }
  result <- c(result, ProbCums[length(ProbCums)])
  # La probabilidad de responder en la categoría más alta es igual a la probabilidad acumulada de esta categoría.
  return(result[cat+1])
}

#######################################################################################################
# Parámetros a recobrar.

# Deltas's item
set.seed(1)
item_delta <- rnorm(15, 0, 1)

print(item_delta)

set.seed(14)
thetas_p <- rnorm(1000, 0, 1)
hist(thetas)
summary(thetas)

set.seed(13)
alphas_i <- rgamma(15, 2, 2)
print(alphas)

set.seed(12)
taus_p <- rlnorm(15, 0, 1)
print(taus)

cate = 4
#######################################################################################################
# Simulación de datos

simulated_data_GRM <- function(thetas, alphas, delta, tau, num_cate){
  data = matrix(data = NA, nrow = 1000, ncol = 15)
  set.seed(10)
  for (item in 1:15) {
    prob_cate <- Prob_SR(theta = thetas_p, alpha = alphas[item], delta = delta[item], tau = tau, num_cate = num_cate)
    for (person in 1:1000) {
      prob = as.vector(prob_cate[person, ])
      resp_item = sample(c(0, 1, 2, 3), size = 1, replace = FALSE, prob = prob)   
      data[person, item] = resp_item
    }   
    
  }
  return(data)
}

generated_data <- simulated_data_GRM(thetas_p, alphas_i, item_delta, taus_p, cate)

#write.csv(item_betas, file = "betas.csv")
#write.csv(alphas, file = "alphas.csv")
#write.csv(thetas, file = "thetas.csv")
#write.csv(generated_data, file = "sim_data.csv")

acum <- same_repa_acum(thetas_p, alphas_i, item_delta, taus_p, cate)

probs <- Prob_SR(thetas_p, alphas_i, item_delta, taus_p, cate)

print(probs)
