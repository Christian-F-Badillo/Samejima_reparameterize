####### Simulación del modelo de Samejima (GRM) ################
####### Y ecuperación de parámetros por Bayes ##################

##########################################################################################
### Modelo modificado de dos parametros


same_repa_acum <- function(theta, alpha, delta, tau, num_cate){
  
  if (length(theta) > 1){
    result <- sapply(theta, FUN = same_repa_acum, alpha = alpha, delta = delta, tau = tau, num_cate)
    return(result)
  }
  # El bloque anterior permite que la función se pueda aplicar a un vector de thetas.
  
  acum <- c()
  
  if (num_cate %% 2 == 0){
    acum[num_cate / 2] <- exp(alpha * (theta - (delta))) / (1 + exp(alpha * (theta - (delta))))
    for (i in 1:(n_cate - 1)) {
      if (i != n_cat/2){
        acum[i] <- exp(alpha * (theta - (delta - tau[i]))) / (1 + exp(alpha * (theta - (delta - tau[i]))))
      }
    }
    return(acum)
  }
  
  if (num_cate %% 2 != 0){
    for (i in 1:(n_cate - 1)) {
      acum[i] <- exp(alpha * (theta - (delta - tau[i]))) / (1 + exp(alpha * (theta - (delta - tau[i]))))
    }
    return(acum)    
  }
  
  
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

#write.csv(item_betas, file = "betas.csv")
#write.csv(alphas, file = "alphas.csv")
#write.csv(thetas, file = "thetas.csv")
#write.csv(generated_data, file = "sim_data.csv")