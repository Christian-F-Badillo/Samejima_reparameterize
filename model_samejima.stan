
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
  }
  
