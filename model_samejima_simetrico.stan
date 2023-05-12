data{
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
}
