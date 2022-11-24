data {
  int<lower=1> D;
  int<lower=1> A;
  int<lower=1> W;
  int<lower=1> M;
  real ymin;
  real ymax;
  int<lower=1> w_vac;
  vector<lower=ymin,upper=ymax>[D] y;
  vector<lower=0,upper=W>[D] tau;
  int<lower=1,upper=W> w[D];
  int<lower=1,upper=M> m[D];
  int<lower=1,upper=A> a[D];
  vector<lower=0,upper=1>[D] V;
  vector<lower=0,upper=1>[M] prop_prev_PCR;
  real<lower=0> r_SV;
}

parameters {
  vector<lower=ymin>[A] s;
  vector<lower=0>[A] h;
  vector<lower=0>[A] alpha;
  vector<lower=0>[A] beta;
  vector<lower=0>[A] lambda;
  
  real<lower=0> r_VP_init;
  real<lower=0> r_VP_end;
  real<lower=w_vac,upper=W> r_VP_switch;
  real<lower=0> psi;
  
  real<lower=0> sigma;
}


transformed parameters {
  real r_VP;
  
  real PI_PpS;
  real PPp_S;
  real PV_S;
  real PV_PpS;
  
  real no_vac;
  real vac;
  
  vector[D] y_pred;
  
  vector[W+1] y_pred_age1;
  vector[W+1] y_pred_age2;
  vector[W+1] y_pred_age3;
  vector[W+1] y_pred_age4;
  vector[W+1] y_pred_age5;
  
  for (i in 1:D) {
    if (w[i] <= r_VP_switch) {
      r_VP = r_VP_init;
    } else {
      r_VP = r_VP_init + (w[i]-r_VP_switch)/(W-r_VP_switch)*(r_VP_end-r_VP_init);
    }
    
    PI_PpS = 1;
    PPp_S = prop_prev_PCR[m[i]];
    PV_S = r_SV*V[i] / (r_SV*V[i] + (1-V[i]));
    PV_PpS = fmin(PV_S*r_VP / (r_VP*PPp_S + (1-PPp_S)), 1);
    
    no_vac = s[a[i]] + PI_PpS * h[a[i]] * gamma_cdf(tau[i], alpha[a[i]], beta[a[i]]) * exp(-lambda[a[i]]*tau[i]);
    vac = no_vac + psi;
    
    y_pred[i] = (1-PV_PpS) * no_vac + PV_PpS * vac;
  }
  
  for (delay_int in 0:W) {
    y_pred_age1[delay_int+1] = s[1] + PI_PpS * h[1] * gamma_cdf(delay_int+.5, alpha[1], beta[1]) * exp(-lambda[1]*(delay_int+.5));
    y_pred_age2[delay_int+1] = s[2] + PI_PpS * h[2] * gamma_cdf(delay_int+.5, alpha[2], beta[2]) * exp(-lambda[2]*(delay_int+.5));
    y_pred_age3[delay_int+1] = s[3] + PI_PpS * h[3] * gamma_cdf(delay_int+.5, alpha[3], beta[3]) * exp(-lambda[3]*(delay_int+.5));
    y_pred_age4[delay_int+1] = s[4] + PI_PpS * h[4] * gamma_cdf(delay_int+.5, alpha[4], beta[4]) * exp(-lambda[4]*(delay_int+.5));
    y_pred_age5[delay_int+1] = s[5] + PI_PpS * h[5] * gamma_cdf(delay_int+.5, alpha[5], beta[5]) * exp(-lambda[5]*(delay_int+.5));
  }
}

model {
  s ~ lognormal(.2,.3);
  h ~ lognormal(-.3,.3);
  alpha ~ lognormal(1.1,.3);
  beta ~ lognormal(.4,.3);
  lambda ~ lognormal(-5,.3);
  
  r_VP_init ~ gamma(3,3);
  r_VP_end ~ gamma(3,3);
  r_VP_switch ~ uniform(w_vac, W);
  psi ~ exponential(2);
  
  sigma ~ exponential(1);
  
  //likelihood
  y ~ normal(y_pred, sigma);
}

