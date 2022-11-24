data {
  int<lower=1> D1;
  int<lower=1> D2;
  int<lower=1> D_grouped;
  int<lower=1> A;
  int<lower=1> W;
  int<lower=1> M;
  real ymin;
  real ymax;
  //vector<lower=0,upper=1>[A] V_max;
  int<lower=1> w_vac;
  vector<lower=ymin,upper=ymax>[D1] y1;
  vector<lower=ymin,upper=ymax>[D2] y2;
  vector<lower=ymin,upper=ymax>[D_grouped] y_mean;
  vector<lower=1,upper=fmax(D1,D2)>[D_grouped] n;
  int<lower=1,upper=D_grouped> row1[D1];
  int<lower=1,upper=D_grouped> row2[D2];
  vector<lower=0,upper=W>[D1] tau1;
  vector<lower=0,upper=W>[D2] tau2;
  int<lower=1,upper=W> w1[D1];
  int<lower=1,upper=W> w2[D2];
  int<lower=1,upper=M> m1[D1];
  int<lower=1,upper=M> m2[D2];
  int<lower=1,upper=A> a1[D1];
  int<lower=1,upper=A> a2[D2];
  vector<lower=0,upper=1>[D1] V1;
  vector<lower=0,upper=1>[D2] V2;
  //vector<lower=0,upper=1>[D1] P1;
  //vector<lower=0,upper=1>[D2] P2;
  vector<lower=0,upper=1>[D2] I;
  int<lower=0,upper=1> result[D2]; //"N"=0, "P"=1
  //vector<lower=0>[D2] r_SI;
  vector<lower=0,upper=1>[M] prop_prev_PCR;
  vector<lower=0,upper=1>[M] D2_prop_pos;
  real<lower=0,upper=1> kS;
  vector<lower=0>[A] s;
  vector<lower=0>[A] h;
  vector<lower=0>[A] alpha;
  vector<lower=0>[A] beta;
  vector<lower=0>[A] lambda;
  real<lower=0> r_SV;
}

parameters {
  //parameters
  //vector<lower=ymin>[A] s;
  //vector<lower=0>[A] h;
  //vector<lower=0>[A] alpha;
  //vector<lower=0>[A] beta;
  //vector<lower=0>[A] lambda;
  
  //real<lower=0,upper=1> factor;
  
  real<lower=0> r_SI;
  //real theta_SI_mid;
  //real theta_SI_end;
  //real<lower=0,upper=W> theta_SI_switch;
  
  //real<lower=0> r_SV;
  real<lower=0> r_VP_init;
  real<lower=0> r_VP_end;
  real<lower=w_vac,upper=W> r_VP_switch;
  real<lower=0> psi;
  //real<lower=0, upper=1> psi;
  
  real<lower=0> sigma;
  //real<lower=0> sigma2;
}

transformed parameters {
  //real theta_SI_init = 1.25;
  //real theta_SI_mid = -.21;
  //real theta_SI_end = 1.37;
  //real theta_SI_switch = 18.7;

  
  //real r_VP_switch=w_vac;
  real r_VP;
  
  real PI_PpS;
  real PPp_S;
  real PV_S;
  real PV_noPpS;
  real PV_PpS;
  
  real PSp_noPpS;
  real PSn_noPpS;
  real PSn_noPpIS;
  real PSp_noPpIS;
  real PI_S;
  real PI_noPpS_prime;
  real PI_noPpS_min;
  real PI_noPpS;
  real PSn_noPpVS;
  real PSp_noPpVS;
  real PI_noPpSp;
  real PI_noPpSn;
  real PV_noPpSp;
  real PV_noPpSn;
  
  real no_vac;
  real vac;
  
  //vector[D1] no_vac1;
  //vector[D1] vac1;
  //vector[D2] no_vac2;
  //vector[D2] vac2;
  
  //vector[D1] Sp1;
  //vector[D2] Sp2;
  
  vector[D1] y_pred1;
  vector[D2] y_pred2;
  
  vector[D_grouped] y_pred_grouped = rep_vector(0, D_grouped);
  
  vector[D1] vac_pred1;
  vector[D2] vac_pred2;
  
  //vector[D1] PV_noPpS1;
  //vector[D2] PV_noPpS2;
  //vector[D1] PV_PpS1;
  //vector[D2] PV_PpS2;
  
  for (i in 1:D1) {
    if (w1[i] <= r_VP_switch) {
      r_VP = r_VP_init;
    } else {
      r_VP = r_VP_init + (w1[i]-r_VP_switch)/(W-r_VP_switch)*(r_VP_end-r_VP_init);
    }
    
    PI_PpS=1;
    PPp_S = prop_prev_PCR[m1[i]];
    PV_S = r_SV*V1[i] / (r_SV*V1[i] + (1-V1[i]));
    //PV_noPpS = fmin(PV_S / (r_VP*PPp_S + (1-PPp_S)), V_max[a1[i]]);
    //PV_PpS = fmin((PV_S - PV_noPpS*(1-PPp_S)) / PPp_S, V_max[a1[i]]);
    PV_noPpS = fmin(PV_S / (r_VP*PPp_S + (1-PPp_S)), 1);
    PV_PpS = fmin(PV_S*r_VP / (r_VP*PPp_S + (1-PPp_S)), 1);
    
    no_vac = s[a1[i]] + PI_PpS * h[a1[i]] * gamma_cdf(tau1[i], alpha[a1[i]], beta[a1[i]]) * exp(-lambda[a1[i]]*tau1[i]);
    vac = no_vac + psi;
    
    y_pred1[i] = (1-PV_PpS) * no_vac + PV_PpS * vac;
    y_pred_grouped[row1[i]] += ((1-PV_PpS) * no_vac + PV_PpS * vac) / n[row1[i]];
    
    //Sp1[i] = PI_PpS;
    
    vac_pred1[i] = PV_PpS;
    
    //PV_noPpS1[i] = PV_noPpS;
    //PV_PpS1[i] = PV_PpS;
  }
  for (i in 1:D2) {
    //if (w2[i] <= theta_SI_switch) {
    //  theta_SI = theta_SI_init + (w2[i]-1)/(theta_SI_switch-1) * (theta_SI_mid-theta_SI_init);
    //} else {
    //  theta_SI = theta_SI_mid + (w2[i]-theta_SI_switch)/(W-theta_SI_switch) * (theta_SI_end-theta_SI_mid);
    //}
    //r_SI = exp(theta_SI);
    
    if (w2[i] <= r_VP_switch) {
      r_VP = r_VP_init;
    } else {
      r_VP = r_VP_init + (w2[i]-r_VP_switch)/(W-r_VP_switch)*(r_VP_end-r_VP_init);
    }
    
    PI_PpS=1;
    PPp_S = prop_prev_PCR[m2[i]];
    PV_S = r_SV*V2[i] / (r_SV*V2[i] + (1-V2[i]));
    //PV_noPpS = fmin(PV_S / (r_VP*PPp_S + (1-PPp_S)), V_max[a1[i]]);
    //PV_PpS = fmin((PV_S - PV_noPpS*(1-PPp_S)) / PPp_S, V_max[a1[i]]);
    PV_noPpS = fmin(PV_S / (r_VP*PPp_S + (1-PPp_S)), 1);
    PV_PpS = fmin(PV_S*r_VP / (r_VP*PPp_S + (1-PPp_S)), 1);
    
    PSp_noPpS = D2_prop_pos[m2[i]];
    PSn_noPpS = 1-PSp_noPpS;
    PSn_noPpIS = PV_noPpS*kS^2 + (1-PV_noPpS)*kS;
    PSp_noPpIS = 1-PSn_noPpIS;
    PI_S = r_SI*I[i] / (r_SI*I[i] + (1-I[i]));
    
    PI_noPpS_prime = fmin(fmax((PI_S-PPp_S)/(1-PPp_S), 0), 1);
    PI_noPpS_min = fmin(fmax((PSp_noPpS-(1-kS)*PV_noPpS)/(PSp_noPpIS+kS*(1-kS)*PV_noPpS), 0), 1);
    PI_noPpS = fmax(PI_noPpS_prime, (PI_noPpS_prime+PI_noPpS_min)/2);
    //PI_noPpS_min = fmax((PSp_noPpS-(1-kS)*PV_noPpS)/(PSp_noPpIS+kS*(1-kS)*PV_noPpS), 0);
    //PI_noPpS = fmax((PI_S-PPp_S*PI_PpS) / (1-PPp_S), 0);
    //if (PI_noPpS < PI_noPpS_min) {
    //  PI_noPpS = fmin((PI_noPpS + PI_noPpS_min) / 2, 1);
    //}
    PSn_noPpVS = PI_noPpS*kS^2 + (1-PI_noPpS)*kS;
    PSp_noPpVS = 1-PSn_noPpVS;
    PI_noPpSp = fmin(PSp_noPpIS * PI_noPpS / PSp_noPpS, 1);
    PI_noPpSn = fmin(PSn_noPpIS * PI_noPpS / PSn_noPpS, 1);
    PV_noPpSp = fmin(PSp_noPpVS * PV_noPpS / PSp_noPpS, 1);
    PV_noPpSn = fmin(PSn_noPpVS * PV_noPpS / PSn_noPpS, 1);
    
    if (result[i]==1) {
      no_vac = s[a2[i]] + PI_noPpSp * h[a2[i]] * gamma_cdf(tau2[i], alpha[a2[i]], beta[a2[i]]) * exp(-lambda[a2[i]]*tau2[i]);
      vac = no_vac + psi;
    
      y_pred2[i] = (1-PV_noPpSp) * no_vac + PV_noPpSp * vac;
      y_pred_grouped[row2[i]] += ((1-PV_noPpSp) * no_vac + PV_noPpSp * vac) / n[row2[i]];
      
      //Sp2[i] = PI_noPpSp;
      
      vac_pred2[i] = PV_noPpSp;
    } else {
      no_vac = s[a2[i]] + PI_noPpSn * h[a2[i]] * gamma_cdf(tau2[i], alpha[a2[i]], beta[a2[i]]) * exp(-lambda[a2[i]]*tau2[i]);
      vac = no_vac + psi;
    
      y_pred2[i] = (1-PV_noPpSn) * no_vac + PV_noPpSn * vac;
      y_pred_grouped[row2[i]] += ((1-PV_noPpSn) * no_vac + PV_noPpSn * vac) / n[row2[i]];
      
      //Sp2[i] = PI_noPpSn;
      
      vac_pred2[i] = PV_noPpSn;
    }
    
    //PV_noPpS2[i] = PV_noPpS;
    //PV_PpS2[i] = PV_PpS;
  }
}


model {
  //factor ~ uniform(0,1);
  
  r_SI ~ gamma(3,3);
  //r_SV ~ gamma(3,3);
  r_VP_init ~ gamma(3,3);
  r_VP_end ~ gamma(3,3);
  r_VP_switch ~ uniform(w_vac, W);
  psi ~ exponential(2);
  //psi ~ uniform(0,1);
  
  sigma ~ exponential(1);
  //sigma2 ~ exponential(1);
  
  //likelihood
  y_pred_grouped ~ normal(y_mean, sigma ./ sqrt(n));
  //y_mean ~ normal(y_pred_grouped, sigma ./ sqrt(n));
  //y1 ~ normal(y_pred1, sigma * sqrt(D1));
  //y2 ~ normal(y_pred2, sigma * sqrt(D2));
}
