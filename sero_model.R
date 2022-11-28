### to be run after thetas.R
require(rstan)
require(rmutil)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(68)


### function from simulated weekly incidence to proportion of age groups infected over time 

#https://covid.cdc.gov/covid-data-tracker/#demographicsovertime
prop_infected_census_data <- group_census_data(prop_infected_age_group_mins)

NJ_incidence_by_age <- read.csv("data/NJ_incidence_by_age.csv")

prop_infected_cap <- prop_vaccinated_cap

prop_infected_by_age <-
  function(I) NJ_incidence_by_age %>% 
  pivot_longer(cols=-date_end, names_to="age_group_min", values_to="incidence") %>% 
  mutate(test_week=get_test_week(date_end),
         age_group_min=to_age_group_min(age_group_min)) %>% 
  left_join(prop_infected_census_data, by="age_group_min") %>%
  mutate(new_infections=incidence/1e5*pop) %>% 
  group_by(test_week) %>% 
  summarise(age_group_min=age_group_min,
            prop_of_new_infections=new_infections/sum(new_infections)) %>%
  ungroup() %>%
  group_by(age_group_min) %>% 
  summarise(test_week=test_week,
            num_infected=cumsum(prop_of_new_infections * I)) %>%
  ungroup() %>%
  left_join(prop_infected_census_data) %>%
  group_by(age_group_min) %>% 
  summarise(test_week=test_week,
            prop_infected=num_infected/pop,
            prop_infected_max=min(max(prop_infected), prop_infected_cap),
            prop_infected=prop_infected/max(prop_infected,1/N)*prop_infected_max) %>%
  ungroup() %>% 
  select(test_week, age_group_min, prop_infected)



## preliminaries for simulating delays for tests with no past PCR positive

# fit generalized gamma distributions to the interval weights for different delays
interval_mins <- c(1:5)
interval_weights <- matrix(0, length(interval_mins), W)

interval_probs <- 
  expand.grid(delay_int=sort(unique(model_sero$delay_int[model_sero$delay>=0])),
              interval=1:length(interval_mins)) %>% 
  left_join(model_sero_prev_PCR %>%
            mutate(interval=findInterval(y_data, interval_mins)) %>% 
            group_by(delay_int) %>%
            summarise(interval=interval,
                       p=1/n()) %>% 
            ungroup() %>% 
            group_by(delay_int,interval) %>%
            summarise(p=sum(p)) %>%
            ungroup()) %>%
  replace(is.na(.), 0) %>% 
  mutate(delay=delay_int+.5)


ggamma_fit <- function(data, par) {
  if (any(par<=0)) {
    return(1e10)
  }
  return(sum((data$p - (pggamma(data$delay_int+1, par[1], par[2], par[3])-sapply(data$delay_int, function(d) ifelse(d==0,0,pggamma(d, par[1], par[2], par[3]))))/par[4])^2))
}

interval_weight_pars <- matrix(nrow=length(interval_mins), ncol=4)
for (i in 1:length(interval_mins)) {
  interval_weight_pars[i,] <- optim(par=c(1,1,1,1), fn=ggamma_fit, data=interval_probs[interval_probs$interval==i,])$par
  interval_weights[i,1:length(min(interval_probs$delay):max(interval_probs$delay))] <-
    dggamma(min(interval_probs$delay):max(interval_probs$delay),
            interval_weight_pars[i,1],
            interval_weight_pars[i,2],
            interval_weight_pars[i,3])
}




## calculate the backwards delay distributions from tests with past PCR positive


# with access to sero_data.csv

# delay_dists <-
#   all_sero %>%
#   filter(prev_PCR) %>% 
#   count(test_week, delay_int) %>% 
#   group_by(test_week) %>% 
#   summarise(delay_int=delay_int,
#             n=n,
#             p=n/sum(n)) %>%
#   ungroup() %>% 
#   pivot_wider(id_cols=c(test_week,delay_int),
#               names_from=delay_int,
#               values_from=p,
#               values_fill=0) %>% 
#   select(-test_week)
# delay_dists <- cbind(delay_dists,matrix(0,nrow(delay_dists),W-ncol(delay_dists)))
# delay_dists <- rbind(do.call(rbind,lapply(1:(W-nrow(delay_dists)), function(i) c(rep(1/i,i),rep(0,W-i)))), as.matrix(delay_dists))
# colnames(delay_dists) <- 0:(W-1)


# without access to sero_data.csv

delay_dists <- as.matrix(read.csv("data/delay_dists.csv"))
colnames(delay_dists) <- 0:(W-1)




### run serology model

stan_summary <- function(s) summary(s)$summary[,"mean"]

reps <- 10
show_progress <- TRUE

I_indices <- sample(1:length(Is[[1]]), reps, replace=F)
best_r <- which.min(abs(sapply(I_indices, function(i) sum(unlist(lapply(Is, function(w) w[i]))))-sum(I_best)))
I_indices <- c(I_indices[-best_r], I_indices[best_r])

## model for main text with r_SV=1
r_SV <- 1
sero_model_summs_normal_r_SV <- vector(mode="list", length=reps)

for (r in 1:reps) {
  if (show_progress) {
    print(r)
  }
  
  model_sero <-
    model_sero %>%
    select(-c(prop_infected)) %>%
    left_join(prop_infected_by_age(unlist(lapply(Is, function(w) w[I_indices[r]]))), by=c("test_week","prop_infected_age_group_min"="age_group_min")) #%>%
  model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
  model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)
  
  # simulate delays for serology tests without past PCR positive
  for (i in 1:W) {
    for (j in 1:length(interval_mins)) {
      ix <- which(model_sero_no_prev_PCR$test_week == i &
                  findInterval(model_sero_no_prev_PCR$y_data, interval_mins) == j)
      model_sero_no_prev_PCR$delay[ix] <- sample(seq(0,W-1), length(ix), replace=T, prob=delay_dists[i,]*interval_weights[j,]) + runif(length(ix))
    }
  }
  model_sero_no_prev_PCR$delay_int <- floor(model_sero_no_prev_PCR$delay)
  
  
  sero_model_prev_PCR_normal_r_SV <- stan(file="sero_model_prev_PCR.stan",
                                          data=list(D=nrow(model_sero_prev_PCR),
                                                    A=length(model_age_group_mins),
                                                    W=W,
                                                    M=max(model_sero$test_month)-min(model_sero$test_month)+1,
                                                    ymin=log(3.8),
                                                    ymax=log(400),
                                                    w_vac=w_vac,
                                                    y=model_sero_prev_PCR$y_data,
                                                    tau=model_sero_prev_PCR$delay,
                                                    w=model_sero_prev_PCR$test_week,
                                                    m=model_sero_prev_PCR$test_month-min(model_sero$test_month)+1,
                                                    a=match(model_sero_prev_PCR$age_group_min, model_age_group_mins),
                                                    V=model_sero_prev_PCR$prop_vaccinated,
                                                    prop_prev_PCR=
                                                      model_sero %>%
                                                      count(test_month, prev_PCR) %>% 
                                                      pivot_wider(id_cols=-n,names_from=prev_PCR,values_from=n) %>%
                                                      mutate(prop_prev_PCR=`TRUE`/(`TRUE`+`FALSE`)) %>%
                                                      pull(prop_prev_PCR),
                                                    r_SV=r_SV),
                                          chains=2,
                                          iter=2000)
  summ_prev_PCR_normal_r_SV <- stan_summary(sero_model_prev_PCR_normal_r_SV)
  
  sero_model_normal_r_SV <- stan(file="sero_model.stan",
                                 data=list(D1=nrow(model_sero_prev_PCR),
                                           D2=nrow(model_sero_no_prev_PCR),
                                           D_grouped=nrow(model_sero_grouped),
                                           A=length(model_age_group_mins),
                                           W=W,
                                           M=max(model_sero$test_month)-min(model_sero$test_month)+1,
                                           ymin=log(3.8),
                                           ymax=log(400),
                                           w_vac=w_vac,
                                           y1=model_sero_prev_PCR$y_data,
                                           y2=model_sero_no_prev_PCR$y_data,
                                           y_mean=model_sero_grouped$y_mean,
                                           n=model_sero_grouped$n,
                                           row1=model_sero_prev_PCR$row,
                                           row2=model_sero_no_prev_PCR$row,
                                           tau1=model_sero_prev_PCR$delay,
                                           tau2=model_sero_no_prev_PCR$delay,
                                           w1=model_sero_prev_PCR$test_week,
                                           w2=model_sero_no_prev_PCR$test_week,
                                           m1=model_sero_prev_PCR$test_month-min(model_sero$test_month)+1,
                                           m2=model_sero_no_prev_PCR$test_month-min(model_sero$test_month)+1,
                                           a1=match(model_sero_prev_PCR$age_group_min, model_age_group_mins),
                                           a2=match(model_sero_no_prev_PCR$age_group_min, model_age_group_mins),
                                           V1=model_sero_prev_PCR$prop_vaccinated,
                                           V2=model_sero_no_prev_PCR$prop_vaccinated,
                                           I=model_sero_no_prev_PCR$prop_infected,
                                           result=1*(model_sero_no_prev_PCR$result=="P"),
                                           prop_prev_PCR=model_sero %>%
                                             count(test_month, prev_PCR) %>% 
                                             pivot_wider(id_cols=-n,names_from=prev_PCR,values_from=n) %>%
                                             mutate(prop_prev_PCR=`TRUE`/(`TRUE`+`FALSE`)) %>%
                                             pull(prop_prev_PCR),
                                           D2_prop_pos=model_sero_no_prev_PCR %>%
                                             count(test_month, result) %>%
                                             pivot_wider(id_cols=-n,names_from="result", values_from=n) %>%
                                             mutate(prop_pos=P/(N+P)) %>%
                                             pull(prop_pos),
                                           kS=kS,
                                           s=summ_prev_PCR_normal_r_SV[c("s[1]","s[2]","s[3]","s[4]","s[5]")],
                                           h=summ_prev_PCR_normal_r_SV[c("h[1]","h[2]","h[3]","h[4]","h[5]")],
                                           alpha=summ_prev_PCR_normal_r_SV[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]")],
                                           beta=summ_prev_PCR_normal_r_SV[c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]")],
                                           lambda=summ_prev_PCR_normal_r_SV[c("lambda[1]","lambda[2]","lambda[3]","lambda[4]","lambda[5]")],
                                           r_SV=r_SV),
                                 chains=2,
                                 iter=1000)
  sero_model_summs_normal_r_SV[[r]] <- stan_summary(sero_model_normal_r_SV)
}

model_sero_normal_r_SV <-
  rbind(model_sero_prev_PCR, model_sero_no_prev_PCR)












########## FOR SUPPLEMENT ONLY

## Supplementary sensitivity analysis with low r_SV
r_SV <- .5
sero_model_summs_low_r_SV <- vector(mode="list", length=reps)

for (r in 1:reps) {
  if (show_progress) {
    print(r)
  }
  
  # simulate delays for serology tests without past PCR positive
  model_sero <-
    model_sero %>%
    select(-c(prop_infected)) %>%
    left_join(prop_infected_by_age(unlist(lapply(Is, function(w) w[I_indices[r]]))), by=c("test_week","prop_infected_age_group_min"="age_group_min")) #%>%
  model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
  model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)
  
  for (i in 1:W) {
    for (j in 1:length(interval_mins)) {
      ix <- which(model_sero_no_prev_PCR$test_week == i &
                    findInterval(model_sero_no_prev_PCR$y_data, interval_mins) == j)
      model_sero_no_prev_PCR$delay[ix] <- sample(seq(0,W-1), length(ix), replace=T, prob=delay_dists[i,]*interval_weights[j,]) + runif(length(ix))
    }
  }
  model_sero_no_prev_PCR$delay_int <- floor(model_sero_no_prev_PCR$delay)
  
  
  sero_model_prev_PCR_low_r_SV <- stan(file="sero_model_prev_PCR.stan",
                                       data=list(D=nrow(model_sero_prev_PCR),
                                                 A=length(model_age_group_mins),
                                                 W=W,
                                                 M=max(model_sero$test_month)-min(model_sero$test_month)+1,
                                                 ymin=log(3.8),
                                                 ymax=log(400),
                                                 w_vac=w_vac,
                                                 y=model_sero_prev_PCR$y_data,
                                                 tau=model_sero_prev_PCR$delay,
                                                 w=model_sero_prev_PCR$test_week,
                                                 m=model_sero_prev_PCR$test_month-min(model_sero$test_month)+1,
                                                 a=match(model_sero_prev_PCR$age_group_min, model_age_group_mins),
                                                 V=model_sero_prev_PCR$prop_vaccinated,
                                                 prop_prev_PCR=
                                                   model_sero %>%
                                                   count(test_month, prev_PCR) %>% 
                                                   pivot_wider(id_cols=-n,names_from=prev_PCR,values_from=n) %>%
                                                   mutate(prop_prev_PCR=`TRUE`/(`TRUE`+`FALSE`)) %>%
                                                   pull(prop_prev_PCR),
                                                 r_SV=r_SV),
                                       chains=2,
                                       iter=2000)
  summ_prev_PCR_low_r_SV <- stan_summary(sero_model_prev_PCR_low_r_SV)
  
  sero_model_low_r_SV <- stan(file="sero_model.stan",
                              data=list(D1=nrow(model_sero_prev_PCR),
                                        D2=nrow(model_sero_no_prev_PCR),
                                        D_grouped=nrow(model_sero_grouped),
                                        A=length(model_age_group_mins),
                                        W=W,
                                        M=max(model_sero$test_month)-min(model_sero$test_month)+1,
                                        ymin=log(3.8),
                                        ymax=log(400),
                                        w_vac=w_vac,
                                        y1=model_sero_prev_PCR$y_data,
                                        y2=model_sero_no_prev_PCR$y_data,
                                        y_mean=model_sero_grouped$y_mean,
                                        n=model_sero_grouped$n,
                                        row1=model_sero_prev_PCR$row,
                                        row2=model_sero_no_prev_PCR$row,
                                        tau1=model_sero_prev_PCR$delay,
                                        tau2=model_sero_no_prev_PCR$delay,
                                        w1=model_sero_prev_PCR$test_week,
                                        w2=model_sero_no_prev_PCR$test_week,
                                        m1=model_sero_prev_PCR$test_month-min(model_sero$test_month)+1,
                                        m2=model_sero_no_prev_PCR$test_month-min(model_sero$test_month)+1,
                                        a1=match(model_sero_prev_PCR$age_group_min, model_age_group_mins),
                                        a2=match(model_sero_no_prev_PCR$age_group_min, model_age_group_mins),
                                        V1=model_sero_prev_PCR$prop_vaccinated,
                                        V2=model_sero_no_prev_PCR$prop_vaccinated,
                                        I=model_sero_no_prev_PCR$prop_infected,
                                        result=1*(model_sero_no_prev_PCR$result=="P"),
                                        prop_prev_PCR=model_sero %>%
                                          count(test_month, prev_PCR) %>% 
                                          pivot_wider(id_cols=-n,names_from=prev_PCR,values_from=n) %>%
                                          mutate(prop_prev_PCR=`TRUE`/(`TRUE`+`FALSE`)) %>%
                                          pull(prop_prev_PCR),
                                        D2_prop_pos=model_sero_no_prev_PCR %>%
                                          count(test_month, result) %>%
                                          pivot_wider(id_cols=-n,names_from="result", values_from=n) %>%
                                          mutate(prop_pos=P/(N+P)) %>%
                                          pull(prop_pos),
                                        kS=kS,
                                        s=summ_prev_PCR_low_r_SV[c("s[1]","s[2]","s[3]","s[4]","s[5]")],
                                        h=summ_prev_PCR_low_r_SV[c("h[1]","h[2]","h[3]","h[4]","h[5]")],
                                        alpha=summ_prev_PCR_low_r_SV[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]")],
                                        beta=summ_prev_PCR_low_r_SV[c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]")],
                                        lambda=summ_prev_PCR_low_r_SV[c("lambda[1]","lambda[2]","lambda[3]","lambda[4]","lambda[5]")],
                                        r_SV=r_SV),
                              chains=2,
                              iter=1000)
  sero_model_summs_low_r_SV[[r]] <- stan_summary(sero_model_low_r_SV)
}

model_sero_low_r_SV <-
  rbind(model_sero_prev_PCR, model_sero_no_prev_PCR)



## supplementary sensitivity analysis with high r_SV
r_SV <- 2
sero_model_summs_high_r_SV <- vector(mode="list", length=reps)

for (r in 1:reps) {
  if (show_progress) {
    print(r)
  }
  
  # simulate delays for serology tests without past PCR positive
  model_sero <-
    model_sero %>%
    select(-c(prop_infected)) %>%
    left_join(prop_infected_by_age(unlist(lapply(Is, function(w) w[I_indices[r]]))), by=c("test_week","prop_infected_age_group_min"="age_group_min")) #%>%
  model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
  model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)
  
  for (i in 1:W) {
    for (j in 1:length(interval_mins)) {
      ix <- which(model_sero_no_prev_PCR$test_week == i &
                    findInterval(model_sero_no_prev_PCR$y_data, interval_mins) == j)
      model_sero_no_prev_PCR$delay[ix] <- sample(seq(0,W-1), length(ix), replace=T, prob=delay_dists[i,]*interval_weights[j,]) + runif(length(ix))
    }
  }
  model_sero_no_prev_PCR$delay_int <- floor(model_sero_no_prev_PCR$delay)
  
  
  sero_model_prev_PCR_high_r_SV <- stan(file="sero_model_prev_PCR.stan",
                                        data=list(D=nrow(model_sero_prev_PCR),
                                                  A=length(model_age_group_mins),
                                                  W=W,
                                                  M=max(model_sero$test_month)-min(model_sero$test_month)+1,
                                                  ymin=log(3.8),
                                                  ymax=log(400),
                                                  w_vac=w_vac,
                                                  y=model_sero_prev_PCR$y_data,
                                                  tau=model_sero_prev_PCR$delay,
                                                  w=model_sero_prev_PCR$test_week,
                                                  m=model_sero_prev_PCR$test_month-min(model_sero$test_month)+1,
                                                  a=match(model_sero_prev_PCR$age_group_min, model_age_group_mins),
                                                  V=model_sero_prev_PCR$prop_vaccinated,
                                                  prop_prev_PCR=
                                                    model_sero %>%
                                                    count(test_month, prev_PCR) %>% 
                                                    pivot_wider(id_cols=-n,names_from=prev_PCR,values_from=n) %>%
                                                    mutate(prop_prev_PCR=`TRUE`/(`TRUE`+`FALSE`)) %>%
                                                    pull(prop_prev_PCR),
                                                  r_SV=r_SV),
                                        chains=2,
                                        iter=2000)
  summ_prev_PCR_high_r_SV <- stan_summary(sero_model_prev_PCR_high_r_SV)
  
  sero_model_high_r_SV <- stan(file="sero_model.stan",
                               data=list(D1=nrow(model_sero_prev_PCR),
                                         D2=nrow(model_sero_no_prev_PCR),
                                         D_grouped=nrow(model_sero_grouped),
                                         A=length(model_age_group_mins),
                                         W=W,
                                         M=max(model_sero$test_month)-min(model_sero$test_month)+1,
                                         ymin=log(3.8),
                                         ymax=log(400),
                                         w_vac=w_vac,
                                         y1=model_sero_prev_PCR$y_data,
                                         y2=model_sero_no_prev_PCR$y_data,
                                         y_mean=model_sero_grouped$y_mean,
                                         n=model_sero_grouped$n,
                                         row1=model_sero_prev_PCR$row,
                                         row2=model_sero_no_prev_PCR$row,
                                         tau1=model_sero_prev_PCR$delay,
                                         tau2=model_sero_no_prev_PCR$delay,
                                         w1=model_sero_prev_PCR$test_week,
                                         w2=model_sero_no_prev_PCR$test_week,
                                         m1=model_sero_prev_PCR$test_month-min(model_sero$test_month)+1,
                                         m2=model_sero_no_prev_PCR$test_month-min(model_sero$test_month)+1,
                                         a1=match(model_sero_prev_PCR$age_group_min, model_age_group_mins),
                                         a2=match(model_sero_no_prev_PCR$age_group_min, model_age_group_mins),
                                         V1=model_sero_prev_PCR$prop_vaccinated,
                                         V2=model_sero_no_prev_PCR$prop_vaccinated,
                                         I=model_sero_no_prev_PCR$prop_infected,
                                         result=1*(model_sero_no_prev_PCR$result=="P"),
                                         prop_prev_PCR=model_sero %>%
                                           count(test_month, prev_PCR) %>% 
                                           pivot_wider(id_cols=-n,names_from=prev_PCR,values_from=n) %>%
                                           mutate(prop_prev_PCR=`TRUE`/(`TRUE`+`FALSE`)) %>%
                                           pull(prop_prev_PCR),
                                         D2_prop_pos=model_sero_no_prev_PCR %>%
                                           count(test_month, result) %>%
                                           pivot_wider(id_cols=-n,names_from="result", values_from=n) %>%
                                           mutate(prop_pos=P/(N+P)) %>%
                                           pull(prop_pos),
                                         kS=kS,
                                         s=summ_prev_PCR_high_r_SV[c("s[1]","s[2]","s[3]","s[4]","s[5]")],
                                         h=summ_prev_PCR_high_r_SV[c("h[1]","h[2]","h[3]","h[4]","h[5]")],
                                         alpha=summ_prev_PCR_high_r_SV[c("alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]")],
                                         beta=summ_prev_PCR_high_r_SV[c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]")],
                                         lambda=summ_prev_PCR_high_r_SV[c("lambda[1]","lambda[2]","lambda[3]","lambda[4]","lambda[5]")],
                                         r_SV=r_SV),
                               chains=2,
                               iter=1000)
  sero_model_summs_high_r_SV[[r]] <- stan_summary(sero_model_high_r_SV)
}

model_sero_high_r_SV <-
  rbind(model_sero_prev_PCR, model_sero_no_prev_PCR)




