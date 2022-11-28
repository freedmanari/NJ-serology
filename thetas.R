### run this file first
require(tidyverse)
require(stringr)
require(lubridate)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(139)


### Upload testing data

#date of first COVID test in NJ
min_date <- as.Date(read.csv("data/q1_pcr.csv")$start[1], "%m/%d/%y")
#converts date to test week, where week 1 is the first epiweek containing any COVID tests in NJ
get_test_week <- function(date) as.numeric(as.Date(date)-floor_date(min_date, unit="week")) %/% 7 + 1
get_test_week_from_epiweek <-
  function(year_week_pairs) unlist(lapply(year_week_pairs, function(ywp) get_test_week(paste(ywp[1],"-01-04",sep="")) + ywp[2] - 1))
get_test_month <- function(date) 12*year(as.Date(date))+month(as.Date(date)) - (12*year(min_date)+month(min_date)) + 1
get_week_start_from_test_week <- function(w) min_date+7*(w-1)

W <- nrow(read.csv("data/q1_pcr.csv")) #number of weeks in testing data = 82
ncols <- ncol(read.csv("data/q1_pcr.csv"))

X1 <- apply(read.csv("data/q1_pcr.csv")[,3:ncols], 1, sum) #first-time positive PCR tests per week
X2 <- apply(read.csv("data/q1_sero.csv")[,3:ncols], 1, sum) #first-time positive serology tests per week

P1 <- apply(read.csv("data/q2_pcr.csv")[,3:ncols], 1, sum) #PCR tests with previous PCR positive per week
P2 <- apply(read.csv("data/q3_pcr.csv")[,3:ncols], 1, sum) #PCR tests without previous PCR positive per week
P <- P1 + P2 #total PCR testing volume per week

S1 <- apply(read.csv("data/q2_sero.csv")[,3:ncols], 1, sum) #serology tests with previous PCR positive per week
S2 <- apply(read.csv("data/q3_sero.csv")[,3:ncols], 1, sum) #serology tests without previous PCR positive per week
S <- S1 + S2 #total serology testing volume per week


Pp <- X1
Pn <- P-Pp
Sp <- X2
Sn <- S-X2

sero_tests <-
  read.csv("data/sero_tests_by_day.csv") %>%
  mutate(test_week=get_test_week(as.Date(Date, "%m/%d/%y"))) %>% 
  filter(!is.na(test_week)) %>% 
  group_by(test_week) %>% 
  summarise(p=sum(Positive.Tests),
            n=sum(Negative.Tests),
            total=sum(Total.Tests)) %>% 
  ungroup() %>%
  mutate(p=round(p/total*S[test_week]), #adjusting to account for discrepancy in total serology testing volume between the two data sets
         n=S[test_week]-p,
         total=S[test_week])
Sp <- rep(0,W)
Sp[sero_tests$test_week] <- sero_tests$p #total serology positives per week
Sn <- rep(0,W)
Sn[sero_tests$test_week] <- sero_tests$n #total serology negatives per week

waiting_times <- #waiting times in weeks from positive PCR to positive serology by PCR positivity week
  read.csv("data/waiting_times.csv")[,3:ncols] %>% 
  apply(1, function(times) paste(times[times != "nnn"], collapse='_')) %>%
  unname %>% 
  lapply(function(times) as.numeric(str_split(times, '_')[[1]])) %>% 
  lapply(function(times) times[times >= 0])
X12 <- sapply(1:W, function(w) sum(sapply(1:w, function(v) sum(waiting_times[[v]]==w-v)))) #positive serology tests with previous positive PCR tests



### NJ demographics

#https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-detail.html
census_data <-
  read.csv("data/sc-est2019-agesex-civ.csv") %>%
  filter(NAME=="New Jersey", SEX==0, AGE<=85) %>% 
  select(age_group_min=AGE, pop=POPEST2019_CIV)

N <- sum(census_data$pop) #New Jersey population size 


# 2021 census population estimates by county, adds up to very different total NJ population size though?
# https://www.census.gov/data/tables/time-series/demo/popest/2020s-counties-total.html
pop_by_county <- read.csv("data/co-est2021-pop-34.csv")




### Mortality data

#Weekly death data
#https://data.cdc.gov/NCHS/Provisional-COVID-19-Death-Counts-by-Week-Ending-D/r8kw-7aab
deaths <- read.csv("data/Provisional_COVID-19_Death_Counts_by_Week_Ending_Date_and_State.csv")
deaths <-
  deaths %>%
  mutate(date = as.Date(End.Date,"%m/%d/%Y")) %>% 
  filter(State == "New Jersey",
         Group == "By Week",
         date > "2020-03-01",
         date < "2021-12-01")
death_dates <- deaths$date
deaths <- deaths %>% pull(COVID.19.Deaths)
deaths[is.na(deaths)] <- 0

W_deaths <- length(deaths)

#Weekly CFR estimates from US
#https://ourworldindata.org/mortality-risk-covid
owid <- read.csv("data/owid-covid-data.csv")
owid <-
  owid %>%
  mutate(date=as.Date(date)) %>% 
  filter(location == "United States",
         date %in% death_dates)
CFR <- owid$total_deaths / owid$total_cases



### Infection fatality ratio (IFR) and time-to-death (TTD) distribution parameters

#Verity et al. (2020)
#https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
IFR_mean <- .00657
IFR_lower <- .00389
IFR_upper <- .0133

conf_diff <- function(sd,level=.95)
  plnorm(IFR_upper,log(IFR_mean^2/sqrt(sd^2+IFR_mean^2)),sqrt(log(sd^2/IFR_mean^2+1)))-
  plnorm(IFR_lower,log(IFR_mean^2/sqrt(sd^2+IFR_mean^2)),sqrt(log(sd^2/IFR_mean^2+1)))-
  level

IFR_sd <- uniroot(conf_diff,c(.001,.005),tol=1e-7)$root
IFR_mu <- log(IFR_mean^2/sqrt(IFR_sd^2+IFR_mean^2))
IFR_sigma <- sqrt(log(IFR_sd^2/IFR_mean^2+1))


#Linton et al. (2020)
#https://www.mdpi.com/2077-0383/9/2/538
TTD_mean <- 14.5/7
TTD_sd <- 6.7/7
TTD_mu <- log(TTD_mean^2/sqrt(TTD_sd^2+TTD_mean^2))
TTD_sigma <- sqrt(log(TTD_sd^2/TTD_mean^2+1))



### age group functions
# takes age group string and gives minimum age of age group
to_age_group_min <- function(age_groups) {
  return(age_groups %>%
           str_extract_all("\\d+|<") %>% 
           lapply(function(x) ifelse(x=="<", 0, suppressWarnings(as.numeric(x)))) %>% 
           lapply(function(x) ifelse(length(x)==0, NA, min(x))) %>% 
           unlist())
}

# gives minimum age of age group corresponding to given age
find_age_group_min <- function(age, age_group_mins) {
  return(age_group_mins[findInterval(age, age_group_mins)])
}

# function to group census data by given age groups

group_census_data <- function(age_group_mins) {
  return(census_data %>%
           group_by(age_group_min=find_age_group_min(age_group_min, age_group_mins)) %>% 
           summarise(pop=sum(pop)) %>%
           ungroup())
}


### age group mins for serology models
model_age_group_mins <- c(0,16,30,50,65)
prop_vaccinated_age_group_mins <-
  read.csv("data/vaccinations_by_age.csv") %>% 
  pull(age_cat) %>%
  to_age_group_min %>% 
  na.omit %>%
  c
prop_infected_age_group_mins <-
  read.csv("data/NJ_incidence_by_age.csv") %>%
  colnames %>%
  grep(pattern="ages", value=TRUE) %>%
  to_age_group_min() %>% 
  sort()

### New vaccinations per week
num_vaccinations_by_age <-
  read.csv("data/vaccinations_by_age.csv") %>%
  pivot_longer(cols=-age_cat, names_to="test_week", values_to="new_vaccinations") %>%
  mutate(test_week=test_week %>% str_remove("X") %>% str_split("_") %>% lapply(as.numeric) %>% get_test_week_from_epiweek,
         age_group_min=to_age_group_min(age_cat)) %>% 
  group_by(test_week) %>%
  summarise(age_group_min=age_group_min,
            new_vaccinations=new_vaccinations,
            week_total=sum(new_vaccinations)) %>%
  ungroup() %>%
  filter(!is.na(age_group_min)) %>% 
  group_by(test_week) %>%
  summarise(age_group_min=age_group_min,
            new_vaccinations=new_vaccinations/sum(new_vaccinations)*week_total)
w_vac <- min(num_vaccinations_by_age$test_week)
V <-
  c(rep(0, w_vac-1),
    num_vaccinations_by_age %>% 
      group_by(test_week) %>%
      summarise(V=sum(new_vaccinations)) %>% 
      pull(V))[1:W]

### Proportion of age groups vaccinated over time
prop_vaccinated_census_data <- group_census_data(prop_vaccinated_age_group_mins)

prop_vaccinated_cap <- .95 #from https://www.cdc.gov/coronavirus/2019-ncov/vaccines/reporting-vaccinations.html

prop_vaccinated_by_age <- 
  num_vaccinations_by_age %>%
  left_join(prop_vaccinated_census_data, by="age_group_min") %>%
  group_by(age_group_min) %>%
  summarise(test_week=test_week,
            prop_vaccinated=cumsum(new_vaccinations)/pop,
            prop_vaccinated_max=min(max(prop_vaccinated), prop_vaccinated_cap),
            prop_vaccinated=prop_vaccinated/max(prop_vaccinated,1/N)*prop_vaccinated_max) %>%
  select(test_week, age_group_min, prop_vaccinated)

### Serology test data for serology model, need here for calculating kS

model_sero <-
  read.csv("data/model_sero.csv") %>%
  mutate(test_date=as.Date(test_date),
         PCR_pos_date=as.Date(PCR_pos_date),
         age_group_min=find_age_group_min(age, model_age_group_mins),
         prop_vaccinated_age_group_min=find_age_group_min(age, prop_vaccinated_age_group_mins),
         prop_infected_age_group_min=find_age_group_min(age, prop_infected_age_group_mins)) %>% 
  left_join(prop_vaccinated_by_age, by=c("test_week","prop_vaccinated_age_group_min"="age_group_min")) %>%
  mutate(across(prop_vaccinated, ~replace_na(., 0)),
         prop_infected = NA,
         y_data = log(numeric)) %>%
  select(test_week,
         test_month,
         test_date,
         PCR_pos_date,
         age,
         age_group_min,
         prop_infected_age_group_min,
         gender,
         result,
         numeric,
         y_data,
         prev_PCR,
         delay,
         delay_int,
         prop_vaccinated,
         prop_infected)


# group serology test results for serology model
model_sero_grouped <-
  model_sero %>%
  group_by(test_week, age_group_min, prev_PCR) %>%
  summarise(y_mean=mean(y_data),
            n=n()) %>%
  ungroup()
model_sero_grouped <- model_sero_grouped %>% mutate(row=1:nrow(model_sero_grouped))
model_sero <-
  model_sero %>%
  left_join(model_sero_grouped, by=c("test_week","age_group_min","prev_PCR"))


model_sero_prev_PCR <- model_sero %>% filter(prev_PCR)
model_sero_no_prev_PCR <- model_sero %>% filter(!prev_PCR)




### Estimating odds ratios (thetas)

Is <- vector(mode="list",length=W)

theta_I_P <- vector(mode="list",length=W)
theta_I_S <- vector(mode="list",length=W)
theta_Pp_P <- vector(mode="list",length=W)
theta_Pp_S <- vector(mode="list",length=W)

reps <- 1000
conf <- .95
show_progress <- TRUE
calculate_odds_ratios <- TRUE

# false negative rates for PCR and serology tests, respectively
kP <- .2 #https://www.acpjournals.org/doi/10.7326/M20-1495
kS <- nrow(filter(model_sero_prev_PCR, test_week<w_vac, result=="N")) /
      nrow(filter(model_sero_prev_PCR, test_week<w_vac))

r <- 1
while (r <= reps) {
  if (show_progress & r %% (reps/100)==0) {
    print(r/(reps/100))
  }

  
  # estimating I as in McCulloh et al. (2020) and Jombart et al. (2020)
  IFR <- rlnorm(1, IFR_mu, IFR_sigma)
  I <- rep(0,W_deaths)
  for (i in 1:W_deaths) {
    n <- floor(deaths[i]/IFR/CFR[i]*mean(CFR))
    death_dist <- sapply(1:i, function(w) plnorm(w, TTD_mu, TTD_sigma) - plnorm(w-1, TTD_mu, TTD_sigma))
    delays <- rmultinom(1,n,prob=death_dist) %>% as.vector()
    I[i:1] <- I[i:1] + delays
  }
  I <- I[1:W]
  Is <- lapply(1:W, function(w) c(Is[[w]], I[w]))
  
  if (calculate_odds_ratios) {
    #entries for contingency table for theta_I_P
    C11_I_P <- round(Pp + kP*I/(kP*I+N-I) * Pn) #X1
    C12_I_P <- round(I - Pp - kP*I/(kP*I+N-I) * Pn) #I-X1
    C21_I_P <- round((N-I)/(kP*I+N-I) * Pn) #P-X1
    C22_I_P <- round(N-I - (N-I)/(kP*I+N-I) * Pn) #N-I-P+X1
    
    #entries for contingency table for theta_I_S
    J <- cumsum(I) + cumsum(V) - cumsum(I)*cumsum(V)/N #MAKE SURE TO CALCULATE V FIRST
    J2 <- cumsum(I)*cumsum(V)/N #MAKE SURE TO CALCULATE V FIRST
    J1 <- J - J2 #MAKE SURE TO CALCULATE V FIRST
    C11_I_S <- round(Sp + (kS*J1+kS^2*J2)/(kS*J1+kS^2*J2+N-J) * Sn) #Sp
    C12_I_S <- round(J - Sp - (kS*J1+kS^2*J2)/(kS*J1+kS^2*J2+N-J) * Sn) #cumsum(I)-Sp
    C21_I_S <- round((N-J)/(kS*J1+kS^2*J2+N-J) * Sn) #Sn
    C22_I_S <- round(N-J - (N-J)/(kS*J1+kS^2*J2+N-J) * Sn) #N-cumsum(I)-Sn
    
    #entries for contingency table for theta_Pp_P
    C11_Pp_P <- P1
    C12_Pp_P <- cumsum(Pp)-P1 #cumsum(X1)-P1
    C21_Pp_P <- P2
    C22_Pp_P <- N-cumsum(Pp)-P2 #N-cumsum(X1)-P2
    
    #entries for contingency table for theta_Pp_S
    C11_Pp_S <- S1
    C12_Pp_S <- cumsum(Pp)-S1 #cumsum(X1)-S1
    C21_Pp_S <- S2
    C22_Pp_S <- N-cumsum(Pp)-S2 #N-cumsum(X1)-S2
    
    for (i in 1:W) {
      C_I_P <- matrix(c(C11_I_P[i],C21_I_P[i],C12_I_P[i],C22_I_P[i]),2)
      C_I_S <- matrix(c(C11_I_S[i],C21_I_S[i],C12_I_S[i],C22_I_S[i]),2)
      C_Pp_P <- matrix(c(C11_Pp_P[i],C12_Pp_P[i],C21_Pp_P[i],C22_Pp_P[i]),2)
      C_Pp_S <- matrix(c(C11_Pp_S[i],C12_Pp_S[i],C21_Pp_S[i],C22_Pp_S[i]),2)
      
      if ((all(C_I_P>0) | length(theta_I_P[[i]])==0) &
          (all(C_I_S>0) | length(theta_I_S[[i]])==0) &
          (all(C_Pp_P>0) | length(theta_Pp_P[[i]])==0) &
          (all(C_Pp_S>0) | length(theta_Pp_S[[i]])==0)) {
        if (all(C_I_P>0)) {
          f <- fisher.test(C_I_P, conf.level=conf)
          mean <- log(f$estimate)
          sd <- (log(f$conf.int[2])-log(f$conf.int[1]))/(2*qnorm((conf+1)/2))
          theta_I_P[[i]] <- c(theta_I_P[[i]], rnorm(1,mean,sd))
        }
        
        if (all(C_I_S>0)) {
          f <- fisher.test(C_I_S, conf.level=conf)
          mean <- log(f$estimate)
          sd <- (log(f$conf.int[2])-log(f$conf.int[1]))/(2*qnorm((conf+1)/2))
          theta_I_S[[i]] <- c(theta_I_S[[i]], rnorm(1,mean,sd))
        }
        
        if (all(C_Pp_P>0)) {
          f <- fisher.test(C_Pp_P, conf.level=conf)
          mean <- log(f$estimate)
          sd <- (log(f$conf.int[2])-log(f$conf.int[1]))/(2*qnorm((conf+1)/2))
          theta_Pp_P[[i]] <- c(theta_Pp_P[[i]], rnorm(1,mean,sd))
        }
        
        if (all(C_Pp_S>0)) {
          f <- fisher.test(C_Pp_S, conf.level=conf)
          mean <- log(f$estimate)
          sd <- (log(f$conf.int[2])-log(f$conf.int[1]))/(2*qnorm((conf+1)/2))
          theta_Pp_S[[i]] <- c(theta_Pp_S[[i]], rnorm(1,mean,sd))
        }
      } else {
        r <- r-1
        
        Is <- lapply(Is, function(w) w[1:r])
        theta_I_P <- lapply(theta_I_P, function(w) w[1:r])
        theta_I_S <- lapply(theta_I_S, function(w) w[1:r])
        theta_Pp_P <- lapply(theta_Pp_P, function(w) w[1:r])
        theta_Pp_S <- lapply(theta_Pp_S, function(w) w[1:r])
        
        break
      }
    }
  }
  
  r <- r+1
}



Is_sorted <- lapply(Is, sort)
I_best <- Is_sorted %>% lapply(function(w) ifelse(is.null(w), NA, mean(w))) %>% unlist
I_lower <- Is_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1-conf)/2)+1])) %>% unlist
I_upper <- Is_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1+conf)/2)])) %>% unlist


if (calculate_odds_ratios) {
  theta_I_P_sorted <- lapply(theta_I_P, sort)
  theta_I_P_best <- theta_I_P_sorted %>% lapply(function(w) ifelse(is.null(w), NA, mean(w))) %>% unlist
  theta_I_P_lower <- theta_I_P_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1-conf)/2)+1])) %>% unlist
  theta_I_P_upper <- theta_I_P_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1+conf)/2)])) %>% unlist
  
  theta_I_S_sorted <- lapply(theta_I_S, sort)
  theta_I_S_best <- theta_I_S_sorted %>% lapply(function(w) ifelse(is.null(w), NA, mean(w))) %>% unlist
  theta_I_S_lower <- theta_I_S_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1-conf)/2)+1])) %>% unlist
  theta_I_S_upper <- theta_I_S_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1+conf)/2)])) %>% unlist
  
  theta_Pp_P_sorted <- lapply(theta_Pp_P, sort)
  theta_Pp_P_best <- theta_Pp_P_sorted %>% lapply(function(w) ifelse(is.null(w), NA, mean(w))) %>% unlist
  theta_Pp_P_lower <- theta_Pp_P_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1-conf)/2)+1])) %>% unlist
  theta_Pp_P_upper <- theta_Pp_P_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1+conf)/2)])) %>% unlist
  
  theta_Pp_S_sorted <- lapply(theta_Pp_S, sort)
  theta_Pp_S_best <- theta_Pp_S_sorted %>% lapply(function(w) ifelse(is.null(w), NA, mean(w))) %>% unlist
  theta_Pp_S_lower <- theta_Pp_S_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1-conf)/2)+1])) %>% unlist
  theta_Pp_S_upper <- theta_Pp_S_sorted %>% lapply(function(w) ifelse(is.null(w), NA, w[round(reps*(1+conf)/2)])) %>% unlist
}


