### to be run after thetas.R
require(ggeasy)
require(bayestestR)
require(tidyverse)
require(colorspace)
require(rstan)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# make nice-looking scientific notation for plots
fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  l <- gsub("\\+","",l)
  parse(text=l)
}


plot_theta <- function(theta_best, theta_lower, theta_upper, zero_line=T, ylab="log odds ratio") {
  data <- data.frame(date=get_week_start_from_test_week(1:W),
                     mean=theta_best,
                     CI_lower=theta_lower,
                     CI_upper=theta_upper)
  if (any(is.na(data))) {
    data <- data[(tail(which(apply(data,1,function(r) any(is.na(r)))),1)+1):W,]
  }
  plot <-
    data %>%
    pivot_longer(cols=-date, names_to="quantile", values_to="theta") %>% 
    mutate(type=ifelse(quantile=="mean","mean","CI")) %>% 
    ggplot() +
    geom_line(aes(date,theta,group=quantile,linetype=type,color=type)) +
    scale_linetype_manual(values=c("mean"=1,"CI"=2)) +
    scale_color_manual(values=c("mean"="black","CI"="steelblue2")) +
    scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=c(0,0)) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.text.x = element_text(angle = 90),
          legend.position = "none",
          axis.line = element_line(colour = "black", size=.3))
    if (zero_line) {
      plot +
        geom_hline(yintercept=0, col="red", size=.3) +
        scale_y_continuous(name=ylab)
    } else {
      plot +
        scale_y_continuous(name=ylab, labels=fancy_scientific, expand=expansion(c(0,.03))) +
        coord_cartesian(ylim = c(0,NA))
    }
}


plot_theta(I_best, I_lower, I_upper, zero_line=F, ylab="new infections")



plot_theta(theta_I_P_best, theta_I_P_lower, theta_I_P_upper)

plot_theta(theta_I_S_best, theta_I_S_lower, theta_I_S_upper)

plot_theta(theta_Pp_P_best, theta_Pp_P_lower, theta_Pp_P_upper)

plot_theta(theta_Pp_S_best, theta_Pp_S_lower, theta_Pp_S_upper)














