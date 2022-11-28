### to be run after thetas.R
require(tidyverse)
require(cowplot)
require(tigris)
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





# population size of each NJ region
north_counties <- c("Bergen","Hudson","Essex","Morris","Passaic","Warren","Sussex")
central_counties <- c("Hunterdon","Somerset","Union","Mercer","Middlesex")
shore_counties <- c("Cape May","Atlantic","Ocean","Monmouth")
south_counties <- c("Cumberland","Salem","Gloucester","Camden","Burlington")

pop_by_region <-
  pop_by_county %>%
  mutate(region =
           case_when(county %in% north_counties ~ "north",
                     county %in% central_counties ~ "central",
                     county %in% shore_counties ~ "shore",
                     county %in% south_counties ~ "south",
                     TRUE ~ "unknown")) %>%
  group_by(region) %>%
  summarise(pop=sum(pop))




# PCR positives broken up by NJ region
X1_region <- read.csv("data/q1_pcr.csv")
X1_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(X1_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(X1_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(X1_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(X1_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))

X1_region <- X1_region %>% pivot_longer(cols=-"date", names_to="region", values_to="X1")
X1_region$region <- factor(X1_region$region, levels=unique(X1_region$region))



# seropositives broken up by NJ region
X2_region <- read.csv("data/q1_sero.csv")
X2_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(X2_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(X2_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(X2_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(X2_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))

X2_region <- X2_region %>% pivot_longer(cols=-"date", names_to="region", values_to="X2")
X2_region$region <- factor(X2_region$region, levels=unique(X2_region$region))





county_list <- read.csv("data/q1_pcr.csv")[,-c(1,2,8)] %>% colnames()
X12_region <-
  cbind(test_week=1:W, read.csv("data/waiting_times.csv")[,-c(1,2,8)]) %>%
  pivot_longer(cols=-"test_week", names_to="county", values_to="delays") %>% 
  filter(delays!="nnn") %>%
  group_by(test_week, county) %>%
  summarise(delay=as.numeric(str_split(delays, '_')[[1]])) %>%
  ungroup() %>%
  filter(delay>=0) %>% 
  group_by(county, date=get_week_start_from_test_week(test_week+delay)) %>%
  summarise(X12=n()) %>%
  ungroup() %>% 
  right_join(expand.grid(date=get_week_start_from_test_week(1:W),county=county_list)) %>% 
  replace_na(list(X12=0)) %>%
  pivot_wider(id_cols=-X12, names_from="county", values_from=X12) %>%
  arrange(date)
X12_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(X12_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(X12_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(X12_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(X12_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))

#seropositives with past PCR positive broken up by NJ region
X12_region <- X12_region %>% pivot_longer(cols=-"date", names_to="region", values_to="X12")
X12_region$region <- factor(X12_region$region, levels=unique(X12_region$region))







## summary plot part a, total cases over study period of all three types, by NJ region
X1_region %>%
  right_join(X2_region) %>%
  right_join(X12_region) %>%
  group_by(region) %>% 
  summarise(X1=sum(X1),
            X2=sum(X2),
            X12=sum(X12)) %>%
  ungroup() %>%
  right_join(pop_by_region) %>%
  mutate(pop=pop/10) %>% 
  pivot_longer(cols=-region, names_to="type", values_to="cases") %>% 
  mutate(type=factor(type, levels=c("X1","X2","X12","pop")),
         region=factor(region, levels=c("north","central","shore","south"))) %>%
  ggplot() +
  geom_col(aes(x=region, y=cases, fill=type, color=type, linetype=type), position="dodge",
           width=.82) +
  scale_fill_manual(values=c("red3","chartreuse4","dodgerblue3","gray97"),
                    labels=c("PCR positives","seropositives","seropositives","")) +   #add "with past PCR positive" to last legend item when editing
  scale_color_manual(values=c("red3","chartreuse4","dodgerblue3","black"),
                     guide="none") +
  scale_linetype_manual(values=c("solid","solid","solid","dashed"),
                     guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,.05)),labels=fancy_scientific,
                     sec.axis = sec_axis(~ ./10, name = "population",labels=fancy_scientific)) +
  ylab("total cases") +
  theme(panel.background=element_blank(),
        axis.line=element_line(size=.4),
        legend.title = element_blank(),
        legend.key.size = unit(17,"pt"),
        legend.text = element_text(size=10),
        legend.position="bottom",
        axis.line.y.right=element_line(linetype="longdash"))






## summary figure part b

# line graph of weekly cases for all three categories
p1 <- data.frame(date=rep(get_week_start_from_test_week(1:(W-1)),3),
           type=rep(c("PCR positives","Seropositives","Seropositives with past PCR positive"),each=(W-1)),
           cases=c(X1[1:(W-1)],X2[1:(W-1)],X12[1:(W-1)])) %>%
      ggplot() +
      geom_line(aes(date,cases,group=type,color=type),size=1) +
      ylab("weekly cases") +
      scale_color_manual(values=c("red3","chartreuse4","dodgerblue3")) +
      scale_x_date(expand=expansion(c(1/(2*(W-2))),1/(2*(W-2)))) +
      scale_y_continuous(lim=c(0,40000), expand=expansion(c(.01,.02)), position="right", labels=fancy_scientific) +
      theme(panel.grid=element_blank(),
            panel.border=element_blank(),
            panel.background=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            axis.line=element_line(size=.3),
            legend.position="none",
            plot.margin = unit(c(2.5,5.5,-2.3,5.5), "pt"),
            legend.key.size = unit(12, "pt"))

# heatmap of new PCR positives
p2 <- X1_region %>% 
      group_by(region) %>%
      summarise(date=date,X1=X1/max(X1)) %>% 
      filter(date < get_week_start_from_test_week(W)) %>% 
      ggplot() +
      geom_tile(aes(x=date, y=region, fill=X1, color=X1), size=1) +
      scale_fill_gradientn(colors = c("white","red3"),
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      scale_color_gradientn(colors = c("white","red3"),
                            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      scale_x_date(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      ylab("") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.x=element_line(size=.3),
            axis.line.y=element_blank(),
            plot.margin = unit(c(0,5.5,-2.3,5.5), "pt"),
            legend.key.height = unit(11, "pt"),
            legend.key.width = unit(14, "pt"),
            legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(0,0,0,0),
            legend.position="none")

# heatmap of new seropositives
p3 <- X2_region %>%
      group_by(region) %>%
      summarise(date=date,X2=X2/max(X2)) %>% 
      filter(date < get_week_start_from_test_week(W)) %>% 
      ggplot() +
      geom_tile(aes(x=date, y=region, fill=X2, color=X2), size=1) +
      scale_fill_gradientn(colors=c("white","chartreuse4"),
                           guide=guide_colorbar(frame.color="black",ticks.color="black")) +
      scale_color_gradientn(colors = c("white","chartreuse4"),
                            guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      scale_x_date(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.line.x=element_line(size=.3),
            axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.margin = unit(c(0,5.5,-2.3,5.5), "pt"),
            legend.key.height = unit(11, "pt"),
            legend.key.width = unit(14, "pt"),
            legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(0,0,0,0),
            legend.position="none")

# heatmap of new seropositives with past PCR positive
p4 <- X12_region %>%
      group_by(region) %>%
      summarise(date=date,X12=X12/max(X12)) %>% 
      filter(date < get_week_start_from_test_week(W)) %>% 
      ggplot() +
      geom_tile(aes(x=date, y=region, fill=X12, color=X12), size=1) +
      scale_fill_gradientn(colors = c("white","dodgerblue3"),
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      scale_color_gradientn(colors = c("white","dodgerblue3"),
                           guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
      scale_x_date(date_labels="%b '%y", date_breaks="1 month", expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      ylab("") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.x=element_line(size=.4),
            axis.line.y=element_blank(),
            plot.margin = unit(c(0,5.5,0,5.5), "pt"),
            legend.key.width = unit(14, "pt"),
            legend.key.height = unit(11, "pt"),
            legend.title = element_blank(),
            legend.margin=margin(0,0,0,0),
            legend.box.margin=margin(0,0,0,0),
            legend.position="none")

# complete summary plot part b
plot_grid(p1, p2, p3, p4, ncol=1, align="v", axis=c("lr"), rel_heights=c(2,1,1,1))

# to get grayscale legend and x axis labels to add into summary plot part b
gray <-
  X1_region %>% 
  group_by(region) %>%
  summarise(date=date,X1=X1/max(X1)) %>% 
  filter(date < get_week_start_from_test_week(W)) %>% 
  ggplot() +
  geom_tile(aes(x=date, y=region, fill=X1)) +
  scale_fill_gradientn(colors = c("white","gray25"), limits=c(0,1),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ylab("region") +
  theme(axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle = 90),
        axis.line.x=element_line(size=.3),
        axis.line.y=element_blank(),
        plot.margin = unit(c(0,5.5,-2.3,5.5), "pt"),
        legend.key.height = unit(11, "pt"),
        legend.key.width = unit(14, "pt"),
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0),
        legend.position="right")
plot_grid(p1, gray, ncol=1, align="v", axis=c("lr"), rel_heights=c(2,1))







## summary plot part c, heatmaps of proportion of seropositive tests with past PCR positive, by region

rainbow_gradient <- c("white", rev(rainbow(7))[-(1:2)])

prop_prev_pcr <-
  left_join(X2_region,X12_region) %>%
  group_by(date) %>%
  summarise(X2=c(X2, sum(X2)),
            X12=c(X12, sum(X12)),
            region=factor(c(as.character(region),"total"),
                          levels=c("south","shore","central","north","total"))) %>%
  ungroup() %>% 
  mutate(prop_prev_pcr=ifelse(X2==0, NA, X12/X2)) %>% 
  select(-c(X2,X12))

prop_prev_pcr %>%
  filter(date < get_week_start_from_test_week(W)) %>% 
  ggplot() +
  geom_tile(aes(x=date, y=region, fill=prop_prev_pcr, color=prop_prev_pcr), size=1) +
  scale_fill_gradientn(colors = rainbow_gradient, na.value = "gray80",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                       name = stringr::str_wrap("proportion of seropositives with past PCR positive",
                                                width = 20)) +
  scale_color_gradientn(colors = rainbow_gradient, na.value = "gray80", guide="none") +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x=element_text(angle = 90),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(color="black", size=.3),
        axis.line.y=element_blank(),
        legend.title=element_text(size=10),
        legend.position="right")




## summary plot part d, map of NJ regions

NJ_counties <-
  counties(state="NJ", cb=T, class="sf") %>% 
  mutate(region =
           case_when(NAME %in% north_counties ~ "north",
                     NAME %in% central_counties ~ "central",
                     NAME %in% shore_counties ~ "shore",
                     NAME %in% south_counties ~ "south",
                     TRUE ~ "unknown"))

NJ_counties_grouped <-
  NJ_counties %>% 
  group_by(region) %>% 
  summarise() %>%
  ungroup()


ggplot() +
  geom_sf(data=NJ_counties, fill="white", color="gray80", size=0.2) +
  geom_sf(data=NJ_counties_grouped, fill="transparent", color="gray20", size=.8) +
  geom_sf_text(data=NJ_counties_grouped,aes(label=region),size=4) +
  theme_void() +
  coord_sf(ndiscr = F)







## summary plots of testing volumes for each all categories

# PCR positives with past PCR positive, broken up by region
P1_region <- read.csv("data/q2_pcr.csv")
P1_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(P1_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(P1_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(P1_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(P1_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))
P1_region <- P1_region %>% pivot_longer(cols=-"date", names_to="region", values_to="P1")
P1_region$region <- factor(P1_region$region, levels=unique(P1_region$region))

# PCR positives without past PCR positive, broken up by region
P2_region <- read.csv("data/q3_pcr.csv")
P2_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(P2_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(P2_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(P2_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(P2_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))
P2_region <- P2_region %>% pivot_longer(cols=-"date", names_to="region", values_to="P2")
P2_region$region <- factor(P2_region$region, levels=unique(P2_region$region))


# seropositives with past PCR positive, broken up by region
S1_region <- read.csv("data/q2_sero.csv")
S1_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(S1_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(S1_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(S1_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(S1_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))
S1_region <- S1_region %>% pivot_longer(cols=-"date", names_to="region", values_to="S1")
S1_region$region <- factor(S1_region$region, levels=unique(S1_region$region))

# seropositives without past PCR positive, broken up by region
S2_region <- read.csv("data/q3_sero.csv")
S2_region <-
  data.frame(date=get_week_start_from_test_week(1:W),
             south=apply(S2_region[,c("CUMBERLAND","SALEM","GLOUCESTER","CAMDEN","BURLINGTON")],1,sum),
             shore=apply(S2_region[,c("CAPE.MAY","ATLANTIC","OCEAN","MONMOUTH")],1,sum),
             central=apply(S2_region[,c("HUNTERDON","SOMERSET","UNION","MERCER","MIDDLESEX")],1,sum),
             north=apply(S2_region[,c("BERGEN","HUDSON","ESSEX","MORRIS","PASSAIC","WARREN","SUSSEX")],1,sum))
S2_region <- S2_region %>% pivot_longer(cols=-"date", names_to="region", values_to="S2")
S2_region$region <- factor(S2_region$region, levels=unique(S2_region$region))


P_region <-
  left_join(P1_region,P2_region) %>% 
  mutate(P=P1+P2)
S_region <-
  left_join(S1_region,S2_region) %>% 
  mutate(S=S1+S2)


# total testing volume for PCR and serology tests by NJ region,
# part a of testing summary summary plot
P_region %>%
  right_join(S_region) %>%
  group_by(region) %>% 
  summarise(P=sum(P),
            S=sum(S)) %>%
  ungroup() %>%
  right_join(pop_by_region) %>%
  pivot_longer(cols=-region, names_to="type", values_to="tests") %>%
  mutate(type=factor(type, levels=c("P","S","pop")),
         region=factor(region, levels=c("north","central","shore","south"))) %>%
  ggplot() +
  geom_col(aes(x=region, y=tests, fill=type, color=type, linetype=type), position="dodge",
           width=.82) +
  scale_fill_manual(values=c("red3","chartreuse4","gray97"),
                    labels=c("PCR tests","serology tests","")) +
  scale_color_manual(values=c("red3","chartreuse4","black"),
                     guide="none") +
  scale_linetype_manual(values=c("solid","solid","dashed"),
                        guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,.05)),labels=fancy_scientific,
                     sec.axis = sec_axis(~ ., name = "population",labels=fancy_scientific)) +
  ylab("total testing volume") +
  theme(panel.background=element_blank(),
        axis.line=element_line(size=.4),
        legend.title = element_blank(),
        legend.key.size = unit(17,"pt"),
        legend.text = element_text(size=10),
        legend.position="bottom",
        axis.line.y.right=element_line(linetype="longdash"))



# PCR and serology testing volume over time by region,
# part b of testing volume summary plot
p1 <-
  data.frame(date=get_week_start_from_test_week(1:(W-1)),
             PCR=P[1:(W-1)],
             serology=S[1:(W-1)]) %>%
    pivot_longer(cols=-date, names_to="type", values_to="tests") %>% 
    ggplot() +
    geom_line(aes(x=date,y=tests,group=type,color=type),size=1) +
    scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=expansion(c(0,.01))) +
    scale_y_continuous(name="weekly testing volume",
                       expand=expansion(c(.005,.03)),
                       position="right",
                       labels=fancy_scientific) +
    scale_color_manual(name="test type",values=c("red3","chartreuse4")) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          panel.border=element_blank(),
          panel.background=element_blank(),
          axis.line=element_line(size=.4),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none",
          plot.margin = unit(c(2.5,5.5,-2.3,5.5), "pt"))

p2 <- 
  P_region %>% 
    group_by(region) %>%
    summarise(date=date,P=P/max(P)) %>% 
    filter(date < get_week_start_from_test_week(W)) %>% 
    ggplot() +
    geom_tile(aes(x=date, y=region, fill=P, color=P), size=1) +
    scale_fill_gradientn(colors = c("white","red3"),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    scale_color_gradientn(colors = c("white","red3"),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    scale_x_date(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ylab("") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x=element_line(size=.3),
          axis.line.y=element_blank(),
          plot.margin = unit(c(0,5.5,-2.3,5.5), "pt"),
          legend.key.height = unit(11, "pt"),
          legend.key.width = unit(14, "pt"),
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.position="none")

p3 <-
  S_region %>% 
    group_by(region) %>%
    summarise(date=date,S=S/max(S)) %>% 
    filter(date < get_week_start_from_test_week(W)) %>% 
    ggplot() +
    geom_tile(aes(x=date, y=region, fill=S, color=S), size=1) +
    scale_fill_gradientn(colors = c("white","chartreuse4"),
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    scale_color_gradientn(colors = c("white","chartreuse4"),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
    scale_x_date(date_labels="%b '%y", date_breaks="1 month", expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    ylab("region") +
    theme(axis.text.x=element_text(angle=90),
          #axis.ticks.x=element_blank(),
          #axis.title.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.x=element_line(size=.4),
          axis.line.y=element_blank(),
          plot.margin = unit(c(0,5.5,-2.3,5.5), "pt"),
          legend.key.height = unit(11, "pt"),
          legend.key.width = unit(14, "pt"),
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.position="none")

aligned <- align_plots(p2, p3, align = "vh", axis=c("lrtb"))
plot_grid(p1, aligned[[1]], aligned[[2]], ncol=1, align="v", axis=c("lr"), rel_heights=c(1,1,1))




# proportion of serology tests each week with past PCR positive by region,
# part c of testing volume summary plot
prop_prev_pcr_all <-
  S_region %>% 
  group_by(date) %>%
  summarise(S=c(S, sum(S)),
            S1=c(S1, sum(S1)),
            region=factor(c(as.character(region),"total"),
                          levels=c("south","shore","central","north","total"))) %>%
  ungroup() %>% 
  mutate(prop_prev_pcr=ifelse(S==0, NA, S1/S)) %>% 
  select(-c(S,S1))

prop_prev_pcr_all %>%
  filter(date < get_week_start_from_test_week(W)) %>% 
  ggplot() +
  geom_tile(aes(x=date, y=region, fill=prop_prev_pcr, color=prop_prev_pcr), size=1) +
  scale_fill_gradientn(colors = rainbow_gradient, na.value = "gray80",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                       name = stringr::str_wrap("proportion of serology tests with past PCR positive",
                                                width = 20)) +
  scale_color_gradientn(colors = rainbow_gradient, na.value = "gray80", guide="none") +
  scale_x_date(date_labels="%b '%y",date_breaks="1 month",expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.x=element_text(angle = 90),
        axis.ticks.y=element_blank(),
        axis.line.x=element_line(color="black", size=.3),
        axis.line.y=element_blank(),
        legend.title=element_text(size=10),
        legend.position="right")











