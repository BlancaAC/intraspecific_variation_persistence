
### Calculate Omega constraining growth rates and using the random sampling method


pb = txtProgressBar(min = 0, max = length(A.list), initial = 0) 

omega.yes.list <- vector("list", length(A.list))
for (i in seq_along(A.list)){
  
  M <- A.list[[i]]
  animal_abund <- abpol.prop 
  n.pol <- length(abpol.prop) 
  n.plant.ind <- nrow(M) - n.pol
  n.plant.sp <- 1
  M <- unname(M)
  
  # these help us determine the type of species, which species each row/column belongs to, and whether some rows/columns correspond to within-species subsets
  types <- c(rep('P', n.plant.ind), rep('A', ncol(M)-n.plant.ind)) 
  species <- c(rep(1, n.plant.ind), seq(2, n.pol +n.plant.sp)) 
  subsets <- c(seq(1, n.plant.ind), rep(1, n.pol))

  # modified M matrix that also includes the constraints on pollinator relative abundances
  M_and_constraints <- Matrix::Matrix(M, sparse=TRUE)
  
  M_and_constraints[which(types=='A')[-1],] <- 0
  M_and_constraints[which(types=='A')[-1],which(types=='A')[1]] <- -animal_abund[-1]
  diag(M_and_constraints)[which(types=='A')[-1]] <- animal_abund[1]

  # we will calculate structural stability via randomly sampled growth rate vectors
  # the number below determines how many samples we draw; the smaller the feasibility domain, the more we need in order to get a good estimate of Omega
  n_samples <- 1E5
  
  # randomly sampled growth rate vector on the unit sphere
  # NOTE: the argument N needs to be 1 less than the required number of growth rates
  growth_rates_yes_constraint <- as.matrix(
    pracma::rands(
      n_samples,
      N=length(unique(species[which(types=='P')]))
    )
  )
  colnames(growth_rates_yes_constraint) <- c(
    paste0('r_P',unique(species[which(types=='P')])),
    paste0('r_A',species[which(types=='A')[1]])
  )
  
  # add constraint columns for the abundances of additional pollinators
  growth_rates_yes_constraint <- Matrix::Matrix(
    cbind(
      growth_rates_yes_constraint,
      matrix(0, n_samples, sum(types=='A')-1)
    ),
    sparse = TRUE
  )
  colnames(growth_rates_yes_constraint)[(length(unique(species[which(types=='P')]))+2):ncol(growth_rates_yes_constraint)] <- paste0('constraint_A',species[which(types=='A')[-1]])
  
  # create a matrix that maps values in the growth_rates_yes_contraint to the M matrix
  Z_mat <- Matrix::Matrix(
    0,
    nrow = ncol(growth_rates_yes_constraint),
    ncol = nrow(M),
    sparse = TRUE
  )
  for(m in 1:ncol(Z_mat)){
    Z_mat[species[m],m] <- 1
  }
  
  # solve for the equilibria
  N_equil <- growth_rates_yes_constraint %*% Z_mat %*% t(solve(M_and_constraints))
  
  # for each sampled growth rate vector determine whether the equilibrium is feasible
  feasible <- colSums(t(N_equil>0)) == ncol(N_equil)
  
  # big Omega is the proportion of vectors with feasible equilibria
  Omega_yes_constraint <- sum(feasible) / n_samples

  # calculate the implied growth rates across all plant individuals and all pollinator species for feasible points
  # NOTE: only the first pollinator respects the "surface of the hypersphere" constraint
  #growth_rates_feasible <- N_equil[feasible,] %*% t(M)
  
  #growth_rates_yes_constraint[feasible, ][,2] # feasible range of growth rates for pollinators
  
  #identify if matrices come from the randomized scenario or the increasing generalization scenario
  scenario <- ifelse(i > length(A.list) - (nrow(int.mat)-1), "gen", "random")
  
  if(per.plant.ind==TRUE){
  omega <- data.frame(N_ind= n.plant.ind, r_constraint="yes",
                      #scenario=scenario,
                      plant_id=rownames(A.list[[i]])[1],
                      min_r_pol=min(growth_rates_yes_constraint[feasible, ][,2]),
                      max_r_pol=max(growth_rates_yes_constraint[feasible, ][,2]),
                      min_r_plant=min(growth_rates_yes_constraint[feasible, ][,1]),
                      max_r_plant=max(growth_rates_yes_constraint[feasible, ][,1]),
                      omega= Omega_yes_constraint)
  }
  else{
    omega <- data.frame(N_ind= n.plant.ind, r_constraint="yes",
                      scenario=scenario,
                      #plant_id=rownames(A.list[[i]])[1],
                      min_r_pol=min(growth_rates_yes_constraint[feasible, ][,2]),
                      max_r_pol=max(growth_rates_yes_constraint[feasible, ][,2]),
                      min_r_plant=min(growth_rates_yes_constraint[feasible, ][,1]),
                      max_r_plant=max(growth_rates_yes_constraint[feasible, ][,1]),
                      omega= Omega_yes_constraint)
    }
  
  
  omega.yes.list[[i]] <- omega
  setTxtProgressBar(pb,i)
}


omega.yes.all <- bind_rows(omega.yes.list) %>% mutate(range_r_pol=max_r_pol-min_r_pol,
                                                      range_r_plant=max_r_plant-min_r_plant)


plot(omega.yes.all$range_r_plant, omega.yes.all$range_r_pol)

#set a different color for each plant species
sp.color <- switch(focal,
                   "CLIB" = "#5DA6A7",
                   "HCOM" = "#F2C57C",
                   "HHAL" = "#E8875A")
sp.label <- switch(focal,
                   "CLIB" = ~italic("Cistus libanotis"),
                   "HCOM" = ~italic("Halimium calycinum"),
                   "HHAL" = ~italic("Halimium halimifolium"))


if(per.plant.ind==FALSE){

# plot results per species
omega_random <- subset(omega.yes.all, scenario == "random") 

omega_random_sum <- as.data.frame(
  as.matrix(aggregate(cbind(omega, range_r_pol, range_r_plant) ~ N_ind + scenario, 
                      omega_random, function(x) c(#mean=mean(x, na.rm=T), 
                                                  #sd=sd(x, na.rm=T), 
                                                  ci=CI(x, ci = 0.95)))))


omega_gen <- subset(omega.yes.all, scenario == "gen") %>% 
  dplyr::rename(omega.ci.mean=omega,
                range_r_pol.ci.mean=range_r_pol,
                range_r_plant.ci.mean=range_r_plant) %>% 
  dplyr::select(-r_constraint)

omega.yes.all.sum <- rbind.fill(omega_random_sum, omega_gen)


# include degree information
k.omega <- omega_gen %>% mutate(degree= int.mat.d[2:nrow(int.mat.d),]$degree) %>%
  mutate(degree = na_if(degree, lag(degree))) %>% na.omit()

# plotting results
if (focal == "CLIB") {
  k.omega %<>% 
    filter(degree %in% c(2, 4, 6, 8, 10, 15))
}

if (focal == "HCOM") {
  k.omega %<>% 
    filter(degree %in% c(2, 4, 6, 8))
}

if (focal == "HHAL") {
  k.omega %<>% 
    filter(degree %in% c(2, 4, 5, 6, 9, 14))
}

glimpse(omega.yes.all.sum)
omega.yes.all.sum <- omega.yes.all.sum[,1:11]
omega.yes.all.sum %<>% mutate(scenario=as.factor(scenario),
                              N_ind=as.integer(N_ind)) %>% mutate_if(is.character,as.numeric) 


# plot size of the feasibility domain
omega.plot <- ggplot() +
  geom_line(data=omega.yes.all.sum, aes(x=N_ind, y=omega.ci.mean, color=scenario)) +
  geom_ribbon(data=omega.yes.all.sum, aes(x=N_ind, y=omega.ci.mean, 
                                          ymin=omega.ci.lower,
                                          ymax=omega.ci.upper, fill=scenario), alpha=0.4) +
  scale_color_manual(values=c(sp.color, "grey70")) + 
  scale_fill_manual(values=c(sp.color, "grey70")) +
  geom_text(parse=T, data=k.omega, color=sp.color, size=4,
            aes(x=N_ind, y=omega.ci.mean, label = paste0("  italic(k)==", degree)), hjust = 0.7,  vjust = -1) +
  geom_point(data=k.omega, color=sp.color, aes(x=N_ind, y=omega.ci.mean)) +
  theme_bw() + ylab("Feasibility \ndomain size") + xlab("Plant population size \n(i.e., number of individuals)") +
  theme(legend.position = "none", text=element_text(size=16),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -0.75),
        plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_log10()


# plot range of feasible pollinator growth rates
r.pol.plot <- ggplot() +
  geom_line(data=omega.yes.all.sum, aes(x=N_ind, y=range_r_pol.ci.mean, color=scenario)) +
  geom_ribbon(data=omega.yes.all.sum, aes(x=N_ind, y=range_r_pol.ci.mean, 
                                          ymin=range_r_pol.ci.lower,
                                          ymax=range_r_pol.ci.upper, fill=scenario), alpha=0.4) +
  scale_color_manual(values=c(sp.color, "grey70")) + 
  scale_fill_manual(values=c(sp.color, "grey70")) +
  geom_text(parse=T, data=k.omega, color=sp.color, size=4,
            aes(x=N_ind, y=range_r_pol.ci.mean, label = paste0("  italic(k)==", degree)), hjust = 0.7,  vjust = -1) +
  geom_point(data=k.omega, color=sp.color, aes(x=N_ind, y=range_r_pol.ci.mean)) +
  theme_bw() + ylab("Range of pollinator \ngrowth rates") + xlab("Plant population size \n(i.e., number of individuals)") +
  theme(legend.position = "none", text=element_text(size=16),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -0.75),
        plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_log10()


# plot range of feasible plant growth rates
r.plant.plot <- ggplot() +
  geom_line(data=omega.yes.all.sum, aes(x=N_ind, y=range_r_plant.ci.mean, color=scenario)) +
  geom_ribbon(data=omega.yes.all.sum, aes(x=N_ind, y=range_r_plant.ci.mean, 
                                          ymin=range_r_plant.ci.lower,
                                          ymax=range_r_plant.ci.upper, fill=scenario), alpha=0.4) +
  scale_color_manual(values=c(sp.color, "grey70")) + 
  scale_fill_manual(values=c(sp.color, "grey70")) +
  geom_text(parse=T, data=k.omega, color=sp.color, size=4,
            aes(x=N_ind, y=range_r_plant.ci.mean, label = paste0("  italic(k)==", degree)), hjust = 0.7,  vjust = -1) +
  geom_point(data=k.omega, color=sp.color, aes(x=N_ind, y=range_r_plant.ci.mean)) +
  theme_bw() + ylab("Range of plant \ngrowth rates") + xlab("Plant population size \n(i.e., number of individuals)") +
  theme(legend.position = "none", text=element_text(size=16),
        axis.title.y = element_text(vjust = +3),
        axis.title.x = element_text(vjust = -0.75),
        plot.margin = unit(c(1,1,1,1), "cm")) + scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_log10()

if (focal == "CLIB") {
plots <- (omega.plot + theme(axis.title.x = element_blank(),
                             plot.title = element_text(size=16, hjust = 0.5),
                             plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")) +
            ggtitle(sp.label)) / 
  (r.pol.plot + theme(axis.title.x = element_blank(),
                      plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))) / 
  (r.plant.plot + theme(axis.title.x = element_blank(),
                        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines")))  + 
  plot_layout(axis_titles="collect") + plot_annotation(tag_levels = 'A')
}

if (focal == "HCOM") {
  plots <- (omega.plot + theme(axis.title.y = element_blank(),
                               plot.title = element_text(size=16, hjust = 0.5),
                               plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))+
              ggtitle(sp.label)) / 
              (r.pol.plot + theme(axis.title.y = element_blank(),
                                  plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))) / 
              (r.plant.plot + theme(axis.title.y = element_blank(),
                                    plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))) + 
    plot_layout(axis_titles="collect")
}

if (focal == "HHAL") {
  plots <- (omega.plot + theme(axis.title.y = element_blank(),
                               axis.title.x = element_blank(),
                               plot.title = element_text(size=16, hjust = 0.5),
                               plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))+
              ggtitle(sp.label)) / 
    (r.pol.plot + theme(axis.title.y = element_blank(),
                        axis.title.x = element_blank(),
                        plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))) / 
    (r.plant.plot + theme(axis.title.y = element_blank(),
                          axis.title.x = element_blank(),
                          plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"))) + 
    plot_layout(axis_titles="collect")
}


}




if(per.plant.ind==TRUE){

# correlation between contribution to feasibility and plant attributes

plant.attr <- read.csv(paste0("Output/fitness_components_", focal, ".csv")) %>% 
  mutate(X=as.factor(X)) %>% dplyr::rename(plant_id=X) %>% arrange(degree)

omega.yes.all <- omega.yes.all[-nrow(omega.yes.all),]
mean.omega <- mean(omega.yes.all$omega)
averaged.ind <- omega.yes.all[nrow(omega.yes.all),]$omega

contr.plot <- ggplot(omega.yes.all, aes(x=omega)) +
  geom_histogram(aes(y=..count..), position="identity", 
                 alpha=0.5, color=sp.color, fill=sp.color) +
  theme_bw() + theme(text=element_text(size = 16),
                     axis.title.x = element_text(
                       margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     plot.title = element_text(size=16, hjust = 0.5)) +
  ylab("Count") + xlab("Feasibility domain size") +
  ggtitle(sp.label) +
  geom_vline(xintercept=mean.omega, linetype="dashed") 
  #geom_vline(xintercept=averaged.ind, linetype="solid")


# for each plant individual

k.omega <- omega.yes.all %>% merge(plant.attr, by="plant_id")

degree.plot <- ggplot(k.omega, aes(x = degree, y = omega)) +
  geom_point(color=sp.color, alpha=0.9) + theme_bw() +
  theme(text=element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  scale_x_log10() + 
  ylab("Feasibility domain size") + xlab("Degree")

k.omega_long <- k.omega %>%
  pivot_longer(cols = c(degree, N_flowers, sum_pol, N_seeds_per_visit_flower), 
               names_to = "predictor", values_to = "pred_value") %>%
  pivot_longer(cols = c(omega, range_r_plant, range_r_pol), 
               names_to = "response", values_to = "resp_value")

# label names for predictor variable
p.labs <- c("Plant degree", "Number of flowers", "Pollinator contribution", "Visitation rate")
names(p.labs) <- c("degree", "N_flowers", "N_seeds_per_visit_flower", "sum_pol")

# label names for predictor variable
r.labs <- c("Feasibility \ndomain size", "Range of plant \ngrowth rates", "Range of pollinator \ngrowth rates")
names(r.labs) <- c("omega", "range_r_plant", "range_r_pol")

plots <- ggplot(k.omega_long, aes(x = pred_value, y = resp_value)) +
  geom_point(color=sp.color, alpha=0.7) +
  facet_grid(response ~ predictor, scales = "free", labeller = labeller(predictor = p.labs, response = r.labs)) +       
  theme_bw() +
  theme(strip.text = element_text(face="bold"),
        text=element_text(size = 13), 
        strip.background = element_blank(),
        axis.title=element_blank()) +
  scale_y_log10() + scale_x_log10()

}

