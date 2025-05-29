
# Perform the feasibility estimation for each focal species
sp.names <- c("CLIB", "HCOM", "HHAL")

# Population-level results (in the absence of heterospecific plant competitors)
per.study.plot <- FALSE
per.plant.ind <- TRUE 
# TRUE when estimating the feasibility of systems composed of a single plant individual type
# FALSE when estimating the feasibility of systems composed of multiple plant individual types
omega.yes.all.list <- list() # save results output for each species
plot.sp.list <- list() # plots for each species
plot.contr.list <- list()
plot.degree.list <- list()

for(name in sp.names){
  focal <- name
  
  # construct interaction matrices
  source('Code/interaction_matrices.R')
  
  # calculate Omega constraining growth rates and using a random sampling approach
  source('Code/population_level_r_constraint_sampling.R')
  
  # save results and plots
  plot.sp.list[[name]] <- plots
  omega.yes.all.list[[name]] <- omega.yes.all
  
  if(per.plant.ind==TRUE){
  plot.contr.list[[name]] <- contr.plot
  plot.degree.list[[name]] <- degree.plot
  }

}

patchwork::wrap_plots(plot.sp.list, nrow = 1, ncol = 3)

if(per.plant.ind==TRUE){
  plot1 <- patchwork::wrap_plots(plot.contr.list, nrow = 1, ncol = 3) + 
    patchwork::plot_layout(axis_titles = "collect")
  plot2 <- patchwork::wrap_plots(plot.degree.list, nrow = 1, ncol = 3) + 
    patchwork::plot_layout(axis_titles = "collect")
  plot1 / plot2
}
  

# Community-level results (with heterospecific plant competitors)
per.study.plot <- FALSE
per.plant.ind <- FALSE
plot.comm.list <- list() # plot for each species
vul.sp.all.list <- list()
metric.all.list <- list()
omega.all.list <- list()

for(name in sp.names){
  focal <- name
  
  # construct interaction matrices per plant species
  source('Code/interaction_matrices.R')
  
  # construct interaction matrix for the entire community and calculate 'vulnerability' of plant species
  source('Code/interaction_matrix_comm.R')

  plot.comm.list[[name]] <- plot.comm
  vul.sp.all.list[[name]] <- vul.sp.all
  metric.all.list[[name]] <- metric.all
  omega.all.list[[name]] <- omega.all
  
}

patchwork::wrap_plots(plot.comm.list, nrow = 1, ncol = 3) + plot_layout(axes="collect_y")


