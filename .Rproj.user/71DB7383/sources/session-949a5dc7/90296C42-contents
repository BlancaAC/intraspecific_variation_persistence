
library(Matrix)
library(pracma)
library(scales)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(patchwork)

get_feasible_growth_rates <- function(M, species, types, subsets, animal_abund, n_samples=1E6, n_minimum = 1E3, sign_constraint=FALSE){
	# modified M matrix that also includes the constraints on pollinator relative abundances
  M_and_constraints <- Matrix::Matrix(M, sparse=TRUE)
  M_and_constraints[which(types=='A')[-1],] <- 0
  M_and_constraints[which(types=='A')[-1],which(types=='A')[1]] <- -animal_abund[-1]
  diag(M_and_constraints)[which(types=='A')[-1]] <- animal_abund[1]
  
  # create a matrix that maps values in the growth_rates_yes_contraint to the M matrix
  Z_mat <- Matrix::Matrix(
    0,
    nrow = (length(animal_abund) + 1),
    ncol = nrow(M),
    sparse = TRUE
  )
  for(m in 1:ncol(Z_mat)){
    Z_mat[species[m],m] <- 1
  }

  all_growth_rates <- matrix(0,0,ncol(M))
  colnames(all_growth_rates) <- paste0("r_",types,species,'_',subsets)

  attempts <- 0
	repeat {
		attempts <- attempts + 1
		# message('trying ',nrow(all_growth_rates))
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
  
	  # solve for the equilibria
	  N_equil <- growth_rates_yes_constraint %*% Z_mat %*% t(solve(M_and_constraints))
	  
	  # for each sampled growth rate vector determine whether the equilibrium is feasible
	  feasible <- colSums(t(N_equil>0)) == ncol(N_equil)

		# recalculate the full vectors of growth rates
		growth_rates <- t(M %*% t(N_equil[feasible,]))
		colnames(growth_rates) <- paste0("r_",types,species,'_',subsets)

		if(sign_constraint){
			meets_sign_constraint <- apply(
				t(growth_rates),
				2,
				function(r, types){
					all(r[types=='P'] < 0) & all(r[types=='A'] > 0)
				},
				types = types
			)
			growth_rates <- growth_rates[meets_sign_constraint,]
		}

		all_growth_rates <- rbind(all_growth_rates, growth_rates)

		if(nrow(all_growth_rates) > n_minimum || attempts > 100) break
	}

	all_growth_rates <- all_growth_rates[order(all_growth_rates[,1]),]
	
	return(list(gr = all_growth_rates, Omega = nrow(all_growth_rates) / (attempts * n_samples)))
}

# just to make sure everyone gets the same results
set.seed(912342)

# indicators for the different components
types <- c('P','P','A','A')
species <- c(1,1,2,3)
subsets <- c(1,2,1,1)

# manually build a mutualistic interaction matrix with:
# two plant species (each with two subset types)
# and two animal pollinator species
M <- diag(4)

# intraspecific coefficients (between two types from the same plant species)
# IMPORTANT NOTE: this is a change that we require in the methods; intraspecific competition must be less than within-subset competition
subset_inter <- 0.9
for(i in which(types == 'P')){
	for(j in which(types == 'P')){
		if(i != j && species[i] == species[j]){
			M[i,j] <- subset_inter
		}
	}
}

# interspecific plant competition
plant_inter <- 0.1
for(i in which(types == 'P')){
	for(j in which(types == 'P')){
		if(species[i] != species[j]){
			M[i,j] <- M[j,i] <- plant_inter
		}
	}
}

# interspecific animal competition
animal_inter <- 0.1
for(i in which(types == 'A')){
	for(j in which(types == 'A')){
		if(species[i] != species[j]){
			M[i,j] <- M[j,i] <- animal_inter
		}
	}
}

# benefits to animals from plants
M[which(types=='A'),which(types=='P')] <- -0.2

# observed pollinator abundances
# NOTE: I again make these random here but for the real data we need the site-specific total visitation rates per pollinator (right?)
animal_abund <- c(0.35, 0.65)
  #runif(sum(types=='A'))
names(animal_abund) <- which(types=='A')

# container to shove the different example networks into
A.list <- list()

# benefits to plants from animals

# two specialized individuals (specialized on different pollinators)
A.list[[1]] <- M
A.list[[1]][which(types=='P'),which(types=='A')] <- c(-0.4,0,0,-0.2)

# one generalized and one specialized individual
A.list[[2]] <- M
pi1 <- c(0.99,0.01)
pi2 <- c(0.0,1.0)
A.list[[2]][which(types=='P'),which(types=='A')] <- -matrix(c(0.4,0.4,0.2,0.2),2,2) * matrix(c(pi1,pi2),2,2)

# one generalized and one specialized individual
A.list[[3]] <- M
pi1 <- c(1.0,0.0)
pi2 <- c(0.01,0.99)
A.list[[3]][which(types=='P'),which(types=='A')] <- -matrix(c(0.4,0.4,0.2,0.2),2,2) * matrix(c(pi1,pi2),2,2)

# two generalized individuals (with not perfect but close visitation overlap)
A.list[[4]] <- M
pi1 <- c(0.5,0.5)
pi2 <- c(0.5,0.5)
A.list[[4]][which(types=='P'),which(types=='A')] <- -matrix(c(0.4,0.4,0.2,0.2),2,2) * matrix(c(pi1,pi2),2,2)


print(A.list)

# # two generalized individuals (and thus no intraspecific variation)
# A.list[[4]] <- M
# A.list[[4]][which(types=='P'),which(types=='A')] <- c(-0.4,-0.4,-0.2,-0.2)

# # one specialized and one generalized individual (cannot be feasible?)
# A.list[[5]] <- M
# A.list[[5]][which(types=='P'),which(types=='A')] <- c(-0.4,-0.4,-0.0,-0.4)

# # one specialized and one generalized individual (cannot be feasible?)
# A.list[[6]] <- M
# A.list[[6]][which(types=='P'),which(types=='A')] <- c(-0.0,-0.8,-0.2,-0.2)

# # run the growth rate function for each input matrix
# gr.list <- lapply(
# 	A.list,
# 	get_feasible_growth_rates,
# 	species = species,
# 	types = types,
# 	subsets = subsets,
# 	animal_abund = animal_abund,
# 	n_samples = 1E3,
# 	n_minimum = 1E3,
# 	sign_constraint = TRUE
# )

# generate a horrid plot of the feasibility domains
par(mfrow=c(1,4))

lapply(A.list,
	function(x){
		# run the growth rate function for each input matrix
		together_gr <- get_feasible_growth_rates(
			x,
			species = species,
			types = types,
			subsets = subsets,
			animal_abund = animal_abund,
			n_samples = 1E3,
			n_minimum = 1E3,
			sign_constraint = FALSE
		)
		plant1_gr <- get_feasible_growth_rates(
			x[c(1,3,4),c(1,3,4)],
			species = species[c(1,3,4)],
			types = types[c(1,3,4)],
			subsets = subsets[c(1,3,4)],
			animal_abund = animal_abund,
			n_samples = 1E3,
			n_minimum = 1E3,
			sign_constraint = FALSE
		)
		plant2_gr <- get_feasible_growth_rates(
			x[c(2,3,4),c(2,3,4)],
			species = species[c(2,3,4)],
			types = types[c(2,3,4)],
			subsets = subsets[c(2,3,4)],
			animal_abund = animal_abund,
			n_samples = 1E3,
			n_minimum = 1E3,
			sign_constraint = FALSE
		)
		together_gr$gr <- rbind(0,together_gr$gr,0)
		plot(
			together_gr$gr[,"r_P1_1"],
			together_gr$gr[,'r_A2_1'],
			type='l',
			xlim = c(-1,1),
			ylim = c(-1,1),
			xlab = 'Plant growth rate',
			ylab = 'Pollinator growth rate'
		)
		plant1_gr$gr <- rbind(0,plant1_gr$gr,0)
		polygon(
			plant1_gr$gr[,"r_P1_1"],
			plant1_gr$gr[,'r_A2_1'],
			col = alpha('blue',0.1),
			border = NA
		)
		plant2_gr$gr <- rbind(0,plant2_gr$gr,0)
		polygon(
			plant2_gr$gr[,"r_P1_2"],
			plant2_gr$gr[,'r_A2_1'],
			col = alpha('red',0.1),
			border = NA
		)
		mtext(paste0('Omega = ', signif(together_gr$Omega, digits=5)))

	}
)


x <- A.list[[1]]
together_gr <- get_feasible_growth_rates(
  x,
  species = species,
  types = types,
  subsets = subsets,
  animal_abund = animal_abund,
  n_samples = 1E3,
  n_minimum = 1E3,
  sign_constraint = TRUE
)
plant1_gr <- get_feasible_growth_rates(
  x[c(1,3,4),c(1,3,4)],
  species = species[c(1,3,4)],
  types = types[c(1,3,4)],
  subsets = subsets[c(1,3,4)],
  animal_abund = animal_abund,
  n_samples = 1E3,
  n_minimum = 1E3,
  sign_constraint = TRUE
)
plant2_gr <- get_feasible_growth_rates(
  x[c(2,3,4),c(2,3,4)],
  species = species[c(2,3,4)],
  types = types[c(2,3,4)],
  subsets = subsets[c(2,3,4)],
  animal_abund = animal_abund,
  n_samples = 1E3,
  n_minimum = 1E3,
  sign_constraint = TRUE
)
together_gr$gr <- rbind(0,together_gr$gr,0)
plot(
  together_gr$gr[,"r_P1_1"],
  together_gr$gr[,'r_A2_1'],
  type='l',
  xlim = c(-1,1),
  ylim = c(-1,1),
  xlab = 'Plant growth rate',
  ylab = 'Pollinator growth rate'
)
plant1_gr$gr <- rbind(0,plant1_gr$gr,0)
polygon(
  plant1_gr$gr[,"r_P1_1"],
  plant1_gr$gr[,'r_A2_1'],
  col = alpha('blue',0.1),
  border = NA
)
plant2_gr$gr <- rbind(0,plant2_gr$gr,0)
polygon(
  plant2_gr$gr[,"r_P1_2"],
  plant2_gr$gr[,'r_A2_1'],
  col = alpha('red',0.1),
  border = NA
)


plot.list <- lapply(A.list,
       function(x){
         # run the growth rate function for each input matrix
         together_gr <- get_feasible_growth_rates(
           x,
           species = species,
           types = types,
           subsets = subsets,
           animal_abund = animal_abund,
           n_samples = 1E3,
           n_minimum = 1E3,
           sign_constraint = FALSE
         )
         plant1_gr <- get_feasible_growth_rates(
           x[c(1,3,4),c(1,3,4)],
           species = species[c(1,3,4)],
           types = types[c(1,3,4)],
           subsets = subsets[c(1,3,4)],
           animal_abund = animal_abund,
           n_samples = 1E3,
           n_minimum = 1E3,
           sign_constraint = FALSE
         )
         plant2_gr <- get_feasible_growth_rates(
           x[c(2,3,4),c(2,3,4)],
           species = species[c(2,3,4)],
           types = types[c(2,3,4)],
           subsets = subsets[c(2,3,4)],
           animal_abund = animal_abund,
           n_samples = 1E3,
           n_minimum = 1E3,
           sign_constraint = FALSE
         )
         
         together_gr$gr <- rbind(0, together_gr$gr, 0)
         plant1_gr$gr <- rbind(0, plant1_gr$gr, 0)
         plant2_gr$gr <- rbind(0, plant2_gr$gr, 0)
         
         together_df <- data.frame(Plant_growth = together_gr$gr[, "r_P1_1"],
                                   Pollinator_growth = together_gr$gr[, "r_A2_1"],
                                   Plant_type = "Together")
         
         plant1_df <- data.frame(Plant_growth = plant1_gr$gr[, "r_P1_1"],
                                 Pollinator_growth = plant1_gr$gr[, "r_A2_1"],
                                 Plant_type = "1")
         
         plant2_df <- data.frame(Plant_growth = plant2_gr$gr[, "r_P1_2"],
                                 Pollinator_growth = plant2_gr$gr[, "r_A2_1"],
                                 Plant_type = "2")
         
         
         polygon_data <- bind_rows(
           plant1_df,
           plant2_df
         )
         
         p <- ggplot() + 
           geom_polygon(data = polygon_data,
                        aes(x = Plant_growth, y = Pollinator_growth,
                            group = Plant_type, fill = Plant_type),
                        alpha = 0.8, color = NA) + 
           geom_polygon_pattern(data = together_df,
                                aes(x = Plant_growth, y = Pollinator_growth),
                                pattern = "stripe",
                                pattern_fill = "#F2C57C",
                                pattern_colour = NA,
                                pattern_density = 0.4,
                                pattern_angle = 45,
                                pattern_spacing=0.04,
                                fill = "#5DA6A7",
                                color = NA, # Remove contour
                                size = 0) +
           #scale_x_continuous(limits = c(-0.5, 1)) +
           scale_y_continuous(limits = c(-0.25, 1), breaks=seq(0,1, by = 0.5)) +
           labs(x = "Plant growth rate", y = "Pollinator growth rate") +
           scale_fill_manual(values = c("#5DA6A7", "#F2C57C")) +
           theme_bw() + theme(axis.title=element_text(size=16), 
                              axis.text=element_text(size=12),
                              plot.title = element_text(size=12),
                              axis.title.y.right = element_text(
                                margin = margin(t = 0, r = 0, b = 0, l = 10)),
                              axis.title.x = element_text(
                                margin = margin(t = 10, r = 0, b = 0, l = 0)),
                              legend.position = "none") 
           return(p)
       }
)

patchwork::wrap_plots(plot.list, nrow = 4, ncol = 1) + 
  plot_layout(guides="collect", axes="collect")


### Create heatmap with different visitation rates

get_feasible_growth_rates <- function(
    M,
    species,
    types,
    subsets,
    animal_abund,
    n_samples=1E6,
    n_minimum = 1E3,
    max_attempts = 10,
    sign_constraint=FALSE
){
  # modified M matrix that also includes the constraints on pollinator relative abundances
  M_and_constraints <- Matrix::Matrix(M, sparse=TRUE)
  M_and_constraints[which(types=='A')[-1],] <- 0
  M_and_constraints[which(types=='A')[-1],which(types=='A')[1]] <- -animal_abund[-1]
  diag(M_and_constraints)[which(types=='A')[-1]] <- animal_abund[1]
  
  # create a matrix that maps values in the growth_rates_yes_contraint to the M matrix
  Z_mat <- Matrix::Matrix(
    0,
    nrow = (length(animal_abund) + 1),
    ncol = nrow(M),
    sparse = TRUE
  )
  for(m in 1:ncol(Z_mat)){
    Z_mat[species[m],m] <- 1
  }
  
  all_growth_rates <- matrix(0,0,ncol(M))
  colnames(all_growth_rates) <- paste0("r_",types,species,'_',subsets)
  
  attempts <- 0
  repeat {
    attempts <- attempts + 1
    # message('trying ',nrow(all_growth_rates))
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
    
    # solve for the equilibria
    N_equil <- growth_rates_yes_constraint %*% Z_mat %*% t(solve(M_and_constraints))
    
    # for each sampled growth rate vector determine whether the equilibrium is feasible
    feasible <- colSums(t(N_equil>0)) == ncol(N_equil)
    
    # recalculate the full vectors of growth rates
    growth_rates <- t(M %*% t(N_equil[feasible,]))
    colnames(growth_rates) <- paste0("r_",types,species,'_',subsets)
    
    if(sign_constraint){
      meets_sign_constraint <- apply(
        t(growth_rates),
        2,
        function(r, types){
          all(r[types=='P'] < 0) & all(r[types=='A'] > 0)
        },
        types = types
      )
      growth_rates <- growth_rates[meets_sign_constraint,]
    }
    
    all_growth_rates <- rbind(all_growth_rates, growth_rates)
    
    if(nrow(all_growth_rates) > n_minimum || attempts >= max_attempts) break
  }
  
  all_growth_rates <- all_growth_rates[order(all_growth_rates[,1]),]
  
  return(list(gr = all_growth_rates, Omega = nrow(all_growth_rates) / (attempts * n_samples)))
}

# just to make sure everyone gets the same results
set.seed(912342)

# indicators for the different components
types <- c('P','P','A','A')
species <- c(1,1,2,3)
subsets <- c(1,2,1,1)

# manually build a mutualistic interaction matrix with:
# two plant species (each with two subset types)
# and two animal pollinator species
M <- diag(4)

# intraspecific coefficients (between two types from the same plant species)
# IMPORTANT NOTE: this is a change that we require in the methods; intraspecific competition must be less than within-subset competition
subset_inter <- 0.9
for(i in which(types == 'P')){
  for(j in which(types == 'P')){
    if(i != j && species[i] == species[j]){
      M[i,j] <- subset_inter
    }
  }
}

# interspecific plant competition
plant_inter <- 0.1
for(i in which(types == 'P')){
  for(j in which(types == 'P')){
    if(species[i] != species[j]){
      M[i,j] <- M[j,i] <- plant_inter
    }
  }
}

# interspecific animal competition
animal_inter <- 0.1
for(i in which(types == 'A')){
  for(j in which(types == 'A')){
    if(species[i] != species[j]){
      M[i,j] <- M[j,i] <- animal_inter
    }
  }
}

# benefits to animals from plants
M[which(types=='A'),which(types=='P')] <- -0.2

# observed pollinator abundances
# NOTE: I again make these random here but for the real data we need the site-specific total visitation rates per pollinator (right?)
animal_abund <- c(0.35, 0.65)
names(animal_abund) <- which(types=='A')

# fix the beta
betas <- matrix(
  c(
    rep(0.4,2),
    rep(0.2,2)
  ),
  2,
  2
)

# adjust the lengths of these vectors to make the grid finer or not
prop1 <- seq(0.0,1.0,length.out = 25)
prop2 <- seq(0.0,1.0,length.out = 25)

omega <- matrix(NA,length(prop1),length(prop2))

for(i in 1:nrow(omega)){
  for(j in 1:ncol(omega)){
    message(i,' ',j)
    A <- M
    pi1 <- c(prop1[i],1-prop1[i])
    pi2 <- c(prop2[j],1-prop2[j])
    visits <- matrix(c(pi1,pi2),2,2)
    A[which(types=='P'),which(types=='A')] <- - betas * visits
    together_gr <- get_feasible_growth_rates(
      A,
      species = species,
      types = types,
      subsets = subsets,
      animal_abund = animal_abund,
      n_samples = 1E4,
      n_minimum = 10,
      max_attempts = 5,
      sign_constraint = FALSE
    )
    omega[i,j] <- together_gr$Omega
    
  }
}

image(
  x = prop1,
  y = prop2,
  z = omega,
  xlab = 'Proportion of visits from "better" pollinator to first plant type',
  ylab = 'Proportion of visits from "worse" pollinator to first plant type'
)
abline(1,-1)



data <- expand.grid(prop1 = prop1, prop2 = prop2)
data$omega <- as.vector(omega)

ggplot(data, aes(x = prop1, y = prop2, fill = omega)) +
  geom_tile(color = "gray90") + 
  scale_fill_gradient(
    low = "grey99",  
    high = "#5DA6A7", limits=c(0,0.4))  +
  labs(
    x = bquote("Proportion of visits from" ~ A[1] ~ "to P"[11]),
    y = bquote("Proportion of visits from" ~ A[2] ~ "to P"[11]),
    fill = "Feasibility \ndomain size"
  )  + 
  #geom_abline(slope = -1, intercept = 1, linetype = "dashed", color = "black") +
  theme_bw() + 
  #scale_y_continuous(position="right") + 
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 10, b = 0, l = 0)
    ),
    axis.title.x = element_text(
      margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.justification=c(1,1), legend.position=c(0.99,0.99)
  ) +
  coord_cartesian(expand = FALSE) 



