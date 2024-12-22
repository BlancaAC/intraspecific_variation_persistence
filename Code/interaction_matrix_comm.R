
scale_to <- function(x, val) x * val/mean(x)


# Read data

# matrix with beta coefficients
filenames.beta <- list.files("Output", pattern="beta_matrix_*", full.names=TRUE)
beta.mat.list <- lapply(filenames.beta, read.csv)
for (i in seq_along(beta.mat.list)){
  beta.mat.list[[i]] <- beta.mat.list[[i]] %>% column_to_rownames("X") 
}

# plant attributes
filenames.attr <- list.files("Output", pattern="plant_attr_*", full.names=TRUE)
attr.list <- lapply(filenames.attr, read.csv)
for (i in seq_along(attr.list)){
  attr.list[[i]] <- attr.list[[i]] %>% column_to_rownames("X") %>%
    dplyr::select(Plant_id, Plant_sp, N_flowers, Flowering_min) %>%
    dplyr::mutate(Prop_N_flowers=N_flowers/sum(N_flowers))
}

attr <- do.call("rbind", attr.list)

# matrix with total number of visits per flower across the season 
filenames.vis <- list.files("Output", pattern="vis_matrix_*", full.names=TRUE)
vis.mat.list <- lapply(filenames.vis, read.csv)
for (i in seq_along(vis.mat.list)){
  temp <- vis.mat.list[[i]] %>% column_to_rownames("X") 
  vis.mat.list[[i]] <- temp
}


# matrix with proportional visitation data 
prop.vis.mat.list <- list()
for (i in seq_along(vis.mat.list)){
  prop.vis.mat.list[[i]] <- vis.mat.list[[i]] %>% 
    dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0)) %>% 
    dplyr::mutate(across(where(is.numeric), ~ ./sum(.))) 
}
#colSums(prop.vis.mat.list[[1]])


# build the matrix of mutualistic benefits (gamma) with proportional visitation data multiplied by beta coefficients and flower number
# (proportion of visits by each pollinator species to each plant individual)
# I need to use this to parameterize the gamma submatrix of A matrix

int.mat.list <- list()
for (i in seq_along(beta.mat.list)){
  temp <- beta.mat.list[[i]] * prop.vis.mat.list[[i]] * attr.list[[i]]$N_flowers
  int.mat.list[[i]] <- temp %>% merge(select(attr.list[[i]], Plant_id, Plant_sp), by="row.names") %>% 
    column_to_rownames("Row.names") 
}


# I need the plant id and the plant species as columns

int.mat <- do.call("rbind.fill", int.mat.list) %>% replace(is.na(.), 0)

int.mat.clean <- int.mat  %>% 
  dplyr::mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
  column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)

# number of pollinator visits
vis.plant.id <- list()
for (j in seq_along(vis.mat.list)){
  temp <- vis.mat.list[[j]] %>% rownames_to_column("Plant_id")
  vis.plant.id[[j]] <- temp %>% merge(select(attr.list[[j]], Plant_id, Plant_sp), by="Plant_id")
}

abpol <- do.call("rbind.fill", vis.plant.id) %>% replace(is.na(.), 0)

# total abundance of pollinator across all plant species
abpol.total <- abpol %>% select(-Plant_id, -Plant_sp) %>% sum()

# proportion of pollinator species across all plant species
abpol.prop <- (colSums(select(abpol, -Plant_id, -Plant_sp))/sum(select(abpol, -Plant_id, -Plant_sp)))
#sum(abpol.prop)


###########

# Generalized-specialized ind mixtures vs specialized individuals

sp <- focal
#rand.output[[39]] # random sample of 40 plant individuals 
#gen.output[[39]] # 40 most specialized plant individuals

glimpse(int.mat)

A.list <- list()
vis.per.sp <- list()

for (k in seq_along(rand.output[[39]])){
  y <- int.mat %>% filter(Plant_sp !=sp | Plant_id %in% rownames(rand.output[[39]][[k]]))
  
  # Interaction matrix A
  
  # for plants as rows
  test <- list()
  for (i in unique(y$Plant_sp)){
    
    con <- y %>% filter(Plant_sp == i) %>% 
      mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
      column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
    
    het <- y %>% filter(Plant_sp != i) %>% 
      mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
      column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
    
    # intraspecific competition submatrix
    intra.plant <- matrix(0.9, nrow=nrow(con), ncol= nrow(con))
    diag(intra.plant) <- 1
    rownames(intra.plant) <- rownames(con)
    colnames(intra.plant) <- rownames(con)
    
    # interspecific competition submatrix
    inter.plant <- matrix(0.1, nrow=nrow(con), ncol= nrow(het))
    rownames(inter.plant) <- rownames(con)
    colnames(inter.plant) <- rownames(het)
    
    # mutualistic benefits received by plants submatrix
    gamma.plant <- y %>% filter(Plant_sp == i) %>% 
      mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
      column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
    
    test[[i]] <- cbind(intra.plant, inter.plant, -gamma.plant)
    
  }
  
  names(test) <- NULL
  plant.submat <- do.call("rbind", test)
  
  # for pollinators as rows
  # interspecific competition
  int.mat.clean <- y  %>% 
    dplyr::mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  inter.pol <- matrix(0.1, nrow=ncol(int.mat.clean), ncol=ncol(int.mat.clean))
  rownames(inter.pol) <- colnames(int.mat.clean)
  colnames(inter.pol) <- colnames(int.mat.clean)
  
  # mutualistic benefits received by pollinators submatrix
  gamma.pol <- matrix(-0.2, nrow=ncol(select(y, -Plant_id, -Plant_sp)), ncol=nrow(y))
  rownames(gamma.pol) <- colnames(int.mat.clean)
  colnames(gamma.pol) <- rownames(int.mat.clean)
  
  pol.submat <- cbind(gamma.pol, inter.pol)
  
  # A matrix
  A <- rbind(plant.submat, pol.submat) %>% as.matrix()
  
  A.list[[k]] <- A
  
  # relative proportion of visitations per plant species
  pol.plant.sp <- abpol %>% filter(Plant_sp != sp | Plant_id %in% rownames(rand.output[[39]][[k]])) %>%
    dplyr::select(-Plant_id) %>% group_by(Plant_sp) %>% summarise_all(list(sum))
  vis.per.sp[[k]] <- rowSums(pol.plant.sp[,-1])/sum(pol.plant.sp[,-1])
  vis.per.sp[[k]] <- scale_to(rowSums(pol.plant.sp[,-1]), 1)
  
}



########## most specialized

y <- int.mat %>% filter(Plant_sp != sp | Plant_id %in% rownames(gen.output[[39]]))
n.focal <- length(rownames(gen.output[[39]]))

# Interaction matrix A

# for plants as rows
test <- list()
for (i in unique(y$Plant_sp)){
  
  con <- y %>% filter(Plant_sp == i) %>% 
    mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  het <- y %>% filter(Plant_sp != i) %>% 
    mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  # intraspecific competition submatrix
  intra.plant <- matrix(0.9, nrow=nrow(con), ncol= nrow(con))
  diag(intra.plant) <- 1
  rownames(intra.plant) <- rownames(con)
  colnames(intra.plant) <- rownames(con)
  
  # interspecific competition submatrix
  inter.plant <- matrix(0.1, nrow=nrow(con), ncol= nrow(het))
  rownames(inter.plant) <- rownames(con)
  colnames(inter.plant) <- rownames(het)
  
  # mutualistic benefits received by plants submatrix
  gamma.plant <- y %>% filter(Plant_sp == i) %>% 
    mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
    column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)
  
  test[[i]] <- cbind(intra.plant, inter.plant, -gamma.plant)
  
}

names(test) <- NULL
plant.submat <- do.call("rbind", test)

# for pollinators as rows
# interspecific competition
int.mat.clean <- y  %>% 
  dplyr::mutate(X=paste0(Plant_sp, "_", Plant_id)) %>% 
  column_to_rownames("X") %>% select(-Plant_id, -Plant_sp)

inter.pol <- matrix(0.1, nrow=ncol(int.mat.clean), ncol=ncol(int.mat.clean))
rownames(inter.pol) <- colnames(int.mat.clean)
colnames(inter.pol) <- colnames(int.mat.clean)

# mutualistic benefits received by pollinators submatrix
gamma.pol <- matrix(-0.2, nrow=ncol(select(y, -Plant_id, -Plant_sp)), ncol=nrow(y))
rownames(gamma.pol) <- colnames(int.mat.clean)
colnames(gamma.pol) <- rownames(int.mat.clean)

pol.submat <- cbind(gamma.pol, inter.pol)

# A matrix

A <- rbind(plant.submat, pol.submat) %>% as.matrix()

# relative proportion of visitations per plant species
pol.plant.sp <- abpol %>% filter(Plant_sp != sp | Plant_id %in% rownames(gen.output[[39]])) %>%
  dplyr::select(-Plant_id) %>% group_by(Plant_sp) %>% summarise_all(list(sum))
vis <- rowSums(pol.plant.sp[,-1])/sum(pol.plant.sp[,-1])
vis <- scale_to(rowSums(pol.plant.sp[,-1]), 1)




# create a list for both scenarios
A.list.all <- c(A.list, list(A))
vis.prop.all <- c(vis.per.sp, list(vis))





###### Using MVN approach to calculate the probability of exclusion of species (just a test, this is not the approach we will use)

#e.prob.list <- list()
#
#for(z in 1:length(A.list.all)){
#  
#  M <- A.list.all[[z]] # this is M
#  animal_abund <- abpol.prop # this is animal_abund
#  n.pol <- length(abpol.prop) # this is animal_abund
#  n.plant.ind <- nrow(M) - n.pol
#  n.plant.sp <- 3
#  
#  n.sp1 <- count <- sum(grepl("^CLIB", rownames(M)))
#  n.sp2 <- count <- sum(grepl("^HCOM", rownames(M)))
#  n.sp3 <- count <- sum(grepl("^HHAL", rownames(M)))
#  
#  M <- unname(M)
#  
#  
#  types <- c(rep('P', n.plant.ind), rep('A', n.pol)) 
#  species <- c(rep(1, n.sp1),
#               rep(2, n.sp2),
#               rep(3, n.sp3),
#               seq(4, n.pol + n.plant.sp)) 
#  subsets <- c(seq(1, n.plant.ind), rep(1, n.pol))
#  
#  names(animal_abund) <- which(types=='A')
#  
#  
#  M_and_constraints <- M
#  
#  for(k in which(types=='A')[-1]){
#    M_and_constraints[k,] <- 0
#    M_and_constraints[k,which(types=='A')[1]] <- -animal_abund[as.character(k)]
#    M_and_constraints[k,k] <- animal_abund[1]
#  }
#  
#  ###
#  # this shows how to impose a constraint that forces growth rates to be equal to each other
#  M_and_constraints2 <- M_and_constraints
#  
#  for(i in 2 : n.sp1){ # first individual of sp1
#    M_and_constraints2[i,] <- M_and_constraints2[i,] - M_and_constraints2[1,] #sp1
#  }
#  
#  for(i in (n.sp1 + 2) : (n.sp1 + n.sp2)){ # first individual of sp2
#    M_and_constraints2[i,] <- M_and_constraints2[i,] - M_and_constraints2[(n.sp1 + 1),] #sp2
#  }
#  
#  for(i in (n.sp1 + n.sp2 + 2) : (n.sp1 + n.sp2 + n.sp3)){ # first individual of sp3
#    M_and_constraints2[i,] <- M_and_constraints2[i,] - M_and_constraints2[(n.sp1 + n.sp2 + 1),] # sp3
#  } 
#  
#  
#  # indices for nodes with different properties
#  free_spp <- c(1, n.sp1+1, n.sp1+n.sp2+1, n.plant.ind+1)
#  abundance_constrained_spp <- c(n.plant.ind+2:n.pol)
#  growth_constrained_spp <- setdiff(1:n.plant.ind, free_spp)
#  
#  # constrained abundance relationship matrix
#  Pi2 <- matrix(0,length(free_spp)+length(abundance_constrained_spp),length(free_spp))
#  Pi2[1,1] <- 1 # plant sp 1
#  Pi2[2,2] <- 1 # plant sp 2
#  Pi2[3,3] <- 1 # plant sp 3
#  Pi2[,4] <- c(rep(0,3), unname(animal_abund))
#    
#  # matrix to map from density of all other species to densities of growth-rate constrained species
#  L <- -1 * M_and_constraints2[growth_constrained_spp,growth_constrained_spp]
#  R <- M_and_constraints2[growth_constrained_spp,-growth_constrained_spp]
#  little_rho <- solve(L) %*% R
#  
#  # the matrix little rho matrix is only section of a new "rho-like" matrix
#  weird_rho <- matrix(0,nrow(M),nrow(Pi2))
#  weird_rho[growth_constrained_spp,] <- little_rho
#  
#  # this maps directly from a 'free' or 'density-constrained' density to all densities in the model
#  
#  weird_rho[1, 1] <- 1 #sp1
#  weird_rho[n.sp1+1, 2] <- 1 #sp2
#  weird_rho[n.sp1+n.sp2+1, 3] <- 1 #sp3
#  
#  # pollinators
#  for (p in (n.plant.ind + 1):length(subsets)) {
#    col_index <- p + 4 - (n.plant.ind + 1)
#        if (col_index >= 4 && col_index <= ncol(weird_rho)) {
#      weird_rho[p, col_index] <- 1
#    }
#  }
#
#  
#  # this is the only portion of the matrix needed for the MVN approach
#  M_sub2 <- (M %*% weird_rho %*% Pi2)[free_spp,]
#  
#  # taking into account that plant species differ in the number of visits received
#  M_sub2[1:3, 4] <- M_sub2[1:3, 4]*vis.prop.all[[z]]
#  
#  number_Omega_replicates <- 200
#  number_boot_replicates <- number_Omega_replicates
#  
#  # use this relevant portion of the matrix for the MVN approach
#  value <- prob_extinction_4_int_matrix(M_sub2, number_Omega_replicates, number_boot_replicates)
#  value.df <- as.data.frame(value) %>% dplyr::select(species, prob_excl_mean) %>% mutate(net=z)
#  e.prob.list[[z]] <- value.df
#  
#}
#  
#
#####
#  
#
#ep.res <- do.call("rbind", e.prob.list)
#
#ep.spec.res <- ep.res %>%
#  filter(net == "11") %>%
#  dplyr::select(species, prob_excl_mean) %>%
#  dplyr::rename(prob_excl_mean_spec = prob_excl_mean)
#
#result <- ep.res %>%
#  filter(net != "11") %>%
#  left_join(ep.spec.res, by = "species") %>%
#  mutate(prob_excl_diff = prob_excl_mean - prob_excl_mean_spec) %>%
#  dplyr::select(species, net, prob_excl_diff)
#
#summ.res <- result %>%
#  group_by(species) %>%
#  dplyr::summarise(
#    mean_prob_excl_diff = mean(prob_excl_diff),
#    sd_prob_excl_diff = sd(prob_excl_diff)
#  )
#
#
#plot.comm <- ggplot(filter(summ.res, species != "sp_4"), aes(species, mean_prob_excl_diff, shape=species,
#                                                             colour=species, group=species)) +  
#  geom_hline(yintercept = 0, linetype="dashed", colour="grey50", size=0.5) +
#  geom_point(alpha=1, size=4, stroke=1.5, position = position_dodge(width = 0.5)) +
#  geom_errorbar(aes(ymin=mean_prob_excl_diff-sd_prob_excl_diff, 
#                    ymax=mean_prob_excl_diff+sd_prob_excl_diff), 
#                width=.1, position = position_dodge(width = 0.5)) +
#  theme_bw() + scale_color_manual(values=c("dodgerblue3","burlywood4", "pink2")) +
#  ylab(bquote("\u0394p"^"E")) + 
#  theme(text=element_text(size=20), axis.title.x = element_blank(), 
#        legend.text = element_text(face="italic", size=13),
#        axis.text.x = element_text(face="italic", size=13),
#        legend.title = element_blank(),
#        legend.spacing.y = unit(3, 'cm'),
#        legend.position = "none") +
#  scale_color_manual(values=c("dodgerblue3","burlywood4", "pink2"), name = "Plant species", 
#                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
#  scale_shape_manual(values=c(16, 15, 17), name = "Plant species", 
#                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
#  scale_x_discrete(labels=c("Cistus \nlibatonis", 
#                            "Halimium \ncalycinum", 
#                            "Halimium \nhalimifolium")) + 
#  ggtitle(paste0("Specialized ", focal))




###### using the random sampling method to calculate a proxy for this measure

#pb = txtProgressBar(min = 0, max = length(A.list.all), initial = 0) 
#
#vul.sp.list <- list()
#
#for(z in 1:length(A.list.all)){
#  
#  M <- A.list.all[[z]] # this is M
#  animal_abund <- abpol.prop # this is animal_abund
#  n.pol <- length(abpol.prop) # this is animal_abund
#  n.plant.ind <- nrow(M) - n.pol
#  n.plant.sp <- 3
#  
#  n.sp1 <- count <- sum(grepl("^CLIB", rownames(M)))
#  n.sp2 <- count <- sum(grepl("^HCOM", rownames(M)))
#  n.sp3 <- count <- sum(grepl("^HHAL", rownames(M)))
#  
#  M <- unname(M)
#  
#  
#  types <- c(rep('P', n.plant.ind), rep('A', n.pol)) 
#  species <- c(rep(1, n.sp1),
#               rep(2, n.sp2),
#               rep(3, n.sp3),
#               seq(4, n.pol + n.plant.sp)) 
#  subsets <- c(seq(1, n.plant.ind), rep(1, n.pol))
#  
#  names(animal_abund) <- which(types=='A')
#  
#  # modified M matrix that also includes the constraints on pollinator relative abundances
#  # modified M matrix that also includes the constraints on pollinator relative abundances
#  M_and_constraints <- Matrix::Matrix(M, sparse=TRUE)
#  
#  M_and_constraints[which(types=='A')[-1],] <- 0
#  M_and_constraints[which(types=='A')[-1],which(types=='A')[1]] <- -animal_abund[-1]
#  diag(M_and_constraints)[which(types=='A')[-1]] <- animal_abund[1]
#  
#  #M_and_constraints <- M
#  #
#  #for(k in which(types=='A')[-1]){
#  #  M_and_constraints[k,] <- 0
#  #  M_and_constraints[k,which(types=='A')[1]] <- -animal_abund[as.character(k)]
#  #  M_and_constraints[k,k] <- animal_abund[1]
#  #}
#  
#  # we will calculate structural stability via randomly sampled growth rate vectors
#  # the number below determines how many samples we draw; the smaller the feasibility domain, the more we need in order to get a good estimate of Omega
#  n_samples <- 1E5
#  
#  ##### random sampling approach
#  
#  # randomly sampled growth rate vector on the unit sphere
#  # NOTE: the argument N needs to be 1 less than the required number of growth rates
#  growth_rates_yes_constraint <- as.matrix(
#    pracma::rands(
#      n_samples,
#      N=length(unique(species[which(types=='P')]))
#    )
#  )
#  colnames(growth_rates_yes_constraint) <- c(
#    paste0('r_P',unique(species[which(types=='P')])),
#    paste0('r_A',species[which(types=='A')[1]])
#  )
#  
#  # add constraint columns for the abundances of additional pollinators
#  growth_rates_yes_constraint <- Matrix::Matrix(
#    cbind(
#      growth_rates_yes_constraint,
#      matrix(0, n_samples, sum(types=='A')-1)
#    ),
#    sparse = TRUE
#  )
#  colnames(growth_rates_yes_constraint)[(length(unique(species[which(types=='P')]))+2):ncol(growth_rates_yes_constraint)] <- paste0('constraint_A',species[which(types=='A')[-1]])
#  
#  # create a matrix that maps values in the growth_rates_yes_contraint to the M matrix
#  Z_mat <- Matrix::Matrix(
#    0,
#    nrow = ncol(growth_rates_yes_constraint),
#    ncol = nrow(M),
#    sparse = TRUE
#  )
#  for(m in 1:ncol(Z_mat)){
#    Z_mat[species[m],m] <- 1
#  }
#  
#  # solve for the equilibria
#  N_equil <- growth_rates_yes_constraint %*% Z_mat %*% t(solve(M_and_constraints))
#  
#  # for each sampled growth rate vector determine whether the equilibrium is feasible
#  feasible <- colSums(t(N_equil>0)) == ncol(N_equil)
#  
#  # to make sure the number of feasible vectors is greater than 1 to calculate the species' normalized densities
#  if (sum(feasible)>100) {
#  # get the density of each species (we need to sum across types) for feasible vectors and normalize
#  feasible_N <- N_equil[feasible,]
#  feasible_N <- feasible_N %*% t(Z_mat)
#  feasible_N <- sweep(feasible_N, 2, colSums(feasible_N), FUN = "/")
#  
#  # assign each vector to the plant species with the lowest density
#  feasible_N <- feasible_N[,1:3]
#  result_matrix <- matrix(0, nrow = nrow(feasible_N), ncol = ncol(feasible_N))
#  for (row in 1:nrow(feasible_N)) {
#    min_index <- which.min(feasible_N[row, ])
#    result_matrix[row, min_index] <- 1
#  }
#  metric <- colSums(result_matrix)/nrow(result_matrix)
#  } else {metric <- rep(NA,3)
#  }
#  
#  # big Omega is the proportion of vectors with feasible equilibria
#  #Omega_yes_constraint <- sum(feasible) / n_samples
#  
#  scenario <- ifelse(z > 100, "gen", "random") # change with the number of randomizations
#  
#  vul.sp <- data.frame(scenario=scenario,
#                      vul_sp_1=metric[1],
#                      vul_sp_2=metric[2],
#                      vul_sp_3=metric[3]
#                      
#  )
#  
#  vul.sp.list[[z]] <- vul.sp
#  setTxtProgressBar(pb,z)
#  
#}



vul.sp.list <- list()
metric.sp.list <- list()

pb <- txtProgressBar(min = 0, max = length(A.list.all), style = 3)

for (z in 1:length(A.list.all)) {
  
  M <- A.list.all[[z]] # this is M
  animal_abund <- abpol.prop # this is animal_abund
  n.pol <- length(abpol.prop) # this is animal_abund
  n.plant.ind <- nrow(M) - n.pol
  n.plant.sp <- 3
  
  n.sp1 <- sum(grepl("^CLIB", rownames(M)))
  n.sp2 <- sum(grepl("^HCOM", rownames(M)))
  n.sp3 <- sum(grepl("^HHAL", rownames(M)))
  
  M <- unname(M)
  
  types <- c(rep('P', n.plant.ind), rep('A', n.pol)) 
  species <- c(rep(1, n.sp1),
               rep(2, n.sp2),
               rep(3, n.sp3),
               seq(4, n.pol + n.plant.sp)) 
  subsets <- c(seq(1, n.plant.ind), rep(1, n.pol))
  
  names(animal_abund) <- which(types == 'A')
  
  # Create a modified M matrix with constraints
  M_and_constraints <- Matrix::Matrix(M, sparse = TRUE)
  M_and_constraints[which(types == 'A')[-1], ] <- 0
  M_and_constraints[which(types == 'A')[-1], which(types == 'A')[1]] <- -animal_abund[-1]
  diag(M_and_constraints)[which(types == 'A')[-1]] <- animal_abund[1]
  
  # we use the inverse of this matrix below; no need to take the inverse multiple times
  M_and_constraints_inverse <- solve(M_and_constraints)

  # map of growth rate and constraint values to the M matrix
  Z_mat <- Matrix::Matrix(
    0,
    nrow = length(unique(species)),
    ncol = nrow(M),
    sparse = TRUE
  )
  for (m in 1:ncol(Z_mat)) {
    Z_mat[species[m], m] <- 1
  }

  # perform this matrix multiplication here since it only needs to be done once
  M_inverse_and_tZ <- M_and_constraints_inverse %*% t(Z_mat)

  n_samples <- 1E5
  cumulative_feasible_N <- NULL
  
  repeat {
    # Randomly sample growth rate vectors
    growth_rates_yes_constraint <- as.matrix(
      pracma::rands(
        n_samples,
        N = length(unique(species[which(types == 'P')]))
      )
    )
    colnames(growth_rates_yes_constraint) <- c(
      paste0('r_P', unique(species[which(types == 'P')])),
      paste0('r_A', species[which(types == 'A')[1]])
    )
    
    growth_rates_yes_constraint <- Matrix::Matrix(
      cbind(
        growth_rates_yes_constraint,
        matrix(0, n_samples, sum(types == 'A') - 1)
      ),
      sparse = TRUE
    )
    colnames(growth_rates_yes_constraint)[
      (length(unique(species[which(types == 'P')])) + 2):ncol(growth_rates_yes_constraint)
    ] <- paste0('constraint_A', species[which(types == 'A')[-1]])
    
    # Solve for equilibria and check if they are feasible
    N_equil <- M_inverse_and_tZ %*% t(growth_rates_yes_constraint)
    feasible <- colSums((N_equil>0)) == nrow(N_equil)
    
    # Accumulate feasible_N across repeated sampling
    if (sum(feasible) > 0) {
      feasible_N <- N_equil[,feasible,drop=FALSE]
      cumulative_feasible_N <- cbind(cumulative_feasible_N, feasible_N)
    }

    # Break the loop if the condition of having at least 100 feasible vectors is satisfied
    if (!is.null(cumulative_feasible_N) && ncol(cumulative_feasible_N) >= 200) break
  }

  # get plant species level densities for the feasible equilibria
  plant_feasible_N <- t((Z_mat %*% cumulative_feasible_N)[1:3 ,])

  # normalize each plant relative to its maximum abundance across equilbria
  plant_normalized_N <- sweep(
    plant_feasible_N,
    2,
    apply(plant_feasible_N, 2, max),
    '/'
  )

  # determine which plant is furthest from its maximum
  plant_most_vulnerable <- t(apply(
    plant_normalized_N,
    1,
    function(x) x == min(x)
  ))
  
  scenario <- ifelse(z > n.random, "gen", "random")
  
  vul.sp <- as.data.frame(as.matrix(plant_normalized_N)) %>% mutate(scenario=scenario, ran=z)
  metric <- colSums(plant_most_vulnerable) / nrow(plant_most_vulnerable)
  
  metric.sp <- as.data.frame(t(metric)) %>% mutate(scenario=scenario)
  
  vul.sp.list[[z]] <- vul.sp
  metric.sp.list[[z]] <- metric.sp
  setTxtProgressBar(pb, z)
}

close(pb)

vul.sp.all <- bind_rows(vul.sp.list)
metric.all <- bind_rows(metric.sp.list)


### just checking the second species
#vul.sp.all <- vul.sp.all.list[[1]]
#metric.all <- metric.all.list[[1]]
#
####
#random_data <- metric.all[metric.all$scenario == "random", ] 
#gen_data <- metric.all[metric.all$scenario == "gen", ]
#
#random_data %<>%
#  pivot_longer(cols = starts_with("V"), names_to = "species", values_to = "value") %>%
#  filter(scenario =="random") %>% dplyr::group_by(species) %>% head(100)
#
#gen_data %<>%
#  pivot_longer(cols = starts_with("V"), names_to = "species", values_to = "value")
#
#ggplot(random_data, aes(species, value)) +
#  geom_violin() + theme_bw() + 
#  geom_point(data=gen_data, aes(species, value), shape=4, color="#5DA6A7", stroke=2)

####

random_data <- metric.all[metric.all$scenario == "random", ] %>% dplyr::select(!scenario)
gen_data <- metric.all[metric.all$scenario == "gen", ] %>% dplyr::select(!scenario)

res <- random_data

for (i in 1:ncol(random_data)) {
  res[, i] <- gen_data[, i] - random_data[, i]
}


long.res <- res %>%
  pivot_longer(cols = starts_with("V"), names_to = "species", values_to = "value")


summ.res <- as.data.frame(as.matrix(aggregate(value ~ species, long.res, 
                                              function(x) c(ci=CI(x, ci = 0.95))))) 
summ.res[, 2:4] <- lapply(summ.res[, 2:4], as.numeric)

#summ.res <- long.res %>% group_by(species) %>% 
#  dplyr::summarise(mean=mean(value), sd=sd(value))

sp.label <- switch(focal,
                   "CLIB" = ~italic("Cistus libanotis"),
                   "HCOM" = ~italic("Halimium calycinum"),
                   "HHAL" = ~italic("Halimium halimifolium"))


plot.comm <- ggplot(summ.res, aes(species, value.ci.mean, shape=species,
                                                             colour=species, group=species)) +  
  geom_hline(yintercept = 0, linetype="dashed", colour="grey50", size=0.5) +
  geom_point(alpha=1, size=4, stroke=2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=value.ci.lower, 
                    ymax=value.ci.upper), 
                width=.1, position = position_dodge(width = 0.5)) +
  theme_bw() + scale_color_manual(values=c("#5DA6A7","#F2C57C", "#E8875A")) +
  ylab("\u0394v") + 
  #ylab(bquote("\u0394p"^"E")) + 
  theme(text=element_text(size=20), axis.title.x = element_blank(), 
        legend.text = element_text(face="italic", size=13),
        axis.text.x = element_text(face="italic", size=13),
        legend.title = element_blank(),
        legend.spacing.y = unit(3, 'cm'),
        legend.position = "none",
        plot.title = element_text(size=16, hjust = 0.5)) +
  scale_color_manual(values=c("#5DA6A7","#F2C57C", "#E8875A"), name = "Plant species", 
                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
  scale_shape_manual(values=c(16, 15, 17), name = "Plant species", 
                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
  scale_x_discrete(labels=c("Cistus \nlibatonis", 
                            "Halimium \ncalycinum", 
                            "Halimium \nhalimifolium")) +
  ggtitle(sp.label)
#+ylim(-0.30, 0.30)




#ran.vul.sp.all <- vul.sp.all %>% filter(scenario=="random") %>% dplyr::select(!c(ran, scenario))
#gen.vul.sp.all <- vul.sp.all %>% filter(scenario=="gen") %>% dplyr::select(!c(ran, scenario))
#
### heterogeneous scenario
#
#cols <- c("V1", "V2", "V3")
#ran_min_col <- apply(ran.vul.sp.all[,c(cols)],1, \(k) cols[which.min(k)])
#ran.vul.sp.all <- cbind(ran.vul.sp.all, ran_min_col)
#
#ran.vul.sp.all <- ran.vul.sp.all %>%
#  slice_sample(n = 1003)
#
#ran.tern.plot <- ggtern(data= ran.vul.sp.all, aes(V1, V2, V3, 
#                                                  color=ran_min_col,
#                                                  fill=ran_min_col)) + 
#  geom_mask() +
#  geom_point(size=1.5, alpha=0.7)  +
#  #stat_mean_ellipse(geom='polygon',steps=100, alpha=0.5) + 
#  theme_bw() +
#  theme_clockwise() +
#  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
#  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
#  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
#  scale_color_manual(values=c("#5DA6A7","#F2C57C", "#E8875A"), name = "Plant species", 
#                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
#  scale_fill_manual(values=c("#5DA6A7","#F2C57C", "#E8875A")) +
#  labs(x = "Clib \nM", y="Hcom \nM", z="Hhal \nM") + theme(legend.position="none") 
#
#
### mostly specialists scenario
#
#gen_min_col <- apply(gen.vul.sp.all[,c(cols)],1, \(k) cols[which.min(k)])
#gen.vul.sp.all <- cbind(gen.vul.sp.all, gen_min_col)
#
#gen.tern.plot <- ggtern(data= gen.vul.sp.all, aes(V1, V2, V3, 
#                                                  color=gen_min_col,
#                                                  fill=gen_min_col)) + 
#  geom_mask() +
#  geom_point(size=1.5, alpha=0.7)  +
#  #stat_mean_ellipse(geom='polygon',steps=100, alpha=0.5) + 
#  theme_bw() +
#  theme_clockwise() +
#  scale_L_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
#  scale_R_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
#  scale_T_continuous(breaks = 0:5 / 5, labels = 0:5/ 5) +
#  scale_color_manual(values=c("#5DA6A7","#F2C57C", "#E8875A"), name = "Plant species", 
#                     labels = c("Cistus \nlibanotis", "Halimium \ncalycinum", "Halimium \nhalimifolium")) +
#  scale_fill_manual(values=c("#5DA6A7","#F2C57C", "#E8875A")) +
#  labs(x = "Clib \nS", y="Hcom \nM", z="Hhal \nM") + theme(legend.position="none") 
#
#
#grid.arrange(arrangeGrob(grobs=list(ran.tern.plot, gen.tern.plot), ncol=2, nrow =1))

