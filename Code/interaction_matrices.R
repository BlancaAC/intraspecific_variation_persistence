

sp <- focal
n.random <- 400

beta.mat <- read.csv(paste0("Output/beta_matrix_", sp, ".csv")) %>% column_to_rownames("X")
vis.mat <- read.csv(paste0("Output/vis_matrix_", sp, ".csv")) %>% column_to_rownames("X") 
attr <- read.csv(paste0("Output/plant_attr_", sp, ".csv")) %>% column_to_rownames("X") 
prop.vis.mat <- vis.mat %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0)) %>% 
  dplyr::mutate(across(where(is.numeric), ~ ./sum(.))) 
#colSums(prop.vis.mat)
int.mat0 <- beta.mat * prop.vis.mat

# relative pollinator abundance
abpol.prop <- (colSums(vis.mat)/sum(vis.mat))
abpol.total <- sum(vis.mat)

# flower abundance 
flower.prop <- attr %>% select(N_flowers) %>% mutate(prop_flowers=N_flowers/sum(N_flowers))
#sum(flower.prop$prop_flowers)



if (per.study.plot == TRUE) { 
## split into different plots
plot.id <- read.csv("Data/all_plant_traits.csv") %>% filter(Plant_sp==sp, Plot==p) 
beta.mat %<>% filter(row.names(beta.mat) %in% plot.id$Plant_id) %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
vis.mat %<>% filter(row.names(vis.mat) %in% plot.id$Plant_id) %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
attr %<>% filter(row.names(attr) %in% plot.id$Plant_id)
prop.vis.mat %<>% filter(row.names(prop.vis.mat) %in% plot.id$Plant_id) %>% 
  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
# relative pollinator abundance
abpol.prop <- (colSums(vis.mat)/sum(vis.mat))
abpol.total <- sum(vis.mat)

int.mat0 <- beta.mat * prop.vis.mat
}
#####


int.mat <- int.mat0 * flower.prop$N_flowers


set.seed(123)
if(per.plant.ind==FALSE){ 
#### Increasing population size by incorporating plant individuals randomly.

rand.output <- list()
for(i in 2:nrow(int.mat)){
  rand.output[[i]] <- lapply(1:n.random, function(x) 
    int.mat[sample(nrow(int.mat), size = i, replace = F),]) 
  for (j in seq_along(rand.output[[i]])) {
    #rand.output[[i]][[j]] %<>% dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))
    rand.output[[i]][[j]]$Plant_id <- NULL
    rand.output[[i]][[j]] <- rand.output[[i]][[j]][order(row.names(rand.output[[i]][[j]])), ] 
    rand.output[[i]][[j]] <- rand.output[[i]][[j]][, order(colnames(rand.output[[i]][[j]]))]
  }
}
rand.output[sapply(rand.output, is.null)] <- NULL # remove elements that only contain 1-4 plant individuals

# mapping of within-species groups to populations
# within-species group X species matrix
Z.plant.ran.list <- rand.output
Z.df.ran.list <- rand.output
Z.ran.list <- rand.output
for (i in seq_along(rand.output)) {
  for (j in seq_along(rand.output[[i]])){
    Z.plant.ran <- matrix(1, nrow(rand.output[[i]][[j]]), 1)
    rownames(Z.plant.ran) <- rownames(rand.output[[i]][[j]])
    colnames(Z.plant.ran) <- sp
    Z.plant.ran.list[[i]][[j]] <- Z.plant.ran # i will need to use this below
    
    Z.pol.ran <- diag(1, ncol(rand.output[[i]][[j]]), ncol(rand.output[[i]][[j]]))
    rownames(Z.pol.ran) <- colnames(rand.output[[i]][[j]])
    colnames(Z.pol.ran) <- colnames(rand.output[[i]][[j]])
    
    ## to collapse only plant individuals
    #Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% 
    #  column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
    #Z.df.list[[i]][[j]] <- Z.df
    #
    #Z <- Z.df %>% as.matrix()
    
    ## to collapse both plant individuals and pollinator species
    Z.df.ran <- cbind(c(rep(1, nrow(rand.output[[i]][[j]])), rep(0, ncol(rand.output[[i]][[j]]))),
                  c(rep(0, nrow(rand.output[[i]][[j]])), rep(1, ncol(rand.output[[i]][[j]]))))
    Z.ran <- Z.df.ran %>% as.matrix()
    
    Z.ran.list[[i]][[j]] <-  Z.ran
  }
}

# partitioning of within-species groups across populations
# within-population group X population matrix

PI.ran.list <- rand.output
for (i in seq_along(Z.df.ran.list)) {
  for (j in seq_along(Z.df.ran.list[[i]])){
    
    # create PI P matrix
    #plant.id <- rownames(Z.df.list[[i]][[j]])
    #plant.if.fl <- filter(attr, rownames(attr) %in% plant.id)
    X.plant.ran <- matrix(1/nrow(Z.df.ran.list[[i]][[j]]), nrow=nrow(Z.df.ran.list[[i]][[j]]))
    #as.matrix(plant.if.fl$N_flowers/sum(plant.if.fl$N_flowers))
    rownames(X.plant.ran) <- rownames(Z.df.ran.list[[i]][[j]])
    colnames(X.plant.ran) <- sp
    X.pol.ran <- diag(1, ncol(int.mat), ncol(int.mat))
    rownames(X.pol.ran) <- colnames(int.mat)
    colnames(X.pol.ran) <- colnames(int.mat)
    
    # PI matrix for plants (relative abundance proportional to flower number)
    PI.ran <- merge(X.plant.ran, X.pol.ran, by = "row.names", all = TRUE) %>% 
      column_to_rownames("Row.names") %>% replace(is.na(.), 0) %>% as.matrix()
    
    # this is to give the same weight to each plant individual within the population
    #PI <- Z.df.list[[i]][[j]] %>% 
    #  mutate(across(c(1), ~ if_else(. ==1, 1/length(Z.plant.list[[i]][[j]]), .))) %>% as.matrix()
    PI.ran.list[[i]][[j]] <- PI.ran
  }
}

# create PI matrix for pollinators (relative abundance proportional to number of visits)
PI.a.ran.list <- rand.output
for (i in seq_along(Z.df.ran.list)) {
  for (j in seq_along(Z.df.ran.list[[i]])){
    PI.a.ran <- as.matrix(cbind(c(1, rep(0, length(abpol.prop))), c(0, (abpol.prop))))
    PI.a.ran.list[[i]][[j]] <- PI.a.ran
  }}


# mutualistic benefits from pollinators to plant individuals

plant.gamma.ran.list <- rand.output
for (i in seq_along(rand.output)) {
  for (j in seq_along(rand.output[[i]])) {
    temp <- as.matrix(rand.output[[i]][[j]])
    # classical approach (assuming a baseline contribution gamma0 and penalizing by plant degree)
    #gamma0 <- coeff.pol # contribution of each pollinator sp (gamma0)
    #d <- rowSums(temp!=0) # degree of plant individuals
    #delta <- 0.2 # mutualistic trade-off
    #plant.gamma <- temp
    #for(row in rownames(temp)) {
    #for(col in colnames(temp)) {
    #    if(temp[row, col]>0){
    #      # mutualistic benefit received by plant individuals from each pollinator sp
    #      plant.gamma[row, col] <- (gamma0[col])/d^delta 
    #    }
    #}
    #}
    #plant.gamma.list[[i]][[j]] <- plant.gamma
    plant.gamma.ran.list[[i]][[j]] <- temp
  }
}

# interaction matrix
# within-species groups X within-species groups matrix

A.ran.list <- rand.output
for (i in seq_along(PI.ran.list)) {
  for (j in seq_along(PI.ran.list[[i]])){
    # setting within-species groups X within-species groups matrix
    A.ran <- diag(0, nrow(PI.ran.list[[i]][[j]]), nrow(PI.ran.list[[i]][[j]])) 
    rownames(A.ran) <- rownames(PI.ran.list[[i]][[j]])
    colnames(A.ran) <- rownames(PI.ran.list[[i]][[j]])
    
    # plant within-species group competition set to 1 above
    nplant <- nrow(Z.plant.ran.list[[i]][[j]])
    ntotal <- ncol(A.ran)
    npol <- ntotal - nplant
    
    #rmat.plant <- abs(matrix(rnorm(nplant*nplant),nrow=nplant))
    #rmat.plant.pol <- abs(matrix(rnorm(nplant*npol),nrow=npol))
    #rmat.pol <- abs(matrix(rnorm(npol*npol),nrow=npol))
    
    A.ran[1:nplant, 1:nplant] <- 0.9
    A.ran[(1+nplant):ntotal, 1:nplant] <- -0.2 # mutualistic benefits of plants given to pollinators (mean-field value)
    A.ran[(1+nplant):ntotal, (1+nplant):ntotal] <- 0.1 # pollinator interspecific competition (mean-field value)
    A.ran[1:nplant, (1+nplant):ntotal] <- -(plant.gamma.ran.list[[i]][[j]]) # mutualistic benefits of pollinators given to plants
    for (row in rownames(A.ran)){
      for (col in colnames(A.ran)){
        if(row==col){
          A.ran[row, col] <- 1 # intraspecific competition set to 1
        }
      }
    }
    
    A.ran.list[[i]][[j]] <- A.ran
  }
}



#### same but unnesting the A.ran.list object
A.list.flat <- flatten(A.ran.list)







#### Increasing population size by incorporating plant individuals with increasing generalization in pollinator use

specif <- specieslevel(vis.mat, index=c("species specificity"), level="lower") %>% 
  rownames_to_column("Plant_id") %>% arrange((-species.specificity.index))

# sorting by degree
gen.output <- list()
int.mat.d <- int.mat %>% mutate(degree= rowSums(.!=0)) %>% 
  rownames_to_column("Plant_id") %>% 
  left_join(specif, by="Plant_id") %>% 
  arrange(degree, -species.specificity.index) %>% column_to_rownames("Plant_id")

for(i in seq(from = 2, to = nrow(int.mat.d), by = 1)) { 
  gen.output[[i]] <- int.mat.d %>% select(-degree, -species.specificity.index) %>% dplyr::slice(1:i) 
  #%>%  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))    
  gen.output[[i]] <- gen.output[[i]][order(row.names(gen.output[[i]])), ] 
  gen.output[[i]] <- gen.output[[i]][, order(colnames(gen.output[[i]]))] 
}

gen.output[sapply(gen.output, is.null)] <- NULL # remove elements that only contain 1-4 plant individuals



# mapping of within-species groups to populations
# within-species group X species matrix

Z.plant.list <- gen.output
Z.df.list <- gen.output
Z.list <- gen.output
for (i in seq_along(gen.output)) {
  Z.plant <- matrix(1, nrow(gen.output[[i]]), 1)
  rownames(Z.plant) <- rownames(gen.output[[i]])
  colnames(Z.plant) <- sp
  Z.plant.list[[i]] <- Z.plant # i will need to use this below
  
  Z.pol <- diag(1, ncol(gen.output[[i]]), ncol(gen.output[[i]]))
  rownames(Z.pol) <- colnames(gen.output[[i]])
  colnames(Z.pol) <- colnames(gen.output[[i]])
  
  # to collapse only plant individuals
  #Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% 
  #  column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
  #Z.df.list[[i]] <- Z.df
  #
  #Z <- Z.df %>% as.matrix()
  
  ## to collapse both plant individuals and pollinator species
  Z.df <- cbind(c(rep(1, nrow(gen.output[[i]])), rep(0, ncol(gen.output[[i]]))),
                c(rep(0, nrow(gen.output[[i]])), rep(1, ncol(gen.output[[i]]))))
  Z <- Z.df %>% as.matrix()
  
  Z.list[[i]] <-  Z
}

# partitioning of within-species groups across populations
# within-population group X population matrix

PI.list <- gen.output
for (i in seq_along(Z.df.list)) {
  
  # create PI P matrix
  #plant.id <- rownames(Z.df.list[[i]])
  #plant.if.fl <- filter(attr, rownames(attr) %in% plant.id)
  #X.plant <- as.matrix(plant.if.fl$N_flowers/sum(plant.if.fl$N_flowers))
  #rownames(X.plant) <- rownames(plant.if.fl)
  X.plant <- matrix(1/nrow(Z.df.list[[i]]), nrow=nrow(Z.df.list[[i]]))
  rownames(X.plant) <- rownames(Z.df.list[[i]])
  colnames(X.plant) <- sp
  X.pol <- diag(1, ncol(int.mat), ncol(int.mat))
  rownames(X.pol) <- colnames(int.mat)
  colnames(X.pol) <- colnames(int.mat)
  
  # PI matrix for plants (relative abundance proportional to flower number)
  PI <- merge(X.plant, X.pol, by = "row.names", all = TRUE) %>% 
    column_to_rownames("Row.names") %>% replace(is.na(.), 0) %>% as.matrix()
  
  # this is to give the same weight to each plant individual within the population
  #PI <- Z.df.list[[i]] %>% 
  #  mutate(across(c(1), ~ if_else(. ==1, 1/length(Z.plant.list[[i]]), .))) %>% as.matrix()
  
  PI.list[[i]] <- PI
}

# create PI matrix for pollinators (relative abundance proportional to number of visits)
PI.a.list <- gen.output
for (i in seq_along(Z.df.list)) {
  PI.a <- as.matrix(cbind(c(1, rep(0, length(abpol.prop))), c(0, (abpol.prop))))
  PI.a.list[[i]] <- PI.a
}

# mutualistic benefits from pollinators to plant individuals

plant.gamma.list <- gen.output
for (i in seq_along(gen.output)) {
  temp <- as.matrix(gen.output[[i]])
  plant.gamma.list[[i]] <- temp
}


#visweb(gen.output[[length(gen.output)]], type="nested")
#hist(total.gamma$degree)

# interaction matrix
# within-species groups X within-species groups matrix

A.gen.list <- gen.output
for (i in seq_along(PI.list)) {
  # setting within-species groups X within-species groups matrix
  A <- diag(0, nrow(PI.list[[i]]), nrow(PI.list[[i]])) 
  rownames(A) <- rownames(PI.list[[i]])
  colnames(A) <- rownames(PI.list[[i]])
  
  # plant within-species group competition set to 1 above
  nplant <- nrow(Z.plant.list[[i]])
  ntotal <- ncol(A)
  npol <- ntotal - nplant
  
  #rmat.plant <- abs(matrix(rnorm(nplant*nplant),nrow=nplant))
  #rmat.plant.pol <- abs(matrix(rnorm(nplant*npol),nrow=npol))
  #rmat.pol <- abs(matrix(rnorm(npol*npol),nrow=npol))
  
  A[1:nplant, 1:nplant] <- 0.9
  A[(1+nplant):ntotal, 1:nplant] <- -0.2 # mutualistic benefits of plants given to pollinators (mean-field value)
  A[(1+nplant):ntotal, (1+nplant):ntotal] <- 0.1 # pollinator interspecific competition (mean-field value)
  A[1:nplant, (1+nplant):ntotal] <- -(plant.gamma.list[[i]]) # mutualistic benefits of pollinators given to plants
  for (row in rownames(A)){
    for (col in colnames(A)){
      if(row==col){
        A[row, col] <- 1 # intraspecific competition set to 1
      }
    }
  }
  
  A.gen.list[[i]] <- A
  
}


A.list <- c(A.list.flat, A.gen.list)

}



##### including plants individually

if(per.plant.ind==TRUE){

gen.output <- list()
int.mat.d <- int.mat %>% mutate(degree= rowSums(.!=0)) %>% 
  rownames_to_column("Plant_id") %>% 
  arrange(degree) %>% column_to_rownames("Plant_id")


for(i in 1:nrow(int.mat.d)) { 
  gen.output[[i]] <- int.mat[i,]
  #%>%  dplyr::select(where( ~ is.numeric(.x) && sum(.x) != 0))    
}



# mapping of within-species groups to populations
# within-species group X species matrix

Z.plant.list <- gen.output
Z.df.list <- gen.output
Z.list <- gen.output
for (i in seq_along(gen.output)) {
  Z.plant <- matrix(1, nrow(gen.output[[i]]), 1)
  rownames(Z.plant) <- rownames(gen.output[[i]])
  colnames(Z.plant) <- sp
  Z.plant.list[[i]] <- Z.plant # i will need to use this below
  
  Z.pol <- diag(1, ncol(gen.output[[i]]), ncol(gen.output[[i]]))
  rownames(Z.pol) <- colnames(gen.output[[i]])
  colnames(Z.pol) <- colnames(gen.output[[i]])
  
  # to collapse only plant individuals
  #Z.df <- merge(Z.plant, Z.pol, by = "row.names", all = TRUE) %>% 
  #  column_to_rownames("Row.names") %>% replace(is.na(.), 0) 
  #Z.df.list[[i]] <- Z.df
  #
  #Z <- Z.df %>% as.matrix()
  
  ## to collapse both plant individuals and pollinator species
  Z.df <- cbind(c(rep(1, nrow(gen.output[[i]])), rep(0, ncol(gen.output[[i]]))),
                c(rep(0, nrow(gen.output[[i]])), rep(1, ncol(gen.output[[i]]))))
  Z <- Z.df %>% as.matrix()
  
  Z.list[[i]] <-  Z
}

# partitioning of within-species groups across populations
# within-population group X population matrix

PI.list <- gen.output
for (i in seq_along(Z.df.list)) {
  
  # create PI P matrix
  #plant.id <- rownames(Z.df.list[[i]])
  #plant.if.fl <- filter(attr, rownames(attr) %in% plant.id)
  #X.plant <- as.matrix(plant.if.fl$N_flowers/sum(plant.if.fl$N_flowers))
  #rownames(X.plant) <- rownames(plant.if.fl)
  X.plant <- matrix(1/nrow(Z.df.list[[i]]), nrow=nrow(Z.df.list[[i]]))
  rownames(X.plant) <- rownames(Z.df.list[[i]])
  colnames(X.plant) <- sp
  X.pol <- diag(1, ncol(int.mat), ncol(int.mat))
  rownames(X.pol) <- colnames(int.mat)
  colnames(X.pol) <- colnames(int.mat)
  
  # PI matrix for plants (relative abundance proportional to flower number)
  PI <- merge(X.plant, X.pol, by = "row.names", all = TRUE) %>% 
    column_to_rownames("Row.names") %>% replace(is.na(.), 0) %>% as.matrix()
  
  # this is to give the same weight to each plant individual within the population
  #PI <- Z.df.list[[i]] %>% 
  #  mutate(across(c(1), ~ if_else(. ==1, 1/length(Z.plant.list[[i]]), .))) %>% as.matrix()
  
  PI.list[[i]] <- PI
}

# create PI matrix for pollinators (relative abundance proportional to number of visits)
PI.a.list <- gen.output
for (i in seq_along(Z.df.list)) {
  PI.a <- as.matrix(cbind(c(1, rep(0, length(abpol.prop))), c(0, (abpol.prop))))
  PI.a.list[[i]] <- PI.a
}

# mutualistic benefits from pollinators to plant individuals

plant.gamma.list <- gen.output
for (i in seq_along(gen.output)) {
  temp <- as.matrix(gen.output[[i]])
  plant.gamma.list[[i]] <- temp
}


# interaction matrix
# within-species groups X within-species groups matrix

A.gen.list <- gen.output
for (i in seq_along(PI.list)) {
  # setting within-species groups X within-species groups matrix
  A <- diag(0, nrow(PI.list[[i]]), nrow(PI.list[[i]])) 
  rownames(A) <- rownames(PI.list[[i]])
  colnames(A) <- rownames(PI.list[[i]])
  
  # plant within-species group competition set to 1 above
  nplant <- nrow(Z.plant.list[[i]])
  ntotal <- ncol(A)
  npol <- ntotal - nplant
  
  #rmat.plant <- abs(matrix(rnorm(nplant*nplant),nrow=nplant))
  #rmat.plant.pol <- abs(matrix(rnorm(nplant*npol),nrow=npol))
  #rmat.pol <- abs(matrix(rnorm(npol*npol),nrow=npol))
  
  A[(1+nplant):ntotal, 1:nplant] <- -0.2 # mutualistic benefits of plants given to pollinators (mean-field value)
  A[(1+nplant):ntotal, (1+nplant):ntotal] <- 0.1 # pollinator interspecific competition (mean-field value)
  A[1:nplant, (1+nplant):ntotal] <- -(plant.gamma.list[[i]]) # mutualistic benefits of pollinators given to plants
  for (row in rownames(A)){
    for (col in colnames(A)){
      if(row==col){
        A[row, col] <- 1 # intraspecific competition set to 1
      }
    }
  }
  
  A.gen.list[[i]] <- A
  
}

A.list <- A.gen.list




# create the interaction matrix with an average individual
Z.plant <- matrix(1, 1, 1)
colnames(Z.plant) <- sp

Z.pol <- diag(1, ncol(gen.output[[1]]), ncol(gen.output[[1]]))
rownames(Z.pol) <- colnames(gen.output[[1]])
colnames(Z.pol) <- colnames(gen.output[[1]])

## to collapse both plant individuals and pollinator species
Z.df <- cbind(c(rep(1, nrow(gen.output[[1]])), rep(0, ncol(gen.output[[1]]))),
              c(rep(0, nrow(gen.output[[1]])), rep(1, ncol(gen.output[[1]]))))
Z <- Z.df %>% as.matrix()

X.plant <- matrix(1/nrow(Z.df.list[[1]]), nrow=nrow(Z.df.list[[1]]))
rownames(X.plant) <- rownames(Z.df.list[[1]])
colnames(X.plant) <- sp
X.pol <- diag(1, ncol(int.mat), ncol(int.mat))
rownames(X.pol) <- colnames(int.mat)
colnames(X.pol) <- colnames(int.mat)

# PI matrix for plants (relative abundance proportional to flower number)
PI <- merge(X.plant, X.pol, by = "row.names", all = TRUE) %>% 
  column_to_rownames("Row.names") %>% replace(is.na(.), 0) %>% as.matrix()

# create PI matrix for pollinators (relative abundance proportional to number of visits)
PI.a <- as.matrix(cbind(c(1, rep(0, length(abpol.prop))), c(0, (abpol.prop))))

# interaction matrix
# within-species groups X within-species groups matrix

A <- diag(0, nrow(PI), nrow(PI)) 
rownames(A) <- rownames(PI)
colnames(A) <- rownames(PI)

# plant within-species group competition set to 1 above
nplant <- nrow(Z.plant)
ntotal <- ncol(A)
npol <- ntotal - nplant


A[(1+nplant):ntotal, 1:nplant] <- -0.2 # mutualistic benefits of plants given to pollinators (mean-field value)
A[(1+nplant):ntotal, (1+nplant):ntotal] <- 0.1 # pollinator interspecific competition (mean-field value)
A[1:nplant, (1+nplant):ntotal] <- -colMeans(int.mat) # mutualistic benefits of pollinators given to plants
for (row in rownames(A)){
  for (col in colnames(A)){
    if(row==col){
      A[row, col] <- 1 # intraspecific competition set to 1
    }
  }
}

A.list <- c(A.list, list(A))

}




#set a different color for each plant species
sp.color <- switch(focal,
                   "CLIB" = "#5DA6A7",
                   "HCOM" = "#F2C57C",
                   "HHAL" = "#E8875A")
sp.label <- switch(focal,
                   "CLIB" = ~italic("Cistus libanotis"),
                   "HCOM" = ~italic("Halimium calycinum"),
                   "HHAL" = ~italic("Halimium halimifolium"))

ggplot(int.mat.d, aes(x=degree)) +
  geom_histogram(binwidth = 1, alpha=0.6, fill=sp.color, color=sp.color)+
  labs(x="Individual plant degree", y = "Count")+
  theme_bw() + theme(axis.text=element_text(size=14),
                     axis.title=element_text(size=18))

