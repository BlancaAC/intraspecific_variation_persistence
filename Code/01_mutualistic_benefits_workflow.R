library(bayesplot)
library(bipartite)
library(dplyr)
library(factoextra)
library(ggdist) #needs version 3.0.0
library(ggExtra)
library(ggplot2) #needs version 3.4.2
library(ggtern)
library(ggpmthemes)
library(magrittr)
library(Matrix)
library(mvtnorm)
library(patchwork) #needs version 1.1.2
library(plyr)
library(posterior)
library(pracma)
library(rlang)
library(Rmisc)
library(scales)
library(textshape)
library(tibble)
library(tidybayes)
library(tidyverse)


df <- list.files(path = "Data/Bernoulli_fruitset", pattern = "nseeds_fruitset_", full.names = T, recursive = TRUE) 
dd.list <- lapply(df, read.csv, sep=' ') 
names(dd.list) <- sapply(
	df,
	function(x){
		nn <- strsplit(x, '_')[[1]]
		nn[4]
	}
)

dd.fl <- read.csv("Data/flowering_weeks_plants.csv", sep=';') %>% select(Plant_id, Flowering_min)

for (i in seq_along(dd.list)){
  temp <- dd.list[[i]] 
  temp <- temp %>% mutate(Plant_id=as.character(Plant_id)) %>% left_join(dd.fl, by="Plant_id")
  
  cn <- colnames(temp)
  pollinators <- cn[which(!cn %in% c('Plant_id','Plant_sp','Response',
                                   'Plot','N_flowers','Type', 'Flowering_min') )]

  temp <- temp %>% mutate(degree= rowSums(select(.,all_of(pollinators))!=0))
  temp <- temp %>% filter(degree !=0)
  dd.list[[i]] <- temp
}

# perform the model fitting process for each focal species to estimate the mutualistic benefits
dd.list <- dd.list[c("CLIB", "HCOM", "HHAL")]
  
for(i in names(dd.list)){
	focal <- i

	# prep the data for model fitting
	source('Code/setup_data.R')

	# fit things with cmdstanr and generate some pretty figures
	source('Code/fit_joint_model_stan.R')
}


