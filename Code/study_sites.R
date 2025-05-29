# check differences in plant individual specialization across sites


sp.names <- c("CLIB", "HCOM", "HHAL")


sp.site.list <- list()

for(name in sp.names){
  focal <- name
  
  #set a different color for each plant species
  sp.color <- switch(focal,
                     "CLIB" = "#5DA6A7",
                     "HCOM" = "#F2C57C",
                     "HHAL" = "#E8875A")
  sp.label <- switch(focal,
                     "CLIB" = ~italic("Cistus libanotis"),
                     "HCOM" = ~italic("Halimium calycinum"),
                     "HHAL" = ~italic("Halimium halimifolium"))
  
  # construct interaction matrices
  source('Code/interaction_matrices.R')
  
  plot.id <- read.csv("Data/all_plant_traits.csv") %>%
    select("Plant_id", "Plot")
  degree.id <- int.mat.d %>% tibble::rownames_to_column("Plant_id") %>%
    merge(plot.id, by="Plant_id") 
  id.all <- degree.id  %>%
    select("Plant_id", "degree", "Plot")
  
  sp.site.plot <- ggplot(id.all, aes(x = Plot, y = degree, fill = Plot)) +
    geom_boxplot(alpha=0.7, color=sp.color, fill=sp.color) +
    theme_bw() +
    ylab("Plant individual degree") + xlab("Study plot") +
    ggtitle(sp.label) + 
    theme(text=element_text(size = 16),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
          plot.title = element_text(size=16, hjust = 0.5)) 
  
  
  anova_result <- aov(degree ~ Plot, data = id.all)
  summary(anova_result)

  # save plots
  sp.site.list[[name]] <- sp.site.plot
  
}

sp.site.list[[1]]  + xlab("") + 
  sp.site.list[[2]]  + ylab("")  + 
  sp.site.list[[3]]  + xlab("") + ylab("") 
  


