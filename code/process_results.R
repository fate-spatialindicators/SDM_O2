library(ggplot2)
library(dplyr)

region=c("goa","wc")[2]
df = readRDS(paste0("output/",region,"/models.RDS"))

for(i in 1:nrow(df)) {
  fname = paste0("output/",region,"/model_",i,".rds")
  if(file.exists(fname)) {
    
  m = readRDS(fname)
  sd_report <- summary(m$sd_report)
  params <- as.data.frame(sd_report[grep("quadratic", row.names(sd_report)), ])
  df$low[i] = params["quadratic_low","Estimate"]
  df$low_se[i] = params["quadratic_low","Std. Error"]
  df$hi[i] = params["quadratic_hi","Estimate"]
  df$hi_se[i] = params["quadratic_hi","Std. Error"]
  df$range[i] = params["quadratic_range","Estimate"]
  df$range_se[i] = params["quadratic_range","Std. Error"]
  df$peak[i] = params["quadratic_peak","Estimate"]
  df$peak_se[i] = params["quadratic_peak","Std. Error"]
  df$reduction[i] = params["quadratic_reduction","Estimate"]
  df$reduction_se[i] = params["quadratic_reduction","Std. Error"]  
  }
}

# remove junk columns
df = dplyr::select(df, -spatial_only, -time_varying)

# save results
write.csv(df, file=paste0("output/",region,"_output.csv"))
#saveRDS(df,file=paste0("output/",region,"_output.rds"))

plot(df$range[which(df$depth_effect==TRUE & df$covariate=="o2")],
  df$range[which(df$depth_effect==FALSE & df$covariate=="o2")],
  xlim=c(0,20),xlab="Depth included",ylab="Depth not included",
  main="o2")

plot(df$range[which(df$depth_effect==TRUE & df$covariate=="temp")],
  df$range[which(df$depth_effect==FALSE & df$covariate=="temp")],
  xlab="Depth included",ylab="Depth not included",
  main="temp")

df = df[complete.cases(df),]

pdf(paste0("plots/",region,"-temp_range.pdf"))
level_order = dplyr::filter(df, !is.na(range), covariate=="temp",
  depth_effect == TRUE, range_se < 5) %>%
  dplyr::arrange(range) %>% select(species)
dplyr::filter(df, !is.na(range), covariate=="temp", 
  range_se < 5,depth_effect == TRUE) %>% 
ggplot(aes(factor(species, level=level_order$species), range)) +
  geom_pointrange(aes(ymin=range-2*range_se, 
    ymax=range+2*range_se)) +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Temperature range - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-temp_range2.pdf"))
level_order = dplyr::filter(df, !is.na(range), covariate=="temp",
  depth_effect == TRUE, range_se < 1) %>%
  dplyr::arrange(range) %>% select(species)
dplyr::filter(df, !is.na(range), covariate=="temp", range_se < 1,
  species %in% level_order$species) %>% 
  ggplot(aes(factor(species, level=level_order$species), range,col=depth_effect)) +
  geom_pointrange(aes(ymin=range-2*range_se, 
    ymax=range+2*range_se)) +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Temperature range - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-temp_low.pdf"))
level_order = dplyr::filter(df, !is.na(low), covariate=="temp",
  depth_effect == TRUE,low_se<5) %>%
  dplyr::arrange(low) %>% select(species)
dplyr::filter(df, !is.na(low), covariate=="temp", depth_effect==TRUE,low_se < 5) %>% 
  ggplot(aes(factor(species, level=level_order$species), low)) +
  geom_pointrange(aes(ymin=low-2*low_se, 
    ymax=low+2*low_se),col="darkblue") +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Temperature low bound - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-temp_hi.pdf"))
level_order = dplyr::filter(df, !is.na(hi), covariate=="temp",
  depth_effect == TRUE, hi_se < 5) %>%
  dplyr::arrange(hi) %>% select(species)
dplyr::filter(df, !is.na(hi), covariate=="temp", depth_effect==TRUE,hi_se < 5) %>% 
  ggplot(aes(factor(species, level=level_order$species), hi)) +
  geom_pointrange(aes(ymin=hi-2*hi_se, 
    ymax=hi+2*hi_se),col="darkblue") +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Temperature upper bound - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-temp_reduction.pdf"))
level_order = dplyr::filter(df, !is.na(reduction), covariate=="temp",
  depth_effect == TRUE, reduction_se < 5) %>%
  dplyr::arrange(reduction) %>% select(species)
dplyr::filter(df, !is.na(reduction), covariate=="temp", depth_effect==TRUE,reduction_se < 5) %>% 
  ggplot(aes(factor(species, level=level_order$species), reduction)) +
  geom_pointrange(aes(ymin=reduction-2*reduction_se, 
    ymax=reduction+2*reduction_se),col="darkblue") +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Reduction from peak - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-temp_reduction2.pdf"))
level_order = dplyr::filter(df, !is.na(reduction), covariate=="temp",
  depth_effect == TRUE, reduction_se < 5) %>%
  dplyr::arrange(reduction) %>% select(species)
dplyr::filter(df, !is.na(reduction), covariate=="temp", reduction_se < 5,
  species %in% level_order$species) %>% 
  ggplot(aes(factor(species, level=level_order$species), reduction,col=depth_effect)) +
  geom_pointrange(aes(ymin=reduction-2*reduction_se, 
    ymax=reduction+2*reduction_se)) +
  coord_flip() + xlab("Species") + ylab("Reduction") + 
  ggtitle(paste0("Temperature range - ",region," survey"))
dev.off()

# look at empirical depth bands versus these estimated ranges.
# for wc species, some spp of rockfish have a big difference between 
# whether depth is included or not


# o2 plots here - wc and bc species
pdf(paste0("plots/",region,"-o2_range.pdf"))
level_order = dplyr::filter(df, !is.na(range), covariate=="o2",
  depth_effect == TRUE, range_se < 1) %>%
  dplyr::arrange(range) %>% select(species)
dplyr::filter(df, !is.na(range), covariate=="o2", range_se < 1,depth_effect == TRUE) %>% 
  ggplot(aes(factor(species, level=level_order$species), range)) +
  geom_pointrange(aes(ymin=range-2*range_se, 
    ymax=range+2*range_se),col="darkblue") + ylim(0,10) +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("o2 range - ",region," survey"))
dev.off()


pdf(paste0("plots/",region,"-o2_range2.pdf"))
level_order = dplyr::filter(df, !is.na(range), covariate=="o2",
  depth_effect == TRUE, range_se < 1) %>%
  dplyr::arrange(range) %>% select(species)
dplyr::filter(df, !is.na(range), covariate=="o2", range_se < 1,
  species %in% level_order$species) %>% 
  ggplot(aes(factor(species, level=level_order$species), range,col=depth_effect)) +
  geom_pointrange(aes(ymin=range-2*range_se, 
    ymax=range+2*range_se)) +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("o2 range - ",region," survey"))
dev.off()


pdf(paste0("plots/",region,"-o2_low.pdf"))
level_order = dplyr::filter(df, !is.na(low), covariate=="o2",
  depth_effect == TRUE, low_se < 1) %>%
  dplyr::arrange(low) %>% select(species)
dplyr::filter(df, !is.na(low), covariate=="o2", low_se < 1,depth_effect == TRUE) %>% 
  ggplot(aes(factor(species, level=level_order$species), low)) +
  geom_pointrange(aes(ymin=low-2*low_se, 
    ymax=low+2*low_se),col="darkblue") +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("o2 low bound - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-o2_hi.pdf"))
level_order = dplyr::filter(df, !is.na(hi), covariate=="o2",
  depth_effect == TRUE, hi_se < 1) %>%
  dplyr::arrange(hi) %>% select(species)
dplyr::filter(df, !is.na(hi), covariate=="o2", hi_se < 1,depth_effect == TRUE) %>% 
  ggplot(aes(factor(species, level=level_order$species), hi)) +
  geom_pointrange(aes(ymin=hi-2*hi_se, 
    ymax=hi+2*hi_se),col="darkblue") +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("o2 upper bound - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-o2_reduction.pdf"))
level_order = dplyr::filter(df, !is.na(hi), covariate=="o2",
  depth_effect == TRUE, reduction_se < 5) %>%
  dplyr::arrange(reduction) %>% select(species)
dplyr::filter(df, !is.na(reduction), covariate=="o2", depth_effect==TRUE,reduction_se < 5) %>% 
  ggplot(aes(factor(species, level=level_order$species), reduction)) +
  geom_pointrange(aes(ymin=reduction-2*reduction_se, 
    ymax=reduction+2*reduction_se),col="darkblue") +
  coord_flip() + xlab("Species") + ylab("Range") + 
  ggtitle(paste0("Reduction from peak - ",region," survey"))
dev.off()

pdf(paste0("plots/",region,"-o2_reduction2.pdf"))
level_order = dplyr::filter(df, !is.na(reduction), covariate=="o2",
  depth_effect == TRUE, reduction_se < 5) %>%
  dplyr::arrange(reduction) %>% select(species)
dplyr::filter(df, !is.na(reduction), covariate=="o2", reduction_se < 5,
  species %in% level_order$species) %>% 
  ggplot(aes(factor(species, level=level_order$species), reduction,col=depth_effect)) +
  geom_pointrange(aes(ymin=reduction-2*reduction_se, 
    ymax=reduction+2*reduction_se)) +
  coord_flip() + xlab("Species") + ylab("Reduction") + 
  ggtitle(paste0("o2 range - ",region," survey"))
dev.off()