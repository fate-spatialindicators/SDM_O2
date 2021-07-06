target_species <- c("dover sole", "sablefish", "lingcod", "petrale sole", "longspine thornyhead","pacific hake")
constraining_species <- c("canary rockfish", "widow rockfish", "darkblotched rockfish", "cowcod")
dat = dplyr::filter(dat, species%in% c(target_species, constraining_species), year%in%seq(2010,2015), 
                    !is.na(temp), !is.na(o2), !is.na(sal),!is.infinite(sal),
                    latitude_dd > min(latitude_dd[which(cpue_kg_km2>0)]),
                    latitude_dd <= max(latitude_dd[which(cpue_kg_km2>0)]),
                    longitude_dd > min(longitude_dd[which(cpue_kg_km2>0)]),
                    longitude_dd < max(longitude_dd[which(cpue_kg_km2>0)]))

wide.dat <- tidyr::spread(data = dat, key = species, value = cpue_kg_km2, fill = 0)

# extract out species
species <- wide.dat[,22:31]

species <- species - colMeans(species)
species <- species / apply(species, MARGIN = 2, sd)

corr.mat <- cor(species, method = "spearman")
