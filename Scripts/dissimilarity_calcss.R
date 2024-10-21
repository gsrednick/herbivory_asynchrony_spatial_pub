# Dissimilarity calculations ####

library(vegan)
library(tidyverse)



# Moorea ####
## MDS ####
MRA_community<-read.csv("./Data/Moorea_data/MRA_PC_analyze.csv")
community<-MRA_community

com_summarized <- community[c(15:36)]

env_summarized <- community[-c(15:36)]

com_summarized_mds <- metaMDS(comm = com_summarized, distance = "bray",
                              trace = FALSE, autotransform = TRUE, na.rm = FALSE)
com_summarized_mds$stress # 0.17


com_summarized_mds$species
com_summarized_mds_points<-data.frame(com_summarized_mds$points)
mds<-merge(env_summarized,com_summarized_mds_points, by="row.names", all.x=TRUE)

ambient_treats_nochange<-read.csv("./Data/Moorea_data/ambient_treatments.csv") # bring in background site treatments

mds_metadata<-merge(mds,ambient_treats_nochange)

# with envfit over the top
ef_sum <- envfit(com_summarized_mds, com_summarized, permu = 999)

# Extract scores
ef_sum.scrs <- as.data.frame(scores(ef_sum, display = "vectors")) # save species coords into dataframe
ef_sum.scrs <- cbind(ef_sum.scrs, species = rownames(ef_sum.scrs)) # add species names to dataframe
ef_sum.scrs <- cbind(ef_sum.scrs, pval = ef_sum$vectors$pvals) # add pvalues to dataframe so you can select species which are significant
ef_sum.scrs <- cbind(ef_sum.scrs, r2 = ef_sum$vectors$r) # add r2 to dataframe

vectors_sum <- subset(ef_sum.scrs, pval<=0.05)


# Write this to csv for plotting in general summary
write.csv(mds_metadata,"./Data/Moorea_data/MRA_mds_points.csv", row.names = F)
write.csv(vectors_sum,"./Data/Moorea_data/MRA_mds_vectors.csv", row.names = F)
write.csv(ef_sum.scrs,"./Data/Moorea_data/MRA_mds_vectors_all.csv", row.names = F)


## BC distance ####
MRA_community<-MRA_exp_data_corrected %>%
  filter(!meta == "CAGE",
         !zone == "B") # remove mid nutrient sites from this analysis

MRA_com_summarized <- MRA_community[c(15:36)]

MRA_env_summarized <- MRA_community[-c(15:36)]

MRA_com_bray<-MRA_com_summarized

MRA_dissim_df<- MRA_community %>%
  mutate(Tile_time = paste0(tile,sep = "_",time_point)) #%>%

row.names(MRA_com_bray)<-MRA_dissim_df$Tile_time

MRA_env_bray <- MRA_community[c(6:12)]

MRA_com_sim <- vegdist(MRA_com_bray, distance = "bray", trace = FALSE, autotransform = FALSE)
### IGNORE WARNING ###

MRA_distmat<-as.matrix(MRA_com_sim,labels=TRUE)

MRA_dist_df <- melt(as.matrix(MRA_distmat), varnames = c("TILE_1", "TILE_2"))
colnames(MRA_dist_df)<-c("TILE_1", "TILE_2","distance")

MRA_dist_df_2<-MRA_dist_df %>% filter(!TILE_1 == TILE_2)

# now I add both to a dataframe with factor labels in the same way I did in the "correlation" df from the MPA paper

MRA_dist_df_2$Tile_1 <- sub("_[^_]*", "", MRA_dist_df_2$TILE_1)
MRA_dist_df_2$Time_point_1 <- sub("[^_]*_", "", MRA_dist_df_2$TILE_1)

MRA_dist_df_2$Tile_2 <- sub("_[^_]*", "", MRA_dist_df_2$TILE_2)
MRA_dist_df_2$Time_point_2 <- sub("[^_]*_", "", MRA_dist_df_2$TILE_2)


MRA_dist_df_2<-MRA_dist_df_2 %>% tibble::rowid_to_column("ID")



# New separate dfs
MRA_Tile1_df<-MRA_dist_df_2[c(1,4,5,6)]
names(MRA_Tile1_df)<-c("ID","value","tile","time_point")

MRA_Tile2_df<-MRA_dist_df_2[c(1,4,7,8)]
names(MRA_Tile2_df)<-c("ID","value","tile","time_point")

# Merge

MRA_Tile_1_merged<-merge(MRA_env_bray,MRA_Tile1_df) # col 12 remains -- this was changed

MRA_Tile_2_merged<-merge(MRA_env_bray,MRA_Tile2_df) # col 12 remains -- this was changed

names(MRA_Tile_1_merged)
# Rename
names(MRA_Tile_1_merged)<-c("tile_1","time_point_1","site_1",
                        "zone_1","meta_1",
                        "com_treat_1","meta_treat_1","ID","value")

names(MRA_Tile_2_merged)<-c("tile_2","time_point_2","site_2",
                        "zone_2","meta_2",
                        "com_treat_2","meta_treat_2","ID","value")

MRA_dist_merged<-merge(MRA_Tile_1_merged,MRA_Tile_2_merged, by=c("ID","value"))


# manually provide zeros for first time point - no differences
MRA_zero_df<-MRA_dist_merged %>%
  filter(time_point_1 == "T2") %>%
  mutate(time_point_1 = "T1",
         time_point_2 = "T1",
         value = 0)

MRA_dist_merged_full <-rbind(MRA_dist_merged,MRA_zero_df)

# Get proper comparisons -- within metacommunity
MRA_dist_complete<-MRA_dist_merged_full %>%
  filter(
    site_1 == site_2,
    zone_1 == zone_2,
    meta_1 == meta_2,
    meta_treat_1 == meta_treat_2,
    time_point_1==time_point_2,
    !tile_1 == tile_2) %>%
  distinct(site_1,zone_1,meta_1,meta_treat_1,time_point_1, .keep_all = TRUE)

dim(MRA_dist_complete)

MRA_dist_complete_clean<-MRA_dist_complete %>%
  select(-contains("_2")) %>%
  rename_with(~str_remove(., "_1"), ends_with("_1")) %>%
  select(-c(com_treat,tile))

#dim(MRA_dist_complete_clean) # for checking
#dim(MRA_dist_complete) # for checking


# Bring in ambient background treatments and merge
MRA_ambient_treats<-read.csv("./Data/Moorea_data/ambient_treatments.csv")

MRA_dist_complete_treats<-merge(MRA_dist_complete_clean,MRA_ambient_treats)

MRA_dist_complete_treats$nutrient_treat_sp = factor(MRA_dist_complete_treats$nutrient_treat_sp, levels=c("low nutrients","mid nutrients","high nutrients"))
MRA_dist_complete_treats$herb_treat_sp = factor(MRA_dist_complete_treats$herb_treat_sp, levels=c("low herbivore","mid herbivore","high herbivore"))

# Make proper dates
MRA_time_points<-read.csv("./Data/Moorea_data/MRA_time_points.csv")
MRA_time_points$time_point_num<-MRA_time_points$Day_count
MRA_time_points<-MRA_time_points %>% filter(!zone == "B")

MRA_dist_complete_treats_time<-merge(MRA_dist_complete_treats,MRA_time_points, all= T)

#dim(MRA_dist_complete_treats_time) # for checking
#dim(MRA_dist_complete_treats) # for checking

MRA_cleaned_merged_dist<-MRA_dist_complete_treats_time

# Write csv for Moorea data
write.csv(MRA_cleaned_merged_dist,"./Data/Moorea_data/MRA_complete.csv",row.names = F)






# Port Phillip ####
## MDS ####
PPB_community<-read.csv('./Data/PPB_data/PPB_PC_analyze.csv')
PPB_com_summarized <- PPB_community[c(18:57)]
PPB_env_summarized <- PPB_community[-c(18:57)]

PPB_com_summarized_mds <- metaMDS(comm = PPB_com_summarized, distance = "bray",
                                  trace = FALSE, autotransform = TRUE, na.rm = FALSE)
PPB_com_summarized_mds$stress # 0.20

PPB_com_summarized_mds$species
PPB_com_summarized_mds_points<-data.frame(PPB_com_summarized_mds$points)
PPB_mds<-merge(PPB_env_summarized,PPB_com_summarized_mds_points, by="row.names", all.x=TRUE)




## plotting
# with envfit over the top
PPB_ef_sum <- envfit(PPB_com_summarized_mds, PPB_com_summarized, permu = 999)

PPB_ef_sum.scrs <- as.data.frame(scores(PPB_ef_sum, display = "vectors")) # save species intrinsic values into dataframe
PPB_ef_sum.scrs <- cbind(PPB_ef_sum.scrs, species = rownames(PPB_ef_sum.scrs)) # add species names to dataframe
PPB_ef_sum.scrs <- cbind(PPB_ef_sum.scrs, pval = PPB_ef_sum$vectors$pvals) # add pvalues to dataframe so you can select species which are significant
PPB_ef_sum.scrs <- cbind(PPB_ef_sum.scrs, r2 = PPB_ef_sum$vectors$r) # add r2 to dataframe
PPB_vectors_sum <- subset(PPB_ef_sum.scrs, pval<=0.05)

PPB_ambient_treats_nochange<-read.csv("./Data/PPB_data/ambient_treatments_PPB.csv")
PPB_mds_metadata<-merge(PPB_mds,PPB_ambient_treats_nochange)

# Write this to csv for assessment - KEEP
write.csv(PPB_mds_metadata,"./Data/PPB_data/PPB_mds_points.csv", row.names = F)
write.csv(PPB_vectors_sum,"./Data/PPB_data/PPB_mds_vectors.csv", row.names = F)
write.csv(PPB_ef_sum.scrs,"./Data/PPB_data/PPB_mds_vectors_all.csv", row.names = F)




## BC distance ####
PPB_com_summarized <- PPB_community[c(18:57)]

PPB_env_summarized <- PPB_community[-c(18:57)]
PPB_com_bray<-PPB_com_summarized

PPB_dissim_df<- PPB_community %>%
  mutate(Tile_time = paste0(tile_number,sep = "_",time_point)) #%>%

#duplicated(dissim_df$Tile_time)
PPB_dups<-PPB_dissim_df %>% group_by(Tile_time) %>% dplyr::summarize(n=n())
PPB_dups %>% filter(n > 1)

row.names(PPB_com_bray)<-PPB_dissim_df$Tile_time

PPB_env_bray <- PPB_community[c(6:17)]

PPB_com_sim <- vegdist(PPB_com_bray, distance = "bray", trace = FALSE, autotransform = FALSE) ### IGNORE WARNING ###

PPB_distmat<-as.matrix(PPB_com_sim,labels=TRUE)

PPB_dist_df <- melt(as.matrix(PPB_distmat), varnames = c("TILE_1", "TILE_2"))
colnames(PPB_dist_df)<-c("TILE_1", "TILE_2","distance")

PPB_dist_df_2<-PPB_dist_df %>% filter(!TILE_1 == TILE_2)


PPB_dist_df_2$Tile_1 <- sub("_[^_]*", "", PPB_dist_df_2$TILE_1)
PPB_dist_df_2$Time_point_1 <- sub("[^_]*_", "", PPB_dist_df_2$TILE_1)

PPB_dist_df_2$Tile_2 <- sub("_[^_]*", "", PPB_dist_df_2$TILE_2)
PPB_dist_df_2$Time_point_2 <- sub("[^_]*_", "", PPB_dist_df_2$TILE_2)


PPB_dist_df_2<-PPB_dist_df_2 %>% tibble::rowid_to_column("ID")



# New separate dfs
PPB_Tile1_df<-PPB_dist_df_2[c(1,4,5,6)]
names(PPB_Tile1_df)<-c("ID","value","tile","time_point")

PPB_Tile2_df<-PPB_dist_df_2[c(1,4,7,8)]
names(PPB_Tile2_df)<-c("ID","value","tile","time_point")

# Merge
PPB_env_bray$tile<-PPB_env_bray$tile_number

PPB_Tile_1_merged<-merge(PPB_env_bray,PPB_Tile1_df) # col 12 remains -- this was changed

PPB_Tile_2_merged<-merge(PPB_env_bray,PPB_Tile2_df) # col 12 remains -- this was changed

PPB_env_bray %>%
  filter(time_point== 1,site == "Jawbone", zone == "Deep", urchin_stat == "Barren") %>%
  #group_by(site,zone,urchin_stat,time_point,meta_treatment) %>%
  filter(metacom_treatment == 12)



PPB_Tile_1_merged$Name<-NULL
PPB_Tile_2_merged$Name<-NULL
PPB_Tile_1_merged$tile_number<-NULL
PPB_Tile_2_merged$tile_number<-NULL
PPB_Tile_1_merged$file<-NULL
PPB_Tile_2_merged$file<-NULL
PPB_Tile_1_merged$site_treatment<-NULL
PPB_Tile_2_merged$site_treatment<-NULL
PPB_Tile_1_merged$site_code<-NULL
PPB_Tile_2_merged$site_code<-NULL

names(PPB_Tile_1_merged)

# Rename
names(PPB_Tile_1_merged)<-c("time_point_1","tile_1","site_1","urchin_stat_1",
                            "zone_1","meta_number_1","com_treat_1","meta_treat_1",
                            "metacom_site_code_1",
                            "ID","value")

names(PPB_Tile_2_merged)<-c("time_point_2","tile_2","site_2","urchin_stat_2",
                            "zone_2","meta_number_2","com_treat_2","meta_treat_2",
                            "metacom_site_code_2",
                            "ID","value")

PPB_dist_merged<-merge(PPB_Tile_1_merged,PPB_Tile_2_merged, by=c("ID","value"))


names(PPB_dist_merged)


PPB_dist_merged_grouped_2<-PPB_dist_merged %>%
  filter(metacom_site_code_1 == metacom_site_code_2,
         meta_number_1 == meta_number_2,
         meta_treat_1 == meta_treat_2,
         time_point_1 == time_point_2) %>%
  distinct(time_point_1,meta_treat_1,meta_treat_1,metacom_site_code_1, .keep_all = TRUE)

dim(PPB_dist_merged_grouped_2)

PPB_dist_check_reps<-PPB_dist_merged_grouped_2 %>%
  group_by(site_1,zone_1,urchin_stat_1,time_point_1,meta_treat_1) %>%
  dplyr::summarize(count = n())

# Add treatments

PPB_ambient_treats_org<-read.csv("./Data/PPB_data/ambient_treatments_PPB.csv")
head(PPB_ambient_treats_org)

head(PPB_dist_merged_grouped_2)

PPB_cleaned_merged_dist<-PPB_dist_merged_grouped_2 %>%
  select(-contains("_2")) %>%
  rename_with(~str_remove(., "_1"), ends_with("_1")) %>%
  select(-c(com_treat,tile))




PPB_exp_data_expanded_meta<-merge(PPB_cleaned_merged_dist,PPB_ambient_treats_org)

PPB_check_reps<-PPB_cleaned_merged_dist %>%
  group_by(site,zone,urchin_stat,time_point,meta_treat) %>%
  dplyr::summarize(count = n())

# add zeros for initial
PPB_zero_df<-PPB_exp_data_expanded_meta %>%
  filter(time_point == "1") %>%
  mutate(time_point = "0",
         value = 0)

PPB_dist_merged_full_pre <-rbind(PPB_exp_data_expanded_meta,PPB_zero_df)

# Make proper dates
PPB_time_points<-read.csv("./Data/PPB_data/PPB_time_points.csv")

PPB_dist_merged_full<-merge(PPB_dist_merged_full_pre,PPB_time_points)
dim(PPB_dist_merged_full_pre)
dim(PPB_dist_merged_full)

# Rename for consistency
PPB_dist_merged_full$time_point_num <- PPB_dist_merged_full$Day_count

PPB_dist_merged_full$nutrient_treat_sp = factor(PPB_dist_merged_full$nutrient_treat_sp, levels=c("low nutrients","mid nutrients","high nutrients"))
PPB_dist_merged_full$herb_treat_sp = factor(PPB_dist_merged_full$herb_treat_sp, levels=c("low herbivore","mid herbivore","high herbivore"))


# Write csv for PPB data
write.csv(PPB_dist_merged_full,"./Data/PPB_data/PPB_complete.csv",row.names = F)




# END #
