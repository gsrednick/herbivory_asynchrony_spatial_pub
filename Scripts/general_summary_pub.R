# General summary script

# In support of:
# Habitat attributes mediate top-down and bottom-up drivers of community development in temperate and tropical algae
# in Ecosphere XXXX
# Authors: Griffin Srednick & Stephen Swearer

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(ggh4x)
library(emmeans)
library(cowplot)
library(patchwork)
library(magick)
library(DHARMa)
library(ggeffects)
library(ggrepel)

# Data import
MRA_community<-read.csv("./Data/Moorea_data/MRA_PC_analyze.csv")
PPB_community<-read.csv("./Data/PPB_data/PPB_PC_analyze.csv")

# Add Icons for figs
FH_NH_image<-magick::image_read("./Figures/inset_text/FH_NH.png") %>%
  magick::image_background("none") %>%
  grid::rasterGrob()

FL_NH_image<-magick::image_read("./Figures/inset_text/FL_NH.png") %>%
  magick::image_background("none") %>%
  grid::rasterGrob()

FH_NL_image<-magick::image_read("./Figures/inset_text/FH_NL.png") %>%
  magick::image_background("none") %>%
  grid::rasterGrob()

FL_NL_image<-magick::image_read("./Figures/inset_text/FL_NL.png") %>%
  magick::image_background("none") %>%
  grid::rasterGrob()



icon_markdown <- paste0("<img src=", list.files("flags/", full.names = TRUE), " width='100'/>")

image_to_grob <- function(image) {
  raster_grob <- grid::rasterGrob(image = image[[1]], interpolate = FALSE)
  return(raster_grob)
}

# Convert images to grobs
FH_NH_grob <- image_to_grob(FH_NH_image)
FL_NH_grob <- image_to_grob(FL_NH_image)
FH_NL_grob <- image_to_grob(FH_NL_image)
FL_NL_grob <- image_to_grob(FL_NL_image)

# Create a list of grobs and corresponding labels
facet_labels <- list(
  "FH_NH" = FH_NH_grob,
  "FL_NH" = FL_NH_grob,
  "FH_NL" = FH_NL_grob,
  "FL_NL" = FL_NL_grob
)

# Create a function to return grob based on label
custom_labeller <- function(variable, value) {
  return(facet_labels[value])
}


# Tropical - Moorea ####
## Cover and diversity ####
MRA_community_sum<-MRA_community %>%
  mutate(total_cover = rowSums(.[16:36]),
         H = diversity(.[16:36]))

ambient_treats_org<-read.csv("./Data/Moorea_data/ambient_treatments.csv")
MRA_community_meta<-merge(MRA_community_sum,ambient_treats_org) %>%
  select(-c(3:7,15:36))


# Add zero time point
MRA_community_meta_zero_point<-MRA_community_meta %>%
  filter(time_point == "T2") %>%
  mutate(time_point = "T1",
         total_cover = 0,
         H = 0)

MRA_merged_full <-rbind(MRA_community_meta,MRA_community_meta_zero_point)

# Add continuous time
MRA_time_points<-read.csv("./Data/Moorea_data/MRA_time_points.csv")
MRA_merged_full<-merge(MRA_merged_full,MRA_time_points)

MRA_merged_full$com_treat<-plyr::revalue(MRA_merged_full$com_treat, c("C" = "Complex",
                                                                      "S"="Simple"))



MRA_merged_full$nutrient_treat_sp = factor(MRA_merged_full$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
MRA_merged_full$herb_treat_sp = factor(MRA_merged_full$herb_treat_sp, levels=c("high herbivore","low herbivore"))


MRA_merged_plot<-MRA_merged_full


## LMER ####
MRA_merged_lmer<-MRA_merged_plot %>% filter(!Day_count =="0") # remove first time point for models -- irrelevant here

tropical_rep<-MRA_merged_lmer %>%
  group_by(com_treat,herb_treat_sp,nutrient_treat_sp,time_point) %>%
  dplyr::summarize(count = n())

# assumptions
logitTransform <- function(p) { log(p/(1.01-p) + 1) } # Transformation for percent cover

MRA_merged_lmer$total_cover_trans<-logitTransform(MRA_merged_lmer$total_cover/100)

qqnorm(MRA_merged_lmer$total_cover_trans) # yes, normal


# Cover - over time
MRA_cover_time_lmer<-lmer(total_cover_trans ~ time_point * com_treat * herb_treat_sp * nutrient_treat_sp + (1|meta/site_treatment),
                          data = MRA_merged_lmer, REML = T)

summary(MRA_cover_time_lmer)
anova(MRA_cover_time_lmer)
ranova(MRA_cover_time_lmer)
emmeans(MRA_cover_time_lmer, list(pairwise ~ time_point * herb_treat_sp * nutrient_treat_sp), adjust = "tukey")

# model assumptions
MRA_cover_time_lmer_fit<-simulateResiduals(fittedModel = MRA_cover_time_lmer, plot = T)
testDispersion(MRA_cover_time_lmer_fit) # good
plot(resid(MRA_cover_time_lmer),MRA_merged_lmer$total_cover_trans)

plot(ggpredict(MRA_cover_time_lmer, terms = c("herb_treat_sp","nutrient_treat_sp")))


# Diversity - over time
MRA_div_time_lmer<-lmer(H ~ time_point * com_treat * herb_treat_sp * nutrient_treat_sp + (1|meta/site_treatment), data = MRA_merged_lmer)
summary(MRA_div_time_lmer)
anova(MRA_div_time_lmer)
ranova(MRA_div_time_lmer)

# model assumptions
MRA_div_time_lmer_fit<-simulateResiduals(fittedModel = MRA_div_time_lmer, plot = F)
plot(MRA_div_time_lmer_fit)
testDispersion(MRA_div_time_lmer_fit) # good

plot(ggpredict(MRA_div_time_lmer, terms = c("time_point","com_treat","herb_treat_sp","nutrient_treat_sp")))


# Plot these
MRA_merged_plot$combined_treat <-paste0(MRA_merged_plot$herb_treat_sp,sep = "_",MRA_merged_plot$nutrient_treat_sp)

MRA_merged_plot$nutrient_treat_sp = factor(MRA_merged_plot$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
MRA_merged_plot$herb_treat_sp = factor(MRA_merged_plot$herb_treat_sp, levels=c("low herbivore","high herbivore"))

# Object from ggpredict
MRA_cover_lmer_predict<-ggpredict(MRA_cover_time_lmer, terms = c("time_point","com_treat","herb_treat_sp","nutrient_treat_sp"))
MRA_cover_lmer_predict_df<-as.data.frame(MRA_cover_lmer_predict)
head(MRA_cover_lmer_predict_df)
names(MRA_cover_lmer_predict_df) <- c("time_point","total_cover","std.error","conf.low","conf.high","com_treat","herb_treat_sp","nutrient_treat_sp")

# Make metadata dataframe
MRA_groups<-MRA_merged_plot %>%
  group_by(com_treat,herb_treat_sp,nutrient_treat_sp,time_point,Day_count) %>%
  reframe()

# Merge ggpredict object with metadata
MRA_cover_lmer_predict_df_complete<-merge(MRA_cover_lmer_predict_df,MRA_groups, all = T) %>%
  mutate(total_cover = if_else(Day_count == 0, 0, total_cover)) # add zero initial point

# Fix ordering
MRA_cover_lmer_predict_df_complete$nutrient_treat_sp = factor(MRA_cover_lmer_predict_df_complete$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
MRA_cover_lmer_predict_df_complete$herb_treat_sp = factor(MRA_cover_lmer_predict_df_complete$herb_treat_sp, levels=c("low herbivore","high herbivore"))

MRA_cover_lmer_predict_df_complete<-MRA_cover_lmer_predict_df_complete %>%
  mutate(com_treat_het = recode(com_treat,
                                        "Simple" = "high access",
                                        "Complex" = "low access"))

MRA_merged_plot<-MRA_merged_plot %>%
  mutate(com_treat_het = recode(com_treat,
                                "Simple" = "high access",
                                "Complex" = "low access"))

MRA_cover_lmer_predict_df_complete$com_treat_het = factor(MRA_cover_lmer_predict_df_complete$com_treat_het, levels=c("high access","low access"))
MRA_merged_plot$com_treat_het = factor(MRA_merged_plot$com_treat_het, levels=c("high access","low access"))


# replicates
MRA_com_reps<-MRA_merged_plot %>%
  group_by(Day_count,com_treat_het,herb_treat_sp,nutrient_treat_sp) %>%
  dplyr::summarize(count = n())

# At what point is cover greater than 95%?
MRA_cover_estimate_standard <- MRA_merged_plot %>%
  group_by(Day_count, com_treat_het, herb_treat_sp, nutrient_treat_sp) %>%
  dplyr::summarize(mean_cover = mean(total_cover)) %>%
  arrange(Day_count) %>%
  filter(mean_cover > 95) %>%
  group_by(com_treat_het, herb_treat_sp, nutrient_treat_sp) %>%
  slice(1)

MRA_cover_estimate<-MRA_cover_estimate_standard
write.csv(MRA_cover_estimate_standard,"./Data/Moorea_data/MRA_sat_point.csv", row.names = F)


MRA_cover_estimate_standard

MRA_cover_plot<-ggplot(MRA_cover_lmer_predict_df_complete, aes(x = Day_count, y = total_cover, color = com_treat_het)) +
  geom_line() +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_color_manual(values = c("black","magenta")) +
  theme_bw()+
  facet_wrap(herb_treat_sp~nutrient_treat_sp, ncol = 4, labeller = as_labeller(facet_labels)) +
  labs(x = "Time (d)", y = "Algal cover \n(logit transformed)", color = "Access")  +
  guides(color = "none", fill = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_text(data = MRA_merged_plot %>%
              group_by(herb_treat_sp,nutrient_treat_sp) %>%
              dplyr::summarise(),
            aes(label = letters[1:4]), x = -Inf, y = Inf, size = 4, color = "black", hjust = -0.5, vjust = 1) +
  geom_text_repel(data = MRA_merged_plot %>%
              group_by(Day_count,com_treat_het,herb_treat_sp,nutrient_treat_sp) %>%
                filter(!Day_count == 0) %>%
              dplyr::summarise(mean_cover = mean(logitTransform(total_cover/100),na.rm =T),
                        count = n()),
            aes(label = count, x = Day_count, y = mean_cover, color = com_treat_het),
            min.segment.length = Inf,
            size = 3) +
  geom_segment(data = MRA_cover_estimate_standard, aes(x = Day_count, y = 0.1, yend = 0.8, color = com_treat_het),
               size = 1,
               arrow = arrow(length = unit(0.2, "cm")),
               position = position_dodge(width = 20),
               lineend = "butt", linejoin = "mitre",
               show_guide = F) +
  geom_point(data = MRA_merged_plot, aes(x = Day_count, y = logitTransform(total_cover/100), color = com_treat_het), alpha = 0.1, inherit.aes = F)



FH_NH_icon<-ggdraw(FH_NH_image)
FL_NH_icon<-ggdraw(FL_NH_image)
FH_NL_icon<-ggdraw(FH_NL_image)
FL_NL_icon<-ggdraw(FL_NL_image)

icon_facet <-FL_NL_icon + FL_NH_icon + FH_NL_icon + FH_NH_icon + plot_layout(ncol = 4)


# Diversity

# Object from ggpredict
MRA_div_lmer_predict<-ggpredict(MRA_div_time_lmer, terms = c("time_point","com_treat","herb_treat_sp","nutrient_treat_sp"))
MRA_div_lmer_predict_df<-as.data.frame(MRA_div_lmer_predict)
head(MRA_div_lmer_predict_df)
names(MRA_div_lmer_predict_df) <- c("time_point","H","std.error","conf.low","conf.high","com_treat","herb_treat_sp","nutrient_treat_sp")

# Merge ggpredict object with metadata
MRA_div_lmer_predict_df_complete<-merge(MRA_div_lmer_predict_df,MRA_groups, all = T) %>%
  mutate(H = if_else(Day_count == 0, 0, H)) # add zero initial point

# Fix ordering
MRA_div_lmer_predict_df_complete$nutrient_treat_sp = factor(MRA_div_lmer_predict_df_complete$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
MRA_div_lmer_predict_df_complete$herb_treat_sp = factor(MRA_div_lmer_predict_df_complete$herb_treat_sp, levels=c("low herbivore","high herbivore"))

MRA_div_lmer_predict_df_complete<-MRA_div_lmer_predict_df_complete %>%
  mutate(com_treat_het = recode(com_treat,
                                "Simple" = "high access",
                                "Complex" = "low access"))

MRA_div_lmer_predict_df_complete$com_treat_het = factor(MRA_div_lmer_predict_df_complete$com_treat_het, levels=c("high access","low access"))

MRA_div_plot<-ggplot(MRA_div_lmer_predict_df_complete, aes(x = Day_count, y = H, color = com_treat_het)) +
  geom_line() +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_color_manual(values = c("black","magenta")) +
  theme_bw()+
  facet_wrap(herb_treat_sp~nutrient_treat_sp, ncol = 4, labeller = as_labeller(facet_labels)) +
  labs(x = "Time (d)", y = "H'", color = "Access")  +
  guides(color = "none", fill = "none") +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_text(data = MRA_merged_plot %>%
              group_by(herb_treat_sp,nutrient_treat_sp) %>%
              dplyr::summarise(),
            aes(label = letters[5:8]), x = -Inf, y = Inf, size = 4, color = "black", hjust = -0.5, vjust = 1) +
  geom_point(data = MRA_merged_plot, aes(x = Day_count, y = H, color = com_treat_het), alpha = 0.1, inherit.aes = F)



# over time
MRA_cover_div_time_plot<-MRA_cover_plot / MRA_div_plot


ggsave(filename = "./Figures/MRA_summary_time_fig.pdf",
       plot = MRA_cover_div_time_plot,
       dpi = 300,
       width = 8,
       height = 4)




# Temperate - Port Phillip ####
PPB_community<-read.csv('./Data/PPB_data/PPB_PC_analyze.csv')

PPB_community_sum<-PPB_community %>%
  #cbind(PPB_com_summarized_filt) %>%
  select(-18) %>% # remove BARE space
  filter(!zone == "Shallow") %>%
  mutate(total_cover = rowSums(.[18:56]),
         H = diversity(.[18:56]))

#PPB_community_sum<-PPB_community %>%
#  select(-18) %>% # remove BARE space
#  filter(!zone == "Shallow") %>%
#  mutate(total_cover = rowSums(.[18:55]),
#         H = diversity(.[18:55]))

PPB_ambient_treats_nochange<-read.csv("./Data/PPB_data/ambient_treatments_PPB.csv")
PPB_time_points<-read.csv("./Data/PPB_data/PPB_time_points.csv")

PPB_community_complete<-merge(PPB_community_sum,PPB_ambient_treats_nochange) %>%
  merge(PPB_time_points)

PPB_community_meta<-PPB_community_complete %>%
  select(-c(5:9,14,18:54))


# check for duplicates?
PP_rep_check<-PPB_community_meta %>%
  filter(site == "Mornington") %>%
  group_by(site,urchin_stat,zone,nutrient_treat_sp,herb_treat_sp,time_point) %>%
  dplyr::summarize(count = n())

MRA_rep_check<-MRA_community_meta %>%
  group_by(site,nutrient_treat_sp,herb_treat_sp,time_point) %>%
  dplyr::summarize(count = n())



# Add zero time point
PPB_community_meta_zero_point<-PPB_community_meta %>%
  filter(time_point == "1") %>%
  mutate(time_point = "0",
         Day_count = 0,
         total_cover = 0,
         H = 0)

PPB_community_full_merge<-PPB_community_complete %>%
  select(all_of(intersect(colnames(PPB_community_complete), colnames(PPB_community_meta_zero_point))))

PPB_merged_full <-rbind(PPB_community_full_merge,PPB_community_meta_zero_point)


PPB_merged_full$nutrient_treat_sp = factor(PPB_merged_full$nutrient_treat_sp, levels=c("low nutrients","mid nutrients","high nutrients"))
PPB_merged_full$herb_treat_sp = factor(PPB_merged_full$herb_treat_sp, levels=c("high herbivore","mid herbivore","low herbivore"))


PPB_merged_plot<-PPB_merged_full %>% filter(!nutrient_treat_sp == "mid nutrients",
                                            !herb_treat_sp == "mid herbivore")





## LMER ####
PPB_merged_lmer<-PPB_merged_plot %>% filter(!Day_count =="0") # There appeared to be an error before but now this is fine

# assumptions
#qqnorm(asinTransform(PPB_merged_lmer$total_cover/100)) # needs work
#qqnorm(PPB_merged_lmer$total_cover_trans) # needs work

qqnorm(PPB_merged_lmer$H) # Acceptable
logitTransform <- function(p) { log(p/(1-p + 0.01) + 1) } # Transformation for percent cover

PPB_merged_lmer$total_cover_trans<-logitTransform(PPB_merged_lmer$total_cover/100)


qqnorm(PPB_merged_lmer$total_cover_trans) # needs work

PPB_merged_lmer %>% group_by(Day_count,time_point) %>% dplyr::summarize(count = n())

PPB_merged_lmer %>%
  filter(site == "Mornington") %>%
  group_by(urchin_stat,zone,nutrient_treat_sp,herb_treat_sp,time_point) %>%
  dplyr::summarize(count = n())

temperate_rep<-PPB_merged_plot %>%
  group_by(complex_treatment,herb_treat_sp,nutrient_treat_sp,Day_count) %>%
  dplyr::summarize(count = n())

temperate_rep<-PPB_merged_plot %>%
  group_by(complex_treatment,herb_treat_sp,nutrient_treat_sp,Day_count,metacom_treatment) %>%
  dplyr::summarize(count = n())

# Cover - basic
PPB_cover_basic_lmer<-lmer(total_cover ~ complex_treatment * herb_treat_sp * nutrient_treat_sp + (1|metacom_treatment/site_treatment), data = PPB_merged_lmer)
summary(PPB_cover_basic_lmer)
anova(PPB_cover_basic_lmer)
emmeans(PPB_cover_basic_lmer, list(pairwise ~ herb_treat_sp * nutrient_treat_sp), adjust = "tukey")

# Cover - over time
#PPB_cover_time_lmer<-lmer(total_cover ~ Day_count * complex_treatment * herb_treat_sp * nutrient_treat_sp + (1|metacom_treatment/site_treatment), data = PPB_merged_lmer)
PPB_cover_time_lmer<-lmer(total_cover_trans ~ as.factor(time_point) * complex_treatment * herb_treat_sp * nutrient_treat_sp + (1|metacom_treatment/site_treatment), data = PPB_merged_lmer)
summary(PPB_cover_time_lmer)
anova(PPB_cover_time_lmer)
plot(ggpredict(PPB_cover_time_lmer, terms = c("time_point","complex_treatment","nutrient_treat_sp","herb_treat_sp")))

# model assumptions
PPB_cover_time_lmer_fit<-simulateResiduals(fittedModel = PPB_cover_time_lmer, plot = T)
testDispersion(PPB_cover_time_lmer_fit) # good


# Diversity - basic
PPB_div_basic_lmer<-lmer(H ~ complex_treatment * herb_treat_sp * nutrient_treat_sp + (1|metacom_treatment), data = PPB_merged_lmer)
summary(PPB_div_basic_lmer)
anova(PPB_div_basic_lmer)
emmeans(PPB_div_basic_lmer, list(pairwise ~ herb_treat_sp * nutrient_treat_sp), adjust = "tukey")

# Diversity - over time
PPB_div_time_lmer<-lmer(H ~ as.factor(time_point) * complex_treatment * herb_treat_sp * nutrient_treat_sp + (1|metacom_treatment/site_treatment), data = PPB_merged_lmer)
summary(PPB_div_time_lmer)
anova(PPB_div_time_lmer)
plot(ggpredict(PPB_div_time_lmer, terms = c("time_point","complex_treatment","nutrient_treat_sp")))

# model assumptions
PPB_div_time_lmer_fit<-simulateResiduals(fittedModel = PPB_div_time_lmer, plot = F)
plot(PPB_div_time_lmer_fit) # heteroscedastic
testDispersion(PPB_div_time_lmer_fit) # good



## Cover plot
# Object from ggpredict
PPB_cover_lmer_predict<-ggpredict(PPB_cover_time_lmer, terms = c("time_point","complex_treatment","herb_treat_sp","nutrient_treat_sp"))
PPB_cover_lmer_predict_df<-as.data.frame(PPB_cover_lmer_predict)
head(PPB_cover_lmer_predict_df)
names(PPB_cover_lmer_predict_df) <- c("time_point","total_cover","std.error","conf.low","conf.high","complex_treatment","herb_treat_sp","nutrient_treat_sp")

# Make metadata dataframe
PPB_groups<-PPB_merged_plot %>%
  group_by(complex_treatment,herb_treat_sp,nutrient_treat_sp,time_point,Day_count) %>%
  reframe()

# Merge ggpredict object with metadata
PPB_cover_lmer_predict_df_complete<-merge(PPB_cover_lmer_predict_df,PPB_groups, all = T) %>%
  mutate(total_cover = if_else(Day_count == 0, 0, total_cover)) # add zero initial point

# Fix ordering
PPB_cover_lmer_predict_df_complete$nutrient_treat_sp = factor(PPB_cover_lmer_predict_df_complete$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
PPB_cover_lmer_predict_df_complete$herb_treat_sp = factor(PPB_cover_lmer_predict_df_complete$herb_treat_sp, levels=c("low herbivore","high herbivore"))

PPB_merged_plot$nutrient_treat_sp = factor(PPB_merged_plot$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
PPB_merged_plot$herb_treat_sp = factor(PPB_merged_plot$herb_treat_sp, levels=c("low herbivore","high herbivore"))

PPB_cover_lmer_predict_df_complete<-PPB_cover_lmer_predict_df_complete %>%
  mutate(complex_treatment_het = recode(complex_treatment,
                                "simple" = "high access",
                                "complex" = "low access"))

PPB_merged_plot<-PPB_merged_plot %>%
  mutate(complex_treatment_het = recode(complex_treatment,
                                        "simple" = "high access",
                                        "complex" = "low access"))

PPB_merged_plot$complex_treatment_het = factor(PPB_merged_plot$complex_treatment_het, levels=c("high access","low access"))
PPB_cover_lmer_predict_df_complete$complex_treatment_het = factor(PPB_cover_lmer_predict_df_complete$complex_treatment_het, levels=c("high access","low access"))

# At what point is cover greater than 95%?
PPB_cover_estimate_standard <- PPB_merged_plot %>%
  group_by(Day_count, complex_treatment_het, herb_treat_sp, nutrient_treat_sp) %>%
  dplyr::summarize(mean_cover = mean(total_cover)) %>%
  arrange(Day_count) %>%
  dplyr::filter(mean_cover > 95) %>%
  group_by(complex_treatment_het, herb_treat_sp, nutrient_treat_sp) %>%
  slice(1)

PPB_cover_estimate<-PPB_cover_estimate_standard
write.csv(PPB_cover_estimate_standard,"./Data/PPB_data/PPB_sat_point.csv", row.names = F)


PPB_cover_plot<-ggplot(PPB_cover_lmer_predict_df_complete, aes(x = Day_count, y = total_cover, color = complex_treatment_het)) +
  geom_line() +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_color_manual(values = c("black","magenta")) +
  theme_bw()+
  facet_wrap(herb_treat_sp~nutrient_treat_sp, ncol = 4, labeller = as_labeller(facet_labels)) +
  labs(x = "Time (d)", y = "Algal cover \n(logit transformed)", color = "Access")  +
  guides(color = "none", fill = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_text(data = PPB_merged_plot %>%
              group_by(herb_treat_sp,nutrient_treat_sp) %>%
              dplyr::summarize(),
            aes(label = letters[9:12]), x = -Inf, y = Inf, size = 4, color = "black", hjust = -0.5, vjust = 1) +
  geom_text_repel(data = PPB_merged_plot %>%
                    group_by(Day_count,complex_treatment_het,herb_treat_sp,nutrient_treat_sp) %>%
                    filter(!Day_count == 0) %>%
                    dplyr::summarise(mean_cover = mean(logitTransform(total_cover/100),na.rm =T),
                              count = n()),
                  aes(label = count, x = Day_count, y = mean_cover, color = complex_treatment_het),
                  min.segment.length = Inf,
                  size = 3) +
  geom_segment(data = PPB_cover_estimate_standard, aes(x = Day_count, y = -0.3, yend = 0.5, color = complex_treatment_het),
               size = 1,
               arrow = arrow(length = unit(0.2, "cm")),
               position = position_dodge(width = 20),
               lineend = "butt", linejoin = "mitre",
               show_guide = F) +
  geom_point(data = PPB_merged_plot, aes(x = Day_count, y = logitTransform(total_cover/100), color = complex_treatment_het), alpha = 0.1, inherit.aes = F)



# Diversity plot

# Object from ggpredict
PPB_div_lmer_predict<-ggpredict(PPB_div_time_lmer, terms = c("time_point","complex_treatment","herb_treat_sp","nutrient_treat_sp"))
PPB_div_lmer_predict_df<-as.data.frame(PPB_div_lmer_predict)
head(PPB_div_lmer_predict_df)
names(PPB_div_lmer_predict_df) <- c("time_point","H","std.error","conf.low","conf.high","complex_treatment","herb_treat_sp","nutrient_treat_sp")

# Merge ggpredict object with metadata
PPB_div_lmer_predict_df_complete<-merge(PPB_div_lmer_predict_df,PPB_groups, all = T) %>%
  mutate(H = if_else(Day_count == 0, 0, H)) # add zero initial point

# Fix ordering
PPB_div_lmer_predict_df_complete$nutrient_treat_sp = factor(PPB_div_lmer_predict_df_complete$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
PPB_div_lmer_predict_df_complete$herb_treat_sp = factor(PPB_div_lmer_predict_df_complete$herb_treat_sp, levels=c("low herbivore","high herbivore"))


PPB_div_lmer_predict_df_complete<-PPB_div_lmer_predict_df_complete %>%
  mutate(complex_treatment_het = recode(complex_treatment,
                                        "simple" = "high access",
                                        "complex" = "low access"))

PPB_div_lmer_predict_df_complete$complex_treatment_het = factor(PPB_div_lmer_predict_df_complete$complex_treatment_het, levels=c("high access","low access"))

PPB_div_plot<-ggplot(PPB_div_lmer_predict_df_complete, aes(x = Day_count, y = H, color = complex_treatment_het)) +
  geom_line() +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  scale_color_manual(values = c("black","magenta")) +
  theme_bw()+
  facet_wrap(herb_treat_sp~nutrient_treat_sp, ncol = 4, labeller = as_labeller(facet_labels)) +
  labs(x = "Time (d)", y = "H'", color = "Access")  +
  guides( fill = "none") +
  theme(#axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  geom_text(data = PPB_merged_plot %>%
              group_by(herb_treat_sp,nutrient_treat_sp) %>%
              dplyr::summarise(),
            aes(label = letters[13:16]), x = -Inf, y = Inf, size = 4, color = "black", hjust = -0.5, vjust = 1) +
  geom_point(data = PPB_merged_plot, aes(x = Day_count, y = H, color = complex_treatment_het), alpha = 0.1, inherit.aes = F)




# over time
PPB_cover_div_time_plot<-PPB_cover_plot / PPB_div_plot


ggsave(filename = "./Figures/PPB_summary_time_fig.pdf",
       plot = PPB_cover_div_time_plot,
       dpi = 300,
       width = 8,
       height = 4)



# Systems combined plot for cover and diversity ####
# add images to top panels

combined_summary_plot<-
  icon_facet/ MRA_cover_plot / MRA_div_plot / PPB_cover_plot / PPB_div_plot +
  plot_layout(heights = unit(c(0.8,2,2,2,2), c('null'))) & theme(legend.position = "bottom")




ggsave(filename = "./Figures/combined_summary_time_fig.pdf",
       plot = combined_summary_plot,
       dpi = 300,
       width = 8,
       height = 7)

ggsave(filename = "./Figures/combined_summary_time_fig.png",
       plot = combined_summary_plot,
       dpi = 300,
       width = 8,
       height = 7)




# LMER tables ####

library(htmltools)
library(kableExtra)
library(modelsummary)

round_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return("; p < 0.001")
  } else {
    return(paste0("; p = ", formatC(p_value, format = "f", digits = 3)))
  }
}

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}


# Tropical ####

# Cover
MRA_cover_time_table<-data.frame(anova(MRA_cover_time_lmer))
MRA_cover_time_table$variable <-c( "time",
                                   "accessibility",
                                   "herbivore",
                                   "nutrient",
                                   "time:access",
                                   "time:herbivore",
                                   "access:herbivore",
                                   "time:nutrient",
                                   "access:nutrient",
                                    "herbivore:nutrient",
                                    "time:access:herbivore",
                                    "time:access:nutrient",
                                    "time:herbivore:nutrient",
                                    "access:herbivore:nutrient",
                                    "time:access:herbivore:nutrient")
rownames(MRA_cover_time_table) <- MRA_cover_time_table$variable


# Add a column with subscript text using the variable
DenDF_MRA_cov_est <- sapply(round(MRA_cover_time_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
MRA_cover_time_table$rounded_p_value <- sapply(MRA_cover_time_table$Pr..F., round_p_value)
MRA_cover_time_table$p_star <- makeStars(MRA_cover_time_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
MRA_cover_time_table$cover_label <- paste0("</p>", DenDF_MRA_cov_est, "</p>", "  =  ",round(MRA_cover_time_table$F.value,2),"<sup>", MRA_cover_time_table$p_star, "</sup>")
MRA_cover_time_table$MRA_cover <-MRA_cover_time_table$cover_label


# Diversity
MRA_div_time_table<-data.frame(anova(MRA_div_time_lmer))
MRA_div_time_table$variable <-c( "time",
                                   "accessibility",
                                   "herbivore",
                                   "nutrient",
                                   "time:access",
                                   "time:herbivore",
                                   "access:herbivore",
                                   "time:nutrient",
                                   "access:nutrient",
                                   "herbivore:nutrient",
                                   "time:access:herbivore",
                                   "time:access:nutrient",
                                   "time:herbivore:nutrient",
                                   "access:herbivore:nutrient",
                                   "time:access:herbivore:nutrient")
rownames(MRA_div_time_table) <- MRA_cover_time_table$variable


# Add a column with subscript text using the variable
DenDF_MRA_div_est <- sapply(round(MRA_div_time_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
MRA_div_time_table$rounded_p_value <- sapply(MRA_div_time_table$Pr..F., round_p_value)
MRA_div_time_table$p_star <- makeStars(MRA_div_time_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
MRA_div_time_table$div_label <- paste0("</p>", DenDF_MRA_div_est, "</p>", "  =  ",round(MRA_div_time_table$F.value,2),"<sup>", MRA_div_time_table$p_star, "</sup>")
MRA_div_time_table$MRA_diversity <-MRA_div_time_table$div_label


# Temperate ####

# Cover
PPB_cover_time_table<-data.frame(anova(PPB_cover_time_lmer))
PPB_cover_time_table$variable <-c( "time",
                                   "accessibility",
                                   "herbivore",
                                   "nutrient",
                                   "time:access",
                                   "time:herbivore",
                                   "access:herbivore",
                                   "time:nutrient",
                                   "access:nutrient",
                                   "herbivore:nutrient",
                                   "time:access:herbivore",
                                   "time:access:nutrient",
                                   "time:herbivore:nutrient",
                                   "access:herbivore:nutrient",
                                   "time:access:herbivore:nutrient")
rownames(PPB_cover_time_table) <- PPB_cover_time_table$variable


# Add a column with subscript text using the variable
DenDF_PPB_cov_est <- sapply(round(PPB_cover_time_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
PPB_cover_time_table$rounded_p_value <- sapply(PPB_cover_time_table$Pr..F., round_p_value)
PPB_cover_time_table$p_star <- makeStars(PPB_cover_time_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
PPB_cover_time_table$cover_label <- paste0("</p>", DenDF_PPB_cov_est, "</p>", "  =  ",round(PPB_cover_time_table$F.value,2),"<sup>", PPB_cover_time_table$p_star, "</sup>")
PPB_cover_time_table$PPB_cover <-PPB_cover_time_table$cover_label


# Diversity
PPB_div_time_table<-data.frame(anova(PPB_div_time_lmer))
PPB_div_time_table$variable <-c( "time",
                                 "accessibility",
                                 "herbivore",
                                 "nutrient",
                                 "time:access",
                                 "time:herbivore",
                                 "access:herbivore",
                                 "time:nutrient",
                                 "access:nutrient",
                                 "herbivore:nutrient",
                                 "time:access:herbivore",
                                 "time:access:nutrient",
                                 "time:herbivore:nutrient",
                                 "access:herbivore:nutrient",
                                 "time:access:herbivore:nutrient")
rownames(PPB_div_time_table) <- PPB_cover_time_table$variable


# Add a column with subscript text using the variable
DenDF_PPB_div_est <- sapply(round(PPB_div_time_table$DenDF,2), function(x) HTML(paste0("<i>F</i><sub>", x, "</sub>")))
PPB_div_time_table$rounded_p_value <- sapply(PPB_div_time_table$Pr..F., round_p_value)
PPB_div_time_table$p_star <- makeStars(PPB_div_time_table$Pr..F.)
#stab_rich_table$label <- paste0("</p>", DenDF_values, "</p>", "  =  ",round(stab_rich_table$F.value,2),stab_rich_table$rounded_p_value)
PPB_div_time_table$div_label <- paste0("</p>", DenDF_PPB_div_est, "</p>", "  =  ",round(PPB_div_time_table$F.value,2),"<sup>", PPB_div_time_table$p_star, "</sup>")
PPB_div_time_table$PPB_diversity <-PPB_div_time_table$div_label
#PPB_div_time_table$PPB_diversity <-paste0('"',PPB_div_time_table$div_label,'"', sep = "")



# Merge
MRA_lmer_results<-data.frame(variable = MRA_cover_time_table$variable,
                             MRA_cover = MRA_cover_time_table$MRA_cover,
                             MRA_diversity = MRA_div_time_table$MRA_diversity)

PPB_lmer_results<-data.frame(variable = PPB_cover_time_table$variable,
                             PPB_cover = PPB_cover_time_table$PPB_cover,
                             PPB_diversity = PPB_div_time_table$PPB_diversity)

# Merge for csv tables
MRA_lmer_results_csv<-data.frame(variable = MRA_cover_time_table$variable,
                             Cover_Ddf = MRA_cover_time_table$DenDF,
                             Cover_Ndf = MRA_cover_time_table$NumDF,
                             Cover_F = MRA_cover_time_table$F.value,
                             Cover_P = MRA_cover_time_table$Pr..F.,
                             Cover_ast = MRA_cover_time_table$p_star,
                             div_Ddf = MRA_div_time_table$DenDF,
                             div_Ddf = MRA_div_time_table$NumDF,
                             div_F = MRA_div_time_table$F.value,
                             div_P = MRA_div_time_table$Pr..F.,
                             div_ast = MRA_div_time_table$p_star)

PPB_lmer_results_csv<-data.frame(variable = PPB_cover_time_table$variable,
                                 Cover_Ddf = PPB_cover_time_table$DenDF,
                                 Cover_Ndf = PPB_cover_time_table$NumDF,
                                 Cover_F = PPB_cover_time_table$F.value,
                                 Cover_P = PPB_cover_time_table$Pr..F.,
                                 Cover_ast = PPB_cover_time_table$p_star,
                                 div_Ddf = PPB_div_time_table$DenDF,
                                 div_Ndf = PPB_div_time_table$NumDF,
                                 div_F = PPB_div_time_table$F.value,
                                 div_P = PPB_div_time_table$Pr..F.,
                                 div_ast = PPB_div_time_table$p_star)

write.csv(MRA_lmer_results_csv,"./Tables/MRA_lmer_table_complete.csv", row.names = F)
write.csv(PPB_lmer_results_csv,"./Tables/PPB_lmer_table_complete.csv", row.names = F)





MRA_lmer_table_complete<-MRA_lmer_results %>%
  kbl(format = "html", escape = F, col.names = c("Fixed effects", "Cover","Diversity")) %>%
  kable_classic(html_font = "Times") %>%
  footnote(general = "Results from linear mixed effects modeling (LME) assessing variation in (A) total algal cover and (B) Shannon-Wiener diversity among tropical marine algae across treatments as fixed effects. Dataset was treated as a random effect. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001",
           escape = T) %>%
  kable_styling()

TPPB_lmer_table_complete<-PPB_lmer_results %>%
  kbl(format = "html", escape = F, col.names = c("Fixed effects", "Cover","Diversity")) %>%
  kable_classic(html_font = "Times") %>%
  footnote(general = "Results from linear mixed effects modeling (LME) assessing variation in (A) total algal cover and (B) Shannon-Wiener diversity among temperate marine algae across treatments as fixed effects. Dataset was treated as a random effect. Subscripts are denominator degrees of freedom from Satterhwaite approximation. Superscripts represent significance based on p-value: ns is p > 0.05; * p < 0.05; ** p < 0.01; *** p < 0.001; **** p < 0.0001",
           escape = T) %>%
  kable_styling()

#rmarkdown::render('LMER_tables.Rmd')


# MDS ####

# function for making smaller legend to fit vectors

addSmallLegend <- function(aPlot, pointSize = 1, textSize = 8, spaceLegend = 0.7) {
  aPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize, linewidth = pointSize, shape = "line")),
           color = guide_legend(override.aes = list(size = pointSize, linewidth = pointSize, shape = "line"))) +
    theme(legend.title = element_text(size = textSize),
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

addLargeLegend <- function(aPlot, pointSize = 2, textSize = 12, spaceLegend = 0.8) {
  aPlot +
    guides(
      shape = guide_legend(override.aes = list(size = pointSize, shape = 21)),
      color = guide_legend(override.aes = list(size = pointSize, shape = 21))
    ) +
    theme(
      legend.title = element_text(size = textSize),
      legend.text  = element_text(size = textSize),
      legend.key.size = unit(spaceLegend, "lines")
    )
}

# For getting SE
se_fn <- function(x) sd(x) / sqrt(length(x))


# Moorea MDS
MRA_mds_df<-read.csv("./Data/Moorea_data/MRA_mds_points.csv")
MRA_vct_df<-read.csv("./Data/Moorea_data/MRA_mds_vectors.csv")

# PPB MDS
PPB_mds_df<-read.csv("./Data/PPB_data/PPB_mds_points.csv")
PPB_vct_df<-read.csv("./Data/PPB_data/PPB_mds_vectors.csv")


# Vector colors
MRA_vect_col<-read.csv("./Data/Moorea_data/MRA_vect_colors.csv")
PPB_vect_col<-read.csv("./Data/PPB_data/PPB_vect_colors.csv")


PPB_vect_col<-PPB_vect_col %>% filter(species %in% PPB_vct_df$species) # filter out old colors


MRA_vct_df$species<-str_to_title(MRA_vct_df$species)
PPB_vct_df$species<-str_to_title(PPB_vct_df$species)

MRA_vect_col$species<-str_to_title(MRA_vect_col$species)
PPB_vect_col$species<-str_to_title(PPB_vect_col$species)

MRA_vct_df_col<-merge(MRA_vct_df,MRA_vect_col)
PPB_vct_df_col<-merge(PPB_vct_df,PPB_vect_col)

# Generate specific colors for vectors within algal functional groups
color_generate_merged_df <- function(main_df, color_list) {
  MRA_colors <- main_df %>%
    group_by(color_agg) %>%
    dplyr::summarize(shades = n())

  MRA_color_list <- lapply(names(color_list), function(color) {
    MRA_color <- MRA_colors %>% filter(color_agg == color)
    if (MRA_color$shades <= 0) {
      warning(paste("No shades found for", color))
      return(NULL)
    }
    data.frame(color_agg = color, shade = grDevices::colorRampPalette(colors = color_list[[color]])(MRA_color$shades))
  }) %>% do.call(rbind, .)

  MRA_color_list <- MRA_color_list[!is.null(MRA_color_list)]

  MRA_colors_merged_df <- merge(main_df, MRA_color_list, by = "color_agg", all.x = TRUE)
  merged_df <- unique(MRA_colors_merged_df)

  main_df <- arrange(main_df, color_agg)
  MRA_color_list <- arrange(MRA_color_list, color_agg)

  MRA_colors_merged_df <- cbind(main_df, shade = MRA_color_list$shade)

  return(MRA_colors_merged_df)
}

# Moorea vectors with unique colors
# These set the ranges for each color group
MRA_vect_colors_merged_df <- color_generate_merged_df(MRA_vct_df_col, list(
  blue = c("blue", "cyan"), # cyanobacteria
  brown = c("burlywood", "chocolate4"), # brown algae; phaeophyta
  red = c("red1", "orange1"), # red algae
  green = c("green", "darkgreen"), # green algae
  pink = c("pink1", "orchid"), # CCA
  grey = c("grey", "black") # additional categories: sediment, bare space
))

# PPB vectors with unique colors
PPB_vect_colors_merged_df <- color_generate_merged_df(PPB_vct_df_col, list(
  blue = c("blue"), # cyanobacteria
  brown = c("burlywood", "chocolate4"), # brown algae; phaeophyta
  red = c("red1", "orange1"), # red algae
  green = c("green", "darkgreen"), # green algae
  pink = c("pink1", "orchid"), # CCA
  grey = c("grey") # additional categories: sediment, bare space
))



# Moorea panels

# Curate
MRA_mds_dftime_correct<-MRA_mds_df %>%
  mutate(time_point_num = recode(time_point,
                                 "T2" = "1",
                                 "T3" = "2",
                                 "T4" = "3",
                                 "T5" = "4"))
# Add zero point
MRA_mds_zero_df<-MRA_mds_dftime_correct %>%
  filter(time_point_num == "1") %>%
  mutate(time_point_num = "0",
         MDS1 = 0,
         MDS2 = 0)

MRA_mds_df_merged <-rbind(MRA_mds_dftime_correct,MRA_mds_zero_df)

MRA_summarized_mds <- MRA_mds_df_merged %>%
  dplyr::group_by(com_treat,time_point_num,nutrient_treat_sp,herb_treat_sp) %>%
  dplyr::summarize(mean_MDS1 = mean(MDS1),
                   mean_MDS2 = mean(MDS2),
                   se_MDS1 = se_fn(MDS1),
                   se_MDS2 = se_fn(MDS2))


MRA_summarized_mds$nutrient_treat_sp = factor(MRA_summarized_mds$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
MRA_summarized_mds$herb_treat_sp = factor(MRA_summarized_mds$herb_treat_sp, levels=c("low herbivore","high herbivore"))

MRA_color_mapping <- setNames(MRA_vect_colors_merged_df$shade, MRA_vect_colors_merged_df$species)

MRA_community_MDS<-ggplot(MRA_summarized_mds[order(MRA_summarized_mds$time_point_num),], aes(x=mean_MDS1,y=mean_MDS2)) +
  geom_hline(yintercept = 0, alpha = 0.1) +
  geom_vline(xintercept = 0, alpha = 0.1) +
  geom_segment(data=MRA_vect_colors_merged_df,aes(x=0,xend=NMDS1/1.5,y=0,yend=NMDS2/1.5, color = species), size = 0.8,arrow = arrow(length = unit(0.3,"cm")),
               key_glyph = "abline") +
  #scale_color_manual(values = unique(MRA_vect_colors_merged_df$shade), labels = unique(MRA_vect_colors_merged_df$species)) +
  scale_color_manual(values = MRA_color_mapping) +
  geom_errorbarh(aes(xmax = mean_MDS1 + se_MDS1, xmin = mean_MDS1 - se_MDS1)) +
  geom_errorbar(aes(ymax = mean_MDS2 + se_MDS2, ymin = mean_MDS2 - se_MDS2)) +
  geom_path(aes(group=com_treat)) +
  geom_point(aes(fill = com_treat), size=5, pch = 21, alpha = 0.8) +
  geom_point(data = . %>% filter(time_point_num == "0"),
             size=6, fill="grey", pch =21) +
  theme_bw() +
  scale_fill_manual(values = c("magenta","black")) +
  geom_text(aes(label = time_point_num)) +
  geom_text(data=MRA_summarized_mds[MRA_summarized_mds$com_treat == "C", ], aes(x=mean_MDS1, y=mean_MDS2, label=time_point_num), color="black") +  # Setting color="white" for the black points
  geom_text(data=MRA_summarized_mds[MRA_summarized_mds$com_treat == "S", ], aes(x=mean_MDS1, y=mean_MDS2, label=time_point_num), color="white") +  # Setting color="white" for the black points
  removeGrid() +
  labs(color = "Tropical algae",fill = "Tile treatment", x = "MDS1",y="Tropical\n\nMDS2") +
  facet_wrap(herb_treat_sp~ nutrient_treat_sp, ncol = 4) +
  theme(text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.spacing = unit(.01,"cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(fill = "none",  color=guide_legend(ncol=3)) +
  coord_cartesian(xlim= c(-0.61,0.61),
                  ylim= c(-0.61,0.61)) +
  geom_text(data = MRA_summarized_mds %>%
              group_by(herb_treat_sp, nutrient_treat_sp) %>%
              dplyr::summarise(),
            aes(label = letters[1:4]), x = -Inf, y = Inf, size = 5, color = "black", hjust = -0.7, vjust = 1.2)



# Vector regions
ggplot(MRA_vct_df,aes(x=NMDS1,y=NMDS2)) +
  geom_text(aes(label = species,color = species)) +
  theme_bw() +
  lims(x = c(-1,1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
  #remove_grid()
# removeGrid()# PPB panels



# PPB ####
# Curate a bit
PPB_time_points<-read.csv("./Data/PPB_data/PPB_time_points.csv")
PPB_mds_df<-merge(PPB_mds_df,PPB_time_points)

# Add zero point
PPB_mds_zero_df<-PPB_mds_df %>%
  filter(time_point == "1") %>%
  mutate(time_point = "0",
         MDS1 = 0,
         MDS2 = 0)

PPB_mds_df_merged <-rbind(PPB_mds_df,PPB_mds_zero_df)



PPB_summarized_mds <- PPB_mds_df_merged %>%
  #dplyr::group_by(complex_treatment,time_point,site,zone) %>%
  dplyr::group_by(complex_treatment,time_point,nutrient_treat_sp,herb_treat_sp) %>%
  dplyr::summarize(mean_MDS1 = mean(MDS1),
                   mean_MDS2 = mean(MDS2),
                   se_MDS1 = se_fn(MDS1),
                   se_MDS2 = se_fn(MDS2))


PPB_summarized_mds$nutrient_treat_sp = factor(PPB_summarized_mds$nutrient_treat_sp, levels=c("low nutrients","high nutrients"))
PPB_summarized_mds$herb_treat_sp = factor(PPB_summarized_mds$herb_treat_sp, levels=c("low herbivore","high herbivore"))

PPB_summarized_mds<-PPB_summarized_mds %>%
  mutate(complex_treatment_het = recode(complex_treatment,
                            "simple" = "high access",
                            "complex" = "low access"))

PPB_color_mapping <- setNames(PPB_vect_colors_merged_df$shade, PPB_vect_colors_merged_df$species)

PPB_community_MDS<-ggplot(PPB_summarized_mds[order(PPB_summarized_mds$time_point),], aes(x=mean_MDS1,y=mean_MDS2)) +
  geom_hline(yintercept = 0, alpha = 0.1) +
  geom_vline(xintercept = 0, alpha = 0.1) +
  geom_segment(data=PPB_vect_colors_merged_df,aes(x=0,xend=NMDS1*1.5,y=0,yend=NMDS2*1.5, color = species), size = 0.8,arrow = arrow(length = unit(0.3,"cm")),
               key_glyph = "abline") +
#  scale_color_manual(values = unique(PPB_vect_colors_merged_df$shade), labels = unique(PPB_vect_colors_merged_df$species)) +
  scale_color_manual(values = PPB_color_mapping) +
  geom_errorbarh(aes(xmax = mean_MDS1 + se_MDS1, xmin = mean_MDS1 - se_MDS1)) +
  geom_errorbar(aes(ymax = mean_MDS2 + se_MDS2, ymin = mean_MDS2 - se_MDS2)) +
  geom_path(aes(group=complex_treatment_het)) +
  geom_point(aes(fill = complex_treatment_het), size=5, pch = 21, alpha = 0.8) +
  geom_point(data = . %>% filter(time_point == "0"),
             size=6, fill="grey", pch =21) +
  theme_bw() +
  scale_fill_manual(values = c("black","magenta")) +
  geom_text(data=PPB_summarized_mds[PPB_summarized_mds$complex_treatment == "complex", ], aes(x=mean_MDS1, y=mean_MDS2, label=time_point), color="black") +  # Setting color="white" for the black points
  geom_text(data=PPB_summarized_mds[PPB_summarized_mds$complex_treatment == "simple", ], aes(x=mean_MDS1, y=mean_MDS2, label=time_point), color="white") +  # Setting color="white" for the black points
  removeGrid() +
  labs(color = "Temperate algae",fill = "Access", x = "MDS1",y="Temperate\n\nMDS2") +
  facet_wrap(herb_treat_sp~nutrient_treat_sp, ncol = 4) +
  theme(text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        panel.spacing = unit(.01,"cm"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.text = element_text(size = 12)) +
  guides(color=guide_legend(ncol=3, order = 1), fill  = guide_legend(order = 1)) +
  coord_cartesian(xlim= c(-1.2,1),
                  ylim= c(-1.2,1)) +
  geom_text(data = PPB_summarized_mds %>%
              group_by(herb_treat_sp, nutrient_treat_sp) %>%
              dplyr::summarise(),
            aes(label = letters[5:8]), x = -Inf, y = Inf, size = 5, color = "black", hjust = -0.7, vjust = 1.2)




MRA_community_MDS_sm<-addLargeLegend(MRA_community_MDS)
PPB_community_MDS_sm<-addLargeLegend(PPB_community_MDS)

# Combine the MDS
combined_mds_plot<-
  icon_facet/MRA_community_MDS_sm / PPB_community_MDS_sm +
  plot_layout(heights = unit(c(0.8,2,2), c('null')),guides = 'collect') &
  theme(legend.position = "right")



ggsave(filename = "./Figures/combined_MDS_new.pdf",
       plot = combined_mds_plot,
       dpi = 600,
       width = 12,
       height = 6.6)

ggsave(filename = "./Figures/combined_MDS_new.png",
       plot = combined_mds_plot,
       dpi = 600,
       width = 12,
       height = 6.6
       )

# Vector regions
ggplot(PPB_vect_colors_merged_df,aes(x=NMDS1,y=NMDS2,label = species,color = species)) +
  #geom_point() +
  geom_text() +
  theme_bw() +
  scale_color_manual(values = PPB_color_mapping) +
  #scale_color_manual(values = unique(PPB_vect_colors_merged_df$shade), labels = unique(PPB_vect_colors_merged_df$species)) +
  #lims(x = c(-1,1)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
#remove_grid()



# Species table for contributions ####
# Moorea
MRA_labels<-read.csv("./Data/Moorea_data/MRA_CN_labels.csv")
MRA_labels$species<-MRA_labels$Short.Code
MRA_species_all<-read.csv("./Data/Moorea_data/MRA_mds_vectors_all.csv")
MRA_species_complete<-merge(MRA_labels,MRA_species_all)


# get percent cover and SD for each species
MRA_long<-MRA_community %>%
  pivot_longer(cols = c(15:36),names_to = "species",values_to = "cover") %>%
  group_by(species) %>%
  dplyr::summarize(mean_cover = mean(cover,na.rm = T),
            SD_cover = sd(cover, na.rm = T))

MRA_long$mean_sd <- sprintf("%.2f ± %.2f", MRA_long$mean_cover, MRA_long$SD_cover)


MRA_species_table<-merge(MRA_species_complete,MRA_long) %>%
  select(species,Name,mean_sd,pval,r2)

write.csv(MRA_species_table,"./Tables/MRA_complete_species_table.csv", row.names = F)


# PPB
PPB_labels<-read.csv("./Data/PPB_data/PPB_CN_labels_updated.csv")
PPB_labels$species<-PPB_labels$Short.Code
PPB_species_all<-read.csv("./Data/PPB_data/PPB_mds_vectors_all.csv")
PPB_species_complete<-merge(PPB_labels,PPB_species_all)


# get percent cover and SD for each species
PPB_long<-PPB_filterd_com_meta %>%
  pivot_longer(cols = c(18:53),names_to = "species",values_to = "cover") %>%
  group_by(species) %>%
  dplyr::summarize(mean_cover = mean(cover,na.rm = T),
            SD_cover = sd(cover, na.rm = T))

PPB_long$mean_sd <- sprintf("%.2f ± %.2f", PPB_long$mean_cover, PPB_long$SD_cover)

PPB_species_table<-merge(PPB_species_complete,PPB_long) %>%
  select(species,Name,mean_sd,pval,r2)

write.csv(PPB_species_table,"./Tables/PPB_complete_species_table_new.csv", row.names = F)




# get all package citations exported for referencing
#knitr::write_bib(c(.packages()), "packages.bib")

# END #
