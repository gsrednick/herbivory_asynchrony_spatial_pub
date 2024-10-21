# GAM models

library(tidyverse)
library(ggeffects)
library(mgcv)
library(ggExtra)
library(ggforce)
library(gratia)
library(ggExtra)
library(patchwork)
library(ggrepel)

# Moorea ####
## Data curation ####
MRA_cleaned_merged_dist<-read.csv("./Data/Moorea_data/MRA_complete.csv")

MRA_cleaned_merged_dist$metacom_unique <- paste(MRA_cleaned_merged_dist$site, MRA_cleaned_merged_dist$zone, MRA_cleaned_merged_dist$meta, sep = "_")
MRA_cleaned_merged_dist$time_point_num<-as.numeric(MRA_cleaned_merged_dist$time_point_num)

# Metadata for merge
MRA_metadata<-MRA_cleaned_merged_dist %>%
  group_by(site,zone,nutrient_treat_sp,herb_treat_sp,metacom_unique,meta_treat) %>%
  reframe()

MRA_metadata$site_zone <- paste(MRA_metadata$site, MRA_metadata$zone, sep = "_")


# ANOCOVA style LMER

MRA_cleaner<-MRA_cleaned_merged_dist %>%
  filter(!time_point_num == 0, !nutrient_treat_sp== "mid nutrients") %>%
  mutate(time_point_num = as.numeric(time_point_num))


MRA_cleaner_W_zero<-MRA_cleaned_merged_dist %>%
  filter(!nutrient_treat_sp== "mid nutrients") %>%
  mutate(time_point_num = as.numeric(time_point_num))


# New try -- try making a new concatenated variable that combines the 3 categorical variables, tests splines, then plot seperately
MRA_cleaner_W_zero_treats<-MRA_cleaner_W_zero %>%
  mutate(treatment_vector = as.factor(paste(meta_treat,nutrient_treat_sp,herb_treat_sp, sep = "_")))



MRA_cleaner_W_zero_treats_summarized<-MRA_cleaner_W_zero_treats %>%
  group_by(meta_treat,nutrient_treat_sp,herb_treat_sp,treatment_vector) %>%
  reframe() %>%
  mutate(group = treatment_vector,
         facet = nutrient_treat_sp,
         panel = herb_treat_sp)

# Points for plots
MRA_cleaner_W_zero_point<-MRA_cleaner_W_zero %>% mutate(x = time_point_num,
                                                        predicted = value,
                                                        group = meta_treat,
                                                        facet = nutrient_treat_sp,
                                                        panel = herb_treat_sp)

# Make factors for analysis
MRA_cleaner_W_zero$meta_treat<-as.factor(MRA_cleaner_W_zero$meta_treat)
MRA_cleaner_W_zero$herb_treat_sp<-as.factor(MRA_cleaner_W_zero$herb_treat_sp)
MRA_cleaner_W_zero$nutrient_treat_sp<-as.factor(MRA_cleaner_W_zero$nutrient_treat_sp)


## GAM ####
kvalue = 3

MRA_meta_rep<-MRA_cleaner_W_zero %>%
  group_by(time_point, meta_treat, nutrient_treat_sp, herb_treat_sp) %>%
  dplyr::summarize(count = n())


gam_V1.5 <- gam(value ~
                  s(time_point_num, by = nutrient_treat_sp, k = kvalue, bs = "fs") +
                  s(time_point_num, by = herb_treat_sp, k = kvalue, bs = "fs") +
                  s(time_point_num, by = meta_treat, k = kvalue, bs = "fs") +
                  t2(time_point_num, meta_treat, nutrient_treat_sp, herb_treat_sp,bs = "fs", m=2),
                method = "REML",
                data = MRA_cleaner_W_zero) # using full dataset --> some replicates missing at T3-5 (crest sites)

summary(gam_V1.5)
anova(gam_V1.5)
AIC(gam_V1.5)
qq_plot(gam_V1.5) # good fit for intercept
plot.gam(gam_V1.5)


# Extract model estimates for figures
MRA_gam_summary <- as.data.frame(summary(gam_V1.5)$s.table)
write.csv(MRA_gam_summary, "./Tables/MRA_gam_summary.csv")

# Annotate summary for figure
MRA_gam_summary$treatment<-rownames(MRA_gam_summary)

MRA_rsq<-summary(gam_V1.5)$r.sq
MRA_int_table<-MRA_gam_summary %>%
  filter(str_detect(treatment,"t2")) %>%
  mutate(rsq = MRA_rsq)


MRA_gam_summary_updated <- MRA_gam_summary %>%
  mutate(NUT_string = ifelse(str_detect(treatment, "nutrient_treat_sp"), str_replace(treatment, ".*nutrient_treat_sp", ""), NA),
         HERB_string = ifelse(str_detect(treatment, "herb_treat_sp"), str_replace(treatment, ".*herb_treat_sp", ""), NA),
         CMP_string = ifelse(str_detect(treatment, "meta_treat"), str_replace(treatment, ".*meta_treat", ""), NA)) %>%
  mutate(facet = str_trim(NUT_string),
         panel = str_trim(HERB_string),
         group = str_trim(CMP_string)) %>%
  select(-NUT_string, -HERB_string, -CMP_string) %>%
  filter(!str_detect(treatment,"t2"))


# Get unique values of facet, panel, and group
unique_facet <- unique(MRA_gam_summary_updated$facet)
unique_panel <- unique(MRA_gam_summary_updated$panel)
unique_group <- unique(MRA_gam_summary_updated$group)

# Create all combinations of unique values
combinations <- expand.grid(facet = unique_facet,
                            panel = unique_panel,
                            group = unique_group)


combinations<-na.omit(combinations)
MRA_gam_summary_updated_merge<-MRA_gam_summary_updated
MRA_gam_summary_updated_merge$treatment<-NULL


expanded_df <- merge(combinations, MRA_gam_summary_updated, by = c("facet"), all.x = T)




# Plotting
gam_V1_predict<-ggeffects::ggpredict(gam_V1.5,terms = c("time_point_num [all]","meta_treat","nutrient_treat_sp","herb_treat_sp"))

gam_V1_ANCOVA<-ggplot() +
  geom_line(data = gam_V1_predict, aes(x, predicted, color = group), size = 1) +
  geom_ribbon(data = gam_V1_predict, aes(x = x, ymin = conf.low, ymax = conf.high, color = group, fill = group), alpha = 0.2, linewidth = 0.01) +
  geom_point(data = MRA_cleaner_W_zero_point, aes(x = x, y = predicted, color = group), alpha = 0.2) +
  facet_grid(panel~facet) +
  theme_bw() +
  guides(fill = "none") +
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_fill_manual(values = c("black", "red", "blue")) +
  labs(x = "Time (days)", y= "Community dissimilarity (BC)", color = "Accessibility") #+
  geom_text(data = MRA_gam_summary_updated, aes(label = edf, color = group), x = 10, y = 0.75, hjust = 0, vjust = 0, size = 3, color = "black", inherit.aes = T)


## Plotting these ####

# Add icons
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



# Modified interactions in GAM -- use this one
sms_2 <- smooth_estimates(gam_V1.5, n = 50) %>%
  filter(.type == "Tensor product (T2)")

est_lim <- c(-1, 1) * max(abs(sms_2[[".estimate"]]), na.rm = TRUE)

color_map <- c("black", "magenta","yellow4")  # Add more colors if needed
names(color_map) <- levels(sms_2$meta_treat)

MRA_effect_plot<-
  sms_2 %>%
  mutate(nutrient_treat_sp = factor(nutrient_treat_sp),
         meta_treat = factor(meta_treat),
         herb_treat_sp = factor(herb_treat_sp)) %>%
  ggplot(aes(x = time_point_num, y = meta_treat, fill = est)) +
  geom_raster(aes(x = time_point_num, y = meta_treat, fill = est)) +
  geom_hline(yintercept = 1.5, color = "black", size = 0.2) +
  geom_hline(yintercept = 2.5, color = "black", size = 0.2) +
  facet_grid(herb_treat_sp~ nutrient_treat_sp) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  expand_limits(fill = est_lim) +
  theme_bw() +
  removeGrid() +
  theme(axis.text.y = element_text(color = color_map)) +
  labs(x = "Time (d)", y = "Accessibility", fill = "Smooth term effect")


SMS_ready<-sms_2 %>%
  mutate(nutrient_treat_sp = factor(nutrient_treat_sp),
         meta_treat = factor(meta_treat),
         herb_treat_sp = factor(herb_treat_sp),
         TDBU_treat = paste(nutrient_treat_sp,herb_treat_sp,sep = "_"))

unique(SMS_ready$TDBU_treat)

SMS_ready$meta_treat = factor(SMS_ready$meta_treat, levels=c("simple","complex","mix"))

SMS_ready<-SMS_ready %>%
  mutate(meta_treat_het = recode(meta_treat,
                            "simple" = "high access",
                            "complex" = "low access",
                            "mix" = "mix"))

# Break into separate ggplot elements for plotting together
p_list_sms_MRA<- lapply(sort(unique(SMS_ready$TDBU_treat)), function(i) {

  ggplot(SMS_ready[SMS_ready$TDBU_treat==i,],aes(x = time_point_num, y = meta_treat_het, fill = .estimate)) +
    geom_raster(aes(x = time_point_num, y = meta_treat_het, fill = .estimate)) +
    scale_fill_distiller(palette = "RdBu", type = "div") +
    expand_limits(fill = est_lim) +
    labs(x = "Time (d)", y = "Accessibility", fill = "Smooth term effect") +
    theme_bw() +
    removeGrid() +
    geom_hline(yintercept = 1.5, color = "black", size = 0.2) +
    geom_hline(yintercept = 2.5, color = "black", size = 0.2) +
    geom_vline(xintercept = c(0, 100, 200, 300), color = "black", size = 0.2) +
    theme(axis.text.y = element_text(color = color_map),
    ) +
    ggtitle(paste(i))


})

MRA_LN_HH_A<-p_list_sms_MRA[[3]]
MRA_HH_HH_B<-p_list_sms_MRA[[1]]
MRA_HN_LH_C<-p_list_sms_MRA[[2]]
MRA_LN_LH_D<-p_list_sms_MRA[[4]]


# ANCOVA split

gam_V2_predict_deco<-as.data.frame(gam_V1_predict)
gam_V2_predict_deco$TDBU_treat <- paste(gam_V2_predict_deco$panel,gam_V2_predict_deco$facet,sep = "_")

MRA_cleaner_W_zero_point_deco<-MRA_cleaner_W_zero_point %>%
  mutate(TDBU_treat = paste(herb_treat_sp,nutrient_treat_sp,sep = "_"))

MRA_cleaner_W_zero_point_deco$meta_treat = factor(MRA_cleaner_W_zero_point_deco$meta_treat, levels=c("simple","complex","mix"))
gam_V2_predict_deco$group = factor(gam_V2_predict_deco$group, levels=c("simple","complex","mix"))

# Rename metacommunity factor levels for plotting
gam_V2_predict_deco<-gam_V2_predict_deco %>%
  mutate(group_het = recode(group,
                                 "simple" = "high access",
                                 "complex" = "low access",
                            "mix" = "mix"))

MRA_cleaner_W_zero_point_deco<-MRA_cleaner_W_zero_point_deco %>%
  mutate(group_het = recode(group,
                            "simple" = "high access",
                            "complex" = "low access",
                            "mix" = "mix"))

MRA_cover_estimate$TDBU_treat <- paste(MRA_cover_estimate$herb_treat_sp,MRA_cover_estimate$nutrient_treat_sp,sep = "_")
MRA_cover_estimate$x = MRA_cover_estimate$Day_count
MRA_cover_estimate$group_het <- MRA_cover_estimate$com_treat_het



p_list_gam_MRA<- lapply(sort(unique(gam_V2_predict_deco$TDBU_treat)), function(i) {

  MRA_filtered_points <- MRA_cleaner_W_zero_point_deco %>%
    filter(TDBU_treat == i)  # Corrected the filtering condition

  MRA_cover_estimate_filt <- MRA_cover_estimate %>%
    filter(TDBU_treat == i)  # Corrected the filtering condition

  ggplot() +
    geom_line(data = gam_V2_predict_deco[gam_V2_predict_deco$TDBU_treat==i,], aes(x, predicted, color = group_het), size = 0.75) +
    geom_ribbon(data = gam_V2_predict_deco[gam_V2_predict_deco$TDBU_treat==i,], aes(x = x, ymin = conf.low, ymax = conf.high, color = group_het, fill = group_het), alpha = 0.1, linewidth = 0.01) +
    geom_point(data = MRA_filtered_points, aes(x = x, y = predicted, color = group_het), alpha = 0.2) +
    geom_segment(data = MRA_cover_estimate_filt, aes(x = x, y = 0.80, yend = 0.7, color = group_het),
                 size = 1,
                 arrow = arrow(length = unit(0.2, "cm")),
                 position = position_dodge(width = 40),
                 lineend = "butt", linejoin = "mitre",
                 show_guide = F) +
    theme_bw() +
    guides(fill = "none") +
    scale_color_manual(values = c("black", "magenta","yellow4")) +
    scale_fill_manual(values = c("black", "magenta","yellow4")) +
    coord_cartesian(ylim=c(0, 0.9)) +
    labs(x = "Time (days)", y= "Community dissimilarity (BC)", color = "Accessibility") +
    geom_text_repel(data = MRA_filtered_points %>%
                      group_by(time_point,group_het) %>%
                      filter(!x == 0) %>%
                      dplyr::summarise(mean_predicted = mean(predicted,na.rm =T),
                                x = max(x,na.rm = T),
                                count = n()),
                    aes(label = count, x = x, y = mean_predicted, color = group_het),
                    min.segment.length = Inf,
                    size = 3) +
    ggtitle(paste(i))



})

MRA_gam_LN_HH_A<-p_list_gam_MRA[[2]]
MRA_gam_HN_HH_B<-p_list_gam_MRA[[1]]
MRA_gam_HN_LH_C<-p_list_gam_MRA[[3]]
MRA_gam_LN_LH_D<-p_list_gam_MRA[[4]]




MRA_gam_LN_HH_A_up<-MRA_gam_LN_HH_A + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            plot.margin = unit(c(0,0,0,0),"cm")) +
  ggplot2::annotation_custom(FH_NL_image,0,100,0.7,0.95) +
  annotate("text", label = "(a)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())


MRA_gam_HN_HH_B_up<-MRA_gam_HN_HH_B + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y = element_blank(),
                                            axis.text.y = element_blank(),
                                            plot.margin = unit(c(0,0,0,0.5),"cm")) +
  ggplot2::annotation_custom(FH_NH_image,0,100,0.7,0.95) +
  annotate("text", label = "(b)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())


MRA_gam_HN_LH_C_up<-MRA_gam_HN_LH_C + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y = element_blank(),
                                            axis.text.y = element_blank(),
                                            plot.margin = unit(c(0,0,0,0.5),"cm")) +
  ggplot2::annotation_custom(FL_NH_image,0,100,0.7,0.95) +
  annotate("text", label = "(d)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())


MRA_gam_LN_LH_D_up<-MRA_gam_LN_LH_D + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            plot.margin = unit(c(0,0,0,0),"cm")) +
  ggplot2::annotation_custom(FL_NL_image,0,100,0.7,0.95) +
  annotate("text", label = "(c)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())


MRA_LN_HH_A_up<-MRA_LN_HH_A + theme(axis.title.y = element_blank(),axis.title.x = element_blank(),
                                    plot.margin = unit(c(0,0,0,0),"cm")) + ggtitle(element_blank())
MRA_HH_HH_B_up<-MRA_HH_HH_B + theme(axis.title.y = element_blank(),axis.title.x = element_blank(), axis.text.y = element_blank(),
                                    plot.margin = unit(c(0,0,0,0.5),"cm")) + ggtitle(element_blank())
MRA_HN_LH_C_up<-MRA_HN_LH_C + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                                    plot.margin = unit(c(0,0,0,0.5),"cm")) + ggtitle(element_blank())
MRA_LN_LH_D_up<-MRA_LN_LH_D + theme(axis.title.y = element_blank(),
                                    plot.margin = unit(c(0,0,0,0),"cm")) + ggtitle(element_blank())



MRA_gam_LN_HH_A_up / MRA_LN_HH_A_up +
  MRA_gam_HN_HH_B_up / MRA_HH_HH_B_up +
  MRA_gam_LN_LH_D_up / MRA_LN_LH_D_up +
  MRA_gam_HN_LH_C_up / MRA_HN_LH_C_up +
  plot_layout()

MRA_panel_A<-MRA_gam_LN_HH_A_up / MRA_LN_HH_A_up +
  plot_layout(heights = c(5, 1))

MRA_panel_B<-MRA_gam_HN_HH_B_up / MRA_HH_HH_B_up +
  plot_layout(heights = c(5, 1))

MRA_panel_D<-MRA_gam_HN_LH_C_up / MRA_HN_LH_C_up +
  plot_layout(heights = c(5, 1))

MRA_panel_C<-MRA_gam_LN_LH_D_up / MRA_LN_LH_D_up +
  plot_layout(heights = c(5, 1))


# Combine the panels
combined_plot <- (MRA_panel_A | MRA_panel_B | MRA_panel_C | MRA_panel_D)

# Set the layout
MRA_complete_gam_plot<-combined_plot + plot_layout(ncol = 2,
                                                   nrow = 2, guides = "collect") +
  plot_annotation(title = paste("tensor product: edf =",round(MRA_int_table$edf,2),
                                ", p < 0.001", # If pvalue is less than 0.001
                                ", r² = ", round(MRA_int_table$rsq, 2)),
                  theme = theme(plot.title = element_text(hjust = 0.4, size = 10)))


# Bring in summary results for top of figure

ggsave(filename = "./Figures/MRA_gam_fig.pdf",
       plot = MRA_complete_gam_plot,
       dpi = 300,
       width = 9)


ggsave(filename = "./Figures/MRA_gam_fig.png",
       plot = MRA_complete_gam_plot,
       dpi = 300,
       width = 9)



# Port Phillip ####
## Data curation
PPB_dist_merged_full<-read.csv("./Data/PPB_data/PPB_complete.csv")

PPB_check_reps<-PPB_dist_merged_full %>%
  group_by(site,zone,urchin_stat,time_point,meta_treat) %>%
  dplyr::summarize(count = n())

PPB_dist_merged_ready<- PPB_dist_merged_full
PPB_dist_merged_ready$metacom_unique <- paste(PPB_dist_merged_ready$metacom_site_code, PPB_dist_merged_ready$meta_number, sep = "_")
PPB_dist_merged_ready$time_point_num<-as.numeric(PPB_dist_merged_ready$time_point_num)
PPB_ready_toanalyze<-PPB_dist_merged_ready

# Metadata for merge
PPB_metadata<-PPB_ready_toanalyze %>%
  group_by(site,zone,nutrient_treat_sp,herb_treat_sp,metacom_unique,meta_treat) %>%
  reframe()

PPB_metadata$site_zone <- paste(PPB_metadata$site, PPB_metadata$zone, sep = "_")


# ANOCOVA style LMER

PPB_cleaner<-PPB_ready_toanalyze %>%
  filter(!time_point_num == 0,
         !nutrient_treat_sp== "mid nutrients",
         !herb_treat_sp== "mid herbivore") %>%
  mutate(time_point_num = as.numeric(time_point_num),
         meta_treat = as.factor(meta_treat))


PPB_cleaner_W_zero<-PPB_ready_toanalyze %>%
  filter(!nutrient_treat_sp== "mid nutrients",
         !herb_treat_sp== "mid herbivore") %>%
  mutate(time_point_num = as.numeric(time_point_num),
         meta_treat = as.factor(meta_treat),
         nutrient_treat_sp = as.factor(nutrient_treat_sp),
         herb_treat_sp = as.factor(herb_treat_sp)
         )

# check replication
ggplot() +
  geom_point(data = PPB_cleaner_W_zero, aes(x = time_point_num, y = value, color = meta_treat)) +
  facet_grid(nutrient_treat_sp~herb_treat_sp) +
  theme_bw()

PPB_cleaner_W_zero_point<-PPB_cleaner_W_zero %>% mutate(x = time_point_num,
                                                        predicted = value,
                                                        group = meta_treat,
                                                        facet = nutrient_treat_sp,
                                                        panel = herb_treat_sp)

# New try -- try making a new concatenated variable that combines the 3 categorical variables, tests splines, then plot seperately
PPB_cleaner_W_zero_treats<-PPB_cleaner_W_zero %>%
  mutate(treatment_vector = as.factor(paste(meta_treat,nutrient_treat_sp,herb_treat_sp, sep = "_")))



PPB_cleaner_W_zero_treats_summarized<-PPB_cleaner_W_zero_treats %>%
  group_by(meta_treat,nutrient_treat_sp,herb_treat_sp,treatment_vector) %>%
  reframe() %>%
  mutate(group = treatment_vector,
         facet = nutrient_treat_sp,
         panel = herb_treat_sp)


ggplot() +
  geom_point(data = PPB_cleaner_W_zero_point, aes(x = x, y = predicted, color = group)) +
  facet_grid(panel~facet) +
  theme_bw()

# GAM ####
PPB_kvalue = 3

str(PPB_cleaner_W_zero)
# modified to get all interactions
PPB_gam_V1.5 <- gam(value ~
                  s(time_point_num, by = nutrient_treat_sp, k = PPB_kvalue, bs = "fs") +
                  s(time_point_num, by = herb_treat_sp, k = PPB_kvalue, bs = "fs") +
                  s(time_point_num, by = meta_treat, k = PPB_kvalue, bs = "fs") +
                  #te(time_point_num,meta_treat, nutrient_treat_sp, herb_treat_sp,bs = "fs", m=2) +
                  #t2(time_point_num,meta_treat, nutrient_treat_sp, bs = "fs") + #+
                  #t2(time_point_num,meta_treat, herb_treat_sp, bs = "fs") + #+
                  #t2(time_point_num,nutrient_treat_sp, herb_treat_sp, bs = "fs") + #+
                  t2(time_point_num,meta_treat, nutrient_treat_sp, herb_treat_sp, bs = "fs"), #+
                method = "REML",
                data = PPB_cleaner_W_zero) # using full dataset --> some replicates missing at T3-5

#coef(PPB_gam_V1.5)
summary(PPB_gam_V1.5)
anova(PPB_gam_V1.5)
AIC(PPB_gam_V1.5)
gam.check(PPB_gam_V1.5)
gam.vcomp(PPB_gam_V1.5)

PPB_gam_predict<-ggeffects::ggpredict(PPB_gam_V1.5,terms = c("time_point_num [all]","meta_treat","nutrient_treat_sp","herb_treat_sp"))
PPB_gam_predict %>% plot()


PPB_gam_ANCOVA<-ggplot() +
  geom_line(data = PPB_gam_predict, aes(x, predicted, color = group), size = 1) +
  geom_ribbon(data = PPB_gam_predict, aes(x = x, ymin = conf.low, ymax = conf.high, color = group, fill = group), alpha = 0.2, linewidth = 0.01) +
  geom_point(data = PPB_cleaner_W_zero_point, aes(x = x, y = predicted, color = group), alpha = 0.2) +
  facet_grid(panel~facet) +
  theme_bw() +
  guides(fill = "none") +
  scale_color_manual(values = c("black", "red", "blue")) +
  scale_fill_manual(values = c("black", "red", "blue")) +
  #removeGrid() +
  labs(x = "Time (days)", y= "Community dissimilarity (BC)", color = "Accessibility")

# Extract model estimates for figures
PPB_gam_summary <- as.data.frame(summary(PPB_gam_V1.5)$s.table)
write.csv(PPB_gam_summary, "./Tables/PPB_gam_summary.csv")

# Annotate summary for figure
PPB_gam_summary$treatment<-rownames(PPB_gam_summary)

PPB_rsq<-summary(PPB_gam_V1.5)$r.sq
PPB_int_table<-PPB_gam_summary %>%
  filter(str_detect(treatment,"t2")) %>%
  mutate(rsq = PPB_rsq)




## Plotting these ####

# Modified interactions in GAM -- use this one
PPB_sms_2 <- smooth_estimates(PPB_gam_V1.5, n = 50) %>%
  #filter(type == "Tensor product int.") %>%
  filter(.type == "Tensor product (T2)",
         .smooth == "t2(time_point_num,meta_treat,nutrient_treat_sp,herb_treat_sp)")

PPB_est_lim <- c(-1, 1) * max(abs(PPB_sms_2[[".estimate"]]), na.rm = TRUE) # old, greater range
#PPB_est_lim <- c(-0.25, 0.25) * max(abs(PPB_sms_2[[".estimate"]]), na.rm = TRUE)

color_map <- c("black", "magenta","yellow4")  # Add more colors if needed
names(color_map) <- levels(PPB_sms_2$meta_treat)

PPB_effect_plot<-
  PPB_sms_2 %>%
  mutate(nutrient_treat_sp = factor(nutrient_treat_sp),
         meta_treat = factor(meta_treat),
         herb_treat_sp = factor(herb_treat_sp)) %>%
  ggplot(aes(x = time_point_num, y = meta_treat, fill = .estimate)) +
  geom_raster(aes(x = time_point_num, y = meta_treat, fill = .estimate)) +
  geom_hline(yintercept = 1.5, color = "black", size = 0.2) +
  geom_hline(yintercept = 2.5, color = "black", size = 0.2) +
  #geom_tile(color = "black", alpha = 0.2) +
  #geom_contour(aes(z = est, y = meta_treat,group = nutrient_treat_sp, fill = NULL), colour = "black") +
  facet_grid(herb_treat_sp~ nutrient_treat_sp) +
  #facet_grid(nutrient_treat_sp~ herb_treat_sp) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  expand_limits(fill = PPB_est_lim) +
  theme_bw() +
  removeGrid() +
  theme(axis.text.y = element_text(color = color_map)) +
  labs(x = "Time (d)", y = "Accessibility", fill = "Smooth term effect")


PPB_SMS_ready<-PPB_sms_2 %>%
  mutate(nutrient_treat_sp = factor(nutrient_treat_sp),
         meta_treat = factor(meta_treat),
         herb_treat_sp = factor(herb_treat_sp),
         TDBU_treat = paste(nutrient_treat_sp,herb_treat_sp,sep = "_"))

unique(PPB_SMS_ready$TDBU_treat)

PPB_SMS_ready$meta_treat = factor(PPB_SMS_ready$meta_treat, levels=c("simple","complex","mix"))

PPB_SMS_ready<-PPB_SMS_ready %>%
  mutate(meta_treat_het = recode(meta_treat,
                                 "simple" = "high access",
                                 "complex" = "low access",
                                 "mix" = "mix"))

# Break into separate ggplot elements for plotting together
p_list_sms_PPB<- lapply(sort(unique(PPB_SMS_ready$TDBU_treat)), function(i) {

  ggplot(PPB_SMS_ready[PPB_SMS_ready$TDBU_treat==i,],aes(x = time_point_num, y = meta_treat_het, fill = .estimate)) +
    geom_raster(aes(x = time_point_num, y = meta_treat_het, fill = .estimate)) +
    #geom_contour(aes(z = est, y = meta_treat,group = nutrient_treat_sp, fill = NULL), colour = "black") +
    #facet_grid(herb_treat_sp~ nutrient_treat_sp) +
    #facet_grid(nutrient_treat_sp~ herb_treat_sp) +
    #geom_tile(color = "black", alpha = 0.1) +
    scale_fill_distiller(palette = "RdBu", type = "div",
                         breaks = c(-0.1, 0, 0.1)) +
    expand_limits(fill = PPB_est_lim) +
    labs(x = "Time (d)", y = "Accessibility", fill = "Smooth term effect") +
    theme_bw() +
    removeGrid() +
    geom_hline(yintercept = 1.5, color = "black", size = 0.2) +
    geom_hline(yintercept = 2.5, color = "black", size = 0.2) +
    geom_vline(xintercept = c(0, 50, 100,150,200,250,300), color = "black", size = 0.2) +
    #geom_vline(xintercept = c(0, 100, 200, 300), color = "black", size = 0.2) +
    theme(axis.text.y = element_text(color = color_map),
          # panel.background = element_rect(fill = NA),
          #  panel.ontop = TRUE,
          # panel.grid.minor = element_line(colour="black", size=0.2),
          #  panel.grid.major = element_blank()
    ) +
    #scale_x_continuous(minor_breaks = c(0, 100, 200)) +
    ggtitle(paste(i))


})

PPB_LN_HH_A<-p_list_sms_PPB[[3]]
PPB_HH_HH_B<-p_list_sms_PPB[[1]]
PPB_HN_LH_C<-p_list_sms_PPB[[2]]
PPB_LN_LH_D<-p_list_sms_PPB[[4]]


# ANCOVA split

PPB_gam_predict_deco<-as.data.frame(PPB_gam_predict)
PPB_gam_predict_deco$TDBU_treat <- paste(PPB_gam_predict_deco$panel,PPB_gam_predict_deco$facet,sep = "_")

PPB_cleaner_W_zero_point_deco<-PPB_cleaner_W_zero_point %>%
  mutate(TDBU_treat = paste(herb_treat_sp,nutrient_treat_sp,sep = "_"))

PPB_gam_predict_deco$group = factor(PPB_gam_predict_deco$group, levels=c("simple","complex","mix"))

# Rename metacommunity factor levels for plotting
PPB_gam_predict_deco<-PPB_gam_predict_deco %>%
  mutate(group_het = recode(group,
                                 "simple" = "high access",
                                 "complex" = "low access",
                                 "mix" = "mix"))

PPB_cleaner_W_zero_point_deco<-PPB_cleaner_W_zero_point_deco %>%
  mutate(group_het = recode(group,
                            "simple" = "high access",
                            "complex" = "low access",
                            "mix" = "mix"))

PPB_cleaner_W_zero_point_deco %>%
  group_by(time_point,group_het) %>%
  filter(!x == 0) %>%
  dplyr::summarise(mean_predicted = mean(predicted,na.rm =T),
                   x = max(x,na.rm = T),
                   count = n())


PPB_cover_estimate$TDBU_treat <- paste(PPB_cover_estimate$herb_treat_sp,PPB_cover_estimate$nutrient_treat_sp,sep = "_")
PPB_cover_estimate$x = PPB_cover_estimate$Day_count
PPB_cover_estimate$group_het <- PPB_cover_estimate$complex_treatment_het

p_list_gam_PPB<- lapply(sort(unique(PPB_gam_predict_deco$TDBU_treat)), function(i) {

  PPB_filtered_points <- PPB_cleaner_W_zero_point_deco %>%
    filter(TDBU_treat == i)  # Corrected the filtering condition

  PPB_cover_estimate_filt <- PPB_cover_estimate %>%
    filter(TDBU_treat == i)  # Corrected the filtering condition


  ggplot() +
    geom_line(data = PPB_gam_predict_deco[PPB_gam_predict_deco$TDBU_treat==i,], aes(x, predicted, color = group_het), size = 0.75) +
    geom_ribbon(data = PPB_gam_predict_deco[PPB_gam_predict_deco$TDBU_treat==i,], aes(x = x, ymin = conf.low, ymax = conf.high, color = group_het, fill = group_het), alpha = 0.1, linewidth = 0.01) +
    #geom_point(data = PPB_cleaner_W_zero_point_deco[PPB_cleaner_W_zero_point_deco$TDBU_treat==i,], aes(x = x, y = predicted, color = group), alpha = 0.2) +
    geom_point(data = PPB_filtered_points, aes(x = x, y = predicted, color = group_het), alpha = 0.2) +
    geom_segment(data = PPB_cover_estimate_filt, aes(x = x, y = 0.80, yend = 0.7, color = group_het),
                 size = 1,
                 arrow = arrow(length = unit(0.2, "cm")),
                 position = position_dodge(width = 20),
                 lineend = "butt", linejoin = "mitre",
                 show_guide = F) +
    #facet_grid(panel~facet) +
    theme_bw() +
    guides(fill = "none") +
    scale_color_manual(values = c("black", "magenta","yellow4")) +
    scale_fill_manual(values = c("black", "magenta","yellow4")) +
    #removeGrid() +
    coord_cartesian(ylim=c(0, 0.9)) +
    labs(x = "Time (days)", y= "Community dissimilarity (BC)", color = "Accessibility") +
    geom_text_repel(data = PPB_filtered_points %>%
                      group_by(time_point,group_het) %>%
                      filter(!x == 0) %>%
                      dplyr::summarise(mean_predicted = mean(predicted,na.rm =T),
                                x = max(x,na.rm = T),
                                count = n()),
                    aes(label = count, x = x, y = 0.03, color = group_het),
                    min.segment.length = Inf,
                    box.padding = 0.2,
                    size = 3) +
    ggtitle(paste(i))


  #gam_V2_ANCOVA

})

PPB_gam_LN_HH_A<-p_list_gam_PPB[[2]]
PPB_gam_HN_HH_B<-p_list_gam_PPB[[1]]
PPB_gam_HN_LH_C<-p_list_gam_PPB[[3]]
PPB_gam_LN_LH_D<-p_list_gam_PPB[[4]]




# Build combined plots
PPB_gam_LN_HH_A_up<-PPB_gam_LN_HH_A + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            # legend.position = "none",
                                            plot.margin = unit(c(0,0,0,0),"cm")) +
  ggplot2::annotation_custom(FH_NL_image,0,120,0.7,0.95) +
  annotate("text", label = "(a)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())


PPB_gam_HN_HH_B_up<-PPB_gam_HN_HH_B + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y = element_blank(),
                                            axis.text.y = element_blank(),
                                            #   legend.position = "none",
                                            plot.margin = unit(c(0,0,0,0.5),"cm")) +
  ggplot2::annotation_custom(FH_NH_image,0,120,0.7,0.95) +
  annotate("text", label = "(b)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())


PPB_gam_HN_LH_C_up<-PPB_gam_HN_LH_C + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            axis.title.y = element_blank(),
                                            axis.text.y = element_blank(),
                                            #  legend.position = "none",
                                            plot.margin = unit(c(0,0,0,0.5),"cm")) +
  ggplot2::annotation_custom(FL_NH_image,0,120,0.7,0.95) +
  annotate("text", label = "(d)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())

PPB_gam_LN_LH_D_up<-PPB_gam_LN_LH_D + theme(axis.text.x = element_blank(),
                                            axis.title.x = element_blank(),
                                            # legend.position = "none",
                                            plot.margin = unit(c(0,0,0,0),"cm")) +
  ggplot2::annotation_custom(FL_NL_image,0,120,0.7,0.95) +
  annotate("text", label = "(c)", x= 0, y = 0.9, size = 5) +
  ggtitle(element_blank())




PPB_LN_HH_A_up<-PPB_LN_HH_A + theme(axis.title.y = element_blank(),axis.title.x = element_blank(), #legend.position = "none",
                                    plot.margin = unit(c(0,0,0,0),"cm")) + ggtitle(element_blank())
PPB_HH_HH_B_up<-PPB_HH_HH_B + theme(axis.title.y = element_blank(),axis.title.x = element_blank(), axis.text.y = element_blank(), #legend.position = "none",
                                    plot.margin = unit(c(0,0,0,0.5),"cm")) + ggtitle(element_blank())
PPB_HN_LH_C_up<-PPB_HN_LH_C + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),#legend.position = "none",
                                    plot.margin = unit(c(0,0,0,0.5),"cm")) + ggtitle(element_blank())
PPB_LN_LH_D_up<-PPB_LN_LH_D + theme(axis.title.y = element_blank(),  #legend.position = "none",
                                    plot.margin = unit(c(0,0,0,0),"cm")) + ggtitle(element_blank())



#PPB_gam_LN_HH_A_up / PPB_LN_HH_A_up +
#  PPB_gam_HN_HH_B_up / PPB_HH_HH_B_up +
#  PPB_gam_LN_LH_D_up / PPB_LN_LH_D_up +
#  PPB_gam_HN_LH_C_up / PPB_HN_LH_C_up +
#  plot_layout()

PPB_panel_A<-PPB_gam_LN_HH_A_up / PPB_LN_HH_A_up +
  #plot_spacer() +
  plot_layout(heights = c(5, 1))

PPB_panel_B<-PPB_gam_HN_HH_B_up / PPB_HH_HH_B_up +
  #plot_spacer() +
  plot_layout(heights = c(5, 1))

PPB_panel_D<-PPB_gam_HN_LH_C_up / PPB_HN_LH_C_up +
  #plot_spacer() +
  plot_layout(heights = c(5, 1))

PPB_panel_C<-PPB_gam_LN_LH_D_up / PPB_LN_LH_D_up +
  #plot_spacer() +
  plot_layout(heights = c(5, 1))

#PPB_panel_A<-ggdraw() +
#  draw_plot(PPB_panel_A) +
#  draw_image(FH_NL_image, x = 0.2, y = 0.7, width = 0.3, height = 0.3)



# Combine the panels
PPB_combined_plot <- (PPB_panel_A | PPB_panel_B | PPB_panel_C | PPB_panel_D)

# Set the layout
PPB_complete_gam_plot<-PPB_combined_plot + plot_layout(ncol = 2,
                                                   nrow = 2, guides = "collect") +
  plot_annotation(title = paste("tensor product: edf =",round(PPB_int_table$edf,2),
                                #", p =",round(PPB_int_table$`p-value`,2),
                                ", p < 0.001", # If pvalue is less than 0.001
                                ", r² = ", round(PPB_int_table$rsq, 2)),
                  theme = theme(plot.title = element_text(hjust = 0.4, size = 10)))


ggsave(filename = "./Figures/PPB_gam_fig.pdf",
       plot = PPB_complete_gam_plot,
       dpi = 300,
       height = 7.16,
       width = 9)


ggsave(filename = "./Figures/PPB_gam_fig.png",
       plot = PPB_complete_gam_plot,
       dpi = 300,
       height = 7.16,
       width = 9)

# END #
