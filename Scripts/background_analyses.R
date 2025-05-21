# Background analyses

# In support of:
# Habitat attributes mediate top-down and bottom-up drivers of community development in temperate and tropical algae
# in Ecosphere XXXX
# Authors: Griffin Srednick & Stephen Swearer


# packages
library(tidyverse)
library(vegan)
library(patchwork)
library(ggExtra)


# Moorea ####
fish<-read.csv("./Data/Moorea_data/Moorea_background/MCR_LTER_Annual_Fish_Survey_20230615.csv") # available at https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-mcr.6.62
CHN<-read.csv("./Data/Moorea_data/Moorea_background/MCR_LTER_Macroalgal_CHN_2005_to_2022_20230713.csv") # available at https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-mcr.20.21

fish_df<-fish %>% filter(Site %in% c("1","2")) %>% select(1,6:8,10,12,14:18)

# Summarize and make years complete; result in wide
fish_wide<-fish_df %>%
  group_by(Taxonomy,Year,Transect,Site,Habitat,Fine_Trophic) %>%
  dplyr::summarize_if(is.numeric,sum) %>%
  select(Taxonomy,Year,Transect,Site,Habitat,Count) %>%
  pivot_wider(id_cols = c(Year,Transect,Site,Habitat),
              names_from = Taxonomy,
              values_from = Count,
              values_fill = 0) %>%
  group_by(Year,Site,Habitat) %>%
  dplyr::summarize_all(sum) %>%
  select(-Transect)

fish_biomass_wide<-fish_df %>%
  group_by(Taxonomy,Year,Transect,Site,Habitat,Fine_Trophic) %>%
  dplyr::summarize_if(is.numeric,sum) %>%
  select(Taxonomy,Year,Transect,Site,Habitat,Biomass) %>%
  pivot_wider(id_cols = c(Year,Transect,Site,Habitat),
              names_from = Taxonomy,
              values_from = Biomass,
              values_fill = 0) %>%
  group_by(Year,Site,Habitat) %>%
  dplyr::summarize_all(sum) %>%
  select(-Transect)

# Make long again for filtering
fish_long<-fish_wide %>%
  pivot_longer(!c(Year,Site,Habitat),
               names_to = "Taxonomy",
               values_to = "Count")

fish_biomass_long<-fish_biomass_wide %>%
  pivot_longer(!c(Year,Site,Habitat),
               names_to = "Taxonomy",
               values_to = "biomass")

fish_long_complete<-merge(fish_long,fish_biomass_long)

herbivores <- fish_df %>%
  #filter(Fine_Trophic == "Herbivore/Detritivore") %>%
  #filter(Coarse_Trophic == "Primary Consumer") %>%
  filter(Fine_Trophic %in% c("Herbivore/Detritivore","Browser")) %>%
  select(Taxonomy,Fine_Trophic)
  #group_by(Year,Site,Habitat,Taxonomy) %>% summarize_at(vars(Count),sum)

herbivores_complete<- fish_long_complete %>%
  filter(Taxonomy %in% herbivores$Taxonomy)


# CHN
CHN_df<-CHN %>% filter(Site %in% c("LTER 1", "LTER 2"))



## Plot herbivorous fish diversity and abundance
herbivores_sum<-herbivores_complete %>%
  group_by(Site,Year,Habitat) %>%
  dplyr::summarise(fish_abund = sum(Count, na.rm = T),
                   fish_biomass = sum(biomass,na.rm = T)) %>% # Sum herbivores
  mutate(Site = recode(Site, "1" = "LTER 1",
                       "2" = "LTER 2"))

herbivores_sum$Habitat<-factor(herbivores_sum$Habitat, levels = c("FO", "BA", "FR"))


herb_plot_df <- herbivores_sum %>%
  filter(Habitat %in% c("BA","FR"))



# calculation for fish per unit area
MRA_fish_area = 5 * 50

herb_plot_df<-herb_plot_df %>%
  mutate(fish_unit_area = fish_abund/MRA_fish_area,
         biomass_unit_area = fish_biomass/MRA_fish_area)


# Filter CHN and herbivores for only the past 2 years
# Then recode names and estimate mean ± SE

# Summary
herbivores_complete %>%
  #filter(Year) %>%
  group_by(Site,Year,Habitat) %>%
  dplyr::summarise(fish_abund = mean(Count),
                   fish_biomass = mean(biomass)) %>% # Sum herbivores
  mutate(Site = recode(Site, "1" = "LTER 1",
                       "2" = "LTER 2")) %>%
  filter(!Habitat == "FO") %>%
  ggplot(aes(x = as.factor(Site), y = fish_biomass, color = Habitat)) +
  geom_boxplot() +
  theme_bw()

herb_plot_summary<-herb_plot_df %>%
  group_by(Site,Habitat) %>%
  dplyr::summarize(
    count_mean = mean(fish_unit_area,na.rm = T),
    count_se = se_fn(fish_unit_area),
    biomass_mean = mean(biomass_unit_area,na.rm = T),
    biomass_se = se_fn(biomass_unit_area)) %>%
  mutate(Habitat_long = recode(Habitat, "BA" = "Back reef",
                       "FR" = "Fringe"),
         site_rev = recode(Site, "LTER 1" = "West",
                           "LTER 2" = "East"))

# NAs ommited and named correctly
herb_plot_df_corrected<-herb_plot_df %>%
  na.omit() %>%
  mutate(Habitat_long = recode(Habitat, "BA" = "Back reef",
                               "FR" = "Fringe"),
         site_rev = recode(Site, "LTER 1" = "West",
                           "LTER 2" = "East")) %>%
  filter(Habitat_long %in% c("Back reef", "Fringe"))

# For ordering correctly
herb_plot_df_corrected$site_rev<-factor(herb_plot_df_corrected$site_rev, levels = c("West","East"))
herb_plot_summary$site_rev<-factor(herb_plot_summary$site_rev, levels = c("West","East"))

# Show temporal replicates
MRA_herb_reps<-herb_plot_df_corrected %>%
  group_by(Habitat_long, site_rev) %>%
  mutate(Year = as.numeric(Year)) %>%
  dplyr::summarize(min_year = min(Year),
                   max_year = max(Year),
                  year_count = max_year - min_year,
                  year_long = paste("n =",year_count))

# Plot with pointrange
MRA_fish_dot_plot<-ggplot(herb_plot_summary,aes(x=Habitat_long, y = count_mean)) +
  geom_jitter(data = herb_plot_df_corrected,aes(x = Habitat_long,y = fish_unit_area)) +
  #geom_point(fill = "red", pch =21, color = "black") +
  geom_pointrange(data = herb_plot_summary, aes(ymin = count_mean - count_se, ymax = count_mean + count_se), color = "red", size = 1, alpha = 0.7) +
  theme_bw() +
  labs(x = "Habitat", y =  expression(paste("herb. fish\nabundance ", m^-2))) +
  facet_grid(~site_rev) +
  removeGrid() +
  geom_text(data=MRA_herb_reps,aes(x=Habitat_long,y=0.22,label=year_long),col="blue") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

MRA_fish_biomass_plot<-ggplot(herb_plot_summary,aes(x=Habitat_long, y = biomass_mean)) +
  geom_jitter(data = herb_plot_df_corrected,aes(x = Habitat_long,y = biomass_unit_area)) +
  geom_boxplot(data = herb_plot_df_corrected,aes(x = Habitat_long,y = biomass_unit_area),color = "red", fill= "NA", outliers = F) +
  #geom_pointrange(data = herb_plot_summary, aes(ymin = biomass_mean - biomass_se, ymax = biomass_mean + biomass_se), color = "red", size = 1) +
  theme_bw() +
  labs(x = "Habitat", y =  expression(paste("herb. fish\nbiomass (g) ", m^-2))) +
  facet_grid(~site_rev) +
  geom_text(data=MRA_herb_reps,aes(x=Habitat_long,y=17,label=year_long),col="blue") +
  removeGrid() +

  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

MRA_fish_dot_plot | MRA_fish_biomass_plot

MRA_fish_lm<-lm(fish_unit_area ~ site_rev*Habitat_long, herb_plot_df_corrected)
anova(MRA_fish_lm)

MRA_biomass_lm<-lm(biomass_unit_area ~ site_rev*Habitat_long, herb_plot_df_corrected)
anova(MRA_biomass_lm)


# Boxplotting
MRA_fish_plot<-ggplot(herb_plot_df_corrected,aes(x = Habitat_long,y = fish_unit_area))+
  geom_jitter(alpha = 0.8) +
  geom_boxplot(color = "red", fill= "NA", outliers = F) +
  theme_bw() +
  labs(x = "Habitat", y =  expression(paste("herb. fish\nabundance ", m^-2))) +
  geom_text(data=MRA_herb_reps,aes(x=Habitat_long,y=0.22,label=year_long),col="blue") +
  facet_grid(~site_rev) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())


herb_plot_df_corrected %>%
  filter(Year > 2011) %>%
  ggplot(aes(x = Habitat_long,y = fish_unit_area))+
  geom_jitter(alpha = 0.8) +
  geom_boxplot(color = "red", fill= "NA", outliers = F) +
  theme_bw() +
  labs(x = "Habitat", y =  expression(paste("herb. fish\nabundance ", m^-2))) +
  #geom_text(data=MRA_herb_reps,aes(x=Habitat_long,y=0.22,label=year_long),col="blue") +
  facet_grid(~site_rev) +
  theme(#axis.title.x = element_blank(),
        #axis.text.x = element_blank()
        )


herb_plot_df_corrected %>%
  filter(Year > 2015) %>%
  ggplot(aes(color = Habitat_long,y = fish_abund, x = Year, group = Habitat_long)) +
  geom_point() +
  geom_path() +
  theme_bw() +
  facet_grid(~site_rev)

herb_plot_df_corrected %>%
  filter(Year > 2015) %>%
  ggplot(aes(color = Habitat_long,y = biomass_unit_area, x = Year, group = Habitat_long)) +
  geom_point() +
  geom_path() +
  theme_bw() +
  facet_grid(~site_rev)

# Ok. Maybe what we do is use the long-term trends to show that there is a general difference between habitats and sides of the pass
# Instead, use paired years to deal with this
herb_plot_df_corrected_filt<-herb_plot_df_corrected %>% filter(Year > 2019)
biomass_lmer<-lmer(biomass_unit_area ~ Habitat_long * site_rev + (1|Year),herb_plot_df_corrected_filt)
summary(biomass_lmer)
anova(biomass_lmer)
ggpredict(biomass_lmer, terms = c("site_rev","Habitat_long")) %>% plot() # so we get a significant interaction here. Filtered from 2012 to 2022
# USE THIS AND MAKE NEW FIG. THEN MOVE ON.


herb_plot_df_corrected_abund_filt<-herb_plot_df_corrected %>% filter(Year > 2019)
abund_lmer<-lmer(fish_unit_area ~ Habitat_long * site_rev + (1|Year),herb_plot_df_corrected_abund_filt)
summary(abund_lmer)
anova(abund_lmer)
ggpredict(abund_lmer, terms = c("Habitat_long","site_rev")) %>% plot() # so we get a significant interaction here. Filtered from 2012 to 2022

abund_lm<-lm(fish_unit_area ~ Habitat_long * site_rev,herb_plot_df_corrected_abund_filt)
summary(abund_lmer)
anova(abund_lmer)

## CHN ####

# Summarize within site reps to size
CHN_df_reduced <- CHN_df %>%
  filter(Genus == "Turbinaria") %>%
  group_by(Year,Site,Habitat) %>%
  dplyr::summarize_all(mean,na.rm = T)


CHN_long<-CHN_df_reduced[-c(4,5,10)] %>%
  gather(variable, amount,-c(1:3)) %>%
  filter(variable == "N")

# reorder for plotting
CHN_long$Habitat<-factor(CHN_long$Habitat, levels = c("Reef Crest", "Back Reef", "Fringe"))

# Assess annual variation; variable pre-2017. Use later -- more recent is more relevant
CHN_long %>%
  ggplot(aes(x=Year, y = amount, color = Habitat)) +
  stat_summary(fun.data = mean_se,geom = "errorbar",width = 0.01) +
  stat_summary(geom = "point",fun = "mean",size = 6) +
  stat_summary(geom = "line",
               fun = "mean")  +
  facet_grid(~Site) +
  labs(x = "Habitat", y = "N (%)") +
  theme_bw() +
  removeGrid()

# Get >2017 and rename; then get annual replicate range;
MRA_CHN_reps<-CHN_long %>%
  filter(Year > 2018) %>%
  mutate(Habitat_long = recode(Habitat, "Back Reef" = "Back reef",
                               "Fringe" = "Fringe"),
         site_rev = recode(Site, "LTER 1" = "West",
                           "LTER 2" = "East")) %>%
  filter(Habitat_long %in% c("Back reef", "Fringe")) %>%
  group_by(site_rev,Habitat_long) %>%
  dplyr::summarize(min_year = min(Year),
            max_year = max(Year),
            max_year - min_year,
            year_count = n(),
            year_long = paste("n =",year_count))



CHN_plot<-CHN_long %>%
  filter(Year > 2018,
         Habitat %in% c("Back Reef", "Fringe")) %>%
  na.omit() %>%
  group_by(Site,Habitat) %>%
  dplyr::summarize(mean = mean(amount,na.rm = T),
            se = se_fn(amount),
            year_count = n()) %>%
  mutate(Habitat_long = recode(Habitat, "Back Reef" = "Back reef",
                               "Fringe" = "Fringe"),
         site_rev = recode(Site, "LTER 1" = "West",
                                      "LTER 2" = "East"))
CHN_corrected<-CHN_long %>%
  filter(Year > 2018) %>%
  na.omit() %>%
  mutate(Habitat_long = recode(Habitat, "Back Reef" = "Back reef",
                               "Fringe" = "Fringe"),
         site_rev = recode(Site, "LTER 1" = "West",
                           "LTER 2" = "East")) %>%
  filter(Habitat_long %in% c("Back reef", "Fringe"))

CHN_corrected$site_rev<-factor(CHN_corrected$site_rev, levels = c("West","East"))
CHN_plot$site_rev<-factor(CHN_plot$site_rev, levels = c("West","East"))

# Plot with pointrange
MRA_CHN_plot<-CHN_plot %>%
  ggplot(aes(x=Habitat_long, y = mean)) +
  geom_jitter(data = CHN_corrected, aes(x = Habitat_long,y = amount)) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se), color = "red", size = 1) +
  theme_bw() +
  labs(x = "Habitat", y = "N (%)") +
  facet_grid(~site_rev) +
  removeGrid() +
  theme(strip.text.x = element_blank())



# boxplot
MRA_CHN_plot<-ggplot(CHN_corrected,aes(x = Habitat_long,y = amount))+
  geom_jitter(alpha = 0.95) +
  geom_boxplot(color = "red", fill= "NA", outliers = F) +
  theme_bw() +
  labs(x = "Habitat", y = "N (%)") +
  geom_text(data=MRA_CHN_reps,aes(x=Habitat_long,y=1,label=year_long),col="blue") +
  ylim(c(0,1)) +
  facet_grid(~site_rev) +
  theme(strip.text.x = element_blank())


#MRA_background_plot <- MRA_fish_plot / MRA_CHN_plot + plot_annotation(tag_levels = "a")
MRA_background_plot <- MRA_fish_dot_plot / MRA_CHN_plot + plot_annotation(tag_levels = "a")
MRA_background_plot <- MRA_fish_biomass_plot / MRA_CHN_plot + plot_annotation(tag_levels = "a")

ggsave(MRA_background_plot,
       filename = "./Figures/MRA_background_updated.jpeg",
       dpi = 300,
       width = 5,
       height = 5)



# Port Phillip ####
# Nutrients
PPB_nutrients<-read.csv("./Data/PPB_data/PPB_background/1984_07-2023_06_Port_Phillip_Bay_Water_quality_data.csv") # Data are available here: https://discover.data.vic.gov.au/dataset/epa-port-phillip-bay-water-quality-data-1984-2024

unique(PPB_nutrients$site_name_short)
PPB_nutrients$year<-format(as.Date(PPB_nutrients$date, format="%m/%d/%y"),"%Y")

# Summarized to mean ± SE for pointrange
PPB_nutrients_summarized<-PPB_nutrients %>%
  filter(site_name_short %in% c("Dromana","Hobsons Bay","Patterson River")) %>%
  select(year,site_name_short,N_TOTAL) %>%
  group_by(site_name_short) %>%
  na.omit() %>%
  dplyr::summarize(mean_N = mean(N_TOTAL,na.rm = T),
                   se_N = se_fn(N_TOTAL))

# Annual duration
PPB_nutrients_year_range<-PPB_nutrients %>%
  filter(site_name_short %in% c("Dromana","Hobsons Bay","Patterson River")) %>%
  select(year,site_name_short) %>%
  group_by(site_name_short) %>%
  mutate(year = as.numeric(year)) %>%
  dplyr::summarize(min_year = min(year),
            max_year = max(year),
            year_count = n(),
            year_long = paste("n =",year_count))

# Number of surveys
PPB_nutrients_survey_count<-PPB_nutrients %>%
  filter(site_name_short %in% c("Dromana","Hobsons Bay","Patterson River")) %>%
  select(year,site_name_short) %>%
  group_by(site_name_short,year) %>%
  mutate(year = as.numeric(year)) %>%
  dplyr::summarize(survey_count = n()) %>%
  group_by(site_name_short) %>%
  dplyr::summarize(min_survey = min(survey_count),
            max_survey = max(survey_count),
            survey_long = paste("n =",max_survey))

# Number of surveys total
PPB_nutrients_survey_count_full<-PPB_nutrients %>%
  filter(site_name_short %in% c("Dromana","Hobsons Bay","Patterson River")) %>%
  group_by(site_name_short) %>%
  dplyr::summarize(survey_count = n()) %>%
  dplyr::mutate(survey_long = paste("n =",survey_count))

# Corrected for boxplotting
PPB_nutrients_corrected<-PPB_nutrients %>%
  filter(site_name_short %in% c("Dromana","Hobsons Bay","Patterson River")) %>%
  select(year,site_name_short,N_TOTAL) %>%
  group_by(site_name_short) %>%
  na.omit()

# ordering
PPB_nutrients_corrected$site_name_short = factor(PPB_nutrients_corrected$site_name_short, levels=c("Hobsons Bay","Patterson River", "Dromana"))
PPB_nutrients_summarized$site_name_short = factor(PPB_nutrients_summarized$site_name_short, levels=c("Hobsons Bay","Patterson River", "Dromana"))

# plot pointrange
PPB_nutrients_plot<-ggplot(PPB_nutrients_summarized,aes(x = site_name_short,y = mean_N))+
  geom_jitter(data = PPB_nutrients_corrected,aes(x = site_name_short,y = N_TOTAL),alpha = 0.4) +
  geom_pointrange(aes(ymin = mean_N - se_N, ymax = mean_N + se_N), color = "red", size = 1) +
  theme_bw() +
  labs(x = "site", y = "ug/L")


# boxplots
PPB_nutrients_plot<-ggplot(PPB_nutrients_corrected,aes(x = site_name_short,y = N_TOTAL)) +
  geom_jitter(alpha = 0.4) +
  geom_boxplot(color = "red", fill= "NA", outliers = F) +
  theme_bw() +
  geom_text(data=PPB_nutrients_survey_count_full,aes(x=site_name_short,y=800,label=survey_long),col="blue") +
  labs(x = "site", y = "ug/L")


## Herbivores ####
PPB_urchin_data<-read.csv("./Data/PPB_data/PPB_background/Stock_assessment/PPB_UrchinStockAssessment_Data_20200402_.csv") # These data are embargoed. Please contact authors for access.

# plot data
PPB_urchins<-PPB_urchin_data %>%
  filter(Region %in% c("Williamstown","Martha Point","St Kilda"),
         Year == 2019) %>%
  select(Year,Region,Large.urchin.density...2.) %>%
  group_by(Region) %>%
  na.omit() %>%
  dplyr::summarize(mean = mean(Large.urchin.density...2.,na.rm = T),
                   se = se_fn(Large.urchin.density...2.))

# Number of surveys
PPB_urchin_range<-PPB_urchin_data %>%
  filter(Region %in% c("Williamstown","Martha Point","St Kilda"),
         Year == 2019) %>%
  na.omit() %>%
  group_by(Region) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::mutate(survey_long = paste("n =",count))


# Ordering
PPB_urchins$Region = factor(PPB_urchins$Region, levels=c("Williamstown","St Kilda", "Martha Point"))

# Pointrange
PPB_urchin_plot<-ggplot(PPB_urchins,aes(x = Region,y = mean))+
  geom_point() +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se)) +
  theme_bw() +
  labs(x = "site", y = expression(paste("urchin abundance ", m^-2)))

# Boxplots instead
PPB_urchin_corrected<-PPB_urchin_data %>%
  filter(Region %in% c("Williamstown","Martha Point","St Kilda"),
         Year == 2019) %>%
  select(Year,Region,Large.urchin.density...2.) %>%
  group_by(Region) %>%
  na.omit()

PPB_urchin_corrected$Region = factor(PPB_urchin_corrected$Region, levels=c("Williamstown","St Kilda", "Martha Point"))

PPB_urchin_plot<-ggplot(PPB_urchin_corrected,aes(x = Region,y = Large.urchin.density...2.)) +
  geom_jitter(alpha = 0.4) +
  geom_boxplot(color = "red", fill= "NA", outliers = F) +
  geom_text(data=PPB_urchin_range,aes(x=Region,y=35,label=survey_long),col="blue") +
  theme_bw() +
  labs(x = "site", y = expression(paste("urchin abundance ", m^-2)))




PPB_background_plot <- PPB_urchin_plot  + PPB_nutrients_plot + plot_annotation(tag_levels = "a")

ggsave(PPB_background_plot,
       filename = "./Figures/PPB_background.jpeg",
       dpi = 300,
       width = 7,
       height = 3.5)




# Urchin culling in PPB ####
urch_cull<-read.csv("./Data/PPB_data/PPB_background/jawbone_urchins_2023-2024.csv") # These data are embargoed. Please contact authors for access.

# Estimate urchin culling at Jawebone
urch_cull_est<-urch_cull %>%
  group_by(time_p,transect,park,area,site) %>%
  dplyr::summarise(mean_dens = mean(count,na.rm = T)) %>%
  filter(site %in% c("barren","cull"),
         transect %in% c("T3","T4"),
         area == "west")

urch_cull_est$time_p = factor(urch_cull_est$time_p, levels=c("precull","postcull_1m","postcull_2m","postcull_3m"))

# What are the urchin densities after cull?
urch_cull_est %>%
  group_by(time_p,park,area,site) %>%
  dplyr::summarise(mean_dens_new = mean(mean_dens,na.rm = T),
                   sd_dens = sd(mean_dens,na.rm =T))






# MRA - herbivory assays ####
MRA_trial_data_ready<-read.csv("./Data/Moorea_data/MRA_trial_data.csv")


MRA_trial_data_ready_cover<-MRA_trial_data_ready %>%
  select(-c(CCA,BEM)) %>%
  mutate(total_cover = rowSums(.[16:34])) %>%
  select(tile,time_point,site,zone,com_treat,period,total_cover)

MRA_trial_data_ready_cover %>%
  mutate(side = ifelse(site %in% c(1,2), "west","east")) %>%
  filter(!is.na(period)) %>%
  ggplot(aes(x = zone, y = total_cover, color = period)) +
  geom_boxplot() +
  facet_grid(~side) +
  theme_bw()


MRA_trial_data_analyze <- MRA_trial_data_ready_cover %>%
  filter(!is.na(period),
         #!(time_point == "E3_1" | tile == "275")
         ) %>%
  pivot_wider(id_cols = c(site,zone,com_treat,time_point,tile),
              names_from = "period",
              values_from = "total_cover",
              values_fn = NULL) %>%
  na.omit() %>%
  filter(!initial == "NULL",
         !final == "NULL") %>%
  mutate(initial = as.numeric(initial),
         final = as.numeric(final))

MRA_trial_data_analyze_act<-MRA_trial_data_analyze %>%
  mutate(log_RR = log(final/initial),
         RR = final/initial,
         side = as.factor(ifelse(site %in% c(1,2), "west","east")),
         habitat = as.factor(ifelse(zone == "C", "Fringe","Backreef"))) %>%
  filter(com_treat =="C") %>%
  filter(!zone == "B")# %>%
  #dplyr::filter(site %in% c("2","4"))


MRA_trial_data_analyze_act %>%
  ggplot(aes(x = habitat, y = log_RR, color = habitat)) +
  geom_boxplot() +
  facet_grid(~side) +
  theme_bw()

trial_data_lm<-lm(RR~habitat*side,MRA_trial_data_analyze_act)
anova(trial_data_lm)
ggpredict(trial_data_lm,terms = c("side","habitat")) %>% plot()


#ggpredict(trial_data_lm,terms = c("side","habitat")) %>% plot()

# So the trials dont show any differences in ambient herbivory across sites
# ultimately this is all about which side of the pass you are on.
#
# END #
