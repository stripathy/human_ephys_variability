library(stringr)
library(lme4)
library(magrittr)
library(dplyr)
library(ggplot2)
library(stats)
library(cowplot)


### get ephys data from aibs

# download summary cell types data from aibs cell types portal
aibs_url = 'http://celltypes.brain-map.org/cell_types_specimen_details.csv'

aibs_ephys_meta = read.csv(url(aibs_url))

# get computed ephys features based on Shreejoy's custom scripts - these are mostly the precomputed features plus a couple extras
# based on the 2018 Fall release of data
aibs_ephys_orig =read.csv(file = 'data/aibs_aggregated_ephys_plus_morpho_v10.csv')

aibs_ephys = aibs_ephys_orig %>% filter(species == 'Homo Sapiens')

# relabel some ephys columns shreejoy and claire used previously for easy manipulation
aibs_ephys %<>% mutate(
  rmp = vrest, 
  rin = ri, 
  apthr = threshold_v_long_square,
  apamp = peak_v_long_square - threshold_v_long_square, 
  ahpamp = threshold_v_long_square - fast_trough_v_long_square,
  aphw = rheo_first_spike_hw * 1000, 
  apvel = upstroke_downstroke_ratio_long_square,
  rheo = threshold_i_long_square,
  maxfreq = max_rate_long_square, 
  adratio = hero_first_mean_adratio, 
  # sag = 1 - sag, # defined by aibs as vm change / max vm deflection
  sagamp.400 = (rmp -vm_for_sag) * sag,
  fislope = f_i_curve_slope, 
  avgisi = avg_isi, 
  cvisi = hero_isi_cv, 
  latency = latency * 1000 # convert to ms
)

# get some extra metadata columns
aibs_ephys_meta_small = aibs_ephys_meta %>% select(contains('donor'), 'specimen__id', 'csl__normalized_depth') %>%
  rename_at(vars(contains('donor')), funs(sub('donor__', '', .))) %>%
  rename(specimen_id = specimen__id, 
         depth_norm = csl__normalized_depth, 
         donor_id = id, donor_name = name)

aibs_ephys_meta_small %<>% filter(species == 'Homo Sapiens')

aibs_ephys_meta_small$age = lapply(aibs_ephys_meta_small$age, function(age_str) strsplit(age_str[[1]] %>% as.character, split = ' ')[[1]][1] %>% as.numeric) %>% unlist

# this is the final data frame to analyze
aibs_human_ephys = merge(aibs_ephys, aibs_ephys_meta_small, by = 'specimen_id')


### generate some basic histograms breaking down cells sampled per subject
basic_cell_histo = aibs_human_ephys %>% filter(dendrite_type %in% c('aspiny', 'spiny')) %>% group_by(donor_id, dendrite_type) %>% tally() %>% 
  ggplot(aes(x = n)) + geom_bar() + facet_wrap(~dendrite_type, scales = "free_y") + 
  xlab('Cells per patient') + ylab('patient count')

spiny_layer_histo = aibs_human_ephys %>% mutate(structure_layer_name_formal = paste('Layer', structure_layer_name)) %>%  filter(dendrite_type %in% c('spiny')) %>% group_by(donor_id, structure_layer_name_formal) %>% tally() %>% 
  ggplot(aes(x = n)) + geom_bar() + facet_wrap(~structure_layer_name_formal, scales = "free_y") + 
  xlab('Spiny cells per patient') + ylab('patient count')

basic_histos = plot_grid(basic_cell_histo, spiny_layer_histo, nrow = 2, rel_heights = c(1, 2))

ggsave('figures/basic_histograms.png', plot = basic_histos, width = 6, height = 8)

### plot ephys measures per subject for just spiny / pyramidal cells
sag_vs_donor_id = aibs_human_ephys %>% filter(dendrite_type %in% c('spiny')) %>% 
  ggplot(aes(x = reorder(donor_name, sag, FUN = median), y = sag)) + 
  geom_boxplot(outlier.alpha = 0) + geom_jitter(alpha = .5, aes(color = structure_layer_name)) + 
  xlab('Donor ID') + ylab('Sag ratio') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  guides(color=guide_legend(title = "Cortical layer"))

apvel_vs_donor_id = aibs_human_ephys %>% filter(dendrite_type %in% c('spiny')) %>% 
  ggplot(aes(x = reorder(donor_name, apvel, FUN = median), y = apvel)) + 
  geom_boxplot(outlier.alpha = 0) + geom_jitter(alpha = .5, aes(color = structure_layer_name)) + 
  xlab('Donor ID') + ylab('AP velocity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  guides(color=guide_legend(title = "Cortical layer"))

ephys_by_donors = plot_grid(sag_vs_donor_id, apvel_vs_donor_id, nrow = 2, align = "v")
ggsave('figures/ephys_by_donors.png', plot = ephys_by_donors, width = 10, height = 10)


### attempt to use mixed effects models to regress out metadata, including demographic and cell type variables
sag_data_frame = aibs_human_ephys %>% filter(dendrite_type %in% c('spiny'))
sag_model = lmer(sag ~ structure_layer_name  + (1|donor_name) + age + sex + disease_state.y, data = sag_data_frame)

sag_data_frame$sag_resids = residuals(sag_model)

sag_resids_vs_donor_id = sag_data_frame %>% filter(dendrite_type %in% c('spiny')) %>% 
  ggplot(aes(x = reorder(donor_name, sag_resids, FUN = median), y = sag_resids)) + 
  geom_boxplot(outlier.alpha = 0) + geom_jitter(alpha = .5, aes(color = structure_layer_name)) + 
  xlab('Donor ID') + ylab('Sag ratio residuals') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  guides(color=guide_legend(title = "Cortical layer"))

apvel_data_frame = aibs_human_ephys %>% filter(dendrite_type %in% c('spiny'))
apvel_model = lmer(apvel ~ structure_layer_name  + (1|donor_name) + age + sex + disease_state.y, data = apvel_data_frame)

apvel_data_frame$apvel_resids = residuals(apvel_model)

apvel_resids_vs_donor_id = apvel_data_frame %>% filter(dendrite_type %in% c('spiny')) %>% 
  ggplot(aes(x = reorder(donor_name, apvel_resids, FUN = median), y = apvel_resids)) + 
  geom_boxplot(outlier.alpha = 0) + geom_jitter(alpha = .5, aes(color = structure_layer_name)) + 
  xlab('Donor ID') + ylab('AP velocity residuals') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  guides(color=guide_legend(title = "Cortical layer"))


ephys_by_donors_resids = plot_grid(sag_resids_vs_donor_id, apvel_resids_vs_donor_id, nrow = 2, align = "v")
ggsave('figures/ephys_by_donors_resids.png', plot = ephys_by_donors_resids, width = 10, height = 10)

