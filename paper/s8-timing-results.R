
rm(list=ls())

library(tidyverse)
library(foreach)
library(cowplot)

# Load all files
df_eales = foreach(ii = c(seq(1, 12), 15, 17, 19), .combine=rbind) %do% {
  if (ii < 19) {
    read.csv(paste0("paper/outputs/s8-timing/eales_R1to", ii, "_warmup200_sampling300_treedepth15.csv"))
  } else {
    read.csv(paste0("paper/outputs/s8-timing/eales_R1to", ii, "_warmup200_sampling300_treedepth16.csv"))
  }
} %>% mutate(Model = "Eales")

df_abbott = foreach(ii = c(seq(1, 13), 15, 17, 19), .combine=rbind) %do% {
  read.csv(paste0("paper/outputs/s8-timing/abbott_R1to", ii, "_warmup200_sampling300.csv"))
} %>% mutate(Model = "Abbott")

df_steyn = foreach(ii = c(seq(1, 13), 15, 17, 19), .combine=rbind) %do% {
  read.csv(paste0("paper/outputs/s8-timing/steyn_R1to", ii, ".csv"))
} %>% mutate(Model = "Steyn")

df_steyn_Rt = foreach(ii = c(seq(1,7), 9, 13, 15, 17, 19), .combine=rbind) %do% {
  if (ii < 13) {
    read.csv(paste0("paper/outputs/s8-timing/steynRt_R1to", ii, ".csv"))
  } else {
    read.csv(paste0("paper/outputs/s8-timing/steynRt_R1to", ii, "_2000.csv"))
  }
} %>% mutate(Model = "Steyn_Rt")

# Temporary fixes
df_steyn$start_round[is.na(df_steyn$start_round)] = 1
df_steyn$end_round[is.na(df_steyn$end_round)] = 2

# Remove all failures after the last success for eales and abbott
abbott_maxiter = df_abbott %>%
  group_by(end_round, success) %>%
  summarise(maxiter = max(iter)) %>%
  filter(success) %>%
  select(-success)

df_abbott = df_abbott %>% left_join(abbott_maxiter, by="end_round") %>% filter(iter<=maxiter)


eales_maxiter = df_eales %>%
  group_by(end_round, success) %>%
  summarise(maxiter = max(iter)) %>%
  filter(success) %>%
  select(-success)

df_eales = df_eales %>% left_join(eales_maxiter, by="end_round") %>% filter(iter<=maxiter)



# Append the number of total successes to abbott and eales
n_success_abbott = df_abbott %>% group_by(end_round) %>% summarise(n_success = sum(success))
df_abbott = left_join(df_abbott, n_success_abbott, by="end_round")

n_success_eales = df_eales %>% group_by(end_round) %>% summarise(n_success = sum(success))
df_eales = left_join(df_eales, n_success_eales, by="end_round")


# Calculate summary statistics
abbott_time_forsuccess = df_abbott %>% group_by(end_round, success, num_obs, num_days) %>% summarise(mean_time = mean(time_taken), lower_time = min(time_taken), upper_time = max(time_taken), n = n()) %>% filter(success)
abbott_time_persuccess = df_abbott %>% group_by(end_round, num_obs, num_days) %>% summarise(total_time = sum(time_taken), n_success = min(n_success)) %>% mutate(mean_time = total_time/n_success)

eales_time_forsuccess = df_eales %>% group_by(end_round, success, num_obs, num_days) %>% summarise(mean_time = mean(time_taken), lower_time = min(time_taken), upper_time = max(time_taken), n = n()) %>% filter(success)
eales_time_persuccess = df_eales %>% group_by(end_round, num_obs, num_days) %>% summarise(total_time = sum(time_taken), n_success = min(n_success)) %>% mutate(mean_time = total_time/n_success)

steyn_time_forsuccess = df_steyn %>% group_by(end_round, num_obs, num_days) %>% summarise(mean_time = mean(time_taken), lower_time = min(time_taken), upper_time = max(time_taken), n = n())

steynRt_time_forsuccess = df_steyn_Rt %>% group_by(end_round, num_obs, num_days) %>% summarise(mean_time = mean(time_taken), lower_time = min(time_taken), upper_time = max(time_taken), n = n())

df_plt_forsuccess = rbind(
  abbott_time_forsuccess %>% mutate(Model = "Abbott"),
  eales_time_forsuccess %>% mutate(Model = "Eales"),
  steyn_time_forsuccess %>% mutate(Model = "SIMPLE"),
  steynRt_time_forsuccess %>% mutate(Model = "SIMPLE (Rt model)")
) %>% mutate(Model=factor(Model, levels=c("SIMPLE", "SIMPLE (Rt model)", "Eales", "Abbott")))

df_plt_persuccess = rbind(
  abbott_time_persuccess %>% mutate(Model = "Abbott"),
  eales_time_persuccess %>% mutate(Model = "Eales"),
  steyn_time_forsuccess %>% mutate(Model = "SIMPLE"),
  steynRt_time_forsuccess %>% mutate(Model = "SIMPLE (Rt model)")
) %>% mutate(Model=factor(Model, levels=c("SIMPLE", "SIMPLE (Rt model)", "Eales", "Abbott")))


custom_theme = theme_bw() + theme(strip.background = element_blank(),
                                  strip.placement="outside",
                                  legend.position = "none",
                                  plot.title = element_text(margin = margin(t = 0, b = 0)))

plt = ggplot(df_plt_forsuccess) +
  geom_point(aes(x=num_days, y=mean_time, color=Model, shape=Model), position=position_dodge(5)) +
  geom_errorbar(aes(x=num_days, ymin=lower_time, ymax=upper_time, color=Model), width=5, position=position_dodge(5)) +
  geom_line(aes(x=num_days, y=mean_time, color=Model, linetype=Model), position=position_dodge(5)) +
  scale_color_manual(values=c("#2d96b3", "#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values=c("solid", "dashed", "solid", "solid")) +
  scale_shape_manual(values=c(16, 4, 16, 16)) +
  xlab("") + ylab("Time for a successful fit\n(seconds)") +
  custom_theme
plt


plt2 = ggplot(df_plt_persuccess) +
  geom_point(aes(x=num_days, y=mean_time, color=Model, shape=Model)) +
  geom_line(aes(x=num_days, y=mean_time, color=Model, linetype=Model)) +
  scale_color_manual(values=c("#2d96b3","#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values=c("solid", "dashed", "solid", "solid")) +
  scale_shape_manual(values=c(16, 4, 16, 16)) +
  xlab("Length of data (total number of days)") + ylab("Average total time per successful fit\n(seconds)") +
  custom_theme
plt2

plt_forlegend = ggplot(df_plt_persuccess) +
  geom_point(aes(x=num_days, y=mean_time, color=Model, shape=Model)) +
  geom_line(aes(x=num_days, y=mean_time, color=Model, linetype=Model)) +
  scale_color_manual(values=c("#2d96b3","#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values=c("solid", "dashed", "solid", "solid")) +
  scale_shape_manual(values=c(16, 4, 16, 16)) +
  theme_bw()

plt_both = plot_grid(plt, plt2, ncol=1, align="v")
plt_both

plt_withlegend = plot_grid(plt_both, get_legend(plt_forlegend), ncol=2, rel_widths=c(1,0.2))
plt_withlegend
ggsave("paper/figures/s8-timing.png", plt_withlegend, width=24, height=14, units="cm", dpi=600)
