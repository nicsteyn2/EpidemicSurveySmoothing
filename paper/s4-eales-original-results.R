
rm(list=ls())

library(tidyverse)
library(cowplot)

df_original = rbind(
  read.csv("paper/outputs/s4-eales-original/original_r1to7_results.csv") %>% mutate(Rounds="1to7", Model="Original"),
  read.csv("paper/outputs/s4-eales-original/original_r8to13_results.csv") %>% mutate(Rounds="8to13", Model="Original"),
  read.csv("paper/outputs/s4-eales-original/original_r14to16_results.csv") %>% mutate(Rounds="14to16", Model="Original"),
  read.csv("paper/outputs/s4-eales-original/original_r17to19_results.csv") %>% mutate(Rounds="17to19", Model="Original")
)

df_modified = rbind(
  read.csv("paper/outputs/s4-eales-original/modified_r1to7_results.csv") %>% mutate(Rounds="1to7", Model="Modified"),
  read.csv("paper/outputs/s4-eales-original/modified_r8to13_results.csv") %>% mutate(Rounds="8to13", Model="Modified"),
  read.csv("paper/outputs/s4-eales-original/modified_r14to16_results.csv") %>% mutate(Rounds="14to16", Model="Modified"),
  read.csv("paper/outputs/s4-eales-original/modified_r17to19_results.csv") %>% mutate(Rounds="17to19", Model="Modified")
)


df_data = read.csv("data/reactdata.csv") %>% mutate(date=as.Date(date))

df = rbind(df_original, df_modified) %>%
  mutate(date=as.Date(date),
         Model = factor(Model, levels=c("Original", "Modified"))) %>%
  left_join(df_data %>% select(date, round, nSamples), by="date") %>%
  group_by(variable, Model) %>%
  complete(date=seq(min(date), max(date), by="day")) %>%
  mutate(obsPosMean = ifelse(variable=="nPos", mean/nSamples, NA),
         obsPosLower = ifelse(variable=="nPos", lower/nSamples, NA),
         obsPosUpper = ifelse(variable=="nPos", upper/nSamples, NA))


colors = c("#c251b5", "#c46e3f")

custom_theme = theme_bw() + theme(strip.background = element_blank(),
                                  strip.placement="outside",
                                  legend.position = "none",
                                  plot.title = element_text(margin = margin(t = 0, b = 0)))


add_rounds_to_plt = function(plt, df, ymin=-10, ymax=100, ytext=0.09) {
  
  df = df %>% filter(!is.na(round))
  unique_rounds = unique(df$round)
  for (current_round in unique_rounds) {
    rnd_st = min(df$date[df$round==current_round])
    rnd_en = max(df$date[df$round==current_round])
    plt = plt + geom_rect(data=df %>% filter(date==rnd_st), xmin = rnd_st, xmax = rnd_en, ymin = ymin, ymax = ymax, fill = "black", alpha = 0.1)
    
    print(current_round)
    print(rnd_st)
    print(rnd_en)
    print("-----")
  }
  return(plt)
}


plt_rt = ggplot(df %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Model)) +
  scale_fill_manual(values=colors) + 
  scale_color_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Growth rate\n(per day)")

plt_rt = add_rounds_to_plt(plt_rt, df %>% filter(variable=="rt", Model=="Original"))

plt_Pt = ggplot(df %>% filter(variable=="Pt")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=100*upper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Model)) +
  scale_fill_manual(values=colors) + 
  scale_color_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Prevalence\n(log scale, %)") + scale_y_log10()

plt_Pt = add_rounds_to_plt(plt_Pt, df %>% filter(variable=="rt", Model=="Original"))


df_obs = df_data %>% mutate(obsPos = nPos/nSamples) %>% filter(obsPos>0)

plt_nPos = ggplot(df %>% filter(variable=="nPos")) +
  geom_ribbon(aes(x=date, ymin=100*obsPosLower, ymax=100*obsPosUpper, fill=Model), alpha=0.4) +
  geom_line(aes(x=date, y=100*obsPosMean, color=Model)) +
  geom_point(aes(x=date, y=100*obsPos), color="black", size=0.3, data=df_obs) +
  # geom_line(aes(x=date, y=100*obsPos), color="black", size=0.2, data=df_obs) +
  scale_fill_manual(values=colors) + 
  scale_color_manual(values=colors) +
  # scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Swab positivity\n(log scale, %)") + scale_y_log10() +
  coord_cartesian(ylim=c(0.003, 10))

plt_nPos = add_rounds_to_plt(plt_nPos, df %>% filter(variable=="rt", Model=="Original"))
plt_nPos

plt_forlegend = ggplot(df %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Model)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme + theme(legend.position="right") + guides(linetype = "none")

plt_nolegend = plot_grid(plt_rt, plt_Pt, plt_nPos, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.1))
plt
ggsave("paper/figures/s4-eales-original.png", plt, width=24, height=14, units="cm", dpi=600)



# Repeat but only for rounds 1-to-7
df2 = df %>% filter(Rounds=="1to7")

plt_rt = ggplot(df2 %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Model)) +
  scale_fill_manual(values=colors) + 
  scale_color_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Growth rate\n(per day)")

plt_rt = add_rounds_to_plt(plt_rt, df %>% filter(variable=="rt", Model=="Original"))

plt_Pt = ggplot(df2 %>% filter(variable=="Pt")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=100*upper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Model)) +
  scale_fill_manual(values=colors) + 
  scale_color_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Prevalence\n(log scale, %)") + scale_y_log10()

plt_Pt = add_rounds_to_plt(plt_Pt, df %>% filter(variable=="rt", Model=="Original"))


df_obs = df_data %>% filter(round <= 7) %>% mutate(obsPos = nPos/nSamples)
plt_nPos = ggplot(df2 %>% filter(variable=="nPos")) +
  geom_ribbon(aes(x=date, ymin=100*obsPosLower, ymax=100*obsPosUpper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=100*obsPosLower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=100*obsPosUpper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=100*obsPosMean, color=Model)) +
  geom_point(aes(x=date, y=100*obsPos), color="black", size=0.3, data=df_obs) +
  # geom_line(aes(x=date, y=100*obsPos), color="black", size=0.2, data=df_obs) +
  scale_fill_manual(values=colors) + 
  scale_color_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Swab positivity\n(log scale, %)") + scale_y_log10() +
  coord_cartesian(ylim=c(0.003, 10))

plt_nPos = add_rounds_to_plt(plt_nPos, df %>% filter(variable=="rt", Model=="Original"))
plt_nPos

plt_forlegend = ggplot(df %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Model), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Model, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Model)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme + theme(legend.position="right") + guides(linetype = "none")

plt_nolegend = plot_grid(plt_rt, plt_Pt, plt_nPos, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.1))
plt
ggsave("paper/figures/s4-eales-original-1to7.png", plt, width=24, height=14, units="cm", dpi=600)



# Calculate coverage

df_coverage = df %>%
  ungroup() %>%
  filter(variable=="nPos") %>%
  select(date, Model, mean, lower, upper, Rounds) %>%
  left_join(df_data %>% select(date, nSamples, nPos), by="date") %>%
  filter(!is.na(nSamples)) %>%
  mutate(inPredInterval = (lower <= nPos) & (upper >= nPos)) %>%
  group_by(Model, Rounds) %>%
  summarise(cov = mean(inPredInterval)) %>%
  pivot_wider(names_from=Rounds, values_from=cov)

