
rm(list=ls())

library(tidyverse)
library(cowplot)
library(clipr)

df_results = rbind(
  read.csv("paper/outputs/s10-variants/variants_r8to13_results.csv") %>% mutate(rounds="8to13"),
  read.csv("paper/outputs/s10-variants/variants_r14to19_results.csv") %>% mutate(rounds="14to19")
)

df_params = rbind(
  read.csv("paper/outputs/s10-variants/variants_r8to13_params.csv") %>% mutate(rounds="8to13"),
  read.csv("paper/outputs/s10-variants/variants_r14to19_params.csv") %>% mutate(rounds="14to19")
)

df_rounds = read.csv("data/reactdata.csv") %>% mutate(date=as.Date(date)) %>% select(date, round)

# Tidy up params
df_params_table = df_params %>%
  mutate(m = ifelse(param=="rho", m*1e4, m),
         l = ifelse(param=="rho", l*1e4, l),
         u = ifelse(param=="rho", u*1e4, u),
         value = sprintf("%.2g (%.2g, %.2g)", m, l, u),
         param = factor(param, levels=c("sigma_WT", "sigma_AL", "sigma_DE", "sigma_OM", "sigma_BA", "rho")),
         rounds = factor(rounds, levels=c("8to13", "14to19"))) %>%
  select(value, rounds, param) %>%
  pivot_wider(names_from = param, values_from = value)

# Copy to clipboard
write_clip(df_params_table)

# Tidy up results
df_results_tidy = df_results %>%
  mutate(Variant = case_when(
    variable == "rt_WT" | variable == "Pt_WT" ~ "Wildtype",
    variable == "rt_AL" | variable == "Pt_AL" ~ "Alpha",
    variable == "rt_DE" | variable == "Pt_DE" ~ "Delta",
    variable == "rt_OM" | variable == "Pt_OM" ~ "Omicron",
    variable == "rt_BA" | variable == "Pt_BA" ~ "Omicron BA.2",
  )) %>%
  mutate(
    variable = str_extract(variable, "^.."),
    Variant = factor(Variant, levels=c("Wildtype", "Alpha", "Delta", "Omicron", "Omicron BA.2")),
    date = as.Date(date)
  ) %>%
  group_by(Variant, variable) %>%
  complete(date=seq(min(date), max(date), by="day")) %>%
  ungroup() %>%
  mutate(drop = (mean==0) & (lower==0) & (upper==0),
         mean = ifelse(drop, NA, mean),
         lower = ifelse(drop, NA, lower),
         upper = ifelse(drop, NA, upper))

# Prepare for plotting
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


plt_rt = ggplot(df_results_tidy %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Variant), alpha=0.3) +
  geom_line(aes(x=date, y=upper, color=Variant, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Variant, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Variant)) +
  custom_theme +
  scale_linetype_manual(values = c("custom" = "11")) +
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("") + ylab("Growth rate\n(per day)")
plt_rt    

plt_rt = add_rounds_to_plt(plt_rt, df_rounds)

plt_Pt = ggplot(df_results_tidy %>% filter(variable=="Pt")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Variant), alpha=0.3) +
  geom_line(aes(x=date, y=100*upper, color=Variant, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Variant, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Variant)) +
  custom_theme +
  scale_linetype_manual(values = c("custom" = "11")) +
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1, 10), labels=c(0.001, 0.01, 0.1, 1,10)) +
  xlab("Date") + ylab("Prevalence\n(log scale, %)") +
  coord_cartesian(ylim=c(0.0006, 10))

plt_Pt = add_rounds_to_plt(plt_Pt, df_rounds)

plt_forlegend = ggplot(df_results_tidy %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Variant), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Variant, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Variant, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Variant)) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme + theme(legend.position="right") + guides(linetype = "none")

plt_nolegend = plot_grid(plt_rt, plt_Pt, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.2))
plt
ggsave("paper/figures/s10-variants.png", plt, width=24, height=15, units="cm", dpi=600)





