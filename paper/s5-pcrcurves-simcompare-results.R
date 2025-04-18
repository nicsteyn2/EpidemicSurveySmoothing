
rm(list=ls())

library(tidyverse)
library(cowplot)

# Load simulated data
df_steyn = read.csv("paper/outputs/2-simcompare/simulations_steyn.csv")
df_abbott = read.csv("paper/outputs/2-simcompare/simulations_abbott.csv")
df_eales = read.csv("paper/outputs/2-simcompare/simulations_eales.csv")

vars = c("rt", "Pt", "nPos")

df_data = rbind(
  df_steyn %>% filter(variable %in% vars) %>% mutate(simulation="SIMPLE"),
  df_abbott %>% filter(variable %in% vars) %>% mutate(simulation="Abbott"),
  df_eales %>% filter(variable %in% vars) %>% mutate(simulation="Eales")
) %>%
  mutate(simulation=factor(simulation, levels=c("SIMPLE", "Eales", "Abbott")))

rm(df_steyn, df_abbott, df_eales)

# Load model results
df_steyn_states = read.csv("paper/outputs/2-simcompare/simcompare_steyn_states.csv") %>% rename(method=model)
df_steyn_params = read.csv("paper/outputs/2-simcompare/simcompare_steyn_params.csv")
df_abbott = read.csv("paper/outputs/2-simcompare/simcompare_abbott.csv")  %>% mutate(method="Abbott (Hellewell et al. (2021))")
df_abbott_binny = read.csv("paper/outputs/s5-pcrcurves/simulations_abbott_binnypcr.csv") %>% mutate(method="Abbott (Binny et al. (2023))")
df_eales = read.csv("paper/outputs/2-simcompare/simcompare_eales_results.csv") %>% mutate(method="Eales")

# Combine model results
df_states = rbind(df_steyn_states %>% select(t, mean, lower, upper, variable, method, simulation, iter),
                  df_abbott %>% select(t, mean, lower, upper, variable, method, simulation, iter),
                  df_abbott_binny %>% select(t, mean, lower, upper, variable, method, simulation, iter),
                  df_eales %>% select(t, mean, lower, upper, variable, method, simulation, iter)) %>%
  filter(variable %in% vars) %>%
  mutate(method = case_when(method=="Steyn" ~ "SIMPLE", TRUE ~ method),
         simulation = case_when(simulation=="Steyn" ~ "SIMPLE", TRUE ~ simulation)) %>%
  left_join(df_data %>% rename(TrueValue=value), by=c("t", "variable", "simulation", "iter"))


# Set prevalence as %
df_states = df_states %>%
  mutate(TrueValue = case_when(variable=="Pt" ~ TrueValue*100, TRUE ~ TrueValue),
         mean = case_when(variable=="Pt" ~ mean*100, TRUE ~ mean),
         lower = case_when(variable=="Pt" ~ lower*100, TRUE ~ lower),
         upper = case_when(variable=="Pt" ~ upper*100, TRUE ~ upper))


# Plot results to check
plt = ggplot(df_states %>% filter(simulation=="Eales", iter==2)) +
  geom_ribbon(aes(x=t, ymin=lower, ymax=upper, fill=method), alpha=0.4) +
  geom_line(aes(x=t, y=mean, color=method)) +
  geom_line(aes(x=t, y=TrueValue), color="black") +
  facet_wrap(~variable, scales="free_y")
plt

# Calculate coverage
df_coverage = df_states %>%
  mutate(inPredInterval = (TrueValue >= lower) & (TrueValue <= upper)) %>%
  group_by(simulation, variable, iter, method) %>%
  summarise(coverage = mean(inPredInterval, na.rm=TRUE))  %>%
  mutate(simulation = factor(simulation, levels=c("SIMPLE", "Eales", "Abbott")),
         method = factor(method, levels=c("SIMPLE", "Eales", "Abbott (Hellewell et al. (2021))", "Abbott (Binny et al. (2023))")),
         variable = case_when(variable == "rt" ~ "Growth rate\n(per day)",
                              variable == "Pt" ~ "Prevalence (%)",
                              variable == "nPos" ~ "Observed\n positive swabs"),
         variable = factor(variable, levels=c("Growth rate\n(per day)", "Prevalence (%)", "Observed\n positive swabs")))

df_width = df_states %>%
  mutate(width = upper - lower) %>%
  group_by(simulation, variable, iter, method) %>%
  summarise(mean_width = mean(width, na.rm=TRUE)) %>%
  mutate(simulation = factor(simulation, levels=c("SIMPLE", "Eales", "Abbott")),
         method = factor(method, levels=c("SIMPLE", "Eales", "Abbott (Hellewell et al. (2021))", "Abbott (Binny et al. (2023))")),
         variable = case_when(variable == "rt" ~ "Growth rate\n(per day)",
                              variable == "Pt" ~ "Prevalence (%)",
                              variable == "nPos" ~ "Observed\n positive swabs"),
         variable = factor(variable, levels=c("Growth rate\n(per day)", "Prevalence (%)", "Observed\n positive swabs")))


custom_theme = theme_bw() + theme(strip.background = element_blank(),
                                  strip.placement="outside",
                                  legend.position = "none",
                                  plot.title = element_text(margin = margin(t = 0, b = 0)))

plt_cov = ggplot(df_coverage) +
  geom_boxplot(aes(x=simulation, y=coverage, color=method), width=0.5, position=position_dodge2(padding = 0.2)) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  facet_wrap(~variable, ncol=3, scales="free_y") +
  ylim(c(0, 1)) +
  xlab("") + ylab("Coverage of 95%\ncredible intervals") +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f", "#274200")) +
  custom_theme

plt_width = ggplot(df_width) +
  geom_boxplot(aes(x=simulation, y=mean_width, color=method), width=0.5, position=position_dodge2(padding = 0.2)) +
  facet_wrap(~variable, ncol=3, scales="free_y") +
  xlab("Simulating model") + ylab("Average width of 95%\ncredible intervals") +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f", "#274200")) +
  custom_theme

plt_forlegend = plt_cov + theme(legend.position="right", legend.direction="horizontal") + labs(color="Method")

plt_nolegend = plot_grid(plt_cov, plt_width, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=1, rel_heights=c(1,0.1))
plt
ggsave("paper/figures/s5-pcrcurves-simcompare.png", plt, width=24, height=14, units="cm", dpi=600)
