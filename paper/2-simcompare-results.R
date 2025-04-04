
rm(list=ls())

library(tidyverse)
library(cowplot)

# Load simulated data
df_steyn = read.csv("paper/outputs/2-simcompare/simulations_steyn.csv")
df_abbott = read.csv("paper/outputs/2-simcompare/simulations_abbott.csv")
df_eales = read.csv("paper/outputs/2-simcompare/simulations_eales.csv")

vars = c("rt", "Pt", "nPos")

df_data = rbind(
  df_steyn %>% filter(variable %in% vars) %>% mutate(simulation="Novel"),
  df_abbott %>% filter(variable %in% vars) %>% mutate(simulation="Abbott"),
  df_eales %>% filter(variable %in% vars) %>% mutate(simulation="Eales")
) %>%
  mutate(simulation=factor(simulation, levels=c("Novel", "Eales", "Abbott")))

rm(df_steyn, df_abbott, df_eales)

# Load model results
df_steyn_states = read.csv("paper/outputs/2-simcompare/simcompare_steyn_states.csv") %>% rename(method=model)
df_steyn_params = read.csv("paper/outputs/2-simcompare/simcompare_steyn_params.csv")
df_abbott = read.csv("paper/outputs/2-simcompare/simcompare_abbott.csv")
df_eales = read.csv("paper/outputs/2-simcompare/simcompare_eales_results.csv") %>% mutate(method="Eales")

# Combine model results
df_states = rbind(df_steyn_states %>% select(t, mean, lower, upper, variable, method, simulation, iter),
                  df_abbott %>% select(t, mean, lower, upper, variable, method, simulation, iter),
                  df_eales %>% select(t, mean, lower, upper, variable, method, simulation, iter)) %>%
  filter(variable %in% vars) %>%
  mutate(method = case_when(method=="Steyn" ~ "Novel", TRUE ~ method),
         simulation = case_when(simulation=="Steyn" ~ "Novel", TRUE ~ simulation)) %>%
left_join(df_data %>% rename(TrueValue=value), by=c("t", "variable", "simulation", "iter"))

# Calculate swab positivity (instead of number of positive swabs)
df_swabpos = df_states %>%
  filter(variable=="nPos") %>%
  mutate(variable="SwabPos", mean=100*mean/5000, lower=100*lower/5000, upper=100*upper/5000, TrueValue=100*TrueValue/5000)

df_states = rbind(df_states, df_swabpos)

# Ensure prevalence is as percentage
df_states = df_states %>%
  mutate(mean = ifelse(variable=="Pt", 100*mean, mean),
         lower = ifelse(variable=="Pt", 100*lower, lower),
         upper = ifelse(variable=="Pt", 100*upper, upper),
         TrueValue = ifelse(variable=="Pt", 100*TrueValue, TrueValue))

# Calculate coverage
df_coverage = df_states %>%
  mutate(inPredInterval = (TrueValue >= lower) & (TrueValue <= upper)) %>%
  group_by(simulation, variable, iter, method) %>%
  summarise(coverage = mean(inPredInterval, na.rm=TRUE))  %>%
  mutate(simulation = factor(simulation, levels=c("Novel", "Eales", "Abbott")),
         method = factor(method, levels=c("Novel", "Eales", "Abbott")),
         variable = case_when(variable == "rt" ~ "Growth rate\n(per day)",
                              variable == "Pt" ~ "Prevalence (%)",
                              variable == "nPos" ~ "Observed\n positive swabs",
                              variable == "SwabPos" ~ "Swab positivity (%)"),
         variable = factor(variable, levels=c("Growth rate\n(per day)", "Prevalence (%)", "Observed\n positive swabs", "Swab positivity (%)")))


df_width = df_states %>%
  # filter(!(simulation=="Eales" & variable=="Growth rate" & t==100)) %>%
  mutate(width = upper - lower) %>%
  group_by(simulation, variable, iter, method) %>%
  summarise(mean_width = mean(width, na.rm=TRUE)) %>%
  mutate(simulation = factor(simulation, levels=c("Novel", "Eales", "Abbott")),
         method = factor(method, levels=c("Novel", "Eales", "Abbott")),
         variable = case_when(variable == "rt" ~ "Growth rate\n(per day)",
                              variable == "Pt" ~ "Prevalence (%)",
                              variable == "nPos" ~ "Observed\n positive swabs",
                              variable == "SwabPos" ~ "Swab positivity (%)"),
         variable = factor(variable, levels=c("Growth rate\n(per day)", "Prevalence (%)", "Observed\n positive swabs", "Swab positivity (%)")))


custom_theme = theme_bw() + theme(strip.background = element_blank(),
                                  strip.placement="outside",
                                  legend.position = "none",
                                  plot.title = element_text(margin = margin(t = 0, b = 0)))

plt_cov = ggplot(df_coverage %>% filter(variable!="Observed\n positive swabs")) +
  geom_boxplot(aes(x=simulation, y=coverage, color=method), width=0.5, position=position_dodge2(padding = 0.2)) +
  geom_hline(yintercept=0.95, linetype="dashed") +
  facet_wrap(~variable, ncol=3, scales="free_y") +
  ylim(c(0.5, 1)) +
  xlab("") + ylab("Coverage of 95%\ncredible intervals") +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  custom_theme

plt_width = ggplot(df_width %>% filter(variable!="Observed\n positive swabs")) +
  geom_boxplot(aes(x=simulation, y=mean_width, color=method), width=0.5, position=position_dodge2(padding = 0.2)) +
  facet_wrap(~variable, ncol=3, scales="free_y") +
  xlab("Simulating model") + ylab("Average width of 95%\ncredible intervals") +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  ggh4x::facetted_pos_scales(y = list(scale_y_continuous(limits=c(0.012, 0.125)),
                                      scale_y_continuous(limits=c(0, 3)),
                                      scale_y_continuous(limits=c(0, 3)))) +
  custom_theme


plt_forlegend = plt_cov + theme(legend.position="right") + labs(color="Method")

plt_nolegend = plot_grid(plt_cov, plt_width, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1,0.1))
plt
ggsave("paper/figures/2-simcompare.png", plt, width=24, height=14, units="cm", dpi=600)
