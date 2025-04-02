
rm(list=ls())

library(tidyverse)
library(foreach)
library(cowplot)

config = read.csv("paper/1-dynamics-configtable.csv") %>%
  mutate(Simulation = case_when(
    xi == 0 & rho == 0 ~ "Basic simulation",
    xi == 0 & rho == 2e-4 ~ "Extra-binomial simulation",
    xi == 0.5 & rho == 0 ~ "Weighted simulation"
  ))

custom_theme = theme_bw() + theme(strip.background = element_blank(),
                                  strip.placement="outside",
                                  legend.position = "none",
                                  plot.title = element_text(margin = margin(t = 0, b = 0)))

# Load data
confignums = seq(1,30)
df_states_in = foreach(ii = seq(1, length(confignums)), .combine=rbind) %do% {
  read.csv(sprintf("paper/outputs/1-dynamics/1-dynamics-config%d-states.csv", confignums[ii])) %>% mutate(confignum=confignums[ii])
}

df_params = foreach(ii = seq(1, length(confignums)), .combine=rbind) %do% {
  read.csv(sprintf("paper/outputs/1-dynamics/1-dynamics-config%d-params.csv", confignums[ii])) %>% mutate(confignum=confignums[ii])
}

df_data_in = foreach(ii = seq(1, length(confignums)), .combine=bind_rows) %do% {
  read.csv(sprintf("paper/outputs/1-dynamics/1-dynamics-config%d-data.csv", confignums[ii])) %>% mutate(confignum=confignums[ii])
}

df_data = df_data_in %>%
  select(-obsmethod) %>%
  pivot_longer(cols=-c(t, iter, confignum), names_to="variable", values_to="TrueValue")

df_states = df_states_in %>% left_join(df_data, by=c("t", "confignum", "iter", "variable")) %>% left_join(config, by="confignum")

rm(df_states_in, df_data_in)


# Calculate coverage
df_plt_byiter = df_states %>%
  mutate(inPredInterval = (lower <= TrueValue) & (upper >= TrueValue),
         width = (upper-lower),
         normalised_width = (upper-lower)/abs(mean)) %>%
  group_by(iter, params, confignum, model, variable, rho, xi, sigma, nSamples, Simulation) %>%
  summarise(cov=mean(inPredInterval),
            width = 100*mean(width),
            widthnrom = mean(normalised_width)) %>%
  ungroup()


df_plt_overall = df_states %>%
  mutate(inPredInterval = (lower <= TrueValue) & (upper >= TrueValue),
         width = (upper-lower),
         normalised_width = (upper-lower)/abs(mean)) %>%
  group_by(params, confignum, model, variable, rho, xi, sigma, nSamples, Simulation) %>% # Not grouping by iter
  summarise(meancov=mean(inPredInterval),
            meanwidth = 100*mean(width),
            meanwidthnorm = mean(normalised_width)) %>%
  ungroup() %>%
  select(params, confignum, model, variable, meancov, meanwidth, meanwidthnorm, Simulation)

df_plt = df_plt_byiter %>%
  left_join(df_plt_overall, by=c("params", "confignum", "model", "variable", "Simulation"))%>%
  filter(params == "estimated", # Only comparing estimated for now
         model %in% c("binomial", "betabinom", "weighted"), # Filter only models we want to compare
         sigma %in% c(0.008, 0.016)) %>% # And fix default parameters
  mutate(model = case_when(
    model == "binomial" ~ "Basic",
    model == "betabinom" ~ "Extra-binomial variation",
    model == "weighted" ~ "Weighted"))


rm(df_plt_byiter, df_plt_overall)



point_alpha = 0.4
point_size = 2
width_lims = c(0, 10)
xbreaks = c(10, 100, 1000, 10000, 100000)
xlabels = parse(text = c("10^1", "10^2", "10^3", "10^4", "10^5"))

plt_cov = ggplot(df_plt %>% filter(variable=="Pt"), aes(x=nSamples, color=model, linetype=factor(sigma))) +
  geom_point(aes(y=cov), alpha=point_alpha, shape="x", size=point_size) +
  geom_point(aes(y=meancov)) +
  geom_line(aes(y=meancov)) +  geom_hline(yintercept=0.95, linetype="dashed") +
  facet_wrap(~Simulation, ncol=3, scales="free_y") +
  # ylim(c(0, 1)) +
  coord_cartesian(ylim=c(0.5, 1)) +
  xlab("") + ylab("Coverage of 95%\ncredible intervals") +
  scale_x_log10(breaks=xbreaks, labels = xlabels) +
  scale_color_manual(values=c("#a12d93", "#2d96b3", "#23851c")) +
  custom_theme
plt_cov

plt_width = ggplot(df_plt %>% filter(variable=="Pt"), aes(x=nSamples, color=model, linetype=factor(sigma))) +
  geom_point(aes(y=width), alpha=point_alpha, shape="x", size=point_size) +
  geom_point(aes(y=meanwidth)) +
  geom_line(aes(y=meanwidth)) + 
  facet_wrap(~Simulation, ncol=3) +
  xlab("Daily sample size (log scale)") + ylab("Average width of 95%\ncredible intervals (log scale)") +
  scale_x_log10(breaks=xbreaks, labels = xlabels) +
  scale_y_log10() +
  scale_color_manual(values=c("#a12d93", "#2d96b3", "#23851c")) +
  custom_theme
plt_width

plt_forlegend = plt_cov + theme(legend.position="right") + labs(color="Observation model") + guides(linetype = "none")

plt_nolegend = plot_grid(plt_cov, plt_width, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1,0.22))
plt
ggsave("paper/figures/1-dynamics.png", plt, width=24, height=14, units="cm", dpi=600)



# Fetch off-plot values
df_offplot = df_plt %>%
  filter(variable=="Pt") %>%
  group_by(params, confignum, model, variable, rho, xi, sigma, nSamples, Simulation) %>%
  summarise(meancov_min = min(meancov), meancov_max=max(meancov))
