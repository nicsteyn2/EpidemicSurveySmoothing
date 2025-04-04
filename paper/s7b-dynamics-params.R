
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

df_params_in = foreach(ii = seq(1, length(confignums)), .combine=rbind) %do% {
  read.csv(sprintf("paper/outputs/1-dynamics/1-dynamics-config%d-params.csv", confignums[ii])) %>% mutate(confignum=confignums[ii])
}

df_params = df_params_in %>%
  left_join(config, by="confignum") %>%
  group_by(param, model, confignum, nSamples, rho, xi, sigma, Simulation) %>%
  summarise(m=mean(m), l=mean(l), u=mean(u)) %>%
  mutate(model = case_when(
    model == "binomial" ~ "Basic",
    model == "betabinom" ~ "Extra-binomial variation",
    model == "weighted" ~ "Weighted")) %>%
  mutate(nSamples = case_when(
    model=="Basic" ~ 0.8*nSamples,
    model=="Extra-binomial variation" ~ nSamples,
    model=="Weighted" ~ 1.2*nSamples))

rm(df_params_in, ii, confignums)

point_alpha = 0.4
point_size = 2
width_lims = c(0, 10)
xbreaks = c(10, 100, 1000, 10000, 100000)
xlabels = parse(text = c("10^1", "10^2", "10^3", "10^4", "10^5"))

chosen_sigma = 0.016

plt_sigma = ggplot(df_params %>% filter(param=="sigma", sigma==chosen_sigma), aes(x=nSamples)) +
  geom_ribbon(aes(ymin=l, ymax=u, fill=model), alpha=0.2) +
  geom_errorbar(aes(ymin=l, ymax=u, color=model), width=0.2) +
  geom_point(aes(y=m, color=model), shape="x", size=3) +
  geom_hline(yintercept=chosen_sigma, linetype="dashed") +
  facet_wrap(~Simulation, ncol=3) +
  scale_x_log10(breaks=xbreaks, labels = xlabels) + scale_y_log10() +
  scale_color_manual(values=c("#a12d93", "#2d96b3", "#23851c")) +
  scale_fill_manual(values=c("#a12d93", "#2d96b3", "#23851c")) +
  xlab("") + ylab("Estimated sigma\n(log scale)") +
  custom_theme

plt_rho = ggplot(df_params %>% filter(param=="rho", sigma==chosen_sigma), aes(x=nSamples)) +
  geom_ribbon(aes(ymin=l, ymax=u), alpha=0.2, fill="#2d96b3") +
  geom_errorbar(aes(ymin=l, ymax=u), width=0.2, color="#2d96b3") +
  geom_point(aes(y=m), shape="x", size=3, color="#2d96b3") +
  geom_hline(aes(yintercept=yint), linetype="dashed", data=data.frame(yint=2e-4, Simulation="Extra-binomial simulation")) +
  facet_wrap(~Simulation, ncol=3) +
  scale_x_log10(breaks=xbreaks, labels = xlabels) + scale_y_log10() +
  xlab("") + ylab("Estimated rho\n(log scale)") +
  custom_theme

plt_c = ggplot(df_params %>% filter(param=="c", sigma==chosen_sigma), aes(x=nSamples)) +
  geom_ribbon(aes(ymin=l, ymax=u), alpha=0.2, fill="#23851c") +
  geom_errorbar(aes(ymin=l, ymax=u), width=0.2, color="#23851c") +
  geom_point(aes(y=m), shape="x", size=3, color="#23851c") +
  geom_hline(aes(yintercept=yint), linetype="dashed", data=data.frame(yint=1, Simulation="Basic simulation")) +
  facet_wrap(~Simulation, ncol=3) +
  scale_x_log10(breaks=xbreaks, labels = xlabels) + scale_y_log10() +
  xlab("Daily sample size (log scale)") + ylab("Estimated c\n(log scale)") +
  custom_theme

legend_ = get_legend(plt_sigma + theme(legend.position="right") + labs(color="Model", fill="Model"))

plt_nolegend = plot_grid(plt_sigma, plt_rho, plt_c, ncol=1, align="vh")
plt = plot_grid(plt_nolegend, legend_, ncol=2, rel_widths=c(1, 0.22))
plt
ggsave("paper/figures/s7b-dynamics-params.png", plt, width=24, height=18, units="cm", dpi=600)

  
  
  