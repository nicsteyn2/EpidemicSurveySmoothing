
rm(list=ls())

library(tidyverse)
library(cowplot)

# Load data
df_data = read.csv("paper/outputs/s1-sensitivity/s1_exampledata.csv")

# Load results
df_r0 = read.csv("paper/outputs/s1-sensitivity/s1_sensitivity_pr0_states.csv") %>% mutate(Prior = "Initial growth rate")
df_r0_params = read.csv("paper/outputs/s1-sensitivity/s1_sensitivity_pr0_params.csv") %>% mutate(Prior = "Initial growth rate")

df_P0 = read.csv("paper/outputs/s1-sensitivity/s1_sensitivity_pP0_states.csv") %>% mutate(Prior = "Initial prevalence")
df_P0_params = read.csv("paper/outputs/s1-sensitivity/s1_sensitivity_pP0_params.csv") %>% mutate(Prior = "Initial prevalence")

df_sig = read.csv("paper/outputs/s1-sensitivity/s1_sensitivity_sigma_states.csv") %>% mutate(Prior = "sigma")
df_sig_params = read.csv("paper/outputs/s1-sensitivity/s1_sensitivity_sigma_params.csv") %>% mutate(Prior = "sigma")

# Set for plotting
df_plt = rbind(df_r0, df_P0, df_sig) %>% filter(variable %in% c("rt", "Pt"))

# Turn prevalence into percentage
df_plt$mean[df_plt$variable=="Pt"] = df_plt$mean[df_plt$variable=="Pt"] * 100
df_plt$lower[df_plt$variable=="Pt"] = df_plt$lower[df_plt$variable=="Pt"] * 100
df_plt$upper[df_plt$variable=="Pt"] = df_plt$upper[df_plt$variable=="Pt"] * 100

custom_theme = theme_bw() + theme(legend.position="none")

# Plot states
plt_r0_rt = ggplot(df_plt %>% filter(Prior == "Initial growth rate", variable=="rt"), aes(x=t)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=prior), alpha=0.3) +
  geom_line(aes(y=mean, color=prior)) +
  custom_theme +
  ylab("Growth rate\n(per day)") + xlab("")

plt_P0_rt = ggplot(df_plt %>% filter(Prior == "Initial prevalence", variable=="rt"), aes(x=t)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=prior), alpha=0.3) +
  geom_line(aes(y=mean, color=prior)) +
  custom_theme +
  ylab("") + xlab("")

plt_sig_rt = ggplot(df_plt %>% filter(Prior == "sigma", variable=="rt"), aes(x=t)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=prior), alpha=0.3) +
  geom_line(aes(y=mean, color=prior)) +
  custom_theme +
  ylab("") + xlab("")

plt_r0_Pt = ggplot(df_plt %>% filter(Prior == "Initial growth rate", variable=="Pt"), aes(x=t)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=prior), alpha=0.3) +
  geom_line(aes(y=mean, color=prior)) +
  custom_theme +
  ylab("Prevalence\n(%)") + xlab("")

plt_P0_Pt = ggplot(df_plt %>% filter(Prior == "Initial prevalence", variable=="Pt"), aes(x=t)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=prior), alpha=0.3) +
  geom_line(aes(y=mean, color=prior)) +
  custom_theme +
  ylab("") + xlab("Time (days)")

plt_sig_Pt = ggplot(df_plt %>% filter(Prior == "sigma", variable=="Pt"), aes(x=t)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=prior), alpha=0.3) +
  geom_line(aes(y=mean, color=prior)) +
  custom_theme +
  ylab("") + xlab("")

plt_states = plot_grid(plt_r0_rt, plt_P0_rt, plt_sig_rt, plt_r0_Pt, plt_P0_Pt, plt_sig_Pt, ncol=3, align="vh")



custom_theme_params = theme_bw() + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank())

plt_r0_param = ggplot(df_r0_params) +
  geom_point(aes(x=prior, y=m, color=prior)) +
  geom_errorbar(aes(x=prior, ymin=l, ymax=u, color=prior), width=0.2) +
  custom_theme_params +
  ylim(c(0, 0.03)) + ylab("sigma") + xlab("")

plt_P0_param = ggplot(df_P0_params) +
  geom_point(aes(x=prior, y=m, color=prior)) +
  geom_errorbar(aes(x=prior, ymin=l, ymax=u, color=prior), width=0.2) +
  geom_hline(yintercept=0.01, linetype="dashed", color="gray") +
  custom_theme_params +
  ylim(c(0, 0.03)) + ylab("sigma") + xlab("")

plt_sig_param = ggplot(df_sig_params) +
  geom_point(aes(x=prior, y=m, color=prior)) +
  geom_errorbar(aes(x=prior, ymin=l, ymax=u, color=prior), width=0.2) +
  geom_hline(yintercept=0.01, linetype="dashed", color="gray") +
  custom_theme_params +
  ylim(c(0, 0.03)) + ylab("sigma") + xlab("")


legend_r0 = get_legend(ggplot(df_r0_params) + geom_point(aes(x=prior, y=m, color=prior)) + theme_bw() + theme(legend.title=element_text(size=12)) + labs(color="Prior dist. on initial growth rate"))
legend_P0 = get_legend(ggplot(df_P0_params) + geom_point(aes(x=prior, y=m, color=prior)) + theme_bw() + theme(legend.title=element_text(size=12)) + labs(color="Prior dist. on initial prevalence"))
legend_sig = get_legend(ggplot(df_sig_params) + geom_point(aes(x=prior, y=m, color=prior)) + theme_bw() + theme(legend.title=element_text(size=12)) + labs(color="Prior dist. on sigma"))

plt_params = plot_grid(plt_r0_param, plt_P0_param, plt_sig_param, ncol=3)
plt_legends = plot_grid(legend_r0, legend_P0, legend_sig, ncol=3)

plt = plot_grid(plt_legends, plt_states, plt_params, ncol=1, rel_heights=c(0.25, 1, 0.4))
plt
ggsave("paper/figures/s1_sensitivity.png", plt, width=24, height=16, dpi=600, units="cm")
