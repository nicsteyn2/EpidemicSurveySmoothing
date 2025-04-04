
rm(list=ls())

source("OtherMethods/Abbott/AbbottMethod.R")

# Set options
start_round = 1
end_round = 19
iter_warmup = 200
iter_sampling = 300

# Fetch start and end dates
dates = read.csv("data/reactdata.csv") %>% mutate(date=as.Date(date)) %>% summarise(st = min(date), en=max(date))

# An Abbott method-specific windin and windout period needs to be added to the data
windin = 50
windout = 21

# Load data and append Abbott method-specific windin and windout periods
df = read.csv("data/reactdata.csv") %>%
  select(date, round, nSamples, nPos) %>%
  filter(round >= start_round & round <= end_round) %>%
  mutate(date=as.Date(date)) %>%
  complete(date=seq(min(date)-windin, max(date)+windout, by="day"), fill=list(round=NA, nSamples=0, nPos=0)) %>%
  mutate(t = as.integer(date - min(date))) %>%
  arrange(t)

# Fit model using Binny et al curves
results_binny = runAbbottMethod(df, inc2prevmethod="Binny")

# Load results from normal model
results_hellewell = read.csv("paper/outputs/3-react/abbott_r1to19_results.csv")

# Combine and tidy
results = rbind(results_binny %>% mutate(Curves = "Binny et al. (2023)"), results_hellewell %>% mutate(Curves = "Hellewell et al. (2021) (default)")) %>%
  filter(variable %in% c("ell", "alpha", "rho", "It", "Pt", "rtinc", "rt", "nPos")) %>%
  mutate(date = date - 1) %>%
  filter(date >= dates$st[1], date <= dates$en[1])

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


# Plot

plt_It = ggplot(results %>% filter(variable=="It")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Curves), alpha=0.2) +
  geom_line(aes(x=date, y=100*upper, color=Curves, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Curves, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Curves)) +
  custom_theme +
  xlab("") + ylab("Infection incidence\n(log scale, %)") +
  scale_y_log10() +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  scale_linetype_manual(values = c("custom" = "11"))
plt_It = add_rounds_to_plt(plt_It, df)

plt_rt = ggplot(results %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Curves), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Curves, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Curves, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Curves)) +
  custom_theme +
  xlab("") + ylab("Growth rate\n(per day)") +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  scale_linetype_manual(values = c("custom" = "11"))
plt_rt = add_rounds_to_plt(plt_rt, df)

plt_Pt = ggplot(results %>% filter(variable=="Pt")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Curves), alpha=0.2) +
  geom_line(aes(x=date, y=100*upper, color=Curves, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Curves, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Curves)) +
  custom_theme +
  xlab("") + ylab("Prevalence\n(log scale, %)") +
  scale_y_log10() +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  scale_linetype_manual(values = c("custom" = "11"))
plt_Pt = add_rounds_to_plt(plt_Pt, df)

plt_nolegend = plot_grid(plt_It, plt_rt, plt_Pt, ncol=1, align="v")

plt_forlegend = plt_Pt = ggplot(results %>% filter(variable=="Pt")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Curves), alpha=0.2) +
  geom_line(aes(x=date, y=100*mean, color=Curves)) +
  labs(color="PCR Positivity Curve", fill="PCR Positivity Curve") +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  scale_linetype_manual(values = c("custom" = "11"))

plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.3))
plt
ggsave("paper/figures/s5-pcrcurves-abbott-react.png", plt, width=24, height=12.5, units="cm", dpi=600)

