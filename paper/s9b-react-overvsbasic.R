
rm(list=ls())

library(tidyverse)
library(cowplot)

# Load raw data
df_data = read.csv("data/reactdata.csv") %>% mutate(date=as.Date(date)) %>% mutate(obsPos = nPos/nSamples)

# Process start and end dates of grouped rounds
df_rounddates = df_data %>%
  mutate(rounds = case_when(
    round <= 7 ~ "1-to-7",
    round <= 13 ~ "8-to-13",
    round <= 16 ~ "14-to-16",
    round <= 19 ~ "17-to-19")) %>%
  group_by(rounds) %>%
  summarise(st_date = min(date), en_date = max(date))


# Load estimates
df_overdisp = rbind(
  read.csv("paper/outputs/3-react/steyn_r1to7_states.csv") %>% mutate(Method="Extra-binomial", rounds="1-to-7"),
  read.csv("paper/outputs/3-react/steyn_r8to13_states.csv") %>% mutate(Method="Extra-binomial", rounds="8-to-13"),
  read.csv("paper/outputs/3-react/steyn_r14to16_states.csv") %>% mutate(Method="Extra-binomial", rounds="14-to-16"),
  read.csv("paper/outputs/3-react/steyn_r17to19_states.csv") %>% mutate(Method="Extra-binomial", rounds="17-to-19")
)

df_basic = rbind(
  read.csv("paper/outputs/s9b-react-basic/steyn_r1to7_states.csv") %>% mutate(Method="Binomial", rounds="1-to-7"),
  read.csv("paper/outputs/s9b-react-basic/steyn_r8to13_states.csv") %>% mutate(Method="Binomial", rounds="8-to-13"),
  read.csv("paper/outputs/s9b-react-basic/steyn_r14to16_states.csv") %>% mutate(Method="Binomial", rounds="14-to-16"),
  read.csv("paper/outputs/s9b-react-basic/steyn_r17to19_states.csv") %>% mutate(Method="Binomial", rounds="17-to-19")
)

df = rbind(df_overdisp, df_basic) %>%
  mutate(date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method, rounds)%>%
  filter(variable %in% c("rt", "Pt", "nPos")) %>%
  group_by(variable, Method) %>%
  complete(date=seq(min(date), max(date), by="day"))%>% 
  left_join(df_data %>% select(date, nSamples, nPos, obsPos, round), by="date")

# Calculate observed swab positivity and join onto main
df_swabpos = df %>%
  filter(variable=="nPos") %>%
  mutate(mean = mean/nSamples, lower=lower/nSamples, upper=upper/nSamples, variable="obsPos")

df = rbind(df, df_swabpos)

# Plot rt
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

plt_rt = ggplot(df %>% filter(variable=="rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Method), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Method)) +
  scale_color_manual(values=c("#a12d93", "#2d96b3")) +
  scale_fill_manual(values=c("#a12d93", "#2d96b3")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Growth rate\n(per day)")

plt_rt = add_rounds_to_plt(plt_rt, df_data)

plt_Pt = ggplot(df %>% filter(variable=="Pt")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Method), alpha=0.3) +
  geom_line(aes(x=date, y=100*upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Method)) +
  scale_color_manual(values=c("#a12d93", "#2d96b3")) +
  scale_fill_manual(values=c("#a12d93", "#2d96b3")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Prevalence\n(log scale, %)") + scale_y_log10() +
  coord_cartesian(ylim=c(0.015, 10))

plt_Pt = add_rounds_to_plt(plt_Pt, df_data)

# df_obs = df_data %>% complete(date=seq(min(df_data$date), max(df_data$date), by="day"))

plt_nPos = ggplot(df %>% filter(variable=="obsPos")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Method), alpha=0.4) +
  geom_line(aes(x=date, y=100*mean, color=Method)) +
  geom_point(aes(x=date, y=100*obsPos), color="black", size=0.3) +
  geom_line(aes(x=date, y=100*obsPos), color="black", size=0.2) +
  scale_color_manual(values=c("#a12d93", "#2d96b3")) +
  scale_fill_manual(values=c("#a12d93", "#2d96b3")) +
  # scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Swab positivity\n(log scale, %)") + scale_y_log10() +
  coord_cartesian(ylim=c(0.015, 10))

plt_nPos = add_rounds_to_plt(plt_nPos, df_data)


plt_forlegend = plt_Pt + theme(legend.position="right") + labs(color="Observation dist.", fill="Observation dist.") + guides(linetype = "none")

plt_nolegend = plot_grid(plt_rt, plt_Pt, plt_nPos, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.2))
plt
ggsave("paper/figures/s9b-react-overvsbasic.png", plt, width=24, height=16, units="cm", dpi=600)
 
