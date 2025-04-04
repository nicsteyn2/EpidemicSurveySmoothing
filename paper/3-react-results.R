
rm(list=ls())

library(tidyverse)
library(Hmisc)
library(cowplot)

# Load raw data
df_data = read.csv("data/reactdata.csv") %>%
  mutate(date=as.Date(date), obsPos = nPos/nSamples)

# Calculate Agresti and Coull confidence intervals
agresticoull = binconf(df_data$nPos, df_data$nSamples, method="exact")
df_data$agrestiLower = agresticoull[,2]
df_data$agrestiUpper = agresticoull[,3]

# Load estimates
df_steyn = read.csv("paper/outputs/3-react/steyn_r1to19_states.csv") %>%
  mutate(Method="Novel", date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method)  %>%
  filter(variable %in% c("rt", "Pt", "nPos"))

df_eales = read.csv("paper/outputs/3-react/eales_r1to19_states.csv") %>%
  mutate(Method="Eales", date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method)  %>%
  filter(variable %in% c("rt", "Pt", "nPos"))

df_abbott = read.csv("paper/outputs/3-react/abbott_r1to19_results.csv") %>%
  mutate(Method="Abbott", date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method)  %>%
  filter(variable %in% c("rt", "Pt", "nPos"))

# Combine
df = rbind(df_steyn, df_eales, df_abbott) %>% 
  left_join(df_data %>% select(date, nSamples, nPos, obsPos, round, agrestiLower, agrestiUpper), by="date") %>%
  mutate(Method = factor(Method, levels=c("Novel", "Eales", "Abbott")))

# Calculate observed swab positivity and join onto main
df_swabpos = df %>%
  filter(variable=="nPos") %>%
  mutate(mean = mean/nSamples, lower=lower/nSamples, upper=upper/nSamples, variable="obsPos")

df = rbind(df, df_swabpos)

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
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Growth rate\n(per day)")

plt_rt = add_rounds_to_plt(plt_rt, df_data)

plt_Pt = ggplot(df %>% filter(variable=="Pt")) +
  geom_errorbar(aes(x=date, ymin=100*agrestiLower, ymax=100*agrestiUpper), color="black", alpha=0.2) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Method), alpha=0.3) +
  geom_line(aes(x=date, y=100*upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=100*mean, color=Method)) +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("Prevalence\n(log scale, %)") + scale_y_log10() +
  coord_cartesian(ylim=c(0.015, 10))

plt_Pt = add_rounds_to_plt(plt_Pt, df_data)

# df_obs = df_data %>% complete(date=seq(min(df_data$date), max(df_data$date), by="day"))

plt_nPos = ggplot(df %>% filter(variable=="obsPos")) +
  geom_ribbon(aes(x=date, ymin=100*lower, ymax=100*upper, fill=Method), alpha=0.3) +
  geom_line(aes(x=date, y=100*mean, color=Method)) +
  geom_point(aes(x=date, y=100*obsPos), color="black", size=0.3, data=.%>%filter(obsPos>0)) +
  # geom_line(aes(x=date, y=100*obsPos), color="black", size=0.2) +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  geom_line(aes(x=date, y=100*upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=100*lower, color=Method, linetype="custom")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("Date") + ylab("Swab positivity\n (log scale, %)") + scale_y_log10() +
  coord_cartesian(ylim=c(0.008, 10))

plt_nPos = add_rounds_to_plt(plt_nPos, df_data)


plt_forlegend = plt_nPos + theme(legend.position="right") + guides(linetype = "none")

plt_nolegend = plot_grid(plt_rt, plt_Pt, plt_nPos, ncol=1, align="v")
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.1))
plt = ggdraw(plt) + draw_line(x=c(0.1, 0.9), y=c(0.355, 0.355), linetype="dashed", color="darkgray")
plt
ggsave("paper/figures/3-react.png", plt, width=24, height=16, units="cm", dpi=600)


