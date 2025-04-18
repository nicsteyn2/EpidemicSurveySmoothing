
rm(list=ls())

library(tidyverse)
library(Hmisc)
library(cowplot)

# Load raw data
df_data = read.csv("data/reactdata.csv") %>%
  mutate(date=as.Date(date), obsPos = nPos/nSamples)

# Load estimates
df_steyn = read.csv("paper/outputs/3-react/steynRt_r1to19_states.csv") %>%
  mutate(Method="SIMPLE", date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method)  %>%
  filter(variable %in% c("Rt","It", "Pt", "nPos"))

df_eales = read.csv("paper/outputs/3-react/eales_r1to19_states.csv") %>%
  mutate(Method="Eales", date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method)  %>%
  filter(variable %in% c("Rt","rt", "Pt", "nPos"))

df_abbott = read.csv("paper/outputs/3-react/abbott_r1to19_results.csv") %>%
  mutate(Method="Abbott", date=as.Date(date)) %>%
  select(variable, date, mean, lower, upper, Method)  %>%
  filter(variable %in% c("Rt","rt", "Pt", "nPos"))

# Create shifted Eales
df_eales_shifted = df_eales %>%
  mutate(date = as.Date(date - 7), Method = "Eales (shifted)")

# Combine
df = rbind(df_steyn, df_eales, df_eales_shifted, df_abbott) %>% 
  left_join(df_data %>% select(date, nSamples, nPos, obsPos, round), by="date") %>%
  mutate(Method = factor(Method, levels=c("SIMPLE", "Eales", "Eales (shifted)", "Abbott")))

# Calculate observed swab positivity and join onto main
df_swabpos = df %>%
  filter(variable=="nPos") %>%
  mutate(mean = mean/nSamples, lower=lower/nSamples, upper=upper/nSamples, variable="obsPos")

df = rbind(df, df_swabpos) #%>% filter()


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


plt_Rta = ggplot(df %>% filter(variable=="Rt", Method!="Eales (shifted)")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Method), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Method)) +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  custom_theme +
  xlab("") + ylab("") +
  geom_rect(aes(xmin=as.Date("2021-09-04"), xmax=as.Date("2022-04-06"), ymin=0.4, ymax=2.2), fill=NA, color="black", linewidth=0.2)

plt_Rtb = ggplot(df %>% filter(variable=="Rt", date>=as.Date("2021-09-09"), Method!="Eales (shifted)")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Method), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Method)) +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#f7a072", "#6aa84f")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  scale_x_date(date_labels = "%b-%Y") +
  custom_theme +
  xlab("") + ylab("Reproduction number")

plt_Rtc = ggplot(df %>% filter(variable=="Rt", date>=as.Date("2021-09-09"), Method!="Eales")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Method), alpha=0.2) +
  geom_line(aes(x=date, y=upper, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=lower, color=Method, linetype="custom")) +
  geom_line(aes(x=date, y=mean, color=Method)) +
  scale_color_manual(values=c("#2d96b3", "#a23e1f", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#a23e1f", "#6aa84f")) +
  scale_linetype_manual(values = c("custom" = "11")) +
  scale_x_date(date_labels = "%b-%Y") +
  custom_theme +
  xlab("Date") + ylab("")


plt_forlegend = ggplot(df %>% filter(variable=="Rt")) +
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, fill=Method), alpha=0.2) +
  geom_line(aes(x=date, y=mean, color=Method)) +
  scale_color_manual(values=c("#2d96b3", "#f7a072", "#a23e1f", "#6aa84f")) +
  scale_fill_manual(values=c("#2d96b3", "#f7a072", "#a23e1f", "#6aa84f"))

plt_Rta = add_rounds_to_plt(plt_Rta, df_data)
plt_Rtb = add_rounds_to_plt(plt_Rtb, df_data)
plt_Rtc = add_rounds_to_plt(plt_Rtc, df_data)

plt_nolegend = plot_grid(plt_Rta, plt_Rtb, plt_Rtc, ncol=1)
plt = plot_grid(plt_nolegend, get_legend(plt_forlegend), ncol=2, rel_widths=c(1, 0.15))
plt
ggsave("paper/figures/3-react-Rt.png", plt, width=24, height=16, units="cm", dpi=600)

