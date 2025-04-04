
rm(list=ls())

# Load data
abbott_raw = read.csv("data/PCR-curves/hellewell.csv")
binny_raw = read.csv("data/PCR-curves/binny.csv")
kucirka_raw = read.csv("data/PCR-curves/kucirka.csv")

# Tidy data
abbott = abbott_raw %>%
  pivot_longer(cols=-sample, values_to="y", names_to="x") %>%
  mutate(t = as.numeric(str_extract(x, "\\d+"))) %>%
  select(-x) %>%
  group_by(t) %>%
  summarise(median=median(y), lower=quantile(y, 0.025), upper=quantile(y, 0.975)) %>%
  mutate(Source="Hellewell et al. (2021)")

binny = binny_raw %>%
  rename(t=days_since_infection, lower=lower_95, upper=upper_95) %>%
  mutate(Source="Binny et al. (2023)")

kucirka <- kucirka_raw %>%
  mutate(median = str_extract(false_neg_rate, "^[0-9.]+") %>% as.numeric(),
         lower = str_extract(false_neg_rate, "(?<=\\()[0-9.]+") %>% as.numeric(),
         upper = str_extract(false_neg_rate, "[0-9.]+(?=%\\))") %>% as.numeric()) %>%
  rename(t=time) %>%
  select(-false_neg_rate) %>%
  mutate(median = 1 - median/100,
         lower = 1 - lower/100,
         upper = 1 - upper/100,
         Source="Kucirka et al. (2020)")

# Combine
df = rbind(abbott, binny, kucirka)

# Plot
plt = ggplot(df, aes(x=t)) +
  geom_ribbon(aes(ymin=100*lower, ymax=100*upper, fill=Source), alpha=0.3) +
  geom_line(aes(y=100*median, color=Source)) +
  theme_bw() + xlab("Time from infection (days)") + ylab("RT-PCR sensitivity (%)") + ylim(c(0, 100))
plt
ggsave("paper/figures/s5-pcrcurves.png", plt, dpi=600, units="cm", width=20, height=10)

