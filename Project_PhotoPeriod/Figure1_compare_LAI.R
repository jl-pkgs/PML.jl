pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  ggplot2, gg.layers
)

gof_glass = fread("./OUTPUT/PMLV2China_flux37_gof.csv")
# df = fread("./OUTPUT/PMLV2China_flux37_OUTPUT.csv")
# df_gof[var == "ET", ] %>% arrange(NSE)

# df_gof_whit = fread("./OUTPUT/PMLV2China_flux37_(LAI_whit)_OUTPUT.csv")
gof_whit = fread("./OUTPUT/PMLV2China_flux37_gof.csv")

df = list(GLASS = gof_glass, WHIT = gof_whit) %>% melt_list("LAI")

dat = df[var == "ET"] %>% select(-n_valid) %>% 
  melt(c("LAI", "site", "var")) %>% 
  dcast(site + var + variable ~ LAI, value.var = "value") %>% 
  .[variable %in% c("NSE", "KGE", "R2"), ]

dat[GLASS <= -2]

p <- ggplot(dat[GLASS > -2], aes(GLASS, WHIT)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  facet_wrap(~variable, scales = "free") 

## 
write_fig(p, 'd:/Rplot.pdf', 10, 5)
# df_gof[var == "ET", ] %>% arrange(NSE)
