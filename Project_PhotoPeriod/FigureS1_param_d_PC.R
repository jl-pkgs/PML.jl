pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  ggplot2, gg.layers, ggrepel
)

d = fread("./OUTPUT/PMLV2China_flux37_LAI_glass,WithPC_par.csv")
hist(d$d_pc) # 他应该是一个不敏感性参数，大胆设置为2

# NonPC = fread("./OUTPUT/PMLV2China_flux37_LAI_glass,NonPC_gof.csv")
# df = fread("./OUTPUT/PMLV2China_flux37_OUTPUT.csv")
# df_gof[var == "ET", ] %>% arrange(NSE)
