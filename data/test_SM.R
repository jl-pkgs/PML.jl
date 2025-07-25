pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  hydroTools, kfold
)

l = read_ufile("data/Daily/CRO_禹城_Day_FluxMet_2003-2010.csv")
df = l$data %>% mutate(ET = W2mm(LE, Ta), .after = LE)
dat = df[, .(ET, Rn, Rs, WS, VPD = cal_es(Ta) - ea, SM1=SM_1, SM2=SM_2)]
cor(dat, use = "pairwise.complete.obs")

X = dat[, -1]
Y = dat[, 1]
r = kfold_rf(X, Y)

## 加入SM不见得变好
X = dat[, -1] %>% select(-SM2, -SM1)
Y = dat[, 1]
r = kfold_rf(X, Y)

fs = dir2("Z:/st261/Flux/ChinaFlux/data/Daily", "EBF")

l = read_ufile(fs[1])
df <- l$data #%>% mutate(ET = W2mm(LE, Ta), .after = LE)
dat <- df[, .(ET, Rn, Rs, WS, VPD = cal_es(Ta) - ea, SM1 = SM_1, SM2 = SM_2)]

dat = df[, -(1:5)]
cor(dat, use = "pairwise.complete.obs")[, 4]

## 一种检测方法
dat = df[, .(
  ET, 
  Rn, 
  # Rs,
  Ta = `Ta_1.5m`, 
  WS = `WS_1.5m`, 
  VPD = cal_es(`Ta_1.5m`) - `ea_1.5m`,
  VPD_canopy = cal_es(Ta_33m) - `ea_33m`,
  SM_5cm, 
  SM_20cm, 
  SM_50cm
)] %>% na.omit()
cor(dat, use = "pairwise.complete.obs")

## 加入SM不见得变好
X <- dat[, -1] %>% select(-SM_50cm) ## 不到20cm
Y <- dat[, 1]
r <- kfold_rf(X, Y)
r
