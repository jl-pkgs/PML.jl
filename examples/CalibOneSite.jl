#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
# # Calibrate PMLV2 model using FLUXNET data
# > US Twt (CRO, 2009)
using PML, Ipaper

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
df_out, df, par = deserialize(file_FLUXNET_CRO_USTwt)
df.GPP_obs = df.GPPobs
df.ET_obs = df.ETobs

# ## 模型参数率定

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
theta, goal, flag = ModelCalib(df, par0)
df_out = PMLV2_sites(df; par=theta2par(theta))
df_out[1:10, :]

# ## 拟合优度

#nb %% A slide [code] {"slideshow": {"slide_type": "fragment"}}
gof = [
  (; var="ET", GOF(df.ET_obs, df_out.ET)...),
  (; var="GPP", GOF(df.GPP_obs, df_out.GPP)...)] |> DataFrame
DataFrame(gof)


# ## 绘图

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
using Plots
gr(framestyle=:box, titlefontsize=12)
t = df.date
inds = 1:46*1

function plot_var(var; label=string(var), title=string(var),
  data=df_out, scale=1, kw...)
  plot(t[inds], data[inds, var] * scale; label, title, kw...)
end
function plot_var!(p, var; label=string(var),
  data=df_out, kw...)
  plot!(p, t[inds], data[inds, var]; label, kw...)
end

p_et = plot(title="ET components (mm/d)")
plot_var!(p_et, :Ec)
plot_var!(p_et, :Es)
plot_var!(p_et, :Ei)
plot_var!(p_et, :ETobs; data=df, label="ET_obs", color=:black)

p_gpp = plot_var(:GPP; title="GPP (gC m-2 d-1)", label="GPP")
plot_var!(p_gpp, :GPPobs; data=df, label="GPP_obs", color=:black)

plot(
  p_et, p_gpp,
  plot_var(:Eeq; title="Eeq (mm/d)", label="Eeq"),
  plot_var(:Gc_w; title="Conductance (m s-1)", label="Gc"),
  plot_var(:fval_soil; title="β_soil", label="β_soil"),
  plot_var(:VPD; data=df, scale=-1, title="-VPD (kPa)", label="-VPD"),
  plot_var(:Pi; title="P - Ei (mm/d)"),
  plot_var(:Es_eq; title="Es_eq (mm/d)"),
  size=(800, 700),
  layout=(4, 2), 
)
