#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
# # Calibrate PMLV2 model using FLUXNET data
# > US Twt (CRO, 2009)
using PML, Ipaper

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
df_out, df, par = deserialize(file_FLUXNET_CRO)
df.GPP_obs = df.GPPobs
df.ET_obs = df.ETobs

# ## 模型参数率定
