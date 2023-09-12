#%%
# Test that bglib works, comparing to original BASGRA_N example run_BASGRA_Saerheim_2000_09_Grindstad.R

import pandas as pd
import numpy as np
import basgra
bg = basgra.bglib

# %%
# Set params and weather based on BASGRA_N/run_BASGRA_Saerheim_2000_09_Grindstad.R

parcol = 14
df_params = pd.read_table("parameters/parameters.txt")
params = df_params.iloc[:, 13].to_numpy()
p = np.zeros(120) #params dim must be 120
p[:len(params)] = params

bg.set_params(p)
assert bg.parameters_plant.claiv == p[9] # Check

#%%
# Weather
weather = pd.read_table("weather/weather_00_Saerheim_format_bioforsk.txt")
year_start = 2000
doy_start = 112
wm = weather.query(f"YR == {year_start}").query(f"doy > {doy_start-1}")

bge = bg.environment

n = wm.shape[0]
NWEATHER = 8
matrix_weather = np.zeros((bge.nmaxdays, NWEATHER))
matrix_weather[0:n,0] = wm.YR
matrix_weather[0:n,1] = wm.doy
matrix_weather[0:n,2] = wm.GR
matrix_weather[0:n,3] = wm["T"]
matrix_weather[0:n,4] = wm["T"]
matrix_weather[0:n,5] = np.exp(17.27*wm["T"]/(wm["T"]+239)) * 0.6108 * wm.RH / 100
matrix_weather[0:n,6] = wm.RAINI
matrix_weather[0:n,7] = wm.WNI

#%%
MATRIX_WEATHER = matrix_weather

bge.yeari  = MATRIX_WEATHER[:,0]
bge.doyi   = MATRIX_WEATHER[:,1]
bge.gri    = MATRIX_WEATHER[:,2]
bge.tmmni  = MATRIX_WEATHER[:,3]
bge.tmmxi  = MATRIX_WEATHER[:,4]
bge.vpi   = MATRIX_WEATHER[:,5]
bge.raini = MATRIX_WEATHER[:,6]
bge.wni   = MATRIX_WEATHER[:,7]

#%% Management from initialisation/initialise_BASGRA_Saerheim_2000_09_Grindstad.R

calendar_fert = np.zeros((100, 3))
calendar_Ndep = np.zeros((100, 3))
calendar_Ndep[0, :] = [1900, 1, 0]
calendar_Ndep[1, :] = [2100, 366, 0]
days_harvest = np.full((100, 2), -1, dtype=int)

import numpy as np
calendar_fert[0, :] = [2000, 115, 140 * 1000 / 10000]  # 140 kg N ha-1 applied on day 115
calendar_fert[1, :] = [2000, 150, 80 * 1000 / 10000]  # 80 kg N ha-1 applied on day 150
calendar_Ndep[0, :] = [1900, 1, 2 * 1000 / (10000 * 365)]  # 2 kg N ha-1 y-1 N-deposition in 1900
calendar_Ndep[1, :] = [1980, 366, 20 * 1000 / (10000 * 365)]  # 20 kg N ha-1 y-1 N-deposition in 1980
calendar_Ndep[2, :] = [2100, 366, 20 * 1000 / (10000 * 365)]  # 20 kg N ha-1 y-1 N-deposition in 2100
days_harvest[0, :] = [2000, 150]
days_harvest[1, :] = [2000, 216]
days_harvest = days_harvest.astype(int)

#%%
bg.bglib.days_harvest = days_harvest
bg.bglib.days_fert = calendar_fert[:,0:2]
bg.bglib.days_ndep = calendar_Ndep[:,0:2]
bg.bglib.nfertv = calendar_fert[:,2] * bg.parameters_site.nfertmult
bg.bglib.ndepv = calendar_Ndep[:, 2]

#%%
bg.bglib.init()

#%%

import copy

out = []
for i in range(254):
    bg.bglib.step()
    #print(bg.bglib.dm)
    #print(bg.bglib.lai)
    d = dict(dm  = bg.bglib.dm,
         davtmp =bg.environment.davtmp,
         lai = bg.bglib.lai,
         lt50 = bg.bglib.lt50,
         day = bg.bglib.day,
         doy = bg.bglib.doy,
         crt = bg.bglib.crt,
         nrt = bg.bglib.nrt,
         nsh = bg.bglib.nsh,
         nshnor = bg.plant.nshnor,
         )
    out.append(copy.deepcopy(d))

# %%
import matplotlib.pyplot as plt
df = pd.DataFrame(out)
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1,4, layout="tight", figsize=(8,4) )
ax1.plot(df["doy"], df["davtmp"])
ax1.set_title("davtmp")
ax2.plot(df["doy"], df["lai"])
ax2.set_title("lai")
ax3.plot(df["doy"], df["nsh"])
ax3.set_title("nsh")
ax4.plot(df["doy"], df["lt50"])
ax4.set_title("lt50")

#ax3.plot(df["doy"], df["dm"])
#ax3.set_title("dm")

plt.show()
# %%
# Compare to R output
# print(f"LAI = {df.lai.max()}, DAVTMP = {df.davtmp.mean()}, NSH = {df.nsh.max()}")

# %%
