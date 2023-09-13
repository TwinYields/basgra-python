#%%
import pandas as pd
import numpy as np
import basgra
import re
import copy

bg = basgra.bglib
bglib = bg.bglib
bge = bg.environment
plant = bg.plant

#%%

def read_met(weatherfile):
    # Header can be of varying length
    # find header line
    with open(weatherfile) as wf:
        lines = wf.readlines()
        row = 1
        for line in lines:
            if line.strip().startswith("year"):
                break
            row += 1

    hdr = re.split("\s+", line.strip())
    wdf = pd.read_fwf(weatherfile, skiprows=row+1,
                names=hdr)
    return wdf

wdf = read_met("weather/maaninka_era5.met")

# %%

# Pad weather arrays
def wpad(x, length = bge.nmaxdays):
    npad = length - len(x)
    return np.pad(x, (0,npad))

class BasgraN(object):

    def __init__(self):
        self.start_date = None
        self.end_date = None
        self.weatherfile = None
        self.set_params()
        self.set_ndep()
        self.Ndays = 0

    def reset_management(self):
        days_harvest = np.full((100, 2), -1, dtype=int)
        bg.bglib.days_harvest = days_harvest.astype(int)
        calendar_fert = np.zeros((100, 3))
        bglib.days_fert = calendar_fert[:,0:2]
        bglib.nfertv = calendar_fert[:,2] * bg.parameters_site.nfertmult

    def set_params(self, parcol = 13):
        #TODO add params file to package
        df_params = pd.read_table("parameters/parameters.txt")
        params = df_params.iloc[:, parcol].to_numpy()
        p = np.zeros(120) #params dim must be 120
        p[:len(params)] = params
        bg.set_params(p)

    def set_ndep(self):
        # Use default N deposition
        calendar_Ndep = np.zeros((100, 3))
        calendar_Ndep[0, :] = [1900, 1, 0]
        calendar_Ndep[1, :] = [2100, 366, 0]
        calendar_Ndep[0, :] = [1900, 1, 2 * 1000 / (10000 * 365)]  # 2 kg N ha-1 y-1 N-deposition in 1900
        calendar_Ndep[1, :] = [1980, 366, 20 * 1000 / (10000 * 365)]  # 20 kg N ha-1 y-1 N-deposition in 1980
        calendar_Ndep[2, :] = [2100, 366, 20 * 1000 / (10000 * 365)]  # 20 kg N ha-1 y-1 N-deposition in 2100
        bglib.days_ndep = calendar_Ndep[:,0:2]
        bglib.ndepv = calendar_Ndep[:, 2]


    def set_weather(self, weatherfile):
        wdf = read_met(weatherfile)
        wdf["T"] = (wdf["maxt"] + wdf["mint"]) / 2.0

        wdf = wdf.query("year >= 2005").query("year <= 2007")
        # TODO cut used weather based on start and end date

        # Zero pad to nmaxdays for Fortran
        self.Ndays = wdf.shape[0]
        npad = bge.nmaxdays - wdf.shape[0]
        padded = np.pad(wdf.to_numpy(), [(0, npad), (0, 0)])
        wdf = pd.DataFrame(padded, columns = wdf.columns)


        bge.yeari = wdf["year"]
        bge.doyi = wdf["day"]
        bge.gri = wdf["radn"]
        bge.tmmni = wdf["mint"]
        bge.tmmxi = wdf["maxt"]
        bge.vpi = np.exp(17.27*wdf["T"]/(wdf["T"]+239)) * 0.6108 * wdf.rh / 100
        bge.raini = wdf["rain"]
        bge.wni = wdf["wind_speed"]

    # TODO set simulation dates, now uses full weatherfile
    def set_dates(self, startdate, enddate):
        pass

    def set_harvest_dates(self, hdates):
        dates = pd.DatetimeIndex(hdates).sort_values()
        days_harvest = np.full((100, 2), -1, dtype=int)
        idx = 0
        for d in dates:
            days_harvest[idx, :] = [d.year, d.day_of_year]
            idx += 1
        bglib.days_harvest = days_harvest.astype(int)

    def set_N(self, Ndata):
        ndf = pd.DataFrame(Ndata, columns=["date", "N"])
        ndf["date"] = pd.DatetimeIndex(ndf["date"])
        ndf = ndf.sort_values("date")
        calendar_fert = np.zeros((100, 3))
        idx = 0
        for d in ndf["date"]:
            calendar_fert[idx, :] = [d.year, d.day_of_year, ndf["N"].iloc[idx]]
            idx += 1
        calendar_fert
        bglib.days_fert = calendar_fert[:,0:2]
        bglib.nfertv = calendar_fert[:,2] * bg.parameters_site.nfertmult

    def init(self):
        bglib.init()

    def step(self):
        bglib.step()

    def run(self):
        out = []
        self.init()

        for i in range(self.Ndays):
            bglib.step()
            d = dict(dm  = bg.bglib.dm,
                    davtmp = bg.environment.davtmp,
                    lai = bglib.lai,
                    lt50 = bglib.lt50,
                    day = bglib.day,
                    dayl = bge.dayl,
                    doy = bglib.doy,
                    crt = bglib.crt,
                    nrt = bglib.nrt,
                    nsh = bglib.nsh,
                    nshnor = plant.nshnor,
                    yield_last = bglib.yield_last
                    )
            out.append(copy.deepcopy(d))

        self.results = pd.DataFrame(out)

model = BasgraN()
model.set_weather("weather/maaninka_era5.met")
model.set_params()
hdates = [f"{y}-06-15" for y in range(2000, 2022)] + [f"{y}-07-30" for y in range(2000, 2022)]
Ndata = [[f"{y}-05-01", 14.0] for y in range(2000, 2023)] +[
         [f"{y}-06-20", 14.0] for y in range(2000, 2023)]

model.reset_management()
#TODO look at harvest code to prevent harvesting 100%
model.set_harvest_dates(hdates)
model.set_N(Ndata)
model.run()

#%%

import matplotlib.pyplot as plt
df = model.results.copy()
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1, layout="tight", figsize=(8,6) )
ax1.plot(df["day"], df["davtmp"])
ax1.set_title("davtmp")
ax2.plot(df["day"], df["lai"])
ax2.set_title("lai")
ax3.plot(df["day"], df["yield_last"])
ax3.set_title("yield_last")
ax4.plot(df["day"], df["dm"])
ax4.set_title("dm")
# %%

#%%
# Maaninka = 11
for p in range(11, 22):
    model.set_params(p)
    model.run()
    df = model.results.copy()
    fig, (ax1, ax2, ax3) = plt.subplots(3,1, layout="tight", figsize=(8,6) )
    ax1.plot(df["day"], df["davtmp"])
    ax1.set_title("davtmp")
    ax2.plot(df["day"], df["lai"])
    ax2.set_title("lai")
    ax3.plot(df["day"], df["dm"])
    ax3.set_title("dm")



# %%
