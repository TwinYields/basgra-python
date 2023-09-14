#%%
import pandas as pd
import numpy as np
import basgra
import re
import copy
from pathlib import Path

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
        #wdf = read_met(weatherfile)
        wdf = pd.read_csv(weatherfile)
        #wdf["T"] = (wdf["maxt"] + wdf["mint"]) / 2.0
        #wdf["T"] = wdf["maxt"]
        wdf["T"] = wdf["avet"]
        wdf = wdf.query("year >= 2015").query("year <= 2017")
        wdf = wdf.iloc[90:, ]
        # TODO cut used weather based on start and end date
        self.weather = wdf.copy()

        # Zero pad to nmaxdays for Fortran
        self.Ndays = wdf.shape[0]
        npad = bge.nmaxdays - wdf.shape[0]
        padded = np.pad(wdf.to_numpy(), [(0, npad), (0, 0)])
        wdf = pd.DataFrame(padded, columns = wdf.columns)

        bge.yeari = wdf["year"]
        bge.doyi = wdf["day"]
        bge.gri = wdf["radn"]
        bge.tmmni = wdf["T"]
        bge.tmmxi = wdf["T"]
        #bge.tmmni = wdf["mint"]
        #bge.tmmxi = wdf["maxt"]

        bge.vpi = np.exp(17.27*wdf["T"]/(wdf["T"]+239)) * 0.6108 * wdf.rh / 100
        bge.raini = wdf["rain"]
        bge.wni = wdf["windspeed"]

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
        ndf["N"] = ndf["N"]/10. # Give input in kg/ha, scale for model
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

    def doy_to_date(self, years, doys):
        import datetime
        years = pd.DatetimeIndex([datetime.date(y, 1,1) for y in years])
        days = pd.to_timedelta(doys -1, unit="days")
        return years + days

    def run(self):
        out = []
        self.init()

        for i in range(self.Ndays):
            bglib.step()
            d = dict(dm  = bg.bglib.dm,
                    davtmp = bg.environment.davtmp,
                    lai = bglib.lai,
                    lt50 = bglib.lt50,
                    year = int(bglib.year),
                    day = int(bglib.day),
                    dayl = bge.dayl,
                    doy = int(bglib.doy),
                    crt = bglib.crt,
                    nrt = bglib.nrt,
                    nsh = bglib.nsh,
                    nshnor = plant.nshnor,
                    yield_last = bglib.yield_last
                    )
            out.append(copy.deepcopy(d))
        df = pd.DataFrame(out)
        df.insert(0, "date", self.doy_to_date(df["year"], df["doy"]))
        self.results = df

model = BasgraN()
model.set_weather(Path.home() / "git/apsim-notebooks/grassmodels/data/maaninka_era5_basgra.csv")
#model.set_weather("weather/maaninka_fmi.csv")
model.set_params(18)
#hdates = [f"{y}-06-15" for y in range(2000, 2022)] + [f"{y}-07-30" for y in range(2000, 2022)]
#hdates = ["2005-06-30", "2006-06-30"]
#Ndata = [[f"{y}-05-01", 14.0] for y in range(2000, 2023)] +[
#         [f"{y}-06-20", 14.0] for y in range(2000, 2023)]

#%%
hdates = pd.read_csv(Path.home() / "git/apsim-notebooks/grassmodels/data/maaninka_harvestdates_grindstad.csv")
Ndates = pd.read_csv(Path.home() /"git/apsim-notebooks/grassmodels/data/maaninka_fertdates_grindstad.csv")
hdates = hdates["date"]
tr6 = Ndates.query("treatment==6")[["date", "N"]]

#%%
model.reset_management()
#TODO look at harvest code to prevent harvesting 100%
model.set_harvest_dates(hdates)
model.set_N(tr6)
#bg.parameters_plant.claiv = 2.0
model.run()


#%%
cutdata = pd.read_csv(Path.home() / "git/apsim-notebooks/grassmodels/data/maaninka_cuts_grindstad.csv")
cutdata["date"] = pd.DatetimeIndex(cutdata["date"])
cut6 = cutdata.query("Nrate == 150")

#%%
import matplotlib.pyplot as plt
df = model.results.copy()
fig, (ax1, ax2, ax3) = plt.subplots(3,1, layout="tight",
                                         sharex=True, figsize=(8,6) )
ax1.plot(df["date"], df["davtmp"])
ax1.set_title("davtmp")
ax1.set_ylabel("Average T")
ax2.plot(df["date"], df["lai"])
ax2.set_ylabel("LAI")
ax2.scatter(cut6["date"], cut6["LAI"], color="red")
ax3.plot(df["date"], df["dm"]*10)
ax3.scatter(cut6["date"], cut6["dm_kg"], color="red")
ax3.set_ylabel("DM (kg)")

#plt.show()
# %%

#%%
# Maaninka = 11
df_params = pd.read_table("parameters/parameters.txt")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates

with PdfPages("sites_avet.pdf") as pdf:
    for p in range(2, 22):
        site = df_params.columns[p]
        sclean = site.replace("/","_")
        print(site, p)
        model.set_params(p)
        model.run()
        df = model.results.copy()

        fig, (ax1, ax2, ax3) = plt.subplots(3,1, layout="tight", figsize=(8,6), sharex = True )
        ax1.plot(df["date"], df["davtmp"])
        ax1.set_title(f"{site}, no = {p}")
        ax1.set_ylabel("Average T")
        ax2.plot(df["date"], df["lai"])
        ax2.scatter(cut6["date"], cut6["LAI"], color="red")
        ax2.set_ylabel("LAI")
        ax3.plot(df["date"], df["dm"]*10)
        ax3.scatter(cut6["date"], cut6["dm_kg"], color="red")
        ax3.set_title("dm")
        ax3.set_ylabel("DM (kg)")

        #ax2.set_title("lai")
        #ax3.plot(df["date"], df["dm"])
        #ax3.set_title("dm")
        #ax3.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
        #ax3.xaxis.set_tick_params(rotation=40)
        pdf.savefig()
        plt.close()

# %%
# %%
