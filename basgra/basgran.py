import pandas as pd
import numpy as np
import basgra
import re
import copy
from pathlib import Path
import os

bg = basgra.bglib
bglib = bg.bglib
bge = bg.environment
plant = bg.plant

class BasgraN(object):

    def __init__(self):
        self.startdate = None
        self.enddate = None
        self.weatherfile = None
        self.weather = None
        self.site = None
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
        df_params = pd.read_table(os.path.dirname(__file__) + "/data/parameters.txt")
        params = df_params.iloc[:, parcol].to_numpy()
        site = df_params.columns[parcol]
        self.site = site.replace("/","_")
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
        wdf = pd.read_csv(weatherfile)
        wdf["T"] = wdf["avet"]

        if self.startdate is not None:
            st = pd.Timestamp(self.startdate)
            et = pd.Timestamp(self.enddate)
            wdf = wdf.query(f"year >= {st.year} and day >= {st.day_of_year}"
                    ).query(f"year <= {et.year} and day <= {et.day_of_year}")

        self.weatherfile = weatherfile
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
        self.startdate = startdate
        self.enddate = enddate
        if self.weatherfile is not None:
            self.set_weather(self.weatherfile)

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
        #calendar_fert
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
            d = dict(dm  = float(bg.bglib.dm),
                    davtmp = bg.environment.davtmp,
                    lai = float(bglib.lai),
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
            bglib.step()
        df = pd.DataFrame(out)
        df.insert(0, "date", self.doy_to_date(df["year"], df["doy"]))
        self.results = df