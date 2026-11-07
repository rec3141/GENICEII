"""Modelled daily light at CMO Churchill (58.76 N, 94.07 W; local CST = UTC-6) for both experiments.
Solar position: NOAA / Spencer (1971) equations, 1-minute integration. Per day: day length (h, sun centre above horizon),
top-of-atmosphere insolation (MJ m-2 d-1; S0 1361 W m-2 with eccentricity correction), clear-sky GHI by Haurwitz (1945):
1098 cosZ exp(-0.057/cosZ) W m-2, PAR = 0.45 GHI, 4.57 umol photons per J -> mol photons m-2 d-1.
Clear-sky and above the ice/roof: an upper bound with no cloud, snow, ice, roof or oil attenuation. No measured PAR exists.
Out: light/out/par_daily.tsv"""
import numpy as np, pandas as pd
LAT, LON, TZ = 58.76, -94.07, -6
def day_values(date):
    d = pd.Timestamp(date); n = d.dayofyear
    t = np.arange(0, 1440) / 60.0                                  # local clock hours, 1-min steps
    g = 2 * np.pi / 365 * (n - 1 + (t - 12) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g)
            - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = t * 60 + eqt + 4 * LON - 60 * TZ                          # true solar time, minutes
    ha = np.radians(tst / 4 - 180)
    la = np.radians(LAT)
    cz = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)
    E0 = 1.000110 + 0.034221 * np.cos(g) + 0.001280 * np.sin(g) + 0.000719 * np.cos(2 * g) + 0.000077 * np.sin(2 * g)
    up = cz > 0
    toa = np.where(up, 1361 * E0 * cz, 0)
    czp = np.where(up, cz, 1)
    ghi = np.where(up, 1098 * czp * np.exp(-0.057 / czp), 0)
    par = (0.45 * ghi * 60).sum() * 4.57e-6                         # J m-2 -> mol photons m-2 d-1
    return dict(date=d.date(), day_length_h=up.sum() / 60, toa_MJ=(toa * 60).sum() / 1e6, ghi_MJ=(ghi * 60).sum() / 1e6,
                par_mol=par, noon_elev_deg=float(np.degrees(np.arcsin(cz.max()))))
if __name__ == '__main__':
    days = list(pd.date_range('2024-11-01', '2025-03-15')) + list(pd.date_range('2025-11-01', '2026-03-15'))
    P = pd.DataFrame([day_values(d) for d in days])
    P.to_csv('/data/dev/rnaseq-exploration/light/out/par_daily.tsv', sep='\t', index=False)
    print(P.iloc[[0, 50, 80, 120]].round(2).to_string())
    print('min/max PAR', P.par_mol.min().round(2), P.par_mol.max().round(2), 'min daylength', P.day_length_h.min().round(2))
