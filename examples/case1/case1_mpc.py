import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logging
import os
import time

from statsmodels.tsa.stattools import acf
from scipy.stats import norm

import mshoot
from fcastgen import error


# Set up logging
logging.basicConfig(filename='mpc_case1.log', filemode='w', level='DEBUG')

# Random seed
# np.random.seed(12345)

# Paths
ms_file = os.path.join('examples', 'case1', 'measurements.csv')
fmu_dir = os.path.join('examples', 'case1', 'models')

# FMU list
fmus = os.listdir(fmu_dir)

# Simulation period
t0 = '2018-04-05 00:00:00'
t1 = '2018-04-08 00:00:00'

# Read measurements
ms = pd.read_csv(ms_file)
ms['datetime'] = pd.to_datetime(ms['datetime'])
ms = ms.set_index('datetime')
ms = ms.loc[t0:t1]

# Resample
ms = ms.resample('1h').mean().ffill().bfill()

# Assign model inputs
inp = ms[['solrad', 'Tout', 'occ', 'dpos', 'vpos']]
inp['time'] = (inp.index - inp.index[0]).total_seconds()
inp = inp.set_index('time')

# Modify inputs
inp['dpos'] = 0

# Initial state
airT0 = 20. + 273.15
x0 = [airT0]

# Cost function
def cfun(xdf, ydf):
    """
    :param ydf: DataFrame, model states
    :param ydf: DataFrame, model outputs
    :return: float
    """
    Qr = ydf['Qr'].iloc[-1]
    return Qr ** 2

# Iterate over FMUs
horizons = [2, 4, 6, 8]
relerr = [0, 10, 20, 30]  # %
runs = 5

for fmu in fmus:
    for hrz in horizons:
        for re in relerr:
            for run in range(1, runs + 1):

                fmu_name = fmu.split('.')[0]
                par_file = os.path.join('examples', 'case1', 'results', 'est',
                                        fmu_name, 'parameters.csv')
                wdir = os.path.join('examples', 'case1', 'results', 'mpc',
                                    fmu_name, "h{}".format(hrz),
                                    "re{}".format(re), "{}".format(run))

                # Skip the case, if results already there
                if os.path.exists(wdir):
                    pass

                else:
                    os.makedirs(wdir)
                    fmu_file = os.path.join(fmu_dir, fmu)

                    # Read parameters and modify heating power (used for heating/cooling in this example)
                    pm = pd.read_csv(par_file)
                    parameters = dict()
                    for p in pm:
                        parameters[p] = pm.iloc[0][p]
                    parameters['maxHeat'] = 6000.  # [W]

                    # Instantiate emulation and control models
                    model_emu = mshoot.SimFMU(
                        fmu_file,
                        outputs=['T', 'Qr', 'vetot'],
                        states=['cair.T'],
                        parameters=parameters,
                        verbose=False)

                    model_ctr = mshoot.SimFMU(
                        fmu_file,
                        outputs=['T', 'Qr', 'vetot'],
                        states=['cair.T'],  # States should be initialized with fixed=False
                        parameters=parameters,
                        verbose=False)

                    # Instantiate MPCEmulation
                    mpc = mshoot.MPCEmulation(model_emu, cfun)
                    step = 1
                    horizon = hrz

                    # Contraints
                    Tmin = np.where((ms.index.hour >= 8) & (ms.index.hour < 17), 21 + 273.15, 19 + 273.15)
                    Tmax = np.where((ms.index.hour >= 8) & (ms.index.hour < 17), 22 + 273.15, 24 + 273.15)

                    constr = pd.DataFrame(data=np.column_stack((Tmin, Tmax)),
                                        columns=['Tmin', 'Tmax'], index=inp.index)
                    constr.to_csv(os.path.join(wdir, 'constr.csv'))

                    # Add error to forecasts (control input): Tout, solrad, occ
                    inp_ctr = inp.copy()
                    inp_emu = inp.copy()
                    n = inp.index.size

                    mae_Tout = re * (inp['Tout'].max() - inp['Tout'].min()) / 100.
                    mae_solrad = re * (inp['solrad'].max() - inp['solrad'].min()) / 100.
                    mae_occ = re * (inp['occ'].max() - inp['occ'].min()) / 100.

                    inp_ctr['Tout'] = inp_ctr['Tout'] + error(0.9, n, mae_Tout, eps=0.01)

                    inp_ctr['solrad'] = np.where(inp_ctr['solrad'] != 0, inp_ctr['solrad'] + error(0.9, n, mae_solrad, eps=0.01), inp_ctr['solrad'])
                    inp_ctr['solrad'] = np.maximum(np.zeros(inp_ctr.index.size), inp_ctr['solrad'])

                    inp_ctr['occ'] = np.where(inp_ctr['occ'] != 0, inp_ctr['occ'] + error(0.9, n, mae_occ, eps=0.01), inp_ctr['occ'])
                    inp_ctr['occ'] = np.maximum(np.zeros(inp_ctr.index.size), inp_ctr['occ'])

                    inp_ctr.to_csv(os.path.join(wdir, 'inp_ctr.csv'))
                    inp_emu.to_csv(os.path.join(wdir, 'inp_emu.csv'))

                    # Run
                    t0 = time.time()
                    u, xctr, xemu, yemu, u_hist = mpc.optimize(
                        model=model_ctr,
                        inp_ctr=inp_ctr,
                        inp_emu=inp_emu,
                        free=['vpos'],
                        ubounds=[(-100., 100.)],
                        xbounds=[(Tmin, Tmax)],
                        x0=x0,
                        maxiter=50,
                        ynominal=[300., 1e7, 2e6],
                        step=step,
                        horizon=horizon
                    )
                    cputime = int(time.time() - t0)

                    # Save results
                    u.to_csv(os.path.join(wdir, 'u.csv'))
                    xctr.to_csv(os.path.join(wdir, 'xctr.csv'))
                    xemu.to_csv(os.path.join(wdir, 'xemu.csv'))
                    yemu.to_csv(os.path.join(wdir, 'yemu.csv'))

                    for i in range(len(u_hist)):
                        u_hist[i].to_csv(os.path.join(wdir, 'u{}.csv'.format(i)))

                    with open(os.path.join(wdir, 'cputime.txt'), 'w') as f:
                        f.write("CPU time: {} s".format(cputime))
