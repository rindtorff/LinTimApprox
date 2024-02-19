import parameters

EXAMPLE_LIST = parameters.EXAMPLE_LIST
T = parameters.T
TIMELIMIT_LIST = parameters.TIMELIMIT_LIST
TIME_INCREASE = parameters.TIME_INCREASE
PATIENCE = parameters.PATIENCE
MAXTIME = parameters.MAXTIME
PLOT_RESULTS = parameters.PLOT_RESULTS

### Run script
BASEPATH = './'
PATH_TO_MODULE = f'{BASEPATH}src/'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append(BASEPATH)
sys.path.append(PATH_TO_MODULE)

from utils import (
    preprocess_data,
    evaluate_algo
)
from solvers import (
    LinTimApprox
)

for TIMELIMIT in TIMELIMIT_LIST:
    for EXAMPLE in EXAMPLE_LIST:

        PATH_TO_EXAMPLE = f'{BASEPATH}example_{EXAMPLE}/data/'

        import logging
        logging.basicConfig(
            filename=f'{BASEPATH}results/log_{EXAMPLE}.log',
            encoding='utf-8',
            level=logging.DEBUG)

        ptn, ean, linepool = preprocess_data(PATH_TO_EXAMPLE)

        try:
            lintimapprox = LinTimApprox(
                ptn=ptn,
                ean=ean,
                linepool=linepool,
                PATIENCE=PATIENCE,
                CSV_NAME=f'{BASEPATH}results/{EXAMPLE}_{TIMELIMIT}_{TIME_INCREASE}.csv',
                TIMELIMIT=TIMELIMIT,
                TIME_INCREASE=TIME_INCREASE,
                MAXTIME=MAXTIME)

        except Exception as e:
            print(e)
            logging.error(f'Error in LinTimApprox for {EXAMPLE} with TIMELIMIT {TIMELIMIT}: {e}')


# Run visualization
if PLOT_RESULTS:
    BASEPATH_mm = f'{BASEPATH}results/'
    for EXAMPLE in EXAMPLE_LIST:
        BASEPATH_sol = f'{BASEPATH}/example_{EXAMPLE}/'
        for TIME in TIMELIMIT_LIST:
            df_mm = pd.read_csv(
                f'{BASEPATH_mm}/{EXAMPLE}_{TIME}_60.csv', index_col=0)
            df_sol = pd.read_csv(f'{BASEPATH_sol}callback_{EXAMPLE}.csv')
            df_sol['Time'] = df_sol['Time'].apply(lambda x: np.round(x))
            df_sol = df_sol.groupby('Time')['Objective Value'].mean().to_frame()
            df_sol.reset_index(inplace=True)
            title_mm_sol = df_mm['nu*(LinTim)'].values[-1]
            title_sol = df_sol['Objective Value'].values[-1]
            gap_percent = np.abs((title_sol - title_mm_sol) / title_sol) * 100  

            if EXAMPLE == 'toy':
                df_sol = pd.read_csv(f'{BASEPATH_sol}callback_{EXAMPLE}.csv')
                title_mm_sol = df_mm['nu*(LinTim)'].values[-1]
                title_sol = df_sol['Objective Value'].values[-1]
                gap_percent = np.abs((title_sol - title_mm_sol) / title_sol) * 100  
                try:
                    evaluate_algo(
                        df_mm=df_mm,
                        df_sol=df_sol,
                        show_A=True,
                        show_sol=True,
                        title=f'Dataset {EXAMPLE.upper()} - Time: {TIME}s \nnu*: {title_mm_sol}, Solution: {title_sol:.2f}, Gap: {gap_percent:.2f}%',
                        grid=True,
                        xlim=(0, 0.05),
                        # ylim=(28000, 60000),
                        save_fig=f'{BASEPATH}plots/{EXAMPLE}_{TIME}')
                except Exception as e:
                    print(e)

            elif EXAMPLE == 'grid':
                try:
                    evaluate_algo(
                        df_mm=df_mm,
                        df_sol=df_sol,
                        show_A=True,
                        show_sol=True,
                        title=f'Dataset {EXAMPLE.upper()} - Time: {TIME}s \nnu*: {title_mm_sol}, Solution: {title_sol:.2f}, Gap: {gap_percent:.2f}%',
                        grid = True,
                        xlim = (0, 35),
                        ylim = (17800, 18200),
                        save_fig=f'{BASEPATH}plots/{EXAMPLE}_{TIME}')
                except Exception as e:
                    print(e)

            elif EXAMPLE == 'ring':
                try:
                    evaluate_algo(
                        df_mm=df_mm,
                        df_sol=df_sol,
                        show_A=True,
                        show_sol=True,
                        title=f'Dataset {EXAMPLE.upper()} - Time: {TIME}s \nnu*: {title_mm_sol}, Solution: {title_sol:.2f}, Gap: {gap_percent:.2f}%',
                        grid=True,
                        xlim=(0, 90),
                        ylim=(28000, 40000),
                        save_fig=f'{BASEPATH}plots/{EXAMPLE}_{TIME}')

                except Exception as e:
                    print(e)
