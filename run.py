### CHANGE PARAMETERS ####

# general parameters
EXAMPLE_LIST = ['toy']
T = 60  # time period

# hyperparameters
TIMELIMIT_LIST = [30, 50, 60, 70, 90]
TIME_INCREASE = 60
PATIENCE = 11
MAXTIME = 120


### Run script
BASEPATH = './'
PATH_TO_MODULE = f'{BASEPATH}src/'

import sys
sys.path.append(BASEPATH)
sys.path.append(PATH_TO_MODULE)

from utils import (
    PTN,
    EAN,
    preprocess_data
)
from solvers import (
    ColppSolver_beta,
    PespSolver,
    LinTimApprox,
    LintimSolver,
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
                MODE_TIME_INC=True,  # mode, that we will incerase time if non-improvement
                TIME_INCREASE=TIME_INCREASE,
                MAXTIME=MAXTIME)
            
        except Exception as e:
            print(e)
            logging.error(f'Error in LinTimApprox for {EXAMPLE} with TIMELIMIT {TIMELIMIT}: {e}')



## Run visualization (todo)