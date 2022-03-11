import numpy as np
import scipy as sp
import scipy.constants as sc
import sys,time,warnings


from numpy              import log, exp
from scipy.optimize     import differential_evolution
from scipy.integrate    import cumtrapz, odeint
from concurrent         import futures
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor
from output_redirect    import *
from compute            import *

np.set_printoptions(precision=4)
warnings.filterwarnings('error','.*Excess.*',)


def main():

    print('Multiple Processing Start')
    with ProcessPoolExecutor(max_workers=8) as ex:
        Pn_list = []
        for mode in (0,1):
            for num in (1,):
                Pn = ex.submit(computing,mode,num)
                Pn_list.append(Pn)
                print(' mode %2d num %2d running\n'%(mode,num))
    print('Multiple Processing Done')
    for Pn_result in futures.as_completed(Pn_list):
        print('Multiple: result: {}\n'.format(Pn.result()))


if __name__ == '__main__':
    main()


