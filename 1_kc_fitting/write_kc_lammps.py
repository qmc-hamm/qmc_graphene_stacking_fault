'''
Generate LAMMPS KC potential files
'''
import os
import numpy as np
import h5py

import kc

if __name__ == '__main__':
    kTs = [f'{kT:.3f}' for kT in np.arange(2, 21)*0.001] + ['inf']
    workdir = os.getcwd()
    for kT in kTs:
        dirname = f'fit/QMC_kT{kT}'
        print(dirname)
        os.chdir(dirname)

        if os.path.isfile('CH_taper.KC'):
            os.remove('CH_taper.KC')
        with h5py.File('result.hdf5', 'r') as f:
            popt = f['popt'][()]

        params = popt[:-1]
        print(params)
        kc.write_kc_potential(*params, kc_filename='CC_QMC.KC', cite=True)
        os.chdir(workdir)
