import h5py
import numpy as np

def readio():
    f = h5py.File("testinput.hdf5", "r")
    dset = f['initial']
    return dset
""""""

def readOut(fname='testout'):
    f = h5py.File("%s.hdf5"%(fname),"r")
    return f
""""""

def writeOut(solutions=None, Jrad=None, Frad=None, vshock=None, fname='testout', ngas=3, ndust=1):
    """
        Write the output
        solutions: x, vgas, tgas, nh, nh2, nhe, ndust, vdust, tdust, adust

    """
    grid = solutions[:,0]
    vgas = solutions[:,1]
    tgas = solutions[:,2]
    ndust = solutions[:,3+ngas]
    vdust = solutions[:,4+ngas]
    tdust = solutions[:,5+ngas]
    adust = solutions[:,6+ngas]
    with h5py.File("%s.hdf5"%(fname), 'w') as fout:
        dset1 = fout.create_dataset("grid", grid.shape, dtype='f')
        dset1[:] = grid[:]
        gas = solutions[:,3:3+ngas]/vgas[:,np.newaxis]
        gas = np.insert(gas, 0, tgas, axis=1)
        gas = np.insert(gas, 0, vgas, axis=1)
        dset = fout.create_dataset("gas", gas.shape, dtype='f')
        dset[:,:] = gas[:,:]
        #
        # write the grid out
        #
        dset.attrs['ngrid'] = grid.shape[0]
        #
        # Dust
        #
        dust = np.array([vdust, tdust, ndust, adust])
        dust = np.transpose(dust)
        dset2 = fout.create_dataset("dust", dust.shape, dtype='f')
        dset2[:,:] = dust[:,:]
        #
        # Other information
        #
        dset.attrs['vshock'] = vshock
        #
        # Radiation
        #
        radiation = np.array([Jrad[0], Jrad[1][1], Frad[:]])
        radiation =  np.transpose(radiation)
        dset3 = fout.create_dataset("radiation", radiation.shape, dtype='f')
        dset3[:,:] = radiation[:,:]
    """"""
    #
    # Done
    #
""""""
