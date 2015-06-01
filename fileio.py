import h5py
import numpy as np

def readio():
    f = h5py.File("testinput.hdf5", "r")
    dset = f['initial']
    return dset
""""""

def readOut(fname='testout'):
    f = h5py.File("%s.hdf5"%(fname),"r")
    dset = f['output']
    return dset
""""""

def writeOut(solutions=None, Jrad=None, vshock=None, fname='testout', ngas=3, ndust=1):
    """
        Write the output
        solutions: x, vgas, tgas, nh, nh2, nhe, ndust, vdust, tdust, adust

    """
    print solutions.shape
    grid = solutions[:,0]
    vgas = solutions[:,1]
    tgas = solutions[:,2]
    ndust = solutions[:,6]
    vdust = solutions[:,7]
    tdust = solutions[:,8]
    adust = solutions[:,9]
    with h5py.File("%s.hdf5"%(fname), 'w') as fout:
        dset = fout.create_dataset("output", solutions[:,3:6].shape, dtype='f')
        dset[:,:] = solutions[:, 3:6]/vgas[:,np.newaxis]
        #
        # write the grid out
        #
        dset.attrs['ngrid'] = grid.shape[0]
        dset.attrs['x'] = grid
        #
        # Dust
        #
        dset.attrs['ndust'] = ndust
        dset.attrs['vdust'] = vdust
        dset.attrs['tdust'] = tdust
        dset.attrs['adust'] = adust
        #
        # Other information
        #
        dset.attrs['vshock'] = vshock
        dset.attrs['tgas'] = tgas
        dset.attrs['vgas'] = vgas
        #
        # Radiation
        #
        dset.attrs['tau'] = Jrad[0]
        dset.attrs['Jrad'] = Jrad[1]
    """"""
    #
    # Done
    #
""""""
