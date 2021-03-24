import numpy as np
import h5py

with h5py.File('jorek_restart.h5', 'r') as hf1:
  with h5py.File('end.h5', 'r') as hf2:
    print(np.amax(np.abs(
      np.array(hf1.get('values')) -
      np.array(hf2.get('values'))
    ), axis=(0,1,2,3)))
