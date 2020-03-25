def get_power(name, z):

   import numpy as np

   k, zs, Pks = get_powers(name)

   i = np.where(zs == z)

   return k, Pks[:,i[0]]

def get_powers(name):

   import numpy as np

   indir = '/Users/Mead/Physics/data/VD20'
   infile = indir+'/'+name+'.dat'
   data = np.loadtxt(infile)

   nz = 15
   nk = 352

   z = np.unique(data[:, 0])
   k = np.unique(data[:, 1])
   Pk = np.swapaxes(data[:,3].reshape(nz, nk), 0, 1)

   z = np.flip(z)

   return k, z, Pk

def get_response(name, z):

   import numpy as np

   k, zs, Rks = get_responses(name)
   
   i = np.where(zs == z)

   return k, Rks[:, i[0]]

def get_responses(name):

   k, z, Pk = get_powers(name)

   name_dmonly = dmonly_counterpart(name)

   k, z, Pk_dmonly = get_powers(name_dmonly)

   return k, z, Pk/Pk_dmonly

def dmonly_counterpart(name):

   match = {
      'BAHAMAS_Theat7.6_nu0_WMAP9': 'DMONLY_2fluid_nu0_WMAP9_L400N1024',
      'BAHAMAS_Theat8.0_nu0_WMAP9': 'DMONLY_2fluid_nu0_WMAP9_L400N1024',
      'BAHAMAS_nu0_WMAP9': 'DMONLY_2fluid_nu0_WMAP9_L400N1024',
      'BAHAMAS_nu0_WMAP9_v2': 'DMONLY_2fluid_nu0_v2_WMAP9_L400N1024',
      'BAHAMAS_nu0_WMAP9_v3': 'DMONLY_2fluid_nu0_v3_WMAP9_L400N1024',
      'BAHAMAS_nu0.06_Planck2015': 'DMONLY_2fluid_nu0.06_Planck2015_L400N1024',
      'BAHAMAS_nu0.06_WMAP9': 'DMONLY_2fluid_nu0.06_WMAP9_L400N1024',
      'BAHAMAS_nu0.12_Planck2015': 'DMONLY_2fluid_nu0.12_Planck2015_L400N1024',
      'BAHAMAS_nu0.12_WMAP9': 'DMONLY_2fluid_nu0.12_WMAP9_L400N1024',
      'BAHAMAS_nu0.24_Planck2015': 'DMONLY_2fluid_nu0.24_Planck2015_L400N1024',
      'BAHAMAS_nu0.24_WMAP9': 'DMONLY_2fluid_nu0.24_WMAP9_L400N1024',
      'BAHAMAS_nu0.48_Planck2015': 'DMONLY_2fluid_nu0.48_Planck2015_L400N1024',
      'BAHAMAS_nu0.48_WMAP9': 'DMONLY_2fluid_nu0.48_WMAP9_L400N1024',
      'BAHAMAS_nu0_BAO_L200N512': 'DMONLY_nu0_BAO_L200N512',
      'BAHAMAS_nu0_Planck2013': 'DMONLY_2fluid_nu0_Planck2013_L400N1024',
      'BAHAMAS_nu0_WMAP9_L100N512': 'DMONLY_2fluid_nu0_WMAP9_L100N512',
      'C-OWLS_AGN_Planck2013': 'DMONLY_Planck2013_L400N1024',
      'C-OWLS_AGN_Theat8.5_Planck2013': 'DMONLY_Planck2013_L400N1024',
      'C-OWLS_AGN_Theat8.7_Planck2013': 'DMONLY_Planck2013_L400N1024',
      'C-OWLS_NOCOOL_UVB_Planck2013': 'DMONLY_Planck2013_L400N1024',
      'C-OWLS_REF_Planck2013': 'DMONLY_Planck2013_L400N1024',
      'C-OWLS_AGN_WMAP7': 'DMONLY_WMAP7_L400N1024',
      'C-OWLS_AGN_Theat8.5_WMAP7': 'DMONLY_WMAP7_L400N1024',
      'C-OWLS_AGN_Theat8.7_WMAP7': 'DMONLY_WMAP7_L400N1024',
      'C-OWLS_NOCOOL_UVB_WMAP7': 'DMONLY_WMAP7_L400N1024',
      'C-OWLS_REF_WMAP7': 'DMONLY_WMAP7_L400N1024',
      }

   return match[name]