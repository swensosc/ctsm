#! /usr/bin/env python
import subprocess
import copy
import numpy as np
import netCDF4 as netcdf4

# creates surface dataset with max patch format

dir_surfdata = '/fs/cgd/csm/inputdata/lnd/clm2/surfdata_esmf/'
dir_output=''

fsurf = f'{dir_surfdata}ctsm5.3.0/surfdata_0.9x1.25_hist_2000_78pfts_c240908.nc'
fsurf_out = f'{dir_output}surfdata_0.9x1.25_hist_2000_ctsm5.3.0_maxpatch_8.nc'
# add new patches
addNewPatches = True
if addNewPatches:
    fsurf_out = f'{dir_output}surfdata_0.9x1.25_hist_2000_ctsm5.3.0_maxpatch_8_new_subtypes.nc'

paramfile = 'ctsm60_params_4_pft_subtypes.nc'

#--  Check whether file exists  ---------------------------------
command=['ls',fsurf_out]
file_exists=subprocess.call(command,stderr=subprocess.PIPE)
if file_exists > 0:
    print('creating new file: ', fsurf_out)
else:
    print('overwriting file: ', fsurf_out)

#--  Open input file
f =  netcdf4.Dataset(fsurf, 'r')
global_attributes  = f.ncattrs()
variables = f.variables
dimensions = f.dimensions
jm, im = len(dimensions['lsmlat']),len(dimensions['lsmlon'])
lon = np.asarray(f.variables['LONGXY'][0,])
lat = np.asarray(f.variables['LATIXY'][:,0])

if 'PFTDATA_MASK' in f.variables.keys():            
    landmask = np.asarray(f.variables['PFTDATA_MASK'][:,])
if 'LANDFRAC_PFT' in f.variables.keys(): # ctsm >= 5.3.0
    landmask = np.where(f.variables['LANDFRAC_PFT'][:,]>0,1,0)

#--  Specify max number of patches per gridcell/landunit
maxpatch = 8

#--  Read in patch to pft mapping from parameter file
f2 =  netcdf4.Dataset(paramfile, 'r')
pft_index   = f2.variables['pftnum'][:,]
patch_index = f2.variables['nat_patch_type'][:,]
f2.close()

#--  Determine pft percentages and indices
pct_nat_patch = np.zeros((maxpatch,jm,im))
nat_patch_subtype = np.zeros((maxpatch,jm,im))
#--  Initialize all points to 100% bare soil, including ocean points
pct_nat_patch[0,] = 100 

for j in range(jm):
    for i in range(im):
        if landmask[j,i] > 0:
            ind = np.where(variables['PCT_NAT_PFT'][:,j,i] > 0)[0]            
            if ind.size > 0:
                sind = np.argsort(variables['PCT_NAT_PFT'][ind,j,i])
                # ensure number of patches <= maxpatch
                if ind.size > maxpatch:
                    ind = ind[sind[::-1]][:maxpatch]
                    # scale so sum = 100
                    sf = 100/np.sum(variables['PCT_NAT_PFT'][ind,j,i])
                    pct_nat_patch[:ind.size,j,i]  = sf*variables['PCT_NAT_PFT'][ind,j,i]
                else:
                    ind = ind[sind[::-1]]
                    pct_nat_patch[:ind.size,j,i]  = variables['PCT_NAT_PFT'][ind,j,i]
                
                # map pft index to new patch (subtype) index
                for k in range(ind.size):
                    indp = np.where(pft_index==ind[k])[0]
                    nat_patch_subtype[k,j,i] = patch_index[indp[0]]
                    # make some changes
                    if addNewPatches:
                        if ind[k]==2 and lon[i]>250: # NEBT
                            nat_patch_subtype[k,j,i] = patch_index[indp[0]]+1
                        if ind[k]==4 and lon[i]>300: # BETT
                            nat_patch_subtype[k,j,i] = patch_index[indp[0]]+1
                        if ind[k]==14 and lon[i]<30: # C4G
                            nat_patch_subtype[k,j,i] = patch_index[indp[0]]+1

#--  Open output file
w = netcdf4.Dataset(fsurf_out, 'w', format='NETCDF3_64BIT')

#--  Set global attributes
for ga in global_attributes:
    setattr(w,ga,f.getncattr(ga))
#--  Set dimensions of output file
for dim in dimensions.keys():
    w.createDimension(dim,len(dimensions[dim]))
# add new dimension
w.createDimension('nmaxpatch',int(maxpatch))

# create new variables (nat_patch_type, nat_patch_subtype)
wvar = w.createVariable('nat_patch_type',np.int32,('nmaxpatch','lsmlat','lsmlon'))
wvar.long_name = 'natural pft subtype index'
wvar.units = 'unitless'
wvar[:,] = nat_patch_subtype
wvar = w.createVariable('pct_nat_patch',float,('nmaxpatch','lsmlat','lsmlon'))
wvar.long_name = 'percent natural pft subtype'
wvar.units = 'percent'
wvar[:,] = pct_nat_patch

for var in variables.keys():
    vdim = f.variables[var].dimensions
    vtype = f.variables[var].datatype
    print(var, vtype, vdim)
    # exclude PCT_NAT_PFT from new surface data file
    if var == 'PCT_NAT_PFT':
        continue
    wvar = w.createVariable(var, vtype, vdim)
    #--  Write variable data to output file
    wvar[:,] = f.variables[var][:,]

    #--  Set attribute values
    att=f.variables[var].ncattrs()
    print(att, '\n')
    km=len(att)
    for attname in att:
        print('name: ',attname,' value: ',f.variables[var].getncattr(attname))
        w.variables[var].setncattr(attname,f.variables[var].getncattr(attname))

#--  Close output file
w.close
print('created ',fsurf_out)

