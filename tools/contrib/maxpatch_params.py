#! /usr/bin/env python
import subprocess
import copy
import numpy as np
import netCDF4 as netcdf4

# creates surface dataset with max patch format

dir_paramdata = '/fs/cgd/csm/inputdata/lnd/clm2/paramdata/'
dir_output    = ''

fparam     = f'{dir_paramdata}ctsm60_params.c241119.nc'
# specify number of pft subtypes via pft_scalar
pft_scalar = 4
fparam_out = f'{dir_output}ctsm60_params_{pft_scalar}_pft_subtypes.nc'

# add new patches
addNewPatches = True
if addNewPatches:
    fparam_out = f'{dir_output}ctsm60_params_{pft_scalar}_pft_subtypes_new_types.nc'

#--  Check whether file exists  ---------------------------------
command=['ls',fparam_out]
file_exists=subprocess.call(command,stderr=subprocess.PIPE)
if file_exists > 0:
    print('creating new file: ', fparam_out)
else:
    print('overwriting file: ', fparam_out)

#--  Open input file
f =  netcdf4.Dataset(fparam, 'r')
global_attributes  = f.ncattrs()
variables = f.variables
dimensions = f.dimensions

nbare = 1
ncft = 64
# set number of pft subtypes via pft_scalar
npft = len(dimensions['pft'])
npft_natveg = npft - ncft - nbare
pft_scalar = 4
npatch = pft_scalar * npft_natveg

new_npft = nbare + npatch + ncft
print(npft_natveg,new_npft)
npft_natveg_plus_bare = npft_natveg + nbare
new_npft_natveg_plus_bare = nbare + npatch

# create map from larger array to smaller array
cftmap = [npft_natveg_plus_bare+i for i in range(ncft)]
pftmap = [1+int(i//pft_scalar) for i in range(npatch)]
mind = [0]+pftmap+cftmap
if len(mind) != new_npft:
    print('array size incorrect')
    print(len(mind),new_npft)
    stop

pftnames = []
x = np.asarray(f.variables['pftname'][:,])
for y in x:
    pftnames.append(''.join([i.decode('utf-8') for i in y]))
pftnames = np.asarray(pftnames)
print(pftnames)

pftnames2 = []
x = np.take(x,np.asarray(mind,dtype=int),axis=0)
for y in x:
    pftnames2.append(''.join([i.decode('utf-8') for i in y]))
pftnames2 = np.asarray(pftnames2)
print(pftnames2[x2])

#--  Open output file
w = netcdf4.Dataset(fparam_out, 'w', format='NETCDF3_64BIT')

#--  Set global attributes
for ga in global_attributes:
    setattr(w,ga,f.getncattr(ga))
#--  Set dimensions of output file
for dim in dimensions.keys():
    if dim == 'pft':
        w.createDimension(dim,new_npft)
    else:
        w.createDimension(dim,len(dimensions[dim]))

# create patch index variable
wvar = w.createVariable('nat_patch_type',np.int32,('pft'))
wvar.long_name = 'natural patch type index'
wvar.units = 'unitless'
wvar[:,] = np.arange(new_npft,dtype=int)

for var in variables.keys():
    vdim = f.variables[var].dimensions
    vtype = f.variables[var].datatype
    print(var, vtype, vdim)
    wvar = w.createVariable(var, vtype, vdim)

    #--  Set attribute values
    att = f.variables[var].ncattrs()
    print(att, '\n')
    km = len(att)
    for attname in att:
        print('name: ',attname,' value: ',f.variables[var].getncattr(attname))
        w.variables[var].setncattr(attname,f.variables[var].getncattr(attname))
    
    if 'pft' in vdim:
        ploc = vdim.index('pft')
        wvar[:,] = np.take(f.variables[var][:,],np.asarray(mind,dtype=int),axis=ploc)
        # mergetoclmpft must be treated separately
        if var=='mergetoclmpft':
            wvar[1:new_npft_natveg_plus_bare] = np.arange(1,new_npft_natveg_plus_bare)
            wvar[new_npft_natveg_plus_bare:] += (new_npft_natveg_plus_bare - npft_natveg_plus_bare)
        # create some new patch values
        if addNewPatches:
            if var=='medlynslope':
                wvar[6] = 1  # low slope for BETT
                wvar[14] = 20 # high slope NEBT
                wvar[54] = 0.1 # high slope C4G
        
    else:
        wvar[:,] = f.variables[var][:,]
        

#--  Close output file
w.close
print('created ',fparam_out)

