; IDL routine to save HDF5 files 

;Input directory
;fname='../Data'
;fname='../testdata/OT-MHD'
;fname='../testdata/OT-MHD-res'
;fname='../testdata/OT-PIP'
;fname='../testdata/RT-MHD-O4'
;fname='../testdata/shockstab-MHD-O4'
;fname='../testdata/shocktube-IRIP'
;fname='../testdata/shocktube-MHD'
;fname='../testdata/shocktube-MHD-O4'
;fname='../testdata/shocktube-PIP'
;fname='/home/snow/datadrive/simdata/MHD-OZ-3D'
fname='/media/snow/Plasma/kink_instability_MHD_0/'
;Define save directory MUST ALREADY EXIST
savdir='/media/snow/store1/kink_instability/MHD_0/'

;Define save directory MUST ALREADY EXIST
;savdir=fname+'/';HDF5/'

;Variables to save
;outvar=['xgrid','ro_p','vx_p','vy_p','vz_p','bx','by','bz','pr_p']
;outvar=['xgrid','ro_p','mx_p','my_p','mz_p','bx','by','bz','en_p']
outvar=['xgrid','ro_p','mx_p','my_p','mz_p','bx','by','bz','en_p']

;Read in the data
for tread=0,20 do begin
;tread=30
rdmpi,ds,datapath=fname,time_step=tread,var=['ro_p','mx_p','my_p','mz_p','bx','by','bz','en_p']

;Number of dimensions. Might be a neater way but this works.
ndim=n_elements(size(ds.ro_p))-3

;Define simulation grid
x=ds.x
if ndim ge 2 then y=ds.y
if ndim ge 3 then z=ds.z
ro=ds.ro_p
;vx=ds.vx_p
;vy=ds.vy_p
;vz=ds.vz_p
bx=ds.bx
by=ds.by
bz=ds.bz
;pr_p=ds.pr_p

if tread eq 0 then begin
if ndim ge 2 then outvar=['ygrid',outvar] 
if ndim ge 3 then outvar=['zgrid',outvar]
endif

    ;; get data type and space, needed to create the dataset

;    file = savdir+outvar(i)+'.h5'
    file = savdir+'t.'+string(tread,FORMAT='(I4.4)')+'.h5'
    fid = H5F_CREATE(file)

    ;datatype_id = H5T_IDL_CREATE(ro)
    ;dataspace_id = H5S_CREATE_SIMPLE(size(ro_p,/DIMENSIONS))

for i=0,n_elements(outvar)-1 do begin
;for i=0,1 do begin

    if outvar(i) eq 'xgrid' then data = x
    if outvar(i) eq 'ygrid' then data = y
    if outvar(i) eq 'zgrid' then data = z
    if outvar(i) eq 'ro_p' then data = ro
    if outvar(i) eq 'vx_p' then data = vx
    if outvar(i) eq 'vy_p' then data = vy
    if outvar(i) eq 'vz_p' then data = vz
    if outvar(i) eq 'bx' then data = bx
    if outvar(i) eq 'by' then data = by
    if outvar(i) eq 'bz' then data = bz
    if outvar(i) eq 'pr_p' then data = pr_p

    if outvar(i) eq 'mx_p' then data = ds.mx_p
    if outvar(i) eq 'my_p' then data = ds.my_p
    if outvar(i) eq 'mz_p' then data = ds.mz_p
    if outvar(i) eq 'en_p' then data = ds.en_p

    datatype_id = H5T_IDL_CREATE(data)
    dataspace_id = H5S_CREATE_SIMPLE(size(data,/DIMENSIONS))

print,outvar(i),datatype_id,dataspace_id
    
    dataset_id = H5D_CREATE(fid,$

    outvar(i),datatype_id,dataspace_id)

    ;; write data to dataset

    H5D_WRITE,dataset_id,data

    ;; close all open identifiers

    H5D_CLOSE,dataset_id

    H5S_CLOSE,dataspace_id

    H5T_CLOSE,datatype_id

endfor

H5F_CLOSE,fid

endfor

END
