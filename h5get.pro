PRO h5get,pv,file,var,grid

tmp=create_struct("info",0)

if grid eq 0 then begin

    print,'Gamma assumed as 5/3'
    gm=5.0/3.0
    ;make array for needed components
    nv=n_elements(var)
    vartarr=[]
    for v=0,nv-1 do begin
	    if (var(v) eq 'bx') or (var(v) eq 'by') or (var(v) eq 'bz') or (var(v) eq 'ro_p') or (var(v) eq 'ro_n') then vartarr=[vartarr,var(v)]
	    if (var(v) eq 'vx_p') then vartarr=[vartarr,'mx_p','ro_p']
	    if (var(v) eq 'vy_p') then vartarr=[vartarr,'my_p','ro_p']
	    if (var(v) eq 'vz_p') then vartarr=[vartarr,'mz_p','ro_p']
	    if (var(v) eq 'pr_p') then vartarr=[vartarr,'en_p','mx_p','my_p','mz_p','ro_p','bx','by','bz']

	    if (var(v) eq 'vx_n') then vartarr=[vartarr,'mx_n','ro_n']
	    if (var(v) eq 'vy_n') then vartarr=[vartarr,'my_n','ro_n']
	    if (var(v) eq 'vz_n') then vartarr=[vartarr,'mz_n','ro_n']
	    if (var(v) eq 'pr_n') then vartarr=[vartarr,'en_n','mx_n','my_n','mz_n','ro_n']
    endfor
    vres=vartarr[uniq(vartarr,sort(vartarr))]
    ;print,vres

    ;Read data
    for i=0,n_elements(vres)-1 do begin
        ;print,var(i)
        file_id = H5F_OPEN(file)

        dataset_id1 = H5D_OPEN(file_id, vres(i))

        ; Read in the actual image data.
        image = H5D_READ(dataset_id1)

        H5D_CLOSE, dataset_id1
        H5F_CLOSE, file_id

        tmp = create_struct(tmp,vres(i),reform(image))

    endfor

    for vloop=0,n_elements(var)-1 do begin
    ;print,var(vloop)

    if var(vloop) eq 'vx_n' then pv=create_struct(pv,var(vloop),reform(tmp.mx_n/tmp.ro_n))
    if var(vloop) eq 'vy_n' then pv=create_struct(pv,var(vloop),reform(tmp.my_n/tmp.ro_n))
    if var(vloop) eq 'vz_n' then pv=create_struct(pv,var(vloop),reform(tmp.mz_n/tmp.ro_n))

    if var(vloop) eq 'vx_p' then pv=create_struct(pv,var(vloop),reform(tmp.mx_p/tmp.ro_p))
    if var(vloop) eq 'vy_p' then pv=create_struct(pv,var(vloop),reform(tmp.my_p/tmp.ro_p))
    if var(vloop) eq 'vz_p' then pv=create_struct(pv,var(vloop),reform(tmp.mz_p/tmp.ro_p))

    if var(vloop) eq 'bx' then pv=create_struct(pv,var(vloop),reform(tmp.bx))
    if var(vloop) eq 'by' then pv=create_struct(pv,var(vloop),reform(tmp.by))
    if var(vloop) eq 'bz' then pv=create_struct(pv,var(vloop),reform(tmp.bz))

    if var(vloop) eq 'ro_n' then pv=create_struct(pv,var(vloop),reform(tmp.ro_n))
    if var(vloop) eq 'ro_p' then pv=create_struct(pv,var(vloop),reform(tmp.ro_p))

    if var(vloop) eq 'pr_n' then pv=create_struct(pv,var(vloop),(gm-1.0)*(reform(tmp.en_n)-0.5*reform(tmp.ro_n)*(reform(tmp.mx_n/tmp.ro_n)^2+reform(tmp.my_n/tmp.ro_n)^2+reform(tmp.mz_n/tmp.ro_n)^2)))
    if var(vloop) eq 'pr_p' then pv=create_struct(pv,var(vloop),(gm-1.0)*(reform(tmp.en_p)-0.5*reform(tmp.ro_p)*(reform(tmp.mx_p/tmp.ro_p)^2+reform(tmp.my_p/tmp.ro_p)^2+reform(tmp.mz_p/tmp.ro_p)^2) -0.5*(reform(tmp.bx)^2+reform(tmp.by)^2+reform(tmp.bz)^2)))
    endfor

endif

if grid eq 1 then begin
    ;Read data
    for i=0,n_elements(var)-1 do begin
        ;print,var(i)
        file_id = H5F_OPEN(file)

        dataset_id1 = H5D_OPEN(file_id, var(i))

        ; Read in the actual image data.
        image = H5D_READ(dataset_id1)

        H5D_CLOSE, dataset_id1
        H5F_CLOSE, file_id

        pv = create_struct(pv,var(i),reform(image))

    endfor
endif

return

END
