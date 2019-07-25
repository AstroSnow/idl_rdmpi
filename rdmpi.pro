pro rdmpi,pv,datapath=datapath,current=current,flag_double=flag_double,$
           time_step=time_step,flag_az=flag_az,flag_te=flag_te
  Compile_Opt DEFINT32  
  if(n_elements(datapath)eq 0 ) then datapath="tmp/" else datapath=datapath+"/"
  if(n_elements(current) eq 0) then current=0
  if(n_elements(flag_double) eq 0) then flag_double=1
  if(n_elements(time_step) eq 0) then time_step=-1
  if(n_elements(flag_az) eq 0) then flag_az = 0 ;; vector potential (for 2D)
  if(n_elements(flag_te) eq 0) then flag_te = 0 ;; temperature
  
  get_param2,datapath,info
  gm=info.gm
  pv=create_struct("info",info)
  ts=info.nt
  if time_step[0] eq -1 then n_read=ts else n_read=n_elements(time_step)

;;;;; test loop to prevent memory overflow
;if n_read GT 151 then n_read=151  

;define margin------------------------------
  margin=info.margin
;-------------------------------------------
  eqs=info.eqs
  flag_mhd=0
  flag_pip=0
  flag_afr=0
  case eqs of 
     'HD':nvar=5
     'MHD':begin
        nvar=8
        flag_mhd=1
     end
     'PIP':begin
        nvar=13
        flag_mhd=1
        flag_pip=1
     end
     "AFR":begin
;        nvar=29
        nvar=13
        flag_mhd=1
        flag_pip=1
        flag_afr=1
     end
  end
  
  pv=create_struct(pv,["eqs","fl_mhd","fl_pip","fl_afr"], $
                   eqs,flag_mhd,flag_pip,flag_afr)
  
  tfile=file_search(datapath+"t.dac.*")
  dacget0s,tfile,t,narg=time_step
  if(eqs eq "AFR") then begin
     gridfile=file_search(datapath+"region.dac.*")
     dacget2s,gridfile,grid,narg=time_step
     pv=create_struct(pv,["grid"],grid)
  endif
  ix=info.ix
  jx=info.jx
  kx=info.kx
  
  if (where(tag_names(info) eq "MPI_DOMAINS"))[0] ne -1 then begin
     mpi_x=info.mpi_domains[0]
     mpi_y=info.mpi_domains[1]
     mpi_z=info.mpi_domains[2]
  endif else begin
     mpi_x=(mpi_y=(mpi_z=1))  
  endelse
  ix_m=(jx_m=(kx_m=1)) 
;  xfile=file_search(datapath+"x.dac.*")
  xfile=file_search(datapath+"x.dac."+string(indgen(mpi_x),form="(i4.4)"))
  ix_m=info.ix
  ix=margin[0]*2+mpi_x*(ix_m-margin[0]*2)
  x=findgen(ix)
  for n=0,mpi_x-1 do begin
     dacget1d,xfile[n],xc
     x0=n*(ix_m-2*margin[0])
     x[x0:x0+ix_m-1]=xc     
  endfor

  pv=create_struct(pv,"t",t,"x",x)
  ndim=pv.info.ndim
  if ndim ge 2 then begin     
     yfile=file_search(datapath+"y.dac."+string(indgen(mpi_y),form="(i4.4)"))
;       yfile=file_search(datapath+"y.dac.*")
     jx_m=info.jx
     jx=margin[1]*2+mpi_y*(jx_m-margin[1]*2)
     y=findgen(jx)
     for n=0,mpi_y-1 do begin
        dacget1d,yfile[n],yc
        y0=n*(jx_m-2*margin[1])
        y[y0:y0+jx_m-1]=yc     
     endfor
     pv=create_struct(pv,"y",y)
     if ndim ge 3 then begin
        zfile=file_search(datapath+"z.dac."+string(indgen(mpi_z),form="(i4.4)"))
        kx_m=info.kx
        kx=margin[2]*2+mpi_z*(kx_m-margin[2]*2)
        z=findgen(kx)
        for n=0,mpi_z-1 do begin
           dacget1d,zfile[n],zc
           z0=n*(kx_m-2*margin[2])
           z[z0:z0+kx_m-1]=zc     
        endfor
        pv=create_struct(pv,"z",z)
     endif
  endif
  pv=create_struct(pv,"ix",n_elements(x))
  if ndim ge 2 then jx=n_elements(y) else jx=1  
  if ndim ge 3 then kx=n_elements(z) else kx=1  
  pv=create_struct(pv,"jx",jx,"kx",kx)  

  n_cpu=mpi_x*mpi_y*mpi_z
  var_names=[]
  if time_step[0] eq -1 then begin
     time_step=indgen(n_read)
  endif
  output=dblarr(ix,jx,kx)


  if(0 eq 1 ) then begin
     add_pv,pv,["ro_amb","en_amb","mx_amb","my_amb","mz_amb",$
                "bx_amb","by_amb","bz_amb"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
     add_pv,pv,["ro_mhd","en_mhd","mx_mhd","my_mhd","mz_mhd",$
                "bx_mhd","by_mhd","bz_mhd"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
  endif  

  if (flag_pip eq 1 or flag_mhd eq 0) then begin
     add_pv,pv,["ro_n","en_n","mx_n","my_n","mz_n"] ,0,var_names,$
            mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath    
  endif  
  if (flag_mhd eq 1) then begin
     add_pv,pv,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],$
            1,var_names,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath
  endif
;===============================================-
; Non-ideal terms
;===============================================
  if info.flag_resi ge 1 or info.flag_artvis eq 2 then begin
     flag_resi=1
     catch,error
     if error ne 0 then begin
        print,"Error loading resistivty"
        flag_resi=0
        catch,/cancel
     endif 
     if flag_resi eq 1 then begin
;        ts=1
        fl_resi=info.flag_resi
        fl_resi=fl_resi mod 10
        if fl_resi eq 1 then begin
           t_resi=1
        endif else begin
           t_resi=ts
        endelse
        et=dblarr(ix,jx,kx,t_resi)
        files=file_search(datapath+"et.dac."+ $
                          string(indgen(n_cpu),form="(i4.4)"))
        mpi_read,et,files,mpi_x,mpi_y,margin,ix_m,jx_m,kx_m,time_step=time_step
        pv=create_struct(pv,["et"], reform(et))
     endif
  endif

  if (info.flag_ir mod 10) ge 1  then begin
     flag_ir=1
;     print, 'hi'
;        files=file_search(datapath+"/"+string(time_step,form="(i4.4)")+"ion.dac.*")
;     stop
;     print,files[0]
;     catch,error
;     if error ne 0 then begin
;        print,"Error loading IR "
;        PRINT, 'Error message: ', !ERROR_STATE.MSG
;        flag_ir=0
;        catch,/cancel
;     endif 
     if flag_ir eq 1 then begin
        ir_flag=info.flag_ir mod 10
        if ir_flag eq 1 then begin
           ion=dblarr(ix,jx,kx,n_read)
           rec=dblarr(ix,jx,kx,n_read)
        endif else begin
           ion=dblarr(ix,jx,kx,n_read)
           rec=dblarr(ix,jx,kx,n_read)
        endelse                   
print, 'ION'
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"ion.dac.*")
;+ $
;                          string(indgen(n_cpu),form="(i4.4)"))

          ;mpi_read2,ion1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,time_step=time_step[np]
	mpi_read,ion1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          ion[*,*,*,np]=ion1
;print,files,ion1,np
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"rec.dac.*")
;+ $
;                          string(indgen(n_cpu),form="(i4.4)"))
          mpi_read,rec1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,time_step=time_step[np]
	mpi_read,ion1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          rec[*,*,*,np]=rec1
        endfor
        pv=create_struct(pv,["ion"], reform(ion))  
        pv=create_struct(pv,["rec"], reform(rec))        
     endif
  endif
  if (info.flag_amb mod 10) ge 1  then begin
     flag_ir=1
     catch,error
     if error ne 0 then begin
        print,"Error loading xin"
        flag_ir=0
        catch,/cancel
     endif 
     if flag_ir eq 1 then begin
        ir_flag=info.flag_amb mod 10
        if ir_flag eq 1 then begin
           xin=dblarr(ix,jx,kx)
        endif else begin
           xin=dblarr(ix,jx,kx,ts)
        endelse                   
        files=file_search(datapath+"xi.dac."+ $
                          string(indgen(n_cpu),form="(i4.4)"))
        mpi_read,xin,files,mpi_x,mpi_y,margin,ix_m,jx_m,kx_m,time_step=time_step
        pv=create_struct(pv,["xi"], reform(xin))                
     endif
  endif

  if info.flag_grav ge 1 then begin
;     dacget3s,datapath+'gr.dac.0000',gr
;     pv=create_struct(pv,["gr"], reform(gr))
  endif

  if flag_az eq 1 then begin
     dx = fltarr(ix)
     dy = fltarr(jx)
     for i=0,ix-2 do begin
        dx[i] = pv.x[i+1]-pv.x[i]
     endfor
     dx[ix-1] = dx[ix-2]
     for j=0,jx-2 do begin
        dy[j] = pv.y[j+1]-pv.y[j]
     endfor
     dy[jx-2] = dy[jx-2]
     cal_az2d,pv.bx,pv.by,dx,dy,az,itype=1,dirct=1,margin=margin[0] $
              ,filename=datapath+'az.sav'
     pv = create_struct(pv,["az"],reform(az))
  endif

  if flag_te eq 1 then begin
     if (flag_pip eq 1 or flag_mhd eq 0) then begin
        te = gm*pv.pr_n/pv.ro_n
        pv = create_struct(pv,["te_n"],reform(te))
     endif
     if (flag_mhd eq 1) then begin
        te = gm*pv.pr_p/pv.ro_p
        pv = create_struct(pv,["te_p"],reform(te))
     endif
  endif
  
  pv=create_struct(pv,["var_names"], var_names)
  if current eq 1 then get_cur,pv
  close,/all
end
