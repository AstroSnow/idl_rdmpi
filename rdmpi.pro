pro rdmpi,pv,datapath=datapath,current=current,flag_double=flag_double,$
           time_step=time_step,flag_az=flag_az,flag_te=flag_te,var=var,$
            h5read=h5read
  Compile_Opt DEFINT32  
  if(n_elements(datapath)eq 0 ) then datapath="tmp/" else datapath=datapath+"/"
  if(n_elements(current) eq 0) then current=0
  if(n_elements(flag_double) eq 0) then flag_double=1
  if(n_elements(time_step) eq 0) then time_step=-1
  if(n_elements(flag_az) eq 0) then flag_az = 0 ;; vector potential (for 2D)
  if(n_elements(flag_te) eq 0) then flag_te = 0 ;; temperature
  if(n_elements(var) eq 0) then var = [] ;; selected variables
  if(n_elements(h5read) eq 0) then h5read=0

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

;Snow
if (n_elements(var) ne 0) then nvar=n_elements(var)
  
  pv=create_struct(pv,["eqs","fl_mhd","fl_pip","fl_afr"], $
                   eqs,flag_mhd,flag_pip,flag_afr)

  if h5read eq 0 then begin  
	tfile=file_search(datapath+"t.dac.*")
	dacget0s,tfile,t,narg=time_step
  endif
  if(eqs eq "AFR") then begin
     gridfile=file_search(datapath+"region.dac.*")
     dacget2s,gridfile,grid,narg=time_step
     pv=create_struct(pv,["grid"],grid)
  endif
  ix=info.ix
  jx=info.jx
  kx=info.kx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Load the grid

;Using dac files
if (h5read eq 0) then begin
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Get time 
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

;Snow
if (n_elements(var) eq 0) then begin
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

endif else begin
     add_pv_var,pv,var ,0,var_names,$
            mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m,$
            time_step,datapath   
;help,pv,/str
endelse

;===============================================-
; Non-ideal terms
;===============================================
if (n_elements(var) eq 0) then begin
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

;restore the alpha coefficient
  if (info.flag_pip ge 1) then begin
	ac=dblarr(ix,jx,kx,n_read)
print, 'ION'
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"ac.dac.*")
	mpi_read,ac1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          ac[*,*,*,np]=ac1
        endfor
        pv=create_struct(pv,["ac"], reform(ac))
  endif

  if (((info.flag_ir mod 10) ge 1) and (info.flag_pip ge 1))  then begin
     flag_ir=1
     if flag_ir eq 1 then begin
        ir_flag=info.flag_ir mod 10
        if ir_flag eq 1 then begin
           ion=dblarr(ix,jx,kx,n_read)
           rec=dblarr(ix,jx,kx,n_read)
        endif else begin
           ion=dblarr(ix,jx,kx,n_read)
           rec=dblarr(ix,jx,kx,n_read)
        endelse                   
;print, 'ION'
	pv=create_struct(pv,["flag_ir"], info.flag_ir)
	pv=create_struct(pv,["flag_ir_type"], info.flag_ir_type) 
	pv=create_struct(pv,["T0"], info.T_norm) 
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"ion.dac.*")
	mpi_read,ion1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          ion[*,*,*,np]=ion1
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"rec.dac.*")
	mpi_read,rec1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          rec[*,*,*,np]=rec1
        endfor
        pv=create_struct(pv,["ion"], reform(ion))  
        pv=create_struct(pv,["rec"], reform(rec))        
     endif
     if (info.flag_ir eq 4)then begin
        Nexcite=dblarr(ix,jx,kx,6,n_read)
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"nexcite1.dac.*")
	mpi_read,nexcite1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          Nexcite[*,*,*,0,np]=(nexcite1)
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"nexcite2.dac.*")
	mpi_read,nexcite1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          Nexcite[*,*,*,1,np]=(nexcite1)
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"nexcite3.dac.*")
	mpi_read,nexcite1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          Nexcite[*,*,*,2,np]=(nexcite1)
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"nexcite4.dac.*")
	mpi_read,nexcite1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          Nexcite[*,*,*,3,np]=(nexcite1)
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"nexcite5.dac.*")
	mpi_read,nexcite1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          Nexcite[*,*,*,4,np]=(nexcite1)
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"nexcite6.dac.*")
	mpi_read,nexcite1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          Nexcite[*,*,*,5,np]=(nexcite1)
        endfor
        pv=create_struct(pv,["nexcite"], reform(Nexcite)) 
     endif
     if (info.flag_ir_type eq 0) then begin
	aheat=dblarr(ix,jx,kx,n_read)
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"aheat.dac.*")
	mpi_read,aheat1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          aheat[*,*,*,np]=aheat1
	endfor
        pv=create_struct(pv,["aheat"], reform(aheat)) 
     endif
  endif

;     if info.flag_visc eq 1 then begin                
;print, 'visc'
;visc=dblarr(ix,jx,kx,3,n_read)
;        for np=0,n_read-1 do begin
;        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"viscx.dac.*")
;	mpi_read,visc1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
;          visc[*,*,*,0,np]=visc1

;        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"viscy.dac.*")
;	mpi_read,visc1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
;          visc[*,*,*,1,np]=visc1


;        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"viscz.dac.*")
;	mpi_read,visc1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
;          visc[*,*,*,2,np]=visc1
;        endfor
;        pv=create_struct(pv,["visc"], reform(visc))        
;     endif


;restore the reference energy cooling term
  if (info.flag_rad ge 1) then begin
	edref=dblarr(ix,jx,kx,n_read)
print, 'Rad_cooling'
        for np=0,n_read-1 do begin
        files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"edref_m.dac.*")
	mpi_read,ac1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
          edref[*,*,*,np]=ac1
        endfor
        pv=create_struct(pv,["edref_p"], reform(edref))
        if (info.flag_pip eq 1) then begin
            for np=0,n_read-1 do begin
            files=file_search(datapath+"/"+string(time_step[np],form="(i4.4)")+"edref_h.dac.*")
	    mpi_read,ac1,files,mpi_x,mpi_y,mpi_z,margin,ix_m,jx_m,kx_m
              edref[*,*,*,np]=ac1
            endfor
            pv=create_struct(pv,["edref_n"], reform(edref))
        endif
    radrho=info.radrhoref
    radt=info.rad_ts
   pv=create_struct(pv,["radt"], reform(radt))
    pv=create_struct(pv,["radrho"], reform(radrho)) 
    rlos_m=(pv.pr_p/(gm-1.0)+(pv.ro_p*pv.vx_p^2+pv.ro_p*pv.vy_p^2+pv.ro_p*pv.vz_p^2)/2.0+$
        (pv.bx^2+pv.by^2+pv.bz^2)/2.0)-pv.edref_p
    rlos_m=rlos_m/(pv.ro_p/pv.radrho(0))^(-1.7)/pv.radt(0)
    pv=create_struct(pv,["rlos_p"], reform(rlos_m)) 
    if (info.flag_pip eq 1) then begin
        rlos_h=pv.pr_n/(gm-1.0)+(pv.ro_n*pv.vx_n^2+pv.ro_n*pv.vy_n^2+pv.ro_n*pv.vz_n^2)/2.0$
            -pv.edref_n
        rlos_h=rlos_h/(pv.ro_n/pv.radrho(0))^(-1.7)/pv.radt(0)
        pv=create_struct(pv,["rlos_n"], reform(rlos_h)) 
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
  
;  pv=create_struct(pv,["var_names"], var_names)
;  if current eq 1 then get_cur,pv

  if current eq 1 then begin
	print,'CURRENT ONLY WORKS IN 2D AT PRESENT'
	j_x=dblarr(ix,jx,kx,n_read)
	j_y=dblarr(ix,jx,kx,n_read)
	j_z=dblarr(ix,jx,kx,n_read)
	dx=abs(pv.x[1]-pv.x[0]) ; grid point size
	dy=abs(pv.y[1]-pv.y[0])
        for np=0,n_read-1 do begin
		j_x[2:(ix-3),2:(jx-3),np]=(-pv.Bz[2:(ix-3),4:(jx-1),np]+8.*pv.Bz[2:(ix-3),3:(jx-2),np]-8.*pv.Bz[2:(ix-3),1:(jx-4),np]+pv.Bz[2:(ix-3),0:(jx-5),np])/(12.*dy)

		j_y[2:(ix-3),2:(jx-3),np]=-(-pv.Bz[4:(ix-1),2:(jx-3),np]+8.*pv.Bz[3:(ix-2),2:(jx-3),np]-8.*pv.Bz[1:(ix-4),2:(jx-3),np]+pv.Bz[0:(ix-5),2:(jx-3),np])/(12.*dx)

		j_z[2:(ix-3),2:(jx-3),np]=(-pv.By[4:(ix-1),2:(jx-3),np]+8.*pv.By[3:(ix-2),2:(jx-3),np]-8.*pv.By[1:(ix-4),2:(jx-3),np]+pv.By[0:(ix-5),2:(jx-3),np])/(12.*dx) - $
		(-pv.Bx[2:(ix-3),4:(jx-1),np]+8.*pv.Bx[2:(ix-3),3:(jx-2),np]-8.*pv.Bx[2:(ix-3),1:(jx-4),np]+pv.Bx[2:(ix-3),0:(jx-5),np])/(12.*dy)
	endfor
        pv = create_struct(pv,["j_x"],reform(j_x))
        pv = create_struct(pv,["j_y"],reform(j_y))
        pv = create_struct(pv,["j_z"],reform(j_z))
  endif

endif

  close,/all

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endif else begin
    RESOLVE_ROUTINE,'h5get'
;Read grid grom h5 files
    ndim=pv.info.ndim
    
    if N_elements(time_step) eq 1 then begin
        fpath=datapath+"/t"+string(string(time_step),form="(i4.4)")+'.h5'
        print,'READING ONE TIME STEP'   
        readvars=['xgrid']
        if ndim ge 2 then readvars=[readvars,'ygrid']
        if ndim ge 3 then readvars=[readvars,'zgrid']
        h5get,pv,fpath,readvars,1

        ;Read conserved variables
        if (n_elements(var) eq 0) then begin
          if (flag_pip eq 1 or flag_mhd eq 0) then begin
    ;         h5get,pv,fpath,["ro_n","en_n","mx_n","my_n","mz_n"],0
             h5get,pv,fpath,["ro_n","pr_n","vx_n","vy_n","vz_n"],0
          endif  
          if (flag_mhd eq 1) then begin
    ;         h5get,pv,fpath,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],0
             h5get,pv,fpath,["ro_p","pr_p","vx_p","vy_p","vz_p","bx","by","bz"],0
          endif

        endif else begin
             h5get,pv,fpath,var,0
        ;help,pv,/str
        endelse

    ;   Non-ideal terms
        if n_elements(var) eq 0 then begin
            if (info.flag_rad ge 1) then begin
            edref=dblarr(ix,jx,kx,n_read)
            print, 'Rad_cooling Loaded'
            h5get,pv,fpath,["edref_m"],1
            radrho=info.radrhoref
            radt=info.rad_ts
            pv=create_struct(pv,["radt"], reform(radt))
            pv=create_struct(pv,["radrho"], reform(radrho)) 
            endif   
	    if (((info.flag_ir mod 10) ge 1) and (info.flag_pip ge 1)) then begin
	     
		    if (info.flag_ir_type eq 0) then begin
		    print, 'Losses loading'
		    h5get,pv,fpath,["aheat"],1
		    h5get,pv,fpath,["ion_loss"],1
		    endif  
		    if (info.flag_ir eq 4) then begin
		    print, 'Hydrogen levels loading'
		    h5get,pv,fpath,["nexcite1","nexcite2","nexcite3","nexcite4","nexcite5","nexcite6"],1
		    endif        
	    endif
        endif
    endif else begin
    
        print,'Reading times'
        for rt = 0,N_elements(time_step)-1 do begin
            print,time_step(rt)
            ;ndim=pv.info.ndim
            fpath=datapath+"/t"+string(string(time_step(rt)),form="(i4.4)")+'.h5'
            
            if rt eq 0 then begin   
                readvars=['xgrid']
                if ndim ge 2 then readvars=[readvars,'ygrid']
                if ndim ge 3 then readvars=[readvars,'zgrid']
                h5get,pv,fpath,readvars,1

                ;Read conserved variables
                if (n_elements(var) eq 0) then begin
                  if (flag_pip eq 1 or flag_mhd eq 0) then begin
            ;         h5get,pv,fpath,["ro_n","en_n","mx_n","my_n","mz_n"],0
                     h5get,pv,fpath,["ro_n","pr_n","vx_n","vy_n","vz_n"],0
                  endif  
                  if (flag_mhd eq 1) then begin
            ;         h5get,pv,fpath,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],0
                     h5get,pv,fpath,["ro_p","pr_p","vx_p","vy_p","vz_p","bx","by","bz"],0
                     pv_full=create_struct("info",0)
                     n_ygrid=1
                     n_zgrid=1
                     n_xgrid=n_elements(pv.xgrid)
                     if ndim ge 2 then n_ygrid=n_elements(pv.ygrid)
                     if ndim ge 3 then n_zgrid=n_elements(pv.zgrid)
                     tmp_ro_p=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_ro_p(*,*,*,0)=pv.ro_p
                     ;tmp_ro_p=reform(tmp_ro_p)
                     tmp_vx_p=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_vy_p=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_vz_p=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_pr_p=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_bx=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_by=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                     tmp_bz=dblarr(n_xgrid,n_ygrid,n_zgrid,n_read)
                  endif

                endif else begin
                     h5get,pv,fpath,var,0
                ;help,pv,/str
                endelse
                
                
            endif else begin
                pv=0
                pv=create_struct("info",0)
                readvars=['xgrid']
                if ndim ge 2 then readvars=[readvars,'ygrid']
                if ndim ge 3 then readvars=[readvars,'zgrid']
                h5get,pv,fpath,readvars,1
            
                ;Read conserved variables
                if (n_elements(var) eq 0) then begin
                  if (flag_pip eq 1 or flag_mhd eq 0) then begin
            ;         h5get,pv,fpath,["ro_n","en_n","mx_n","my_n","mz_n"],0
                     h5get,pv,fpath,["ro_n","pr_n","vx_n","vy_n","vz_n"],0
                  endif  
                  if (flag_mhd eq 1) then begin
            ;         h5get,pv,fpath,["ro_p","en_p","mx_p","my_p","mz_p","bx","by","bz"],0
                     h5get,pv,fpath,["ro_p","pr_p","vx_p","vy_p","vz_p","bx","by","bz"],0
                     if ndim eq 1 then tmp_ro_p(*,0,0,rt)=pv.ro_p
                     if ndim eq 2 then tmp_ro_p(*,*,0,rt)=pv.ro_p
                     if ndim eq 3 then tmp_ro_p(*,*,*,rt)=pv.ro_p
                     
                     if ndim eq 1 then tmp_vx_p(*,0,0,rt)=pv.vx_p
                     if ndim eq 2 then tmp_vx_p(*,*,0,rt)=pv.vx_p
                     if ndim eq 3 then tmp_vx_p(*,*,*,rt)=pv.vx_p
                     
                     if ndim eq 1 then tmp_vy_p(*,0,0,rt)=pv.vy_p
                     if ndim eq 2 then tmp_vy_p(*,*,0,rt)=pv.vy_p
                     if ndim eq 3 then tmp_vy_p(*,*,*,rt)=pv.vy_p
                     
                     if ndim eq 1 then tmp_vz_p(*,0,0,rt)=pv.vz_p
                     if ndim eq 2 then tmp_vz_p(*,*,0,rt)=pv.vz_p
                     if ndim eq 3 then tmp_vz_p(*,*,*,rt)=pv.vz_p
                     
                     if ndim eq 1 then tmp_pr_p(*,0,0,rt)=pv.pr_p
                     if ndim eq 2 then tmp_pr_p(*,*,0,rt)=pv.pr_p
                     if ndim eq 3 then tmp_pr_p(*,*,*,rt)=pv.pr_p
                     
                     if ndim eq 1 then tmp_bx(*,0,0,rt)=pv.bx
                     if ndim eq 2 then tmp_bx(*,*,0,rt)=pv.bx
                     if ndim eq 3 then tmp_bx(*,*,*,rt)=pv.bx
                     
                     if ndim eq 1 then tmp_by(*,0,0,rt)=pv.by
                     if ndim eq 2 then tmp_by(*,*,0,rt)=pv.by
                     if ndim eq 3 then tmp_by(*,*,*,rt)=pv.by
                     
                     if ndim eq 1 then tmp_bz(*,0,0,rt)=pv.bz
                     if ndim eq 2 then tmp_bz(*,*,0,rt)=pv.bz
                     if ndim eq 3 then tmp_bz(*,*,*,rt)=pv.bz
                     
                     
                  endif

                endif else begin
                     h5get,pv,fpath,var,0
                ;help,pv,/str
                endelse
            
            endelse
        endfor
        tmp=pv
        pv=create_struct("info",0)
        pv=create_struct(pv,'xgrid',tmp.xgrid)
        if ndim ge 2 then pv=create_struct(pv,'ygrid',tmp.ygrid)
        if ndim ge 3 then pv=create_struct(pv,'zgrid',tmp.zgrid)
        pv=create_struct(pv,'ro_p',reform(tmp_ro_p))
        pv=create_struct(pv,'vx_p',reform(tmp_vx_p))
        pv=create_struct(pv,'vy_p',reform(tmp_vy_p))
        pv=create_struct(pv,'vz_p',reform(tmp_vz_p))
        pv=create_struct(pv,'pr_p',reform(tmp_pr_p))
        pv=create_struct(pv,'bx',reform(tmp_bx))
        pv=create_struct(pv,'by',reform(tmp_by))
        pv=create_struct(pv,'bz',reform(tmp_bz))
    endelse

endelse
end
