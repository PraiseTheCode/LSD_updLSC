subroutine model_spectrum(x,Jacobian,n,mode) ! calculate model spectra on a grid of observed wavelengths
use common_vars
real(8) Jacobian(n_points,n), vp, prof, p, dv
character(4) mode
real(8) x(n)

 model = 0.d0 ! initialize model spectrum array

 if(mode == 'grad') Jacobian = 0.d0 ! initialize Jacobian if needed

 !dv = vLSD(2) - vLSD(1) ! step width in RV
 
 do imask = 1, n_masks	! loop over the masks
	mask = masks(imask)

	 do i = 1, mask%n_points ! loop over the line mask entries

	  do k = 1, n_points ! loop over the observed wavelength
	   vp = (obs%wave(k) - mask%wave(i))/mask%wave(i)*vc ! observed wavelength (with respect to the wavelength of a given line in mask) in RV space
	   
	   if(vp < vLSD(1)) cycle	! check if particular entry from the line mask contributes to the observed wavelength
	   if(vp > vLSD(nRV)) cycle
	
	   if(nLSD > 1) then
		call interp_LSD(imask,x, n, mask%R(i), kk)	! interpolate LSD profile for a given line strength
		id = map1(vLSD,flux_lsd,nRV,vp,prof,1)	! interpolate LSD profile for a given velocity point
	   else
		id = map1(vLSD,x,nRV,vp,prof,1) 
	   endif

	   model(k) = model(k) + prof*mask%R(i) ! model spectrum through a convolution of the LSD profile with the line mask

	   if(mode == 'grad') then     ! calculate Jacobian as function was linearly interpolated

		ki = 0
		do l = 1, nRV
		 if(vp <= vLSD(l)) then
		  ki = l; exit
		 endif
		enddo
		
		if(ki == 1) ki = 2; if(ki == 0) ki = nRV
		dv = vLSD(ki) - vLSD(ki-1)
		p = (vp - vLSD(ki-1))/dv

		if(nLSD > 1) then
			Jacobian(k,nRV*(imask-1)*nLSD+nRV*(kk-1)+ki-1) = Jacobian(k,nRV*(imask-1)*nLSD+nRV*(kk-1)+ki-1) + mask%R(i) * (1.d0 - p)
			Jacobian(k,nRV*(imask-1)*nLSD+nRV*(kk-1)+ki)   = Jacobian(k,nRV*(imask-1)*nLSD+nRV*(kk-1)+ki) + mask%R(i) * p
		else
			Jacobian(k,ki-1) = Jacobian(k,ki-1) + mask%R(i)*(1.d0-p)
			Jacobian(k,ki) = Jacobian(k,ki) + mask%R(i)*p
		endif

	   endif

	  enddo

	 enddo
	enddo

end


subroutine interp_LSD(icomp,x,nx,flux_line,kk)	! interpolate LSD profile for a given line strength
	use common_vars
	implicit real(8) (a-h,o-z)
	real(8) x(nx)

	lsd_temp = 0.d0; flux_lsd = 0.d0
	do iv = 1, nRV
		do ilsd = 1, nLSD
			lsd_temp(iv,ilsd) = x((nRV*nLSD*(icomp-1)+nRV*(ilsd-1))+iv)
		enddo
	enddo

	do iv = 1, nRV
		ios = map1(LS_mean,lsd_temp(iv,1:nLSD),nLSD,flux_line,flux_lsd(iv),1)
	enddo

	if(flux_line <= LS_limits(1)) then
		kk = 1
	elseif(flux_line > LS_limits(nLSD-1)) then
		kk = nLSD
	else
		do ilimit = 1, nLSD-2
			if(flux_line > LS_limits(ilimit) .and. flux_line <= LS_limits(ilimit+1)) then
				kk = ilimit + 1
				exit
			endif
		enddo
	endif

end subroutine interp_LSD