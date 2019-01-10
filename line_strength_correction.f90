subroutine LSC
	use common_vars 
	implicit real(8) (a-h,o-z)

	print *, 'Line strength correction...'
	
	do imask = 1, n_masks
		do ilsd = 1, nLSD
			do iv = 1,nRV
				rLSD(imask, ilsd, iv) = finLSD(nRV*nLSD*(imask-1)+nRV*(ilsd-1)+iv) 
			enddo
		enddo
		allocate(newmasks(imask)%line_index_optim(newmasks(imask)%n_points), stat=ios)
		if(ios /= 0) stop 'Memory allocation failed (line_index_optim in subroutine LSC)'
		newmasks(imask)%line_index_optim = 0
	enddo
	if(alloc_flag == 0) allocate(ilin_left(n_masks),ilin_right(n_masks),flux_temp_lsd(nRV), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (ilins/flux_temp_lsd in subroutine LSC)'
	deallocate(flux_temp_lsd)
	allocate(flux_temp_lsd(nRV))
	!velocity_step = vLSD(2) - vLSD(1)

	call LSD_center
	iter = 0; call spectrum_calc(1,iter); tol = 0.01d0
	do iter = 1, n_iter
		iiter = iter
		do imask = 1, n_masks
			if (iter == 1) then
				ilin = 0
				newmasks(imask)%initial_corrs = 0.0d0
				newmasks(imask)%el_corrs = 0.d0
				newmasks(imask)%el_corrs_counters = 0
				do
					ilin = ilin + 1; if(ilin > newmasks(imask)%n_points) exit
					if(imask > 1) then
						if(newmasks(imask)%line_index_optim(ilin) == 1) cycle 
					endif
					wl_lin_cent = (LSD_cent(imask,nLSD) + vc)*newmasks(imask)%wave(ilin)/vc		! wavelength of the line in the observed spectrum
					call define_wave_range(ilin,imask)
					
					if (masks_LSCb(imask)%R(ilin) > 0.2d0) then
						if (nblends == 1) then
							x2 = 0.d0
							do k = 1, nblends
								x2 = x2 + masks_LSCb(imask)%R(ilin)	! for very close lines, sum of their contributions will be fitted. So, compute the sum
							enddo
							x1 = x2/10.d0; x3 = x2*10.d0		! the optimization ranges
							call contribution(x2,nblends,1,imask)		! compute individual contributions in per cent
							ff = golden(x1,x2,x3,tol,xmini)		! optimize the sum of individual contribution
							call contribution(xmini,nblends,2,imask)		! compute individual line strengths by keeping percentage contributions
							corr_b = (newmasks(imask)%R(ilin)-masks_LSCb(imask)%R(ilin))/masks_LSCb(imask)%R(ilin)
							if (abs(corr_b) < 2.d0) then
								do iel = 1, size(newmasks(imask)%elements)
									if (newmasks(imask)%el_ID(ilin) == newmasks(imask)%elements(iel)) then
										newmasks(imask)%el_corrs(iel) = newmasks(imask)%el_corrs(iel) + corr_b
										newmasks(imask)%el_corrs_counters(iel) = newmasks(imask)%el_corrs_counters(iel) + 1
										exit
									endif
								enddo
							endif
						endif
					endif
					deallocate(line_index)
				enddo
				do iel = 1, size(newmasks(imask)%elements)
					if (newmasks(imask)%el_corrs_counters(iel) /= 0) newmasks(imask)%el_corrs(iel) = newmasks(imask)%el_corrs(iel) / newmasks(imask)%el_corrs_counters(iel) 
					!print *, newmasks(imask)%el_corrs(iel)
				enddo
				do ilin = 1, newmasks(imask)%n_points
					do iel = 1, size(newmasks(imask)%elements)
						if (newmasks(imask)%el_ID(ilin) == newmasks(imask)%elements(iel)) then
							newmasks(imask)%initial_corrs(ilin) = newmasks(imask)%el_corrs(iel)
						endif
					enddo
				enddo
				
				ilin = 0
				do
					ilin = ilin + 1; if(ilin > newmasks(imask)%n_points) exit
					if(imask > 1) then
						if(newmasks(imask)%line_index_optim(ilin) == 1) cycle 
					endif
					wl_lin_cent = (LSD_cent(imask,nLSD) + vc)*newmasks(imask)%wave(ilin)/vc		! wavelength of the line in the observed spectrum
					call define_wave_range(ilin,imask)

					if (nblends /= 1) then
						newmasks(imask)%R(ilin) = masks_LSCb(imask)%R(ilin)*(1+newmasks(imask)%initial_corrs(ilin))
					endif

					deallocate(line_index)
				enddo
			endif

			newmasks(imask)%line_index_optim = 0
			ilin = 0
			do 
				ilin = ilin + 1; if(ilin > newmasks(imask)%n_points) exit
				if(imask > 1) then
					if(newmasks(imask)%line_index_optim(ilin) == 1) cycle 
				endif
				wl_lin_cent = (LSD_cent(imask,nLSD) + vc)*newmasks(imask)%wave(ilin)/vc		! wavelength of the line in the observed spectrum
				call define_wave_range(ilin,imask)
				x2 = 0.d0; ilin = ilin + nblends_imask - 1
				do k = 1, nblends
					x2 = x2 + newmasks(imask)%R(line_index(k))	! for very close lines, sum of their contributions will be fitted. So, compute the sum
				enddo
				x1 = x2/10.d0; x3 = x2*10.d0		! the optimization ranges
				call contribution(x2,nblends,1,imask)		! compute individual contributions in per cent
				ff = golden(x1,x2,x3,tol,xmini)		! optimize the sum of individual contribution
				call contribution(xmini,nblends,2,imask)		! compute individual line strengths by keeping percentage contributions
				deallocate(line_index)
			enddo
			
			masks_LSCb(imask) = newmasks(imask)
		enddo
		if(iter == n_iter) then
			call spectrum_calc(2,iter)
		else
			call spectrum_calc(1,iter)
		endif

	enddo

end subroutine LSC


subroutine LSD_center		! compute center of gravity using five points around the center of the LSD-profile
use common_vars
	implicit real(8) (a-h,o-z)
	integer imax(1)

	if (.not.allocated(LSD_cent)) allocate(LSD_cent(n_masks,nLSD), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (LSD_cent in subroutine LSD_center)'
	do imask = 1, n_masks
		do ilsd = 1, nLSD
			imax = maxloc(rLSD(imask,ilsd,1:nRV))

			svi = 0.d0; sv = 0.d0

			do j = imax(1)-6, imax(1)+6, 1
				svi = svi + vLSD(j)*rLSD(imask,ilsd,j)
				sv  = sv + rLSD(imask,ilsd,j)
			enddo
			LSD_cent(imask,ilsd) = svi/sv		! center of gravity
		enddo
	enddo

end subroutine LSD_center


subroutine contribution(x2,j,switch,imask)
use common_vars
implicit real(8) (a-h,o-z)
integer switch, ifobs
real(8), allocatable :: arR(:)

if(switch == 1) then
 allocate(contrib(j))
 iiii = 0
 do i = 1, j
  contrib(i) = newmasks(imask)%R(line_index(i))*100.d0/x2		! define percentage contribution
 enddo
else
 do i = 1, j
  newmasks(imask)%R(line_index(i)) = x2*contrib(i)/100.d0		! compute individual line strengths by keeping percentage contributions
 enddo
 deallocate(contrib)
endif

end subroutine contribution


subroutine define_wave_range(ilin,imask)
use common_vars
implicit real(8) (a-h,o-z)

index_left = 0; index_right = 0; ilin_left = 1000000; ilin_right = -1; nmod = 0
wl_lin_left = (vLSD(2) + vc)*newmasks(imask)%wave(ilin)/vc		! left boundary
wl_lin_right = (vLSD(nRV-1) + vc)*newmasks(imask)%wave(ilin)/vc		! right boundary

do iwl = 1, obs%n_points-1
 if(obs%wave(iwl) <= wl_lin_left .and. obs%wave(iwl+1) > wl_lin_left) index_left = iwl
 if(obs%wave(iwl) <= wl_lin_right .and. obs%wave(iwl+1) > wl_lin_right) index_right = iwl
 if(index_right /= 0) exit
enddo
if(index_left == 0) index_left = 1; if(index_right == 0) index_right = obs%n_points; nmod = index_right - index_left + 1
do iwl = index_left, index_right
 do j = 1, n_masks
  do i = 1, newmasks(j)%n_points
   vv = (obs%wave(iwl) - newmasks(j)%wave(i))/newmasks(j)%wave(i)*vc		! velocity space
   if(vv < vLSD(nRV-1) .and. vv > vLSD(1)) then		! check if the lines contributes to a given observed wavelength
    if(i < ilin_left(j)) ilin_left(j) = i
    if(i > ilin_right(j)) ilin_right(j) = i
   endif
  enddo
 enddo
enddo
do k = 1, 2
 nblends_imask = 1; nblends = 1; if(k == 2) line_index(nblends) = ilin
 do i = ilin_left(imask), ilin_right(imask)
  if(i <= ilin) cycle
  deltav = (newmasks(imask)%wave(i) - newmasks(imask)%wave(ilin))/(newmasks(imask)%wave(i) + newmasks(imask)%wave(ilin))*2.d0*vc		! check if the next to the close neihbour another close line exists
  if(abs(deltav) <= velocity_step) then
   nblends = nblends + 1 
   if(k == 2) then
    line_index(nblends) = i
	nblends_imask = nblends
   endif
  endif
 enddo
 if(imask /= n_masks) then
  do i = imask+1, n_masks
   do ii = ilin_left(i), ilin_right(i)
    wl_lin_cent1 = newmasks(i)%wave(ii)*LSD_cent(i,nLSD)/vc; wl_lin_cent1 = newmasks(i)%wave(ii) + wl_lin_cent1
    deltav = (wl_lin_cent1 - wl_lin_cent)/(wl_lin_cent1 + wl_lin_cent)*2.d0*vc
    if(abs(deltav) <= velocity_step) then
     nblends = nblends + 1 
	 if(k == 2) then
	  line_index(nblends) = ii
	  newmasks(i)%line_index_optim(ii) = 1
	 endif
    endif
   enddo
  enddo
 endif
 if(k == 1) then
	allocate(line_index(nblends), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (line_index in subroutine define_wave_range)'
 endif
enddo

end subroutine define_wave_range


subroutine spectrum_calc(icheck,iter)		! compute model spectrum
use common_vars
implicit real(8) (a-h,o-z)

do iwl = 1, obs%n_points
 model(iwl) = 0.d0
 do imask = 1, n_masks
  do ilin = 1, newmasks(imask)%n_points
   vv = (obs%wave(iwl) - newmasks(imask)%wave(ilin))/newmasks(imask)%wave(ilin)*vc
   if(vv < vLSD(nRV) .and. vv > vLSD(1)) then
    if(nLSD > 1) then
	 !if (alloc_flag == 1) call interp_LSD_test(newmasks(imask)%R(ilin),imask)
     call interp_LSD_n(newmasks(imask)%R(ilin),imask)
     ios = map1(vLSD,flux_temp_lsd,nRV,vv,Zi,1)		! interpolation in velocity space
	else
     ios = map1(vLSD,rLSD(imask,1,1:nRV),nRV,vv,Zi,1)		! interpolation in velocity space
    endif
    model(iwl) = model(iwl) + newmasks(imask)%R(ilin) * Zi
   endif
  enddo
 enddo
enddo

call correction
call compute_st_deviation(iter)

end subroutine spectrum_calc


subroutine compute_st_deviation(iter)
use common_vars
implicit real(8) (a-h,o-z)

dev = 0.d0
do iwl = 1, obs%n_points
  dev = dev + (obs%R(iwl) - model(iwl))**2.0		! sum of the O-C squared
enddo
dev = dev/real(obs%n_points,8); dev = sqrt(dev)		! 1-sigma level
if(iter == 0) then
 write(*,"('Initial model evaluation: Standart dev = ', f8.6,' (',f6.3,'%)')") dev, dev*100.d0
else
 write(*,"('Iteration ',i2,': Standart dev = ', f8.6,' (',f6.3,'%)')") iter, dev, dev*100.d0
endif

end subroutine compute_st_deviation