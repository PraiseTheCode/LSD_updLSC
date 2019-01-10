subroutine obs_LSD	
	use common_vars
	external LSD_fun_obs
	real(8), allocatable :: x(:), absc(:)
	real(8) regpar(2)

	n_points = obs%n_points
	
	call initial_guess_gauss	! get velocities limits and first approximation for LSD
	!print *, gauss_pars(1), gauss_pars(2), gauss_pars(3)

	allocate(lsd_temp(nRV,nLSD),flux_lsd(nRV),x(nRV*nLSD*n_masks),absc(nRV*nLSD*n_masks),stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (lst_temp/flux_lsd/x in subroutine obs_LSD)'
	do imask = 1, n_masks
		do ilsd = 1, nLSD
			do iv = 1,nRV
				x(nRV*nLSD*(imask-1)+nRV*(ilsd-1)+iv) = rLSD(imask, ilsd, iv)	! fill in vector of unknowns
				absc(nRV*nLSD*(imask-1)+nRV*(ilsd-1)+iv) = vLSD(iv)
			enddo
		enddo
	enddo

	print *, "Trying with LSD..."

	if (reguli == 1) then		! Tikhonov`s regularization
		regpar(1) = regul_par1
		regpar(2) = regul_par2
		call mini_reg(x,absc,obs%R,weights,regpar,nRV*nLSD*n_masks,n_points,LSD_fun_obs)
	else						! no regularization
		call mini(x,obs%R,weights,nRV*nLSD*n_masks,n_points,LSD_fun_obs) ! call solver
	endif

	if (alloc_flag == 1) deallocate(finLSD)
	allocate(finLSD(n_masks*nRV*nLSD), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (finLSD in subroutine obs_LSD)'
	finLSD = x

	!call approx

	deallocate(lsd_temp, flux_lsd, x, absc, stat=ios)
	if(ios /= 0) stop 'Memory deallocation failed (lst_temp/flux_lsd/x in subroutine obs_LSD)'

end subroutine obs_LSD


subroutine LSD_fun_obs(x,f,Jacobian,n,nwl,mode) ! calculate model spectrum and Jacobian at each iteration
	use common_vars
	real(8) x(n), f(n_points), Jacobian(n_points,n)
	character(4) mode

	call model_spectrum(x,Jacobian,n,mode) 
	call correction		! take into account correction based on synth
	f = model

end subroutine LSD_fun_obs
