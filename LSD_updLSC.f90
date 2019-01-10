program LSD_reg

	use common_vars

	call input	! load masks and observed spectra

	do i = 1, n_spectra		!	loop over the observed spectra
		masks = masks_buff	!	for the next spectrum, return the initial mask because the LSC of the previous spectrum has changed the mask
		obs = obs_spectra(i)
		print *, 'obs ',i,': ',obs%name
		allocate(model(obs%n_points), weights(obs%n_points), stat=ios)
		if(ios /= 0) stop 'Memory allocation failed (model/weights)'
		weights = 1.d0
		alloc_flag = 0
		do j = 1, n_main_iters
			call obs_LSD	!	compute LSD profile of observed spectrum and model spectrum with correction
			call LSC	!	line strengths correction
			alloc_flag = 1
			masks = newmasks
			masks_LSCb = newmasks
		enddo
		call output
	enddo

end program
