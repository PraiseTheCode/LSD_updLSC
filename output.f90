subroutine output
	use common_vars

	k = index(obs%name,'.',back=.true.)
	open(100,file = obs%name(1:k)//'oc')
	write(100, *) gauss_pars(2)
	do i=1,obs%n_points
		write(100, *) obs%wave(i), 1.d0-obs%R(i), 1.d0-model(i)
	enddo 
	close(100) 

	print *, "Output file produced: ", obs%name(1:k)//'oc'

	k = index(obs%name,'.',back=.true.)
	open(100,file = obs%name(1:k)//'lsd')
	do iv=1,nRV
		write(100, '(1x,<nLSD*n_masks+1>f10.5)') vLSD(iv), ((1.d0-finLSD(nRV*(imask-1)*nLSD+nRV*(i-1)+iv), i = 1, nLSD), imask = 1, n_masks)
	enddo  
	close(100) 

	print *, "Output file produced: ", obs%name(1:k)//'lsd'

	k = index(masks(1)%name,'.',back=.true.)
	open(100,file = masks(1)%name(1:k)//'after')
	do i=1,newmasks(1)%n_points
		write(100, '(f10.4,1x,a4,f8.4)') newmasks(1)%wave(i), newmasks(1)%el_ID(i), 1.d0-newmasks(1)%R(i)
	enddo  
	close(100) 

	print *, "Output file produced: ", masks(1)%name(1:k)//'after'

	deallocate(model, weights, vLSD, rLSD, finLSD, ilin_left, ilin_right, flux_temp_lsd, LSD_cent)

	gauss_pars(1) = 1.d0
	gauss_pars(2) = 1.d0
	gauss_pars(3) = 3.d0


end subroutine
