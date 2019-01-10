subroutine input

	use common_vars
	use dflib
	character(100) config_name
	character(100) rcorr_name, mask_ext, synth_name, obs_ext
	character(100), allocatable :: obs_names(:), mask_names(:)
	real(8) lower, upper, Res
	integer number_regions
	integer(kind=int_ptr_kind()) handle
	integer(4) ress
	type (FILE$INFO) info

	config_name = "LSD_updLSC.conf"
	open(1,file=trim(adjustl(config_name)),status='old',iostat=ios)
	if(ios /= 0) stop 'Configuration file does not exist'

	read(1,'((a))') obs_ext; k = 0; k = index(obs_ext,'!',back=.true.); if(k /= 0) obs_ext(k:) = ' '
	read(1, *) Res	!	spectral resolution
	read(1,'((a))') mask_ext; k = 0; k = index(mask_ext,'!',back=.true.); if(k /= 0) mask_ext(k:) = ' '
	read(1,'((a))') synth_name; k = 0; k = index(synth_name,'!',back=.true.); if(k /= 0) synth_name(k:) = ' '
	read(1,'((a))') rcorr_name; k = 0; k = index(rcorr_name,'!',back=.true.); if(k /= 0) rcorr_name(k:) = ' '

	read(1, *) nLSD
	if(nLSD > 1) then
		allocate(LS_limits(nLSD-1), LS_mean(nLSD))
		read(1,*) (LS_limits(i), i = 1, nLSD-1)
	else
		read(1,*)
	endif

	read(1,*) number_regions
	allocate(regionss(number_regions), stat = ios)
	if(ios /= 0) stop 'Memory allocation failed (regions)'
	do i = 1,number_regions
		read(1,*) lower, upper
		regionss(i) = new_region(lower, upper)
	enddo

	read(1,*) n_iter	! number of iterations of the line strength correction

	read(1,*) regul_par1, regul_par2
	if (regul_par1 == 0 .and. regul_par2 == 0) then
		reguli = 0
	else if (regul_par1 == 0) then
		stop 'Wrong regularization parameters'
	else 
		reguli = 1
	endif

	read(1,*) n_main_iters	! number of iterations of the LSD-LSC cycle

	close(1)
	
	k = index(obs_ext,'!'); if(k /= 0) obs_ext = obs_ext(1:k-1) 

	handle = FILE$FIRST; n_spectra = 0
	do
		ress = getfileinfoqq(trim(adjustl(obs_ext)),info,handle)
		if(handle == FILE$LAST) exit
		n_spectra = n_spectra + 1 
	enddo
	print *, 'N obs spec: ', n_spectra

	allocate(obs_names(n_spectra))
	allocate(obs_spectra(n_spectra))

	handle = FILE$FIRST; i = 0
	do
		ress = getfileinfoqq(trim(adjustl(obs_ext)),info,handle)
		if(handle == FILE$LAST) exit
		i = i + 1; obs_names(i) = trim(adjustl(info%name)) 
		obs_spectra(i) = read_spectrum(obs_names(i), number_regions, regionss, "spectrum")
		obs_spectra(i)%Res = Res
	enddo

	k = index(mask_ext,'!'); if(k /= 0) mask_ext = mask_ext(1:k-1) 

	handle = FILE$FIRST; n_masks = 0
	do
		ress = getfileinfoqq(trim(adjustl(mask_ext)),info,handle)
		if(handle == FILE$LAST) exit
		n_masks = n_masks + 1 
	enddo
	print *, 'N masks: ', n_masks

	allocate(mask_names(n_masks))
	allocate(masks(n_masks), newmasks(n_masks), masks_buff(n_masks), masks_LSCb(n_masks))

	handle = FILE$FIRST; i = 0
	do
		ress = getfileinfoqq(trim(adjustl(mask_ext)),info,handle)
		if(handle == FILE$LAST) exit
		i = i + 1; mask_names(i) = trim(adjustl(info%name)) 
		masks(i) = read_spectrum(mask_names(i), number_regions, regionss, "mask")
		newmasks(i) = read_spectrum(mask_names(i), number_regions, regionss, "mask")
		masks_buff(i) = read_spectrum(mask_names(i), number_regions, regionss, "mask")
		masks_LSCb(i) = read_spectrum(mask_names(i), number_regions, regionss, "mask")
	enddo

	rcorr = read_spectrum(rcorr_name, number_regions, regionss, "spectrum")

	synth = read_spectrum(synth_name, number_regions, regionss, "spectrum")

	
	k = index(masks(1)%name,'.',back=.true.)
	open(100,file = masks(1)%name(1:k)//'before')
	do i=1,masks(1)%n_points
		write(100, '(f10.4,1x,a4,f8.4)') masks(1)%wave(i), masks(1)%el_ID(i), 1.d0-masks(1)%R(i)
	enddo  
	close(100) 

	if(nLSD > 1) then
		call group_lines	! check if there are enough lines in groups and move line strength limits if not

		do ilsd = 1, nLSD
			if(ilsd == 1) then
				LS_mean(ilsd) = (0.d0 + LS_limits(ilsd))*0.5d0
			elseif(ilsd == nLSD) then
				LS_mean(ilsd) = (1.d0 + LS_limits(ilsd-1))*0.5d0
			else
				LS_mean(ilsd) = (LS_limits(ilsd-1) + LS_limits(ilsd))*0.5d0
			endif
		enddo
	endif

	print *, "Output file produced: ", masks(1)%name(1:k)//'before'

end subroutine
