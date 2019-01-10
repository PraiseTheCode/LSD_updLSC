module regions

	type region	!	"class" for specifying working ranges of spectra
		real(8) left, right
	end type

	contains

	function new_region(left, right) !	constructor
		type(region) new_region
		real(8), intent(in) :: left, right

		new_region%left = left
		new_region%right = right
	end function

	function point_inside(point, reg) ! check if given point enters the working range
		logical point_inside
		type(region), intent(in) :: reg
		real(8), intent(in) :: point

		if (point >= reg%left .and. point <= reg%right) then
			point_inside = .true.
		else 
			point_inside = .false.
		endif
	end function

end module regions


module spectra 
	use regions

	type spectrum !	derived type that contains information about spectra/masks
		real(8), allocatable :: wave(:), R(:)
		real(8), allocatable :: line_index_optim(:), initial_corrs(:), el_corrs(:) ! arrays for LSC
		character(4), allocatable :: el_ID(:), elements(:), els(:)
		integer, allocatable :: el_num(:), el_corrs_counters(:)
		real(8) wavelength_step, Res
		type(region), dimension(:), allocatable :: regionss
		character(100) name
		integer n_regions, n_points
	end type

	contains

	function read_spectrum(name, n_reg, regionss, mode)	!	initialize type fields - read spectrum or mask and save the parts included in the working ranges
		type(spectrum) read_spectrum
		character(100), intent(in) :: name
		integer, intent(in) :: n_reg
		type(region), dimension(:), intent(in) :: regionss
		character(8), intent(in) :: mode
		character(8) :: ss = "spectrum"
		character(4) :: mm = "mask"
		real(8) ww, rr
		character(4) el
		integer counter, n_points, i_count_eq, i_count_filled

		read_spectrum%name = name
		read_spectrum%n_regions = n_reg
		allocate(read_spectrum%regionss(n_reg), stat=ios)
		if(ios /= 0) stop 'Memory allocation failed (regions in module spectra)'
		do i = 1, n_reg
			read_spectrum%regionss(i) = regionss(i)
		enddo

		open(10,file=trim(adjustl(name)),status='old',iostat=ios)
		if(ios /= 0) then
			write(*,"('File ',a,' does not exist')") trim(adjustl(name))
			stop
		endif

		n_points = 0 
		if (mode == ss) then
			do          
				read(10,*,iostat=ios) ww
				if (ios /= 0) exit
				do i = 1, n_reg
					if (point_inside(ww, read_spectrum%regionss(i))) then
						n_points = n_points + 1
					endif
				enddo
			enddo  
		else if (mode(1:4) == mm) then
			do          
				read(10,'(f10.4,1x,a4,f8.4)',iostat=ios) ww, el, rr
				if (ios /= 0) exit
				do i = 1, n_reg
					if (point_inside(ww, read_spectrum%regionss(i)).and.rr<0.99d0) then
						n_points = n_points + 1
					endif
				enddo
			enddo
		endif  
		read_spectrum%n_points = n_points

		allocate(read_spectrum%wave(n_points),read_spectrum%R(n_points), stat=ios)
		if(ios /= 0) stop 'Memory allocation failed (wave/R in module spectra)'
		rewind(10)

		counter = 0
		if (mode == ss) then
			do while (counter < n_points)
				read(10,*,iostat=ios) ww, rr
				do i = 1, n_reg
					if (point_inside(ww, read_spectrum%regionss(i))) then
						counter = counter + 1
						read_spectrum%wave(counter) = ww; read_spectrum%R(counter) = 1.d0 - rr
					endif
				enddo
			enddo
		else if (mode(1:4) == mm) then
			allocate(read_spectrum%el_ID(n_points), stat = ios)
			if(ios /= 0) stop 'Memory allocation failed (el_ID in module spectra)'
			do while (counter < n_points)
				read(10,'(f10.4,1x,a4,f8.4)',iostat=ios) ww, el, rr
				do i = 1, n_reg
					if (point_inside(ww, read_spectrum%regionss(i)).and.rr<0.99d0) then
						counter = counter + 1
						read_spectrum%wave(counter) = ww; read_spectrum%R(counter) = 1.d0 - rr; read_spectrum%el_ID(counter) = el
					endif
				enddo
			enddo
			allocate(read_spectrum%initial_corrs(n_points))
			allocate(read_spectrum%els(1000))
			read_spectrum%els = "null"
			i_count_filled = 0
			do i = 1, n_points
				i_count_eq = 0
				do j = 1, size(read_spectrum%els)
					if (read_spectrum%el_ID(i) == read_spectrum%els(j)) then
						i_count_eq = 1
						exit
					endif
				enddo
				if (i_count_eq == 0) then
					i_count_filled = i_count_filled + 1
					read_spectrum%els(i_count_filled) = read_spectrum%el_ID(i)
				endif
			enddo
			allocate(read_spectrum%elements(i_count_filled))
			read_spectrum%elements = read_spectrum%els(1:i_count_filled)
			deallocate(read_spectrum%els)
			allocate(read_spectrum%el_corrs(i_count_filled))
			allocate(read_spectrum%el_corrs_counters(i_count_filled))
			allocate(read_spectrum%el_num(i_count_filled))
			do j = 1, i_count_filled
				i_count_eq = 0
				do i = 1, n_points
					if (read_spectrum%el_ID(i) == read_spectrum%elements(j)) i_count_eq = i_count_eq + 1
				enddo
				read_spectrum%el_num(j) = i_count_eq
			enddo

		else 
			print *, "Wrong mode in read_spectrum: ", mode
			stop
		endif

		close(10)
		read_spectrum%wavelength_step = read_spectrum%wave(2) - read_spectrum%wave(1)

	end function

end module spectra


module gaussians

	type gaussian !	"class" for specifying initial gaussian with function for redetermination velocity grid limits during the iterative fitting with Gaussians
		real(8) gx(3)
		real(8) dV
		real(8), allocatable :: V(:), R(:)
		integer nV, lower_index, upper_index
	end type 

	contains

	function new_gaussian(gx, nV, dV) 
		type(gaussian) new_gaussian
		real(8), intent(in) :: gx(3), dV
		integer, intent(in) :: nV
		real(8) max
		integer count
	
		new_gaussian%dV = dV
		!new_gaussian%dV = 3.d0

		allocate(new_gaussian%V(nV+1), stat=ios)
		if(ios /= 0) stop 'Memory allocation failed (V in module gaussians)'
		new_gaussian%V(nV/2 + 1) = 0
		do i=1,nV/2
			new_gaussian%V(nV/2+1-i) = - i*new_gaussian%dV
			new_gaussian%V(nV/2+1+i) = i*new_gaussian%dV
		enddo

		allocate(new_gaussian%R(nV+1), stat=ios)
		if(ios /= 0) stop 'Memory allocation failed (R in module gaussians)'
		count = 0
		do i = 1, nV+1
			new_gaussian%R(i) = gx(1) * exp(-(new_gaussian%V(i)-gx(2))**2/gx(3)**2)
		enddo
		new_gaussian%nV = nV

		max = maxval(new_gaussian%R)
		do i = 1, new_gaussian%nV
			if (new_gaussian%R(i) == max) count_l = i
		enddo
		count_u = count_l	
		do while (new_gaussian%R(count_l) > 1.d-5 .or. new_gaussian%R(count_u) > 1.d-5)
			new_gaussian%upper_index = count_u
			new_gaussian%lower_index = count_l
			count_u = count_u + 1
			count_l = count_l - 1
		enddo

	end function

	function get_lims(gaus, gx)	!	find new velocity grid limits for updated Gaussian
		real(8) get_lims(2)
		type(gaussian) gau
		type(gaussian), intent(in) :: gaus
		real(8), intent(in) :: gx(3)
		integer count_l, count_u
		real(8) max
		gau = gaus
		
		max = maxval(gau%R)
		do i = 1, gau%nV
			if (gau%R(i) == max) count_l = i
		enddo
		count_u = count_l
		gau%R(count_l) = gx(1) * exp(-(gau%V(count_l)-gx(2))**2/gx(3)**2)
		do while (gau%R(count_l) > 1.d-5 .or. gau%R(count_u) > 1.d-5)
			gau%R(count_u+1) = gx(1) * exp(-(gau%V(count_u+1)-gx(2))**2/gx(3)**2)
			gau%R(count_l-1) = gx(1) * exp(-(gau%V(count_l-1)-gx(2))**2/gx(3)**2)
			count_u = count_u + 1
			count_l = count_l - 1
		enddo
		get_lims(1) = gau%V(count_l)
		get_lims(2) = gau%V(count_u)
	end function

end module gaussians

