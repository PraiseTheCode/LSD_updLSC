subroutine initial_guess_gauss
	use common_vars
	external gauss_fun
	real(8) lims(2), gx(3), x0, x1,lin,func, x2
	real(8), allocatable :: oldV(:), newV(:)
	integer flag

	call get_Vrad
	print *, 'initial Vrad = ', Vrad
	gauss_pars(2) = Vrad
	velocity_step = vc/obs%res*velstep_coef	!	initial velocity step determined by spectral resolution
	print *, 'dV = ', velocity_step
	gauss = new_gaussian(gauss_pars, int(400.d0/velocity_step), velocity_step)	!	initialize gaussian fields
	!lims = get_lims(gauss, gauss_pars)	!	find limits of velocity grid
	!lower_lim = lims(1)
	!upper_lim = lims(2)

	print *, "Trying with gaussian..."
	call mini(gauss_pars, obs%R, weights, 3, n_points, gauss_fun)	!	minimize residuals of model and synth
	
	print *, 'Vrad = ', gauss_pars(2)
	gauss = new_gaussian(gauss_pars, int(400.d0/velocity_step), velocity_step)	!	update gaussian fields
	nRV = size(gauss%V(gauss%lower_index:gauss%upper_index)) + 1
	lims = get_lims(gauss, gauss_pars)	!	find limits of velocity grid
	lower_lim = lims(1)
	upper_lim = lims(2)

					
	allocate(oldV(nRV), newV(10*nRV))	!	specify the velocity grid with variable step
	oldV = gauss%V(gauss%lower_index-1:gauss%upper_index) 
	newV = 0
	i = 3
	new = 0
	x2 = oldV(3)
	x1 = oldV(2)
	x0 = oldV(1)
	gx = gauss_pars
	flag = 0
	do	!	loop over the segments with initial step
		if (x2>oldV(nRV)) exit
		lin =  G(x0,gx) + ( G(x1,gx) - G(x0,gx) ) / (x1-x0) * ((x1+x0)/2 - x0)
		func =  G((x0+x1)/2 , gx)
		if (abs(lin-func) < 0.001) then	!	if there is no significant difference between line approximation and real value then go to the next segmetn
			new = new + 1
			newV(new) = x1
			x0 = x1
			x1 = x2
			x2 = x1 + velocity_step
			!print *, x0, x1, x2
			flag = 0
		else	!	if there is significant difference then divide the segment in half
			x2 = x1
			x1 = (x0+x1)/2	
		endif
	enddo

	do while (newV(new) < oldV(nRV))
		new = new + 1
		newV(new) = newV(new-1)+velocity_step
	enddo

	nRV = new

	if (alloc_flag == 1) deallocate(vLSD, rLSD)
	allocate(vLSD(nRV), rLSD(n_masks,nLSD,nRV), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (vLSD/rLSD in subroutine initial_guess_gauss)'
	!vLSD = gauss%V(gauss%lower_index:gauss%upper_index) ! set velocity grid and first approximation for LSD profile
	vLSD = newV(1:new)
	do imask = 1, n_masks
		do i = 1, nLSD	
			do iv = 1, nRV
				rLSD(imask, i, iv) = G(vLSD(iv),gx)
			enddo
		enddo
	enddo

	deallocate(oldV, newV)

end subroutine initial_guess_gauss


subroutine gauss_fun(gx,f,Jacobian,n, nwl, modeg)	! calculate model spectrum and Jacobian at each iteration
	use common_vars

	real(8) gx(3), f(n_points), Jacobian(n_points, 3)
	character(4) modeg
	real(8) lims(2)
	
	lims = get_lims(gauss, gx)	!	update limits of velocity grid
	lower_lim = lims(1)
	upper_lim = lims(2)
	!print *, lower_lim, upper_lim

	call model_gauss(Jacobian, gx, modeg)
	f = model

end subroutine gauss_fun


subroutine model_gauss(Jacobian, gx, modeg)
	use common_vars
	real(8) vp, prof 
	real(8) gx(3), Jacobian(n_points, 3),const,const2,const3
	character(4) modeg
	if (modeg == "grad") Jacobian = 0.d0
	model = 0.d0

	do imask = 1, n_masks
		mask = masks(imask)
		do i = 1, mask%n_points 
			do k = 1, n_points 

				vp = (obs%wave(k) - mask%wave(i))/mask%wave(i)*vc

				if(vp < lower_lim) cycle
				if(vp > upper_lim) cycle

				model(k)=model(k)+gx(1)*exp(-((vp-gx(2))/gx(3))**2)*mask%R(i)

				if (modeg == "grad") then
					const2=(vp-gx(2))/((gx(3)))**2
					const=exp(-const2*(vp-gx(2)))
					const3=2.d0*gx(1)*const*const2
					Jacobian(k,1) = Jacobian(k,1)+ const
					Jacobian(k,2) = Jacobian(k,2)+ const3
					Jacobian(k,3) = Jacobian(k,3)+ const3*(vp-gx(2))
				endif

			enddo
		enddo
	enddo

end subroutine model_gauss

function G(x, gx)
	real(8) G, gx(3), x

	G = gx(1) * exp(-(x-gx(2))**2/gx(3)**2)
end function


subroutine get_Vrad
	use common_vars
	real(8) CCF1, CCF2, Rv1, Rv2, i
	external CCF_max
	!tolerant=1.d-5; ! Code should correct radial velocity by search maximum CCF in region Rv_old +- dv

	!Rv1=fmin(-200.d0,0.d0,CCF_max,tolerant,CCF1, 1)
	!Rv2=fmin(0.d0,200.d0,CCF_max,tolerant,CCF2, 1)
	!Vrad=Rv1
	!if (CCF1 > CCF2) Vrad=Rv2

	CCF1 = CCF_max(0.d0, 1)
	Vrad = 0.d0
	do i = -200.d0, 200.d0, 1.d0
		CCF2 = CCF_max(i, 1)
		if (CCF1 > CCF2) then
			Vrad=i
			CCF1 = CCF2
		endif
	enddo

end subroutine get_Vrad


real(8) function CCF_max(v, n)
	use common_vars
	implicit real(8) (a-h,o-z)

	allocate(xint(obs%n_points), yint(obs%n_points), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (xint/yint in subroutine correction)'
	xint = obs%wave*(1-v/vc)
	id=map1(synth%wave,synth%R,synth%n_points,xint,yint,obs%n_points)

	sx=sum(obs%R**2); sy=sum(yint**2); sxy=sum(obs%R*yint)

	if(sx <= 0.or.sy <= 0.or.sxy <= 0)then
		CCF_max=100
	else  
		CCF_max=sqrt(sx)*sqrt(sy)/sxy ! really 1/CCF - fmin search minimum of function, so CCF_max=1/CCF
	endif
	deallocate (xint,yint)
end