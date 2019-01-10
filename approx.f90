subroutine approx
	use common_vars
	external sym_fun
	external asym_fun
	real(8) weis(nRV), YY(nRV)
	real(8) sym_gauss_pars(3)
	weis = 1.d0
	
	if (alloc_flag == 0) allocate(approxFS(nRV),approxFA(nRV))
	!gauss_pars(2) = Vrad
	YY = rLSD(1,1,1:nRV)
	sym_gauss_pars(1) = 1.d0
	sym_gauss_pars(2) = gauss_pars(2)
	sym_gauss_pars(3) = 3.d0
	call mini(sym_gauss_pars, YY, weis, 3, nRV, sym_fun)
	print *, gauss_pars
	print *, '!!!'
	!print *, approxFS

	!open(666,file = 'RZPsc_weak_sym.gau')
	!do i = 1, nRV
	!	write(666,*) vLSD(i), 1.d0-approxFS(i)
	!enddo
	!close(666)

	
	asym_gauss_pars(1) = 1.d0
	asym_gauss_pars(2) = gauss_pars(2)
	asym_gauss_pars(3) = 3.d0
	asym_gauss_pars(4) = 0.1d0
	call mini(asym_gauss_pars, YY, weis, 4, nRV, asym_fun)

	
	k = index(obs%name,'.',back=.true.)
	open(666,file = obs%name(1:k)//'approx1')
	write(666,'(f9.4,"    ",e11.5)') asym_gauss_pars(2), asym_gauss_pars(4)
	do i = 1, nRV
		write(666,*) vLSD(i), 1.d0-approxFS(i), 1.d0-approxFA(i)
	enddo
	close(666)

	!YY = rLSD(1,2,1:nRV)
	!sym_gauss_pars(1) = 1.d0
	!sym_gauss_pars(2) = gauss_pars(2)
	!sym_gauss_pars(3) = 3.d0
	!call mini(sym_gauss_pars, YY, weis, 3, nRV, sym_fun)
	!!print *, gauss_pars
	!!print *, '!!!'
	!!print *, approxFS

	!!open(666,file = 'RZPsc_weak_sym.gau')
	!!do i = 1, nRV
	!!	write(666,*) vLSD(i), 1.d0-approxFS(i)
	!!enddo
	!!close(666)

	
	!asym_gauss_pars(1) = 1.d0
	!asym_gauss_pars(2) = gauss_pars(2)
	!asym_gauss_pars(3) = 3.d0
	!asym_gauss_pars(4) = 0.1d0
	!call mini(asym_gauss_pars, YY, weis, 4, nRV, asym_fun)

	
	!k = index(obs%name,'.',back=.true.)
	!open(666,file = obs%name(1:k)//'approx2')
	!write(666,'(f9.4,"    ",e11.5)') asym_gauss_pars(2), asym_gauss_pars(4)
	!do i = 1, nRV
	!	write(666,*) vLSD(i), 1.d0-approxFS(i), 1.d0-approxFA(i)
	!enddo
	!close(666)

	!YY = rLSD(1,3,1:nRV)
	!sym_gauss_pars(1) = 1.d0
	!sym_gauss_pars(2) = gauss_pars(2)
	!sym_gauss_pars(3) = 3.d0
	!call mini(sym_gauss_pars, YY, weis, 3, nRV, sym_fun)
	!print *, gauss_pars
	!print *, '!!!'
	!!print *, approxFS

	!!open(666,file = 'RZPsc_weak_sym.gau')
	!!do i = 1, nRV
	!!	write(666,*) vLSD(i), 1.d0-approxFS(i)
	!!enddo
	!!close(666)

	
	!asym_gauss_pars(1) = 1.d0
	!asym_gauss_pars(2) = gauss_pars(2)
	!asym_gauss_pars(3) = 3.d0
	!asym_gauss_pars(4) = 0.1d0
	!call mini(asym_gauss_pars, YY, weis, 4, nRV, asym_fun)

	
	!k = index(obs%name,'.',back=.true.)
	!open(666,file = obs%name(1:k)//'approx3')
	!write(666,'(f9.4,"    ",e11.5)') asym_gauss_pars(2), asym_gauss_pars(4)
	!do i = 1, nRV
	!	write(666,*) vLSD(i), 1.d0-approxFS(i), 1.d0-approxFA(i)
	!enddo
	!close(666)


end subroutine approx


subroutine sym_fun(gx,f,Jacobian,n, nwl, modeg)

	use common_vars

	real(8) gx(3), f(nwl), Jacobian(nwl, 3)
	character(4) modeg

	call function_gauss(nwl, 3, gx)
	if (modeg == 'grad') call jacobian_gauss(nwl, 3, gx, Jacobian,nwl)

	f = approxFS

end subroutine sym_fun



subroutine asym_fun(gx,f,Jacobian,n, nwl, modeg)

	use common_vars

	real(8) gx(4), f(nwl), Jacobian(nwl, 4)
	character(4) modeg

	!print *, 1

	call function_asym_gauss(nwl, 4, gx)
	if (modeg == 'grad') call jacobian_asym_gauss(nwl, 4, gx, Jacobian,nwl)
		
	!print *, asym_gauss_pars

	f = approxFA

end subroutine asym_fun


subroutine function_gauss(nfunc,nunknowns,xi) ! calculate approximation function
 use common_vars
 real(8) xi(nunknowns),sign,s,x(nfunc)
 x = vLSD

 do i=1,nfunc
   approxFS(i)=xi(1)*exp(-((x(i)-xi(2))/(xi(3)))**2)
 enddo

!f=f-y
!s=0.d0
!do i=1,nfunc
! s=s+f(i)**2
!enddo
!s=sqrt(s)/nfunc
!write(4,*)s
end

subroutine jacobian_gauss(nfunc,nunknowns,xi,jacobian,ldfjac) ! calculate Jacobian
 use common_vars
real(8) jacobian(ldfjac,nunknowns),xi(nunknowns),sign,const,const2,const3,x(nfunc)
 x = vLSD

 do i=1,nfunc
	const2=(x(i)-xi(2))/((xi(3)))**2
	const=exp(-const2*(x(i)-xi(2)))
	const3=2.d0*xi(1)*const*const2
	Jacobian(i,1) = const
	Jacobian(i,2) = const3
	Jacobian(i,3) = const3*(x(i)-xi(2))
 enddo
end subroutine jacobian_gauss


subroutine function_asym_gauss(nfunc,nunknowns,xi) ! calculate approximation function
 use common_vars
 real(8) xi(nunknowns),sign,s,x(nfunc)
  
 x = vLSD

 !print *, 7
 do i=1,nfunc
  sign=1.d0
   if(x(i) < xi(2))sign=-1.d0
   approxFA(i)=xi(1)*exp(-((x(i)-xi(2))/(xi(3)*(1.d0+sign*xi(4))))**2)
 enddo

 !print *, 8
 !print *, approxFA
 !pause

end subroutine function_asym_gauss

subroutine jacobian_asym_gauss(nfunc,nunknowns,xi,jacobian,ldfjac) ! calculate Jacobian
 use common_vars
 real(8) jacobian(ldfjac,nunknowns),xi(nunknowns),sign,const,const2,const3,x(nfunc)
 x = vLSD
	!print *, 9

 do i=1,nfunc
  sign=1.d0
   if(x(i) < xi(2))sign=-1.d0
  const2=(x(i)-xi(2))/((xi(3)*(1.d0+sign*xi(4))))**2
  const=exp(-const2*(x(i)-xi(2)))
  const3=2.d0*xi(1)*const*const2
  jacobian(i,1)=const
  jacobian(i,2)=const3
  jacobian(i,3)=const3*(x(i)-xi(2))/xi(3)
  jacobian(i,4)=sign*const3*((x(i)-xi(2))/(1.d0+sign*xi(4)))
 enddo

	!print *, 10
end subroutine jacobian_asym_gauss