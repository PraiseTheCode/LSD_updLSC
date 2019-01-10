real(8) function f(x)
use common_vars
implicit real(8) (a-h,o-z)
real(8) temp(nblends), flux_mod(nmod)
!real(8), allocatable :: wl_temp(:), flux_temp_obs(:), flux_temp_mod(:)

do i = 1, nblends		! loop over the number of lines contributing to the sum
 temp(i) = x*contrib(i)/100.d0		! their line strengths computed according to the percentage contribution
! if(temp(i) > 1.d0) temp(i) = 0.9999d0
! if(temp(i) <= 0.d0) temp(i) = 1.d-4
enddo

flux_mod = 0.d0; i = 0

do iwl = index_left, index_right		! loop over the number of observed wavelentghs
	i = i + 1
	do imask = 1, n_masks
		do ilin = ilin_left(imask), ilin_right(imask)		! loop over the number of theoretical lines contributing to an observational block
			vv = (obs%wave(iwl) - newmasks(imask)%wave(ilin))/newmasks(imask)%wave(ilin)*vc		! velocity space
			if(vv < vLSD(nRV) .and. vv > vLSD(1)) then		! check if the lines contributes to a given observed wavelength
				do j = 1, nblends 
					if(ilin == line_index(j)) exit
				enddo
				if(j <= nblends) then
					if(nLSD > 1) then	
						call interp_LSD_n(temp(j),imask)
						ios = map1(vLSD,flux_temp_lsd,nRV,vv,Zi,1)		! interpolation in velocity space	
					else
						ios = map1(vLSD,rLSD(imask,1,1:nRV),nRV,vv,Zi,1)		! interpolation in velocity space	 
					endif
					flux_mod(i) = flux_mod(i) + temp(j) * Zi		! for the line we are optimizing
				else
					if(nLSD > 1) then
						call interp_LSD_n(newmasks(imask)%R(ilin),imask)
						ios = map1(vLSD,flux_temp_lsd,nRV,vv,Zi,1)		! interpolation in velocity space	 
					else
						ios = map1(vLSD,rLSD(imask,1,1:nRV),nRV,vv,Zi,1)		! interpolation in velocity space	 
					endif
					flux_mod(i) = flux_mod(i) + newmasks(imask)%R(ilin) * Zi		! all other lines
				endif
			endif
		enddo
	enddo
enddo


allocate (xint(synth%n_points),yint(nmod), stat=ios)		! correction 
if(ios /= 0) stop 'Memory allocation failed (xint/yint in function f)'
xint=rcorr%wave*(1.d0+Vrad/2.997925d5)
id=map1(xint,rcorr%R,synth%n_points,obs%wave(index_left:index_right),yint,nmod)
do i=1,nmod
	dr=(1.d0-yint(i))*(1.d0-flux_mod(i))
	flux_mod(i)=flux_mod(i)-dr
enddo
deallocate (xint,yint)

f = 0.d0; i = 0
do iwl = index_left, index_right
	i = i + 1
	f = f + (obs%R(iwl) - flux_mod(i))**2.d0
enddo

end


subroutine interp_LSD_test(x,imask)
	use common_vars
	implicit real(8)(a-h,o-z)

	flux_temp_lsd = 0.d0

	do iv = 1, nRV
	 print *, iv
	 print *, imask
	 print *, rLSD(imask,1:nLSD,iv)
	 print *, size(flux_temp_lsd)
	 print *, flux_temp_lsd(iv)
	 ios = map1(LS_mean,rLSD(imask,1:nLSD,iv),nLSD,x,flux_temp_lsd(iv),1)
	enddo
	print *, 'exit'
end

subroutine interp_LSD_n(x,imask)
	use common_vars
	implicit real(8)(a-h,o-z)

	flux_temp_lsd = 0.d0

	do iv = 1, nRV
	 ios = map1(LS_mean,rLSD(imask,1:nLSD,iv),nLSD,x,flux_temp_lsd(iv),1)
	enddo
end


REAL*8 FUNCTION GOLDEN(AX,BX,CX,TOL,XMIN) ! Golden search algorithm
!Given a function F, and given a bracketing triplet of abscissas 
!AX, BX, CX (such that BX is between AX and CX, and F(BX) is less 
!than both F(AX) and F(CX)), this routine performs a golden section 
!search for the minimum, isolating it to a fractional precision of 
!about TOL. The abscissa of the minimum is returned as XMIN, and the minimum
!function value is returned as GOLDEN, the returned function value.
PARAMETER(R=.61803399,C=1.-R)  !Golden ratios
IMPLICIT REAL(8) (A-H,O-Z)
  if(AX > 1.d0) AX = 1.d0; if(AX < 0.d0) AX = 1.d-10
  if(CX > 1.d0) CX = 1.d0; if(CX < 0.d0) CX = 1.d-10
  X0=AX  !At any given time we will keep trace of 4 points: 
  X3=CX  !X0,X1,X2,X3. 
  IF(ABS(CX-BX).GT.ABS(BX-AX)) THEN
    X1=BX; X2=BX+C*(CX-BX)
  ELSE
    X2=BX; X1=BX-C*(BX-AX)
  ENDIF
  F1=F(X1); F2=F(X2)  !Initial function evaluations
1 IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2))) THEN
    IF(F2.LT.F1) THEN
	  X0=X1; X1=X2
	  X2=R*X1+C*X3
	  F0=F1; F1=F2
      F2=F(X2)
    ELSE
	  X3=X2; X2=X1
	  X1=R*X2+C*X0
	  F3=F2; F2=F1
      F1=F(X1)
    ENDIF
	GOTO 1
  ENDIF
  IF(F1.LT.F2) THEN
    GOLDEN=F1
	XMIN=X1
  ELSE
    GOLDEN=F2
	XMIN=X2
  ENDIF
  RETURN
END