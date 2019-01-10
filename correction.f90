subroutine correction ! correct model spectrum corresponding differencies between LSD synthetic model and synthetic spectrum
	use common_vars

	allocate (xint(synth%n_points),yint(obs%n_points), stat=ios)
	if(ios /= 0) stop 'Memory allocation failed (xint/yint in subroutine correction)'

	xint=rcorr%wave*(1.d0+Vrad/2.997925d5)
	id=map1(xint,rcorr%R,synth%n_points,obs%wave,yint,obs%n_points)
	do i=1,obs%n_points
		dr=(1.d0-yint(i))*(1.d0-model(i))
		model(i)=model(i)-dr
	enddo
	deallocate (xint,yint)
end
