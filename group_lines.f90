subroutine group_lines	!	check if line groups contain enough lines
	use common_vars

	integer count, n
	character(1) iii, jjj

	do j = 1, n_masks	!	open files for each group of every mask
		k = index(masks(j)%name,'.',back=.true.)
		do i = 1, nLSD
			write(iii, '(i1)') i
			write(jjj, '(i1)') j
			n = (j-1)*nLSD+i
			open(100+n,file = masks(j)%name(1:k)//'m'//jjj//'gr'//iii)
		enddo
	enddo

	!	for the most strong group check number of lines
	do j = 1, n_masks
		do 	
			count = 0
			do k = 1, masks(j)%n_points	
				if (masks(j)%R(k) > LS_limits(nLSD-1)) count = count + 1
			enddo
		
			print *, nLSD, 'group lines: ', count
			if (count > 10) then	!	if there are more than 10 lines in group then ok, write the line list for this group end exit
				n = (j-1)*nLSD+nLSD
				do k = 1, masks(j)%n_points	
					if (masks(j)%R(k) > LS_limits(nLSD-1)) then	
						write(100+n, '(f10.4,1x,a4,f8.4)') masks(j)%wave(k), masks(j)%el_ID(k), 1.d0-masks(j)%R(k)
					endif
				enddo
				close(100+n)
				exit
			endif
			LS_limits(nLSD-1) = LS_limits(nLSD-1) - 0.01d0	!	else decrease lower limit
			print *, LS_limits(nLSD-1)
		enddo
	enddo


	!	for inner groups
	do  j = 1, n_masks
		do i = nLSD-1, 2, -1
			do			
				count = 0
				do k = 1, masks(j)%n_points
					if (masks(j)%R(k) > LS_limits(i-1) .and. masks(j)%R(k) < LS_limits(i)) count = count + 1
				enddo
				print *, i, ' group lines: ', count
				if (count > 10) then
					n = (j-1)*nLSD+i
					do k = 1, masks(j)%n_points	
						if (masks(j)%R(k) > LS_limits(i-1) .and. masks(j)%R(k) < LS_limits(i)) then	
							write(100+n, '(f10.4,1x,a4,f8.4)') masks(j)%wave(k), masks(j)%el_ID(k), 1.d0-masks(j)%R(k)
						endif
					enddo
					close(100+n)
					exit
				endif
				print *, LS_limits(i-1)
				LS_limits(i-1) = LS_limits(i-1) - 0.01d0
			enddo
		enddo
	enddo

	!	for the most weak group
	do j = 1, n_masks
		do 
			count = 0
			do k = 1, masks(j)%n_points
				if (masks(j)%R(k) < LS_limits(1)) count = count + 1
			enddo
			
			print *, '1 group lines: ', count
			if (count > 10) then
				n = (j-1)*nLSD+1
				do k = 1, masks(j)%n_points	
					if (masks(j)%R(k) < LS_limits(1)) then	
						write(100+n, '(f10.4,1x,a4,f8.4)') masks(j)%wave(k), masks(j)%el_ID(k), 1.d0-masks(j)%R(k)
					endif
				enddo
				close(100+n)
				exit
			endif
			LS_limits(1) = LS_limits(1) + 0.01d0
		enddo
	enddo

	print *, 'Line strength limits for LSDs: ', LS_limits


end subroutine group_lines