module common_vars
	use spectra
	use regions
	use gaussians

	type(spectrum) rcorr, synth, obs, mask
	type(spectrum), allocatable :: obs_spectra(:), masks(:), newmasks(:), masks_buff(:), masks_LSCb(:)
	type(region), allocatable :: regionss(:)
	type(gaussian) gauss

	real(8), parameter :: vc = 2.997925d5, velstep_coef = 0.5

	real(8), allocatable :: model(:), finLSD(:),approxFS(:),approxFA(:)
	real(8), allocatable :: vLSD(:), rLSD(:,:,:), weights(:), lsd_temp(:,:), flux_lsd(:)
	real(8), allocatable :: LS_limits(:), LS_mean(:)
	real(8) :: gauss_pars(3) = (/1.d0, 1.d0, 3.d0/)
	real(8) :: asym_gauss_pars(4) = (/1.d0,1.d0,3.d0,0.d0/)
	real(8) upper_lim, lower_lim, Vrad
	integer :: nRV, n_points, n_spectra, nLSD, n_masks, n_iter
	real(8), allocatable :: xint(:), yint(:)
	real(8), allocatable :: LSD_cent(:,:), contrib(:), flux_temp_lsd(:)
	integer, allocatable :: line_index(:), ilin_right(:), ilin_left(:)
	integer nblends, nblends_imask, index_left, index_right, nmod
	real(8) velocity_step, xmini, wl_lin_cent
	integer alloc_flag, reguli, n_main_iters, iiter
	real(8) regul_par1, regul_par2

end module common_vars
