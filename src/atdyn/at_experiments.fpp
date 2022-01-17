!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_experiments_mod
!> @brief   experimental data
!! @authors Osamu Miyashita (OM), Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
  
module at_experiments_mod

  use molecules_str_mod
  use at_restraints_str_mod
  use at_experiments_str_mod
  use at_enefunc_str_mod
  use fileio_sit_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use string_mod
  use timers_mod
  use fileio_pdb_mod
#ifdef MPI
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_exp_info
    logical                         :: emfit           = .false.
    character(MaxFilename)          :: emfit_target    = ''
    real(wp)                        :: emfit_sigma     = 2.0_wp
    real(wp)                        :: emfit_tolerance = 0.001_wp
    integer                         :: emfit_period    = 1
  end type s_exp_info

  type(s_experiments), target, save :: experiments
  real(wp),       allocatable, save :: dot_exp_drhodx(:), dot_exp_drhody(:), dot_exp_drhodz(:)
  real(wp),       allocatable, save :: erfa(:,:), expa(:,:)
  real(wp),       allocatable, save :: derfa_x(:,:), derfa_y(:,:), derfa_z(:,:)
  real(wp),       allocatable, save :: dexpa_x(:,:), dexpa_y(:,:), dexpa_z(:,:)
  real(wp),       allocatable, save :: dot_sim_drhodx(:), dot_sim_drhody(:), dot_sim_drhodz(:)
  real(wp),       allocatable, save :: emfit_force(:,:)

  real(wp),       allocatable, save :: exp_image(:,:) ! input exp image
  real(wp),       allocatable, save :: sim_image(:,:) ! input exp image

  real(wp),       			   save :: emfit_roll_ang, emfit_tilt_ang, emfit_yaw_ang,emfit_shift_x, emfit_shift_y
  character(MaxFilename),      save :: emfit_exp_image,emfit_target_pdb
  character(MaxFilename),      save :: emfit_mode, emfit_image_type, emfit_im_gen_mode
  integer,       			   save :: image_size, emfit_image_size
  logical, 					   save :: image_out
  real(wp),       			   save :: emfit_pixel_size, pixel_size
  logical, 					   save :: optimizer_mode, EMfit_image_mode
  real*4, dimension(256),	   save :: file_header
  real(wp),       			   save :: gradient_optimizer, image_gen_optimizer 
  character(MaxFilename),      save :: MPI_mode
  !define for image emfit
  real(wp),       allocatable, save :: I_sim(:,:), I_exp(:,:)
  ! subroutines
  public  	:: show_ctrl_experiments
  public  	:: read_ctrl_experiments
  public  	:: setup_experiments
  private 	:: setup_experiments_emfit
  public  	:: compute_energy_experimental_restraint
  private 	:: compute_energy_experimental_restraint_emfit_with_image
  public  	:: Euler_angle2Matrix
  public  	:: Euler_angle2Matrix_Inv
  public	:: writeSPIDER
  public 	:: readSPIDER
  public	:: readSPIDERheader
  public	:: creatSpiderHeader
contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_experiments
  !> @brief        show EXPERIMENTS section usage
  !! @authors      TM
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "remd", "min", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_experiments(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md', 'remd')

        write(MsgOut,'(A)') '[EXPERIMENTS]'
        write(MsgOut,'(A)') 'emfit           = NO          # EM fit'
        write(MsgOut,'(A)') '# emfit_target    = sample.sit  # EM data file'
        write(MsgOut,'(A)') '# emfit_sigma     = 2.0         # resolution parameter of the simulated map'
        write(MsgOut,'(A)') '# emfit_tolerance = 0.001       # Tolerance for error'
        write(MsgOut,'(A)') '# emfit_period    = 1           # emfit force update period'
        write(MsgOut,'(A)') ''

      end select

    end if

    return

  end subroutine show_ctrl_experiments
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      read_ctrl_experiments
  !> @brief        read EXPERIMENTS section in the control file
  !! @authors      TM
  !! @param[in]    handle   : unit number of control files
  !! @param[out]   exp_info : EXPERIMENTS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_experiments(handle, exp_info) 

    ! parameters
    character(*),            parameter     :: Section = 'Experiments'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_exp_info),        intent(inout) :: exp_info 

    ! local variables


    ! read parameters from control file
    ! 
	emfit_exp_image    = ''
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_logical(handle, Section, 'emfit',           &
                               exp_info%emfit)

    call read_ctrlfile_real   (handle, Section, 'emfit_sigma',     &
                              exp_info%emfit_sigma)

    call read_ctrlfile_real   (handle, Section, 'emfit_tolerance', &
                              exp_info%emfit_tolerance)

    call read_ctrlfile_integer(handle, Section, 'emfit_period',    &
                              exp_info%emfit_period)

    call read_ctrlfile_string (handle, Section, 'emfit_exp_image',    &
                              emfit_exp_image)

	call read_ctrlfile_real   (handle, Section, 'emfit_roll_angle', &
                              emfit_roll_ang)

	call read_ctrlfile_real   (handle, Section, 'emfit_tilt_angle', &
                              emfit_tilt_ang)

	call read_ctrlfile_real   (handle, Section, 'emfit_yaw_angle',	&
                              emfit_yaw_ang)

	call read_ctrlfile_real   (handle, Section, 'emfit_shift_x',	&
								emfit_shift_x)

	call read_ctrlfile_real   (handle, Section, 'emfit_shift_y',	&
								emfit_shift_y)

    call read_ctrlfile_integer(handle, Section, 'emfit_image_size',		&
                              emfit_image_size)

	call read_ctrlfile_real(handle, Section, 'emfit_pixel_size',		&
				emfit_pixel_size)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then

      if (exp_info%emfit) then

        write(MsgOut,'(a)') 'Read_Ctrl_Experiments > Parameters for experimental data fitting'
        write(MsgOut,'(a20,F10.4,a20,F10.4)')                      &
              '  emfit_sigma     = ', exp_info%emfit_sigma,        &
              '  emfit_tolerance = ', exp_info%emfit_tolerance
        write(MsgOut,'(a20,I10)')                                  &
              '  emfit_period    = ', exp_info%emfit_period
		write(MsgOut,'(a20,F10.4)')                      		   &
              '  emfit_roll_angle     = ', emfit_roll_ang 
		write(MsgOut,'(a20,F10.4)')                      		   &
              '  emfit_tilt_angle     = ', emfit_tilt_ang 
		write(MsgOut,'(a20,F10.4)')                      		   &
              '  emfit_yaw_angle     = ', emfit_yaw_ang
		write(MsgOut,'(a20,I10)')                                  &
              '  image_size    = ', emfit_image_size       			   
        write(MsgOut,'(a)') ''

      end if

    end if
	
    return

  end subroutine read_ctrl_experiments
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments
  !> @brief        setup experiments information
  !! @authors      TM
  !! @param[in]    exp_info    : EXPERIMENTS section control parameters
  !! @param[inout] experiments : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments(exp_info, molecule, restraints, enefunc)

    ! formal arguments
    type(s_exp_info),        intent(in)    :: exp_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_restraints),      intent(in)    :: restraints
    type(s_enefunc),         intent(inout) :: enefunc

    ! local variables
    integer    :: i
    logical    :: do_emfit_res

    experiments%do_emfit = .false.

    if (exp_info%emfit) then

      do_emfit_res = .false.
      do i = 1, restraints%nfunctions
        if (enefunc%restraint_kind(i) == RestraintsFuncEM) do_emfit_res = .true.
      end do

      if (.not. do_emfit_res) then
        call error_msg('Setup_Experiments> EM is not defined in [RESTRAINTS]')
      end if

      experiments%do_emfit = .true.
      call setup_experiments_emfit(exp_info, molecule)
	
	  EMfit_image_mode = .true.

    end if

    return

  end subroutine setup_experiments

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      setup_experiments_emfit
  !> @brief        setup experiments information
  !! @authors      OM, TM
  !! @param[in]    exp_info    : EXPERIMENTS section control parameter
  !! @param[inout] experiments : experiments information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_experiments_emfit(exp_info, molecule)

    ! formal arguments
    type(s_exp_info),        intent(in)    	:: exp_info
    type(s_molecule),        intent(in)    	:: molecule

    ! local variables
    type(s_sit)                 			:: sit
    type(s_image)               			:: image_exp ! Alex
    integer                     			:: i, j, k, n(3), max_nxyz, num_atoms
    real(wp)                    			:: pi, r, y, sum1
    real(wp)                   				:: A1(4,4)
    real(wp),dimension(:,:),allocatable		:: coord_pdb
    real(wp),dimension(:,:),allocatable		:: tmp_img
    integer                     			:: size_at
    integer                    				:: close_status

    !read image and storge in the experiment Alex
    ! image generation variables
    !integer                     			:: image_size_image
    character(MaxFilename)           		:: filename1_pdb
    type(s_pdb)             				:: pdb1
    integer									:: alloc_stat1
	
    !spider image read and write parameters
    integer 								:: cmd_status
    character(128) 							:: cmd_msg
    character(128) 							:: command0
    logical 								:: file_exists

    !###########***initialization of fitting with 2D image***##########
	! clean all the previous temporary files
    call execute_command_line('rm -f run.*') 	!this line is to remove the result of the previous experiment
 	call execute_command_line('rm -f images/sim/* images/exp/* images/dif/*') 	!image generate by the previous run will be deleted
    !image_size = 128 						! xmipp value for this parameter is 84
	
	!initialization of INP parameters
	image_out		 	= .FALSE.					!TRUE will allow to write experimental and simulated image in a folder
	pixel_size 	=emfit_pixel_size
	image_size = emfit_image_size

	experiments%emfit%sigma        = exp_info%emfit_sigma
    experiments%emfit%tolerance    = exp_info%emfit_tolerance
    experiments%emfit%emfit_period = exp_info%emfit_period


    allocate(exp_image(image_size,image_size))
    allocate(sim_image(image_size,image_size))
    allocate(tmp_img(image_size,image_size))

	call readSPIDER(emfit_exp_image, tmp_img)
	print*, "READING IS OK "

	! print*, "/////////////////////////////////////////////////////////////////////////////////////////"
	! do i=1, image_size
	! 	do j=1, image_size
	! 		exp_image(i,j) = tmp_img(image_size-i +1, image_size-j+1)
	! 	end do
	! end do

	exp_image=tmp_img

	open (unit=66, file="img.txt")
	write(66,*) exp_image
	close(66)


    allocate(emfit_force(3,molecule%num_atoms))
	deallocate(tmp_img)

    return

  end subroutine setup_experiments_emfit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint
  !> @brief        calculate restraint energy from experimental data
  !! @authors      TM
  !! @param[in]    enefunc : potential energy functions information
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint(enefunc, coord, inum, &
                                           calc_force, force, virial, eexp, cv)

    ! formal arguments
    type(s_enefunc), target, intent(inout) :: enefunc
    real(wp),                intent(in)    :: coord(:,:)
    integer,                 intent(in)    :: inum
    logical,                 intent(in)    :: calc_force
    real(wp),                intent(inout) :: force(:,:)
    real(wp),                intent(inout) :: virial(3,3)
    real(wp),                intent(inout) :: eexp
    real(wp),                intent(inout) :: cv
    
    logical                                :: control_para
	if(EMfit_image_mode) then
    	control_para = .TRUE. !.TRUE. is used to perfom MDfit by images
	else
		control_para = .FALSE.
	endif
    if (experiments%do_emfit .and. control_para) then
		call compute_energy_experimental_restraint_emfit_with_image &
             (enefunc, coord, inum, calc_force, force, virial, eexp, cv)
    else
		write(*,*)"this GENESIS package allows to run EMfit with image"
		write(*,*)"set parameter emfit to YES and re-run the program"
		stop
    endif

    return

  end subroutine compute_energy_experimental_restraint

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_energy_experimental_restraint_emfit_with_image
  !> @brief        calculate restraint energy from experimental data, the data includes two set of images
  !! @authors      Rémi Vuillemot
  !! @param[in]    enefunc : potential energy functions information: this potential engergy function is extracted from CHARMM forcefield
  !! @param[in]    coord   : coordinates of target systems
  !! @param[in]    inum    : pointer for restraint function
  !! @param[in]    const   : force constants
  !! @param[in]    ref     : reference value
  !! @param[inout] force   : forces of target systems
  !! @param[inout] virial  : virial of target systems
  !! @param[out]   eexp    : restraint energy
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_energy_experimental_restraint_emfit_with_image(enefunc, coord, inum, &
	calc_force, force, virial, eexp, cv)
	! formal arguments
	type(s_enefunc), target, intent(inout) 	:: enefunc
	real(wp),                intent(in)    	:: coord(:,:)
	integer,                 intent(in)    	:: inum
	logical,                 intent(in)    	:: calc_force
	real(wp),                intent(inout) 	:: force(:,:)
	real(wp),                intent(inout) 	:: virial(3,3)
	real(wp),                intent(inout) 	:: eexp
	real(wp),                intent(inout) 	:: cv

	integer                   :: i, j, a, n_atoms, n_pix, group_id, n, n_atoms_group, cutoff_pixel

	real(wp) :: mu1,mu2, gaussian, sigma, norm, cutoff, threshold 	! Generate Image
	real(wp) :: sum_sim2, sum_exp2, sum_simpexp, force_constant, dcc		! CC computation
	real(wp) :: dpsim1, dpsim2,  dcc1, dcc2, const1, const2					! gradient computation
    real(wp) :: rot_matrix(4,4),  inv_rot_matrix(4,4)
	real(wp),dimension(:,:),allocatable	:: rot_coord
	integer,dimension(:,:),allocatable	:: pixels
	real(wp),dimension(:,:,:),allocatable	:: gaussians_saved

	integer,  pointer ::numatoms(:), atom_id(:,:)

	n_atoms = size(coord,2)
	atom_id        => enefunc%restraint_atomlist
	group_id       =  enefunc%restraint_grouplist(1,inum)
	numatoms       => enefunc%restraint_numatoms
	n_atoms_group = numatoms(group_id)


	! ------------------------------------------------------------------------
	! ALLOCATE
	! ------------------------------------------------------------------------
	allocate(rot_coord(3, n_atoms_group))
	allocate(pixels(2, n_atoms_group))
	allocate(gaussians_saved(image_size, image_size,n_atoms_group ))

	! ------------------------------------------------------------------------
	! ROTATE PDB
	! ------------------------------------------------------------------------
	call Euler_angle2Matrix(emfit_roll_ang, emfit_tilt_ang, emfit_yaw_ang, rot_matrix)
	do a=1, n_atoms_group
		n=  atom_id(a,group_id)
		rot_coord(1:3,a)=matmul(rot_matrix(1:3,1:3),coord(1:3,n))
		rot_coord(1,a) = rot_coord(1,a) - (emfit_shift_x * pixel_size)
		rot_coord(2,a) = rot_coord(2,a) - (emfit_shift_y * pixel_size)
	end do
	inv_rot_matrix = transpose(rot_matrix)

	
	! ------------------------------------------------------------------------
	! SELECT PIXELS TO INTEGRATE
	! ------------------------------------------------------------------------
	sigma = experiments%emfit%sigma

	! distance in Angstrom where the gaussian is truncated
	cutoff =  sqrt((-(2*(sigma**2))) * log(experiments%emfit%tolerance))  
	cutoff_pixel = 	ceiling(cutoff/pixel_size)

	! number of pixels considered per atom
    n_pix = int(cutoff_pixel*2)  											

    do a = 1, n_atoms_group
        pixels(1,a) = 1+nint(rot_coord(1,a)/pixel_size + (image_size/2))-cutoff_pixel ! starting pixel
        pixels(2,a) = 1+nint(rot_coord(2,a)/pixel_size + (image_size/2))-cutoff_pixel

		! Check if the atom is extending outside the box
		if (((pixels(1,a) <= 0) .or. (pixels(2,a) <= 0)) .or. &
		 (((pixels(1,a)+n_pix) > image_size) .or. ((pixels(2,a)+n_pix) > image_size))) then
		 write(MsgOut,'(A)') 'Gaussian kernel is extending outside the map box for atom '
		 write(MsgOut,'(A,I10)') '   pix1 : ', pixels(1,a)
		 write(MsgOut,'(A,I10)') '   pix2 : ', pixels(2,a)
		 call error_msg('Gaussian kernel is extending outside the map box ')
		endif
	end do

	! ------------------------------------------------------------------------
	! GENERATE SIM IMAGE
	! ------------------------------------------------------------------------

	sim_image(:,:) = 0

	!$omp parallel do default(none)                                 &
	!$omp private(a,i,j , mu1, mu2, norm, gaussian)                   &
	!$omp shared(pixels, n_pix,image_size, pixel_size, &
	!$omp 	rot_coord, sigma, sim_image, gaussians_saved, atom_id, group_id)
	!

	do a=1, n_atoms_group
		do i=pixels(1,a), pixels(1,a) + n_pix
			do j=pixels(2,a), pixels(2,a) + n_pix

				mu1 = (i -1- (image_size / 2)) * pixel_size
				mu2 = (j -1- (image_size / 2)) * pixel_size

				norm = (rot_coord(1,a)-mu1)**2 + (rot_coord(2,a)-mu2)**2
				gaussian = exp(-norm / (2 * (sigma ** 2)))
				
				sim_image(i,j) = sim_image(i,j) + gaussian

				! save the gaussian to avoid computing them for the gradient
				gaussians_saved(i,j,a) = gaussian    
			end do
		end do
	end do

	!$omp end parallel do

	! ------------------------------------------------------------------------
	! COMPUTE CC
	! ------------------------------------------------------------------------

	sum_sim2 = 0.0
	sum_exp2 = 0.0
	sum_simpexp = 0.0

	do i=1, image_size
		do j=1, image_size
			sum_sim2 = sum_sim2 + sim_image(i,j) ** 2
			sum_exp2 = sum_exp2 + exp_image(i,j) ** 2
			sum_simpexp = sum_simpexp + exp_image(i,j) * sim_image(i,j)
		end do
	end do

	cv= (sum_simpexp / sqrt(sum_sim2 * sum_exp2))
	force_constant = enefunc%restraint_const(1,inum)
	eexp = force_constant * (1.0_wp - cv)

	! ------------------------------------------------------------------------
	! GRADIENT COMPUTATION
	! ------------------------------------------------------------------------
	if (calc_force) then

		const1 = 1/ sqrt(sum_sim2 * sum_exp2)
		const2 = sum_simpexp / (sqrt(sum_exp2) * sum_sim2**(1.5_wp))
		emfit_force(:,:) = 0.0_wp

		!$omp parallel do default(none)                                 &
		!$omp private(a,n, i, j,mu1, mu2, dpsim1, dpsim2,dcc1, dcc2)     &
		!$omp shared(pixels, n_pix, image_size, pixel_size, &
		!$omp 	rot_coord, sigma, sim_image, gaussians_saved, atom_id, group_id,&
		!$omp 	exp_image, const1, const2, force_constant,&
		!$omp 	emfit_force,inv_rot_matrix, force)
		!

		do a= 1, n_atoms_group

			n=  atom_id(a,group_id)
			dcc1 = 0.0_wp
			dcc2 = 0.0_wp

			do i=pixels(1,a), pixels(1,a) + n_pix
				do j=pixels(2,a), pixels(2,a) + n_pix

					mu1 = (i -1- (image_size / 2)) * pixel_size
					mu2 = (j -1- (image_size / 2)) * pixel_size

					dpsim1 = -(rot_coord(1,a) - mu1) * gaussians_saved(i,j,a) / (sigma ** 2)
					dpsim2 = -(rot_coord(2,a) - mu2) * gaussians_saved(i,j,a) / (sigma ** 2)

					dcc1 = dcc1 + ((exp_image(i,j) * dpsim1 * const1) - ( sim_image(i,j) * dpsim1 * const2))
					dcc2 = dcc2 + ((exp_image(i,j) * dpsim2 * const1) - ( sim_image(i,j) * dpsim2 * const2))
				end do
			enddo

			emfit_force(1,n) = dcc1 * force_constant
			emfit_force(2,n) = dcc2 * force_constant

			force(1:3,n) = force(1:3,n) + matmul(inv_rot_matrix(1:3,1:3),emfit_force(1:3,n))
		end do

		!$omp end parallel do
	end if


	! ------------------------------------------------------------------------
	! DEALLOCATE
	! ------------------------------------------------------------------------

	deallocate(rot_coord)
	deallocate(pixels)
	deallocate(gaussians_saved)

	! open (unit=66, file="force.txt")
	! write(66,*) force
	! close(66)
	! open (unit=66, file=trim(emfit_exp_image)//"_img.txt")
	! write(66,*) sim_image
	! close(66)
	! call error_msg('Stop')


    return
  end subroutine compute_energy_experimental_restraint_emfit_with_image
  !

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         Euler_angle2Matrix
  !> @brief             convert three rotation angle roll, title and yaw in dgree into a roation matrix
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         roll : angle in degree
  !! @param[in]         pitch: angle in degree
  !! @param[in]          yew : angle in degree
  !! @param[inout]       A   : 4*4 matrix for rotation and translation
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine Euler_angle2Matrix(roll, tilt, yaw, A)
 	real(wp),                intent(in)    :: roll
    real(wp),                intent(in)    :: tilt
    real(wp),                intent(in)    :: yaw
    real(wp),                intent(inout) :: A(4,4)

	!local variables
	real(wp) :: pi
	!real(wp) :: a, b, g
	real(wp) :: a1, b1, g1, sa, ca, sb, cb, sg, cg, cc, cs, sc, ss

    !a = 45.0_wp 
    !b = 90.0_wp
	!g = 45.0_wp

	pi = acos(-1.0_wp)
	a1=pi*roll/180.0_wp
	b1=pi*tilt/180.0_wp
	g1=pi*yaw/180.0_wp

	sa= sin(a1)
	ca= cos(a1)

	sb= sin(b1)
	cb= cos(b1)

	sg= sin(g1)
	cg= cos(g1)

	cc = cb * ca
	cs = cb * sa
	sc = sb * ca
	ss = sb * sa
	! finally in this section we will calculate the rotation matrix which will
	! be saved in A. 
	A(4,1:3)=0
	A(1:3,4)=0
	A(4,4)= 1
	A(1, 1) =  cg * cc - sg * sa
	A(1, 2) =  cg * cs + sg * ca
	A(1, 3) = -cg * sb
	A(2, 1) = -sg * cc - cg * sa
	A(2, 2) = -sg * cs + cg * ca
	A(2, 3) =  sg * sb
	A(3, 1) =  sc
	A(3, 2) =  ss
	A(3, 3) =  cb
	!write(*,*)"transation matrix:"
	!write(*,*) A(1,1),A(1,2),A(1,3)
	!write(*,*) A(2,1),A(2,2),A(2,3)
	!write(*,*) A(3,1),A(3,2),A(3,3)
  end subroutine Euler_angle2Matrix
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         Euler_angle2Matrix_Inv
  !> @brief             convert three rotation angle roll, title and yaw in dgree into a roation matrix
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]          B   : 4*4 input matrix for rotation and translation 
  !! @param[inout]       A   : 4*4 invers matrix
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine Euler_angle2Matrix_Inv(B,outflag, A )
 		real(wp),                intent(in)    :: B(4,4)
		logical,                 intent(in)    :: outflag
    	real(wp),                intent(inout) :: A(4,4)

	!local variables
	real(wp)					:: det1

	if(outflag .eqv. .TRUE.) then
		write(*,*)"transation matrix:"
		write(*,*) B(1,1),B(1,2),B(1,3)
		write(*,*) B(2,1),B(2,2),B(2,3)
		write(*,*) B(3,1),B(3,2),B(3,3)
	end if
	det1 = B(2,2)*B(3,3)*B(1,1)+B(1,2)*B(2,3)*B(3,1)+B(2,1)*B(3,2)*B(1,3)&
			-B(1,3)*B(2,2)*B(3,1)-B(1,2)*B(3,3)*B(2,1)-B(2,3)*B(3,2)*B(1,1)
	A(4,1:3)= 	0
	A(1:3,4)= 	0
	A(4,4)  = 	1
	A(1, 1) =  	(1.0_wp/det1)*(B(2,2)*B(3,3)-B(2,3)*B(3,2))
	A(2, 1) =  	(1.0_wp/det1)*(B(3,1)*B(2,3)-B(2,1)*B(3,3))
	A(3, 1) = 	(1.0_wp/det1)*(B(2,1)*B(3,2)-B(2,2)*B(3,1))
	A(1, 2) = 	(1.0_wp/det1)*(B(1,3)*B(3,2)-B(1,2)*B(3,3))
	A(2, 2) = 	(1.0_wp/det1)*(B(1,1)*B(3,3)-B(1,3)*B(3,1))
	A(3, 2) =  	(1.0_wp/det1)*(B(3,1)*B(1,2)-B(1,1)*B(3,2))
	A(1, 3) = 	(1.0_wp/det1)*(B(1,2)*B(2,3)-B(2,2)*B(1,3))
	A(2, 3) =  	(1.0_wp/det1)*(B(2,1)*B(1,3)-B(1,1)*B(2,3))
	A(3, 3) =  	(1.0_wp/det1)*(B(2,2)*B(1,1)-B(1,2)*B(2,1))
	if(outflag .eqv. .TRUE.) then
		write(*,*)"invers transation matrix:"
		write(*,*) det1
		write(*,*) A(1,1),A(1,2),A(1,3)
		write(*,*) A(2,1),A(2,2),A(2,3)
		write(*,*) A(3,1),A(3,2),A(3,3)
	end if 
  end subroutine Euler_angle2Matrix_Inv


   !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         readSPIDER
  !> @brief		takes filename as an input and read the image data from an image file in SPIDER format and return Image in the image plane as an output.
  !!			the spider image (ex. Image.spi) can be displayed by showSPIDER function
  !!			2D image data which will be returned as a SPIDER image shall be formated I(i,j,In) where i and j is pixel coordiante and In is the pixel intensity.
  !!			3D data may be written as either a SPIDER volume (default) or as a stack of 2D images (requires optional 3rd argument='stack'). 
  !!                    in this case we  are not dealing with 3D images. it is better to use situs for 3D purpose
  !!			SPIDER files do have a specific file extension, *.spi.
  !!			Examples
  !!			readSPIDERfile('img001.dat', I)
  !!			version 1.0 (Mar 2021) A. Mirzaei
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         datasize        : size of the image data
  !! @param[in]         stackarg        : stacking defined by the write function and imposed by the data type
  !! @param[inout]      header          : spi file header generated for spi write function
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine readSPIDER(filenameSR, I)
	character(MaxFilename),         intent(in)    	 :: filenameSR
	real(wp),						intent(inout)    :: I(:,:)

	
	!local variable
	real*4								:: nx
	!real*4, dimension(256)				:: file_header
	real*4, dimension(-int(image_size/2):int(image_size/2)-1,-int(image_size/2):int(image_size/2)-1)	:: image_out,image_in
	character(MaxFilename)				:: filenameW
	integer								:: integer_var,i1,j1,k1,l1
	integer 							:: ios
	integer								:: im_size, stop_condition
	real(wp)							:: norm_image
	write(*,*)"input experiment image:",filenameSR
	open(UNIT=10, FILE=filenameSR, FORM="UNFORMATTED", access="STREAM", convert='little_endian')
	integer_var = 1
	i1 = 1
	l1=1
	im_size =image_size
	stop_condition = im_size*im_size+256
	norm_image = 0_wp
	do while (integer_var .GE. 0)
		read(10,IOSTAT=integer_var) nx
		if(i1 .LE. 256)then
			file_header(i1)=nx
			!print *, nx
		else 
			if(i1 .EQ. 257)then
				write(*,*)"file header reading is finished"	
				j1 = -int(im_size/2)-1
				k1 = -int(im_size/2)
			end if
			!print *, nx, integer_var
			if(j1 .LT. int(im_size/2)-1)then
				j1= j1+1
			else
				j1=-int(im_size/2)
				k1=k1+1
			end if
			if((k1 .LT. int(im_size/2)).and.(j1 .LT. int(im_size/2)))then
				if(i1 .LT. stop_condition)then
					image_out(j1,k1)= nx
					norm_image= norm_image + nx*nx
					l1=l1+1
					!print *,i1,j1,k1, nx
				end if
				if(j1 == int(image_size/2)-2)then
					!write(*,*)image_out(j1,k1)
				end if
			end if
			if(integer_var .LT. 0)then
				write(*,*)"end of file reached"
				print *, i1, l1
				exit
			end if
		end if
		i1=i1+1
	end do
	close(10)
	print *, i1, l1
	
	write(*,*)"pixel inside function (61,61)=", image_out(int(image_size/2)-2,int(image_size/2)-2)
	write(*,*) "Spider file was successfully read and reported"
	write(*,*) "norm of image inside read funciton: ",sqrt(norm_image)
	!open(UNIT=20, FILE="exp000010.spi", FORM="UNFORMATTED", access="STREAM", status ='new', convert='little_endian',iostat=ios)
	!print *, ios
	!file_header(i)
	!write(20)file_header
	!write(20)image_out
	!close(20)
	I = real(image_out,8)
	write(*,*)"########### image data is successfully returned at the end of read function ###############"
  end subroutine readSPIDER

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         readSPIDERheader
  !> @brief		takes filename as an input and read the image header of an image file in SPIDER format and return header 
  !!			the spider header can be used by SPIDERread and SPIDERwrite functions
  !!			header data of a SPIDER image will be returned as header(256) array 
  !!			Examples
  !!			readSPIDERfile('img001.dat', header)
  !!			version 1.0 (Mar 2021) A. Mirzaei
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         filenameSR      : SPIDER image name and directory
  !! @param[inout]      header          : spi file header read out of SPIDER image
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine readSPIDERheader(filenameSR, header)
	character(MaxFilename),         intent(in)    	 :: filenameSR
	real*4,						    intent(inout)    :: header(:)

	
	!local variable
	real*4								:: nx
	real*4,dimension(1:256)				:: Spider_header
	character(MaxFilename)				:: filenameW
	integer								:: integer_var,i1,j1,k1,l1
	integer 							:: ios, none_zero_mem
	write(*,*)"input spider image name and directory:",filenameSR
	open(UNIT=10, FILE=filenameSR, FORM="UNFORMATTED", access="STREAM", convert='little_endian')
	integer_var = 1
	i1 = 1

	do while (integer_var .GE. 0)
		read(10,IOSTAT=integer_var) nx
		if(i1 .LE. 256)then
			 Spider_header(i1)=nx
			!print *,i1, nx
		else 
			if(i1 .EQ. 257)then
				write(*,*)"file header is finished"	
				none_zero_mem = count( Spider_header/=0)
			end if
		end if
		i1=i1+1
	end do
	close(10)
	print *, i1,none_zero_mem
	
	if(none_zero_mem /= 0) then
		write(*,*) "the header of the Spider image was successfully read and reported"
	else
		write(*,*) "the header of the Spider image was not successful, check your image file"
		return
	end if

	header = Spider_header

  end subroutine readSPIDERheader
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         writeSPIDER
  !> @brief             takes filename, and Image in the image plane as an input and write the image data into an image file in SPIDER format.
  !!			the spider image (ex. Image.spi) can be displayed by showSPIDER function
  !!			2D image data which will be written as a SPIDER image shall be formated I(i,j,In) where i and j is pixel coordiante and In is the pixel intensity.
  !!			3D data may be written as either a SPIDER volume (default) or as a stack of 2D images (requires optional 3rd argument='stack'). 
  !!                    in this case we  are not dealing with 3D images.it is better to use situs for 3D purpose
  !!			SPIDER files do have a specific file extension, *.spi.
  !!			Examples
  !!			writeSPIDER('img001.spi', I)
  !!			version 1.0 (Mar 2021) A. Mirzaei
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         filename        : filename for instance 'image.spi' (do not forget quotation around filename: '')
  !! @param[in]         I               : Image data (image peixel in i and j coordiante)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine writeSPIDER(filenameSR, I)
	character(24),     	intent(in)    	:: filenameSR
	real*4,				intent(inout)   :: I(:,:)

	! local variable 
	integer								:: i1,j1, x, y
	real*4, dimension(-int(image_size/2):int(image_size/2)-1,-int(image_size/2):int(image_size/2)-1)			:: image_out
	integer								:: im_size
	integer 							:: ios


	im_size = image_size
	call creatSpiderHeader(I,file_header)
	open(UNIT=20, FILE=filenameSR, FORM="UNFORMATTED", access="STREAM", status ='new', convert='little_endian',iostat=ios)
	print *, ios
	write(20)file_header
	write(20)I

	close(20)
  end subroutine writeSPIDER
	!======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         creatSpiderHeader
  !> @brief             takes filename, and Image in the image plane as an input and write the image data into an image file in SPIDER format.
  !!			the spider image (ex. Image.spi) can be displayed by showSPIDER function
  !!			2D image data which will be written as a SPIDER image shall be formated I(i,j,In) where i and j is pixel coordiante and In is the pixel intensity.
  !!			3D data may be written as either a SPIDER volume (default) or as a stack of 2D images (requires optional 3rd argument='stack'). 
  !!                    in this case we  are not dealing with 3D images.it is better to use situs for 3D purpose
  !!			SPIDER files do have a specific file extension, *.spi.
  !!			Examples
  !!			creatSpiderHeader(I, header)
  !!			version 1.0 (Mar 2021) A. Mirzaei
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         filename        : filename for instance 'image.spi' (do not forget quotation around filename: '')
  !! @param[in]         I               : Image data (image peixel in i and j coordiante)
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine creatSpiderHeader(I, header)
	real*4,				intent(in)   	:: I(:,:)
	real*4,				intent(inout)   :: header(:)

	! local variable 
	integer								:: i1,j1, x, y
	real*4, dimension(256)				:: generated_header
	integer								:: nsam,nrow,nslice
	integer 							:: ios
	real*4								:: labrec1,labbyt1,lenbyt1
	
	! header parameter initialization
	header = 0.0
	nsam = size(I,1)
	nrow = size(I,2)
	 
	lenbyt1 = nsam*4
	labrec1 = ceiling(1024/lenbyt1)
	if(MOD(1024.0,lenbyt1).NE. 0)then
		labrec1 = labrec1 + 1
	end if
	labbyt1 = labrec1 * lenbyt1 
	nslice = 1.0
	if(nslice == 1)then
		write(*,*)"writing an image"
		generated_header(5)		=1.0
		generated_header(1)		=nslice
		generated_header(2)		=nrow
		generated_header(12)	=nsam
		generated_header(13)	=labrec1
		generated_header(22)	=labbyt1
		generated_header(23) 	=lenbyt1
	else
		write(*,*)"writing a volume"
		generated_header(5)		=3.0
		generated_header(1)		=nslice
		generated_header(2)		=nrow
		generated_header(12)	=nsam
		generated_header(13)	=labrec1
		generated_header(22)	=labbyt1
		generated_header(23) 	=lenbyt1
	end if

	!call creatSpiderHeader(file_header)
	header = generated_header
  end subroutine creatSpiderHeader

end module at_experiments_mod
