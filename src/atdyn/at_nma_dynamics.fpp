!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_md_vverlet_mod
!> @brief   perform molecular dynamics simulation with velocity verlet algorithm
!! @authors Jaewoon Jung (JJ), Takaharu Mori (TM), Chigusa Kobayashi (CK),
!           Tadashi Ando (TA), Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_nma_dynamics_mod

  use fileio_control_mod !EDIT REMI

  implicit none
  private

  ! subroutines
  public  :: vverlet_dynamics
  private :: initial_vverlet
  private :: control_temp_pres_vver1
  private :: control_temp_pres_vver2
  private :: evans_thermostat_vverlet 
  private :: berendsen_thermostat_vverlet
  private :: andersen_thermostat_vverlet
  private :: nosehoover_thermostat_vverlet
  private :: langevin_thermostat_vv1
  private :: langevin_thermostat_vv2
  private :: langevin_barostat_vv1
  private :: langevin_barostat_vv2
  private :: update_barostat
  private :: bussi_thermostat_vverlet
  private :: bussi_barostat_vv1
  private :: bussi_barostat_vv2
  private :: simulated_annealing_vverlet


contains  
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    vverlet_dynamics
  !> @brief        velocity verlet integrator
  !! @authors      JJ, TM, CK, TA, KY
  !! @param[inout] output      : output information
  !! @param[inout] molecule    : molecule information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynamic variables information
  !! @param[inout] dynamics    : dynamics information
  !! @param[inout] pairlist    : non-bond pair list information
  !! @param[inout] boundary    : boundary conditions information
  !! @param[inout] constraints : bond constraints information
  !! @param[inout] ensemble    : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine vverlet_dynamics(output, molecule, enefunc, dynvars, dynamics, &
    pairlist, boundary, constraints, ensemble)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_dynamics), target, intent(inout) :: dynamics
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary
    type(s_constraints),      intent(inout) :: constraints
    type(s_ensemble),         intent(inout) :: ensemble

    ! local variables
    real(wp)                  :: simtim, dt
    integer                   :: i, j, natom
    integer                   :: nsteps, istart, iend

    real(wp),         pointer :: coord(:,:), coord_ref(:,:)
    real(wp),         pointer :: vel(:,:), vel_ref(:,:)
    real(wp),         pointer :: force(:,:), mass(:), inv_mass(:)
    character(256)            :: folder,basename
    character                 :: num*9
    logical                   :: savefile

    !EDIT REMI
    integer :: nmodes, start_modes
    integer :: atm
    real(wp) :: global_dt
    real(wp), dimension (:,:,:), allocatable :: normalModeVec 
    real(wp), dimension (:), allocatable :: global 
    real(wp), dimension (:,:), allocatable :: local
    real(wp), dimension (:,:), allocatable :: coord0
    real(wp), dimension (:,:), allocatable :: global_coord
    real(wp), dimension (:), allocatable :: global_vel
    real(wp), dimension (:), allocatable :: global_force
    character(len=:), allocatable :: fprefix
    character(256)            :: fmodes
    logical :: global_fit

    integer :: icount

    mass      => molecule%mass
    inv_mass  => molecule%inv_mass
    coord     => dynvars%coord
    coord_ref => dynvars%coord_ref
    force     => dynvars%force
    vel       => dynvars%velocity
    vel_ref   => dynvars%velocity_ref

    natom     =  molecule%num_atoms
    istart    =  dynamics%istart_step
    iend      =  dynamics%iend_step
    nsteps    =  dynamics%nsteps
    dt        =  dynamics%timestep/AKMA_PS
    simtim    =  dynamics%initial_time

    ! <EDIT REMI>
    print*, 'RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR'
    call open_ctrlfile
    call read_ctrlfile_string(handle, 'Experiments', 'emfit_nma',  &
    fmodes)
    print*, fmodes
    global_fit = .false.
    global_dt = 10.0_wp
    fprefix = '/home/guest/PycharmProjects/bayesian-md-nma/data/AK/modes_psf/vec.'
    nmodes = 3
    start_modes= 6

    global_dt = global_dt*dt    
    allocate ( normalModeVec(3,natom,nmodes) ) 
    allocate ( global(nmodes) ) 
    allocate ( local(3, natom) ) 
    allocate ( coord0(3, natom) ) 
    allocate ( global_coord(3, natom) ) 
    allocate ( global_vel(nmodes) ) 
    allocate ( global_force(nmodes) ) 
    do i = 1,nmodes
    if (i+ start_modes<10) then
    write(fmodes, '(A,I1)') fprefix,  i + start_modes
    else
    write(fmodes, '(A,I2)') fprefix,  (i + start_modes)
    endif
    open(unit=66, file=fmodes)
    read(66,*)  normalModeVec(:,:, i)
    close(66)
    end do
    print*, 'GLOBAL DT='
    print*, global_dt
    ! <\EDIT REMI>

    if (dynamics%target_md) enefunc%rmsd_force = 1.0_wp / (dt*dt)
    if (abs(dynamics%initial_value) .lt. 0.001_wp)  then
    if (dynamics%target_md) &
    dynamics%initial_value =  &
    dynvars%energy%restraint_cv(enefunc%target_function)
    if (dynamics%steered_md) &
    dynamics%initial_value =  &
    dynvars%energy%restraint_cv(enefunc%steered_function)
    endif
    enefunc%target_value = dynamics%initial_value 

    ! first-step MD
    ! 
    if (.not. dynamics%restart) then
    call initial_vverlet(output, molecule, enefunc, dynamics, pairlist, &
    boundary, ensemble, constraints, dynvars)
    end if

    !
    ! stop NPAT & NPgT
    ! 
    if (ensemble%ensemble == EnsembleNPAT .or.  &
    ensemble%ensemble == EnsembleNPT  .or.  &
    ensemble%ensemble == EnsembleNPgT)      & 
    call error_msg('Vverlet_dynamics> Barostats are not allowed in ATDYN')

    !EDIT REMI init velocities
    print*, '//////////////// Enetering Remi Code ///////////////////'
    do j =1, nmodes
    global_vel(j) = 0.0_wp
    end do
    do i =1, natom
    do j =1, nmodes
    global_vel(j) =global_vel(j) + normalModeVec(1,i,j)* vel(1, i)
    global_vel(j) =global_vel(j) + normalModeVec(2,i,j)* vel(2, i) 
    global_vel(j) =global_vel(j) + normalModeVec(3,i,j)* vel(3, i)
    end do
    end do
    print*, 'GLOBAL VEL'
    print*, global_vel
    do i = 1,nmodes
    global(i) =0.0_wp
    end do
    do i =1, natom
    do j =1, 3
    local(j, i) =0.0_wp
    coord0(j, i) = coord(j, i)
    global_coord(j,i) = 0.0_wp
    end do
    end do

    #ifdef KCOMP
    ! Start performance check on K computer
    !
    call fpcoll_start()
    #endif

    ! Reset qm_count
    if(enefunc%qmmm%do_qmmm) enefunc%qmmm%qm_count = istart

    ! Main MD loop
    ! at this point
    !   coord, vel, and force are at t = 0   , if restart off
    !   coord, vel, and force are at t = t+dt, if restart on
    !

    do i = istart, iend

    call timer(TimerIntegrator, TimerOn)

    simtim = simtim + dt*AKMA_PS
    dynvars%time = simtim
    dynvars%step = i

    if (dynamics%target_md .or. dynamics%steered_md) &
    enefunc%target_value = dynamics%initial_value &
      + (dynamics%final_value-dynamics%initial_value) &
      *real(dynvars%step,wp)/real(nsteps,wp)

    enefunc%rpath_sum_mf_flag = enefunc%rpath_flag

    ! Save coordinates(t) and velocities(t)
    !
    do j = 1, natom
    coord_ref(1,j) = coord(1,j)
    coord_ref(2,j) = coord(2,j)
    coord_ref(3,j) = coord(3,j)
    vel_ref(1,j) = vel(1,j)
    vel_ref(2,j) = vel(2,j)
    vel_ref(3,j) = vel(3,j)
    end do

    ! Compute velocities(t + 1/2 dt) and coordinates(t + dt)
    !
    if (ensemble%tpcontrol /= TpcontrolNoseHoover .and. &
    ensemble%tpcontrol /= TpcontrolLangevin   .and. &
    ensemble%tpcontrol /= TpcontrolBussi      ) then

    ! Newtonian dynamics
    !   v(t+1/2dt) = v(t) + 0.5dt*F(t)/m
    !   r(t+dt)    = r(t) + dt*v(t+1/2dt)
    !
    do j = 1, natom
    vel(1,j) = vel(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
    vel(2,j) = vel(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
    vel(3,j) = vel(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)

    ! EDIT REMI
    local(1,j) = local(1,j) + dt*vel(1,j)
    local(2,j) = local(2,j) + dt*vel(2,j)
    local(3,j) = local(3,j) + dt*vel(3,j)
    end do

    !EDIT REMI
    if (global_fit) then
    do j =1, nmodes
    global_force(j) = 0.0_wp
    end do
    do j = 1, natom
    global_coord(1,j)=0.0_wp
    global_coord(2,j)=0.0_wp
    global_coord(3,j)=0.0_wp
    end do
    do atm = 1, natom
    do j =1, nmodes
    global_force(j) = global_force(j) + normalModeVec(1, atm,j)* force(1,atm) * inv_mass(j)
    global_force(j) = global_force(j) + normalModeVec(2, atm,j)* force(2,atm) * inv_mass(j)
    global_force(j) = global_force(j) + normalModeVec(3, atm,j)* force(3,atm) * inv_mass(j)
    end do
    end do
    do j =1, nmodes
    global_vel(j) = global_vel(j) + 0.5_wp*global_dt*global_force(j)
    global(j) = global(j) + global_dt*global_vel(j)
    end do

    do atm = 1, natom
    do j =1, nmodes
    global_coord(1,atm) = global_coord(1,atm) + (global(j) * normalModeVec(1,atm,j))
    global_coord(2,atm) = global_coord(2,atm) + (global(j) * normalModeVec(2,atm,j))
    global_coord(3,atm) = global_coord(3,atm) + (global(j) * normalModeVec(3,atm,j))
    end do
    end do
    do j = 1, natom
    coord(1,j) = coord0(1,j) + local(1,j) + global_coord(1,j) 
    coord(2,j) = coord0(2,j) + local(2,j) + global_coord(2,j) 
    coord(3,j) = coord0(3,j) + local(3,j) + global_coord(3,j) 
    end do
    else
    do j = 1, natom
    coord(1,j) = coord0(1,j) + local(1,j)
    coord(2,j) = coord0(2,j) + local(2,j)
    coord(3,j) = coord0(3,j) + local(3,j)
    end do
    endif
    ! Coordinate constraint (RATTLE VV1)
    !
    if (constraints%rigid_bond) then
    call compute_constraints(ConstraintModeVVER1, dt, molecule, &
            dynvars, constraints)
    dynvars%virial_const(1:3,1:3) = 2.0_wp * dynvars%virial_const(1:3,1:3)

    end if

    else

    call control_temp_pres_vver1(molecule, dynamics, ensemble,  &
              constraints, boundary, dynvars)

    end if

    call timer(TimerIntegrator, TimerOff)

    ! calculate potential energy(t + dt), force(t + dt), and virial(t + dt)
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, &
    mod(i,dynamics%eneout_period) == 0,    &
    enefunc%nonb_limiter,  &
    dynvars%coord,         &
    dynvars%energy,        &
    dynvars%temporary,     &
    dynvars%force,         &
    dynvars%virial,        &
    dynvars%virial_extern)


    call timer(TimerIntegrator, TimerOn)

    if (ensemble%tpcontrol == TpcontrolLangevin) then

    if (ensemble%ensemble == EnsembleNVT) then

    call langevin_thermostat_vv2(molecule, dynamics, ensemble, &
                boundary, constraints, dynvars)

    else

    call langevin_barostat_vv2  (molecule, dynamics, ensemble, &
                constraints, boundary, dynvars)

    end if

    else

    ! Compute velocities(t + dt)
    !   v(t+dt) = v(t+1/2dt) + 0.5dt*F(t+dt)/m
    !
    do j = 1, natom
    vel(1,j) = vel(1,j) + 0.5_wp*dt*force(1,j)*inv_mass(j)
    vel(2,j) = vel(2,j) + 0.5_wp*dt*force(2,j)*inv_mass(j)
    vel(3,j) = vel(3,j) + 0.5_wp*dt*force(3,j)*inv_mass(j)
    end do

    !EDIT REMI
    ! do j =1, nmodes
    !   global_vel(j) = global_vel(j) + 0.5_wp*global_dt*global_force(j)
    ! end do

    ! Velocty constraint (RATTLE VV2)
    !   coord     is at t + dt, vel (unconstrained) is at t + dt
    !   coord_ref is at t,      vel_ref             is at t
    !   compute constrained velocities(t+dt)
    !   add constraint virial(t+dt)
    !   Note that constraint virial(t) is added to virial(t+dt) in leapfrog
    !             which is correct?
    !
    if (constraints%rigid_bond) then

    call compute_constraints(ConstraintModeVVER2, &
            dt, molecule, dynvars, constraints) 

    endif

    ! Control temperature VV2
    !
    if (ensemble%ensemble /= EnsembleNVE) then

    call control_temp_pres_vver2(molecule, dynamics, ensemble, dynvars)

    end if
    end if

    call timer(TimerIntegrator, TimerOff)

    ! Remove translational and rotational motion about COM(t + dt)
    !
    if (dynamics%stoptr_period > 0) then

    if (mod(i,dynamics%stoptr_period) == 0) then
    call stop_trans_rotation(molecule%num_atoms, molecule%mass, &
            dynamics%stop_com_translation,     &
            dynamics%stop_com_rotation,        &
            dynvars%coord,                     &
            boundary%fixatm,                   &
            dynvars%velocity)
    end if

    end if


    ! Update nonbond pairlist for coordinates(t + 2dt)
    !
    if (dynamics%nbupdate_period > 0) then
    if (mod(i,dynamics%nbupdate_period) == 0 .and. real_calc) then

    ! Update boundary if pressure is controlled
    !
    if (ensemble%use_barostat) &
    call update_boundary(enefunc%table%table,  &
          enefunc%pairlistdist, &
          boundary)
    call update_pairlist(enefunc, boundary, coord, pairlist)

    !if(enefunc%spot_use) &
    !  call update_spot_atomlist(molecule%num_atoms, coord, boundary)

    end if
    end if


    ! Simulated annealing
    !
    if (dynamics%anneal_period > 0) then

    if (mod(i,dynamics%anneal_period) == 0) then
    call simulated_annealing_vverlet(dynamics, enefunc, ensemble)
    end if

    end if

    ! Update GAMD
    !
    if (enefunc%gamd%update_period > 0) then
    if (mod(i,enefunc%gamd%update_period) == 0) then
    call update_gamd_vverlet(output, enefunc, dynvars)
    end if
    end if

    ! Output energy(t + dt) and dynamical variables(t + dt)
    !
    call output_md(output, molecule, enefunc, dynamics, boundary, &
    ensemble, dynvars)

    end do

    !EDIT REMI
    deallocate(normalModeVec) 
    deallocate(global) 
    deallocate(global_force) 
    deallocate(global_vel)
    deallocate(local) 
    deallocate(coord0) 
    deallocate(global_coord)

    #ifdef KCOMP
    ! Stop performance check on K computer
    !
    call fpcoll_stop()
    #endif


    return

    end subroutine vverlet_dynamics