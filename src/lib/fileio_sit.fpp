!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_sit_mod
!> @brief   read situs density map file
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_sit_mod

  use fileio_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_sit
    real(wp), allocatable  :: map_data(:,:,:)
    real(wp)               :: dx
    real(wp)               :: dy
    real(wp)               :: dz
    real(wp)               :: x0
    real(wp)               :: y0
    real(wp)               :: z0
    integer                :: nx
    integer                :: ny
    integer                :: nz
  end type s_sit

  type, public :: s_image
    real(wp), allocatable  :: image_data(:,:,:)
    real(wp)               :: xs0
    real(wp)               :: ys0
    integer                :: nsx
    integer                :: nsy
    real(wp)               :: Vs0
  end type s_image
  ! subroutines
  public  :: input_sit
  public  :: alloc_sit
  public  :: dealloc_sit
  public  :: read_sit
  public  :: read_image
  public  :: read_pdb

  public  :: ch_cap
  public  :: pdb_check
  public  :: get_unit
  public  :: pdb_init
  public  :: pdb_print_coord
  public  :: pdb_summary
  public  :: pdb_write_atom
  public  :: s_eqi
  public  :: timestamp

  public  :: writeSPIDERfile
  public  :: makeSPIDERheader
contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_sit
  !> @brief        open, read, and close sitfile 
  !! @authors      TM
  !! @param[in]    sit_filename : filename of sitfile
  !! @param[out]   sit          : EM data information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_sit(sit_filename, sit)

    ! formal arguments
    character(*),            intent(in)    :: sit_filename
    type(s_sit),             intent(inout) :: sit

    ! local variables
    integer                  :: unit_no

    ! open sitfile
    !
    call open_file(unit_no, sit_filename, IOFileInput)

    ! read sitfile
    !
    call read_sit(unit_no, sit)

    ! close sitfile
    !
    call close_file(unit_no)

    ! write sumary
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Input_Emap> Summary of sitfile'
      write(MsgOut,'(A20,F10.3)') '  voxel size x    = ', sit%dx
      write(MsgOut,'(A20,F10.3)') '  voxel size y    = ', sit%dy
      write(MsgOut,'(A20,F10.3)') '  voxel size z    = ', sit%dz
      write(MsgOut,'(A20,I10  )') '  num x increments= ', sit%nx
      write(MsgOut,'(A20,I10  )') '  num y increments= ', sit%ny
      write(MsgOut,'(A20,I10  )') '  num z increments= ', sit%nz
      write(MsgOut,'(A20,F10.3)') '  first voxel xcrd= ', sit%x0
      write(MsgOut,'(A20,F10.3)') '  first voxel ycrd= ', sit%y0
      write(MsgOut,'(A20,F10.3)') '  first voxel zcrd= ', sit%z0
      write(MsgOut,'(A)') ''
    end if

    return

  end subroutine input_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_sit
  !> @brief        allocate EM data information
  !! @authors      TM
  !! @param[inout] sit     : EM data information
  !! @param[in]    var_size : allocation size
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine alloc_sit(sit, var_size1, var_size2, var_size3)

    ! formal arguments
    type(s_sit),             intent(inout) :: sit
    integer,                 intent(in)    :: var_size1
    integer,                 intent(in)    :: var_size2
    integer,                 intent(in)    :: var_size3
    
    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat

    
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables
    !
    if (allocated(sit%map_data)) then
      if (size(sit%map_data(:,0,0)) == var_size1) &
        return
      deallocate(sit%map_data,      &
                 stat = dealloc_stat)
    end if

    allocate(sit%map_data(0:var_size1-1, 0:var_size2-1, 0:var_size3-1), &
             stat = alloc_stat)

    sit%map_data(0:var_size1-1,0:var_size2-1,0:var_size3-1) = 0.0_wp

    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_sit
  !> @brief        deallocate EM data information
  !! @authors      TM
  !! @param[inout] sit     : EM data information
  !! @param[in]    variable : an variable to be allocated 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_sit(sit)

    ! formal arguments
    type(s_sit),             intent(inout) :: sit

    ! local variables
    integer                  :: dealloc_stat
    

    dealloc_stat = 0

    if (allocated(sit%map_data)) then
      deallocate(sit%map_data,       &
                 stat = dealloc_stat)
    end if

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_sit
  !> @brief        read data from sitfile
  !! @authors      TM
  !! @param[in]    unit_no : unit number of sitfile
  !! @param[out]   sit    : EM data inforation
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_sit(unit_no, sit)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_sit),             intent(inout) :: sit

    ! local variables
    real(wp)       :: dx, x0, y0, z0
    integer        :: nx, ny, nz

    read(unit_no, *) dx, x0, y0, z0, nx, ny, nz

    sit%dx = dx
    sit%dy = dx
    sit%dz = dx
    sit%x0 = x0
    sit%y0 = y0
    sit%z0 = z0
    sit%nx = nx
    sit%ny = ny
    sit%nz = nz

    call alloc_sit(sit, nx, ny, nz)

    read(unit_no, *) sit%map_data

    return

  end subroutine read_sit

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_image
  !> @brief        open, read, and close image file 
  !! @authors      Alex Mirzaei
  !! @param[in]    image_filename : filename of pseudoatom
  !! @param[out]   image          : image information including pixel value
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine read_image(image_filename, image)
    ! formal arguments
    !character(*),            intent(in)    :: Pseudo_filename
    integer,                 intent(in)    :: image_filename
    type(s_image),           intent(inout) :: image
    
    ! local variables
    CHARACTER(128) 	    :: buffer
    integer                 :: strlen, rows, cols , io, i
    real, dimension(:,:), allocatable :: xcoord
    
    ! local variables
    real(wp)       ::  xs0, ys0
    integer        :: nsx, nsy
    write(MsgOut,'(A)') 'Input_image> Summary of image detail'

    ! read image file and show the coordinates and values
    write(*,*) "*************---------*************"
    write(*,*) "reading input image for energy minimization written by Alex"
    write(*,*) "*************---------*************"
    
    
    OPEN (1, file = 'xmipp_image.txt', status='old', action='read')

    OPEN (2, file = 'genesis_log.txt', status='old', action='write')

    !Count the number of columns

     read(1,'(a)') buffer !read first line WITH SPACES INCLUDED
     REWIND(1) !Get back to the file beginning

     strlen = len(buffer) !Find the REAL length of a string read
     do while (buffer(strlen:strlen) == ' ') 
  	strlen = strlen - 1 
     enddo

     cols=0 !Count the number of spaces in the first line
     do i=0,strlen
       if (buffer(i:i) == ' ') then
         cols=cols+1
       endif
     enddo

     cols = cols+1

     !Count the number of rows

     rows = 0 !Count the number of lines in a file
     DO
       READ(1,*,iostat=io)
       IF (io/=0) EXIT
       rows = rows + 1
     END DO
     REWIND(1)

     print*, 'Number of rows:', rows
     print*, 'Number of columns:', cols
     allocate(xcoord(rows,cols))

     do I=1,rows,1
       read(1,*) xcoord(I,:)
       write(*,*) xcoord(I,:)
       write(2,*) xcoord(I,:)
    enddo
    write(*,*) "*************---------*************"
    write(*,*) "end of reading part written by Alex"
    write(*,*) "*************---------*************"
    
    return
  end subroutine read_image
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_pdb structure
  !> @brief        open, read, and close standard pdb file 
  !! @authors      Alex Mirzaei
  !! @param[in]    pdb_filename : filename of pdb structure
  !! @param[out]   pdb_data     : pdb data information
  !! @param[out]   coord        : coordinate of pdb atoms
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine read_pdb(pdb_filename, pdb_data, coord_pdb)

 
    !*****************************************************************************80
!
!! MAIN is the main program for PDB_READ_TEST.
!
!  Discussion:
!
!    PDB_READ_TEST tests the PDB_READ library.
!
!    PDB_READ_TEST can be invoked from the command line, with the file
!    to be read included as a command line argument:
!
!      pdb_read_TEST file.pdb > report.txt
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2006
!
!  Author:
!
!    John Burkardt
!
   implicit none 

  ! formal argument
   integer,                 intent(in)    :: pdb_filename
   type(s_image),           intent(inout) :: pdb_data
   real(wp), dimension(:,:),allocatable,intent(inout) :: coord_pdb

  ! parameters
  integer ( kind = 4 ), parameter :: atom_max = 1656
  integer ( kind = 4 ), parameter :: maxres = 1000
  integer ( kind = 4 ), parameter :: mxpatm = 18

  !integer              pdb_filename
  !type(s_image)           pdb_data

  integer ( kind = 4 ) atom_num,I
  real ( kind = 8 ) coord(atom_max,3)
  !real ( kind = 8 ) coord_pdb(atom_max,3)
  character ( len = 80 ) filepdb
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ilen
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipxfargc
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) numarg
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  integer ( kind = 4 ) prtatm(maxres,mxpatm)
  character ( len = 3 ) resnam(atom_max)
  integer ( kind = 4 ) resnum(atom_max)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin

   
  !atom_coord(atom_max,3) = 0.0_wp
  
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_READ_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the PDB_READ library.'
!
!  Get the number of command line arguments.
!
!  Old style:
!
  numarg = iargc ( )
!
!  New style:
!
! numarg = ipxfargc ( )

 ! if ( 1 <= numarg ) then

 !   iarg = 1
!
!  Old style:
!
!    call getarg ( iarg, filepdb )
!
!  New style:
!
!   call pxfgetarg ( iarg, filepdb, ilen, ierror )
!
!   if ( ierror /= 0 ) then
!     write ( *, '(a)' ) ' '
!     write ( *, '(a)' ) 'PDB_READ_TEST - Fatal error!'
!     write ( *, '(a)' ) '  Could not read a command line argument.'
!     stop
!   end if

 ! else

  !  write ( *, '(a)' ) '  Enter the PDB file name:'
  !  read ( *, '(a)', iostat = ios ) filepdb
  !  if ( ios /= 0 ) then
  !    write ( *, '(a)' ) ' '
   !   write ( *, '(a)' ) 'PDB_READ_TEST - Fatal error!'
    !  write ( *, '(a)' ) '  Could not read the PDB file name.'
  !    stop
  !  end if

  !end if
  !filepdb ="2FW0.pdb"
  filepdb ="AK.pdb"
  write ( *, '(a)') 'Reading PDB file : ' // trim ( filepdb )
!
!  Open the PDB file.
!
  call get_unit ( input_unit )

  open ( unit = input_unit, file = filepdb, status = 'old', iostat = ios )

  if ( ios /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_READ_TEST - Fatal error!'
    write ( *, '(a)' ) '  Could not open the PDB file.'
    stop
  end if
!

!  Read the ATOM information from the file.
!
  allocate(coord_pdb(atom_max,3))
  call pdb_read_atom ( coord, hiatom, hires, input_unit, atom_max, maxres, &
    mxpatm, atom_num, numchain, numline, numlost, numres, prtatm, &
    resnam, resnum, xmax, xmin, ymax, ymin, zmax, zmin )
!
 ! write(*,*)"########the pdb file is read here and coordinates are displayed########"
  
  do I=1,atom_max
    coord_pdb(I,:) = coord(I,:)
  enddo
 ! write(*,*)coord_pdb
!  Close the PDB file.
! coord_pdb = coord
  close ( unit = input_unit )
!
!  Print a summary of the information.
!
  call pdb_summary ( hiatom, hires, atom_max, maxres, mxpatm, &
    atom_num, numchain, numline, numlost, numres, xmax, xmin, &
    ymax, ymin, zmax, zmin )
!
!  Reset any data that is out of bounds.
!
  call pdb_check ( hires, hiatom, atom_max, maxres, atom_num, numres )
!
!  Print out the PRTATM array.
!
  call pdb_print_prtatm ( maxres, mxpatm, numres, prtatm )
!
!  Print out the COORD array.
!
  call pdb_print_coord ( coord, atom_max, atom_num, resnam, resnum )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_READ_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  !stop 0
  end subroutine read_pdb


  !================ this file just added=======================================8
  subroutine ch_cap ( c )

!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none

  character c
  integer ( kind = 4 ) itemp

  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT.
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5 and 6).
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
subroutine pdb_check ( hires, hiatom, atom_max, maxres, atom_num, numres )

!*****************************************************************************80
!
!! PDB_CHECK truncates the PDB size data to internal maxima, if neccessary.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) HIATOM, the maximum atom index 
!    encountered.
!
!    Input/output, integer ( kind = 4 ) HIRES, the maximum residue index 
!    encountered.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input/output, integer ( kind = 4 ) ATOM_NUM, the number of atoms read in.
!
!    Input/output, integer ( kind = 4 ) NUMRES, the number of residuals read. 
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) atom_num
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) numres

  if ( maxres < hires ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The highest residual index exceeds the maximum.'
    hires = maxres
  end if

  if ( maxres < numres ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The residual total exceeds the maximum.'
    numres = maxres
  end if

  if ( atom_max < hiatom ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The highest atom index exceeds the maximum.'
    hiatom = atom_max
  end if

  if ( atom_max < atom_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PDB_CHECK: Warning!'
    write ( *, '(a)' ) '  The atom total index exceeds the maximum.'
    atom_num = atom_max
  end if

  return
end
subroutine pdb_init ( coord, hiatom, hires, atom_max, maxres, mxpatm, &
  atom_num, numchain, numline, numlost, numres, prtatm, resnam, &
  resnum, xmax, xmin, ymax, ymin, zmax, zmin )

!*****************************************************************************80
!
!! PDB_INIT initializes data associated with a given PDB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COORD(ATOM_MAX,3) the coordinates of atoms.
!
!    Output, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Output, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Output, integer ( kind = 4 ) ATOM_NUM, the number of atoms read in.
!
!    Output, integer ( kind = 4 ) NUMCHAIN, the number of chains.
!
!    Output, integer ( kind = 4 ) NUMLINE, the number of lines of text in 
!    the file.
!
!    Output, integer ( kind = 4 ) NUMLOST, the number of atom items lost.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of residuals read in. 
!
!    Output, integer ( kind = 4 ) PRTATM(MAXRES,MXPATM), contains the indices
!    of the atoms that make up each residue.  For the I-th residue, 
!    the following entries are significant:
!      PRTATM(I,1) - N, the nitrogen;
!      PRTATM(I,2) - CA, the C-Alpha;
!      PRTATM(I,3) - C, the carbon;
!      PRTATM(I,4) - O, the oxygen.
!      PRTATM(I,5) - CB, another carbon?
!      PRTATM(I,5+) - the other atoms.
!
!    Output, character ( len = 3 ) RESNAM(ATOM_MAX), the residue name for each
!    entry.
!
!    Output, integer ( kind = 4 ) RESNUM(ATOM_MAX), the index of the residue 
!    to which the atom belongs.
!
!    Output, real ( kind = 8 ) XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, the 
!    maximum and minimum X, Y and Z atomic coordinates.
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm

  integer ( kind = 4 ) atom_num
  real ( kind = 8 ) coord(atom_max,3)
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  integer ( kind = 4 ) prtatm(maxres,mxpatm)
  character ( len = 3 ) resnam(atom_max)
  integer ( kind = 4 ) resnum(atom_max)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin

  atom_num = 0
  coord(1:atom_max,1:3) = 0.0D+00
  hiatom = 0
  hires = 0
  numchain = 0
  numline = 0
  numlost = 0
  numres = 0
  prtatm(1:maxres,1:mxpatm) = 0
  resnam(1:atom_max) = ' '
  resnum(1:atom_max) = 0
  xmax = 0.0D+00
  xmin = 0.0D+00
  ymax = 0.0D+00
  ymin = 0.0D+00
  zmax = 0.0D+00
  zmin = 0.0D+00

  return
end
subroutine pdb_print_coord ( coord, atom_max, atom_num, resnam, resnum )

!*****************************************************************************80
!
!! PDB_PRINT_COORD prints out the COORD array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) COORD(ATOM_MAX,3) the coordinates of atoms.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) ATOM_NUM, the maximum atom serial number 
!    read in.
!
!    Input, character ( len = 3 ) RESNAM(ATOM_MAX), the residue name for
!    each entry.
!
!    Input, integer ( kind = 4 ) RESNUM(ATOM_MAX), the index of the residue 
!    to which the atom belongs.
!
  implicit none

  integer ( kind = 4 ) atom_max

  integer ( kind = 4 ) atom
  integer ( kind = 4 ) atom_num
  real ( kind = 8 ) coord(atom_max,3)
  integer ( kind = 4 ) j
  character ( len = 3 ) resnam(atom_max)
  integer ( kind = 4 ) resnum(atom_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Atom  ResName   Residue   Atomic XYZ coordinates:'
  write ( *, '(a)' ) ' '
 
  do atom = 1, atom_num
    write ( *, '(1x,i4,6x,a3,4x,i4,2x,3g14.6)' ) &
      atom, resnam(atom), resnum(atom), coord(atom,1:3)
  end do

  return
end
subroutine pdb_print_prtatm ( maxres, mxpatm, numres, prtatm )

!*****************************************************************************80
!
!! PDB_PRINT_PRTATM prints out the PRTATM array.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Input, integer ( kind = 4 ) NUMRES, the number of residuals read in.  
!
!    Input, integer ( kind = 4 ) PRTATM(MAXRES,MXPATM), contains the indices of 
!    the atoms that make up each residue.  In particular, for the I-th residue, 
!    the following entries are significant:
!      PRTATM(I,1) - N, the nitrogen;
!      PRTATM(I,2) - CA, the C-Alpha;
!      PRTATM(I,3) - C, the carbon;
!      PRTATM(I,4) - O, the oxygen.
!      PRTATM(I,5) - CB, another carbon?
!      PRTATM(I,5+) - the other atoms.
!
  implicit none

  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm

  integer ( kind = 4 ) ires
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) numres
  integer ( kind = 4 ) prtatm(maxres,mxpatm)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRTATM entries:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Residue, PRTATM(I,*)'
  write ( *, '(a)' ) ' '

  do ires = 1, numres
    jhi = 0
    do j = 1, mxpatm
      if ( prtatm(ires,j) /= 0 ) then
        jhi = j
      end if
    end do
    write ( *, '(1x,19i5)' ) ires, prtatm(ires,1:jhi)
  end do

  return
end
subroutine pdb_read_atom ( coord, hiatom, hires, input_unit, atom_max, &
  maxres, mxpatm, atom_num, numchain, numline, numlost, numres, &
  prtatm, resnam, resnum, xmax, xmin, ymax, ymin, zmax, zmin )

!*****************************************************************************80
!
!! PDB_READ_ATOM reads the ATOM records in a PDB file.
!
!  Discussion:
!
!    The PDB file is presumed to have been opened by the user.
!
!    PDB_READ builds the COORD and PRTATM arrays from the ATOM data.
!
!  Format:
!
!    COLUMNS  DATA TYPE     FIELD       DEFINITION
!    --------------------------------------------------------------------------
!     1 -  6  Record name   "ATOM  "
!     7 - 11  Integer       serial      Atom serial number.
!    13 - 16  Atom          name        Atom name.
!    17       Character     altLoc      Alternate location indicator.
!    18 - 20  Residue name  resName     Residue name.
!    22       Character     chainID     Chain identifier.
!    23 - 26  Integer       resSeq      Residue sequence number.
!    27       AChar         iCode       Code for insertion of residues.
!    31 - 38  Real(8.3)     x           Orthogonal coordinates for X, Angstroms.
!    39 - 46  Real(8.3)     y           Orthogonal coordinates for Y, Angstroms.
!    47 - 54  Real(8.3)     z           Orthogonal coordinates for Z, Angstroms.
!    55 - 60  Real(6.2)     occupancy   Occupancy.
!    61 - 66  Real(6.2)     tempFactor  Temperature factor.
!    73 - 76  LString(4)    segID       Segment identifier, left-justified.
!    77 - 78  LString(2)    element     Element symbol, right-justified.
!    79 - 80  LString(2)    charge      Charge on the atom.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 June 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COORD(ATOM_MAX,3) the coordinates of atoms.
!
!    Output, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Output, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) INPUT_UNIT, the FORTRAN unit number associated
!    with the file.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Output, integer ( kind = 4 ) ATOM_NUM, the maximum atom serial number 
!    read in.
!
!    Output, integer ( kind = 4 ) NUMCHAIN, the number of chains.
!
!    Output, integer ( kind = 4 ) NUMLINE, the number of lines of text in 
!    the file.
!
!    Output, integer ( kind = 4 ) NUMLOST, the number of atom items lost.
!
!    Output, integer ( kind = 4 ) NUMRES, the number of residuals read in.  
!
!    Output, integer ( kind = 4 ) PRTATM(MAXRES,MXPATM), contains the indices 
!    of the atoms that make up each residue.  For the I-th residue, 
!    the following entries are significant:
!      PRTATM(I,1) - N, the nitrogen;
!      PRTATM(I,2) - CA, the C-Alpha;
!      PRTATM(I,3) - C, the carbon;
!      PRTATM(I,4) - O, the oxygen.
!      PRTATM(I,5) - CB, another carbon?
!      PRTATM(I,5+) - the other atoms.
!
!    Output, character ( len = 3 ) RESNAM(ATOM_MAX), the residue name for each
!    entry.
!
!    Output, integer ( kind = 4 ) RESNUM(ATOM_MAX), the index of the residue 
!    to which the atom belongs.
!
!    Output, real ( kind = 8 ) XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, the 
!    maximum and minimum X, Y and Z atomic coordinates.
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm

  character altloc
  integer ( kind = 4 ) atom
  character ( len = 4 ) atom_name
  integer ( kind = 4 ) atom_num
  character chains
  character ( len = 2 ) charge
  real ( kind = 8 ) coord(atom_max,3)
  character ( len = 2 ) element
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) ibase
  character icode
  integer ( kind = 4 ) input_unit
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iset
  integer ( kind = 4 ) j
  logical s_eqi1
  logical str_check1
  logical str_check2
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  real ( kind = 8 ) occ
  integer ( kind = 4 ) prtatm(maxres,mxpatm)
  character prvchn
  integer ( kind = 4 ) prvnum
  character ( len = 3 ) prvres
  real ( kind = 8 ) r1
  character ( len = 3 ) resnam(atom_max)
  character ( len = 3 ) resname
  integer ( kind = 4 ) resno
  integer ( kind = 4 ) resnum(atom_max)
  real ( kind = 8 ) rfactr
  character ( len = 4 ) segid
  character ( len = 128 ) string
  real ( kind = 8 ) temp
  character ( len = 4 ) w1
  real ( kind = 8 ) x
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) y
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) z
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin
!
!  Initialize PDB information.
!
  call pdb_init ( coord, hiatom, hires, atom_max, maxres, mxpatm, &
    atom_num, numchain, numline, numlost, numres, prtatm, resnam, &
    resnum, xmax, xmin, ymax, ymin, zmax, zmin )

  ibase = 0
  prvchn = ' '
  prvnum = 0
  prvres = ' '

  iset = 0

  do

    read ( input_unit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      exit
    end if

    numline = numline + 1
    !s_eqi (str_check1, string(1:6), 'ENDMDL' )  
 	!s_eqi (str_check2, string(1:4), 'ATOM' ) 
    if (string(1:6).eq. 'ENDMDL'  ) then

      exit

    else if (string(1:4).eq. 'ATOM' ) then

      read ( string, &
        '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)', &
        iostat = ios ) &
        w1, atom, atom_name, altloc, resname, chains, resno, &
        icode, x, y, z, occ, temp, segid, element, charge

      if ( ios /= 0 ) then
        exit
      end if

      atom_num = atom_num + 1
 
      hiatom = max ( hiatom, atom )
!
!  Remove a possible initial blank in ATOM_NAME or RESNAME.
!
      if ( atom_name(1:1) == ' ' ) then
        atom_name = atom_name(2:)
      end if
 
      if ( resname(1:1) == ' ' ) then
        resname = resname(2:)
      end if
!
!  If necessary, increment the number of residues read.
!
      if ( resno /= prvnum .or. resname /= prvres .or. chains /= prvchn ) then

        prvnum = resno

        prvres = resname

        if ( chains /= prvchn .and. 1 < atom ) then
          ibase = numres
          numchain = numchain + 1
        end if

        prvchn = chains

        numres = numres + 1

      end if
!
!  Correct the residue index by accounting for the chain it is on.
!
      resno = resno + ibase
!
!  Keep track of the highest residue index.
!
      hires = max ( hires, resno )
!
!  For each atom, store the atomic coordinates, and the name and number 
!  of the residue to which the atom belongs.
!
      if ( 1 <= atom .and. atom <= atom_max ) then

        coord(atom,1) = x
        coord(atom,2) = y
        coord(atom,3) = z

        resnam(atom) = resname
        resnum(atom) = numres
		
		write(*,*)coord(atom,:)

      end if

      if ( iset == 0 ) then
        xmax = x
        xmin = x
        ymax = y
        ymin = y
        zmax = z
        zmin = z
        iset = 1
      else
        xmax = max ( xmax, x )
        xmin = min ( xmin, x )
        ymax = max ( ymax, y )
        ymin = min ( ymin, y )
        zmax = max ( zmax, z )
        zmin = min ( zmin, z )
      end if
!
!  For each residue, add the atomic index to the PRTATM database.
!
      if ( 1 <= numres .and. numres <= maxres ) then
 
        if ( atom_name == 'C' ) then

          prtatm(numres,3) = atom

        else if ( atom_name == 'CA ' ) then

          prtatm(numres,2) = atom

        else if ( atom_name == 'CB' ) then

          prtatm(numres,5) = atom

        else if ( atom_name == 'N' ) then

          prtatm(numres,1) = atom

        else if ( atom_name == 'O' ) then

          prtatm(numres,4) = atom

        else

          numlost = numlost + 1

          do j = 6, mxpatm
            if ( prtatm(numres,j) == 0 ) then
              prtatm(numres,j) = atom
              numlost = numlost - 1
              exit
            end if
          end do

        end if

      end if
 
    end if
  
  end do
 
  return
end
subroutine pdb_summary ( hiatom, hires, atom_max, maxres, mxpatm, atom_num, &
  numchain, numline, numlost, numres, xmax, xmin, ymax, ymin, zmax, zmin )

!*****************************************************************************80
!
!! PDB_SUMMARY prints a summary of the data read from a PDB file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) HIATOM, the maximum atom index encountered.
!
!    Input, integer ( kind = 4 ) HIRES, the maximum residue index encountered.
!
!    Input, integer ( kind = 4 ) ATOM_MAX, the maximum number of atoms.
!
!    Input, integer ( kind = 4 ) MAXRES, the maximum number of residues.
!
!    Input, integer ( kind = 4 ) MXPATM, the maximum number of atoms that can be
!    associated with a residue.
!
!    Input, integer ( kind = 4 ) ATOM_NUM, the number of atoms read in.
!
!    Input, integer ( kind = 4 ) NUMCHAIN, the number of chains.
!
!    Input, integer ( kind = 4 ) NUMLINE, the number of lines of text in 
!    the file.
!
!    Input, integer ( kind = 4 ) NUMLOST, the number of atom items lost.
!
!    Input, integer ( kind = 4 ) NUMRES, the number of residuals read in.  
!
!    Input, real ( kind = 8 ) XMAX, XMIN, YMAX, YMIN, ZMAX, ZMIN, the 
!    maximum and minimum X, Y and Z atomic coordinates.
!
  implicit none

  integer ( kind = 4 ) atom_max
  integer ( kind = 4 ) atom_num
  integer ( kind = 4 ) hiatom
  integer ( kind = 4 ) hires
  integer ( kind = 4 ) maxres
  integer ( kind = 4 ) mxpatm
  integer ( kind = 4 ) numchain
  integer ( kind = 4 ) numline
  integer ( kind = 4 ) numlost
  integer ( kind = 4 ) numres
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) zmax
  real ( kind = 8 ) zmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PDB_SUMMARY:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Storage limits:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '    Maximum number of atoms         ', atom_max
  write ( *, '(a,i6)' ) '    Maximum number of residues      ', maxres
  write ( *, '(a,i6)' ) '    Maximum atoms per residue       ', mxpatm
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data from the file:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '    Number of text lines in file    ', numline
  write ( *, '(a,i6)' ) '    Highest residue index           ', hires
  write ( *, '(a,i6)' ) '    Number of residues              ', numres
  write ( *, '(a,i6)' ) '    Highest atom index              ', hiatom
  write ( *, '(a,i6)' ) '    Number of atoms                 ', atom_num
  write ( *, '(a,i6)' ) '    Number of chains                ', numchain
  write ( *, '(a,i6)' ) '    Number of atom items lost       ', numlost
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Atomic coordinate ranges, in Angstroms:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '    X: ', xmin, xmax
  write ( *, '(a,2g14.6)' ) '    Y: ', ymin, ymax
  write ( *, '(a,2g14.6)' ) '    Z: ', zmin, zmax

  return
end
subroutine pdb_write_atom ( atom_num, atom_xyz, pdb_file_unit )

!*****************************************************************************80
!
!! PDB_WRITE_ATOM writes the ATOM records in a PDB file.
!
!  Discussion:
!
!    This routine accepts a set of 3D coordinates, and writes them
!    out as plausible ATOM records of a PDB file.  
!
!    The reason for doing this is to allow the data to be viewed
!    by programs which can display configurations of atoms described by
!    a PDB file.
!
!    In order to write plausible ATOM records, all the fields must
!    be filled in.  This routine makes up plausible values for these
!    other fields.
!
!  PDB ATOM Record Format:
!
!    COLUMNS  DATA TYPE     FIELD       DEFINITION
!    --------------------------------------------------------------------------
!     1 -  6  Record name   "ATOM  "
!     7 - 11  Integer       serial      Atom serial number.
!    13 - 16  Atom          name        Atom name.
!    17       Character     altLoc      Alternate location indicator.
!    18 - 20  Residue name  resName     Residue name.
!    22       Character     chainID     Chain identifier.
!    23 - 26  Integer       resSeq      Residue sequence number.
!    27       AChar         iCode       Code for insertion of residues.
!    31 - 38  Real(8.3)     x           Orthogonal coordinates for X, Angstroms.
!    39 - 46  Real(8.3)     y           Orthogonal coordinates for Y, Angstroms.
!    47 - 54  Real(8.3)     z           Orthogonal coordinates for Z, Angstroms.
!    55 - 60  Real(6.2)     occupancy   Occupancy.
!    61 - 66  Real(6.2)     tempFactor  Temperature factor.
!    73 - 76  LString(4)    segID       Segment identifier, left-justified.
!    77 - 78  LString(2)    element     Element symbol, right-justified.
!    79 - 80  LString(2)    charge      Charge on the atom.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 January 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ATOM_NUM, the maximum atom serial number read.
!
!    Input, real ( kind = 8 ) ATOM_XYZ(3,ATOM_NUM) the coordinates of atoms.
!
!    Input, integer ( kind = 4 ) PDB_FILE_UNIT, the FORTRAN unit number 
!    associated with the file.
!
  implicit none

  integer ( kind = 4 ) atom_num

  character altloc
  integer ( kind = 4 ) atom
  character ( len = 4 ) atom_name
  real ( kind = 8 ) atom_xyz(3,atom_num)
  character chains
  character ( len = 2 ) charge
  character ( len = 2 ) element
  character icode
  real ( kind = 8 ) occ
  integer ( kind = 4 ) pdb_file_unit
  character ( len = 3 ) resname
  integer ( kind = 4 ) resno
  character ( len = 4 ) segid
  real ( kind = 8 ) temp

  atom_name = 'H   '
  altloc = ' '
  resname = '   '
  chains = ' '
  resno = 0
  icode = ' '
  occ = 1.0D+00
  temp = 1.0D+00
  segid = '    '
  element = '  '
  charge = '  '

  do atom = 1, atom_num

    write ( pdb_file_unit, &
      '(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)' ), &
      'ATOM  ', atom, atom_name, altloc, resname, chains, resno, &
      icode, atom_xyz(1:3,atom), occ, temp, segid, element, charge

  end do

  return
end
subroutine s_eqi (s_eqi1, strng1, strng2 )

!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) STRNG1, STRNG2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi1
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2

  len1 = len ( strng1 )
  len2 = len ( strng2 )
  lenc = min ( len1, len2 )

  s_eqi1 = .false.

  do i = 1, lenc

    s1 = strng1(i:i)
    s2 = strng2(i:i)
    call ch_cap ( s1 )
    call ch_cap ( s2 )

    if ( s1 /= s2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( strng1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( strng2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi1 = .true.

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 ) time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         writeSPIDER
  !> @brief             takes filename, and Image in the image plane and write the spider image 
  !!			the spider image (ex. Image.spi) will be displayed by Xmipp tool     	
  !!			writeSPIDER writes data to a file in SPIDER format.
  !!			2D data are written as a SPIDER image.
  !!			3D data may be written as either a SPIDER volume (default) or
  !!			as a stack of 2D images (requires optional 3rd argument='stack').
  !!			SPIDER files do not have a specific file extension.
  !!			Examples
  !!			writeSPIDERfile('img001.dat', I)
  !!			writeSPIDERfile('abc001.stk', V, 'stack')
  !!			version 1.0 (Mar 2021) A. Mirzaei
  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         filename        : filename for instance image.spi
  !! @param[in]         I               : Image data (image peixel in i and j coordiante)
  !! @param[in]      	varargin        : additional agrument for stocking
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine writeSPIDERfile(filenameS, I, varargin)
 	character(128),         intent(in)    :: filenameS
	character(128),		intent(in)    :: varargin
	real(wp),		intent(in)    :: I(:,:)

	integer                               :: stackarg
	integer                               :: write_stack
	integer                               :: datasize(3)
	integer                               :: stat,fp,ndim,nimages,nsam,nrow,i5
	integer                               :: header(27),hdr(27)
	integer                               :: imgsize(3)
	real(wp),pointer                      :: imgdata(:,:)

	if (varargin =='stack') then
	    stackarg = 0;      ! make a stack 
	    write_stack = 1;
	else
 	   stackarg = -1;     ! not a stack
	   write_stack = 0;
	endif
	datasize(1) = size(I,1);
	datasize(2) = size(I,2);
	datasize(3) = 1;

	call makeSPIDERheader(datasize, stackarg, header)
	if (header(1) < 0) then
	    write(*,*)"Unable to create header"
	    return
	endif
	fp = 27
	OPEN (fp, file = filenameS, status='old', action='write',iostat=stat)
	if (stat .ne. 0)then
	  write(*,*)"File cannot be opened"
	  return
	endif
	ndim = datasize(3);
	if (ndim == 1) then
	    write(fp,*)header
	    write(fp,*)I
	elseif (ndim == 3) then
	   ! write the volume header or overall stack header
	   write(fp,*)header
   
	   nimages = datasize(3)
	   nsam = datasize(1)
	   nrow = datasize(2)
	   allocate(imgdata(nsam, nrow))
	   imgsize = [nsam, nrow, 1]
   
	   do i5 = 1,nimages
	       if (write_stack==1) then
	          call makeSPIDERheader(imgsize, i5,hdr)
	          write(fp,*)hdr
	       endif
	       imgdata = I(:,:)
	       write(fp,*)imgdata
	   enddo 
	elseif (ndim > 3 .or. ndim < 1) then
	    write(*,*)"output data must be 2 or 3 dimensions"
	endif
	close(fp)	
	
  end subroutine writeSPIDERfile
   !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine         makeSPIDERheader
  !> @brief             takes datasize and stacking as input and generate header file for spi image 
  !!			the spider image (ex. Image.spi) is written using writeSPIDER subroutine
  !!			main use: Create a header for a SPIDER file.
  !!			Stack files present a special case. They have an overall
  !!			header, plus each image in the stack has its own header.
  !!			In the overall header, istack > 0, maxim = total # of images. 
  !!			In individual stack images, imgnum > 0
  !!			stackarg: < 0 : not a stack
  !!			          = 0 : make overall stack header
  !!			          > 0 : imgnum in individual stack image
  !!			version 1.0 (Mar 2021) A. Mirzaei

  !! @authors           Alex Mirzaei, IMPMC pars
  !! @param[in]         datasize        : size of the image data
  !! @param[in]         stackarg        : stacking defined by the write function and imposed by the data type
  !! @param[inout]      header          : spi file header generated for spi write function
  ! 
  !======1=========2=========3=========4=========5=========6=========7=========8
  subroutine makeSPIDERheader(datasizeS, stackargS, headerS)
 	integer,        	intent(in)    	:: datasizeS(3)
	integer,		intent(in)    	:: stackargS
	integer,		intent(inout)   :: headerS(27)

	integer                                 ::n,nsam, nrow,nslice, iform,maxim, imgnum,istack
	integer                                 ::labrec,lenbyt,nvalues,labbyt 
	n = datasizeS(3)
	write(*,*)"image type and data size",n
	nsam = datasizeS(2)   ! switch due to row, column indexing
	nrow = datasizeS(1)
	if (n == 1 .and. stackargS < 0)then  ! 2D image
	    nslice = 1
	    iform = 1
	elseif (n == 3 .and. stackargS < 0)then ! 3D volume
	    nslice = datasizeS(3)
	    iform = 3
	elseif (stackargS == 0)  then  ! overall stack file header
	    nslice = 1
	    iform = 1
	    istack = 2
	    maxim = datasizeS(3)
 	   imgnum = 0
	elseif (stackargS > 0) then   ! image within a stack file
	    nslice = 1
	    iform = 1
	    istack = 0
	    maxim = 0
	    imgnum = stackargS
	endif
	lenbyt = nsam * 4
	labrec = INT(1024 / lenbyt)
	if (mod(1024,lenbyt) /= 0) then
	    labrec = labrec + 1
	endif
	labbyt = labrec * lenbyt
	nvalues = INT(labbyt/4)
	if (nvalues < 256) then
	    headerS(1) = -1
	    return
	endif
	headerS = 0;
	headerS(1)  = nslice 	! nslice (=1 for an image) 
	headerS(2)  = nrow    	! number of rows per slice
	headerS(5)  = iform  	! iform: 1 for 2D image, 3 for volume
	headerS(12) = nsam    	! number of pixels per line
	headerS(13) = labrec  	! number of records in file header
	headerS(22) = labbyt  	! total number of bytes in header
	headerS(23) = lenbyt  	! record length in bytes
	if (stackargS >= 0) then
	   headerS(24) = istack  ! > 0 in overall stack file header
	   headerS(26) = maxim   ! total number of images in a stack file
	   headerS(27) = imgnum  ! current image in a stack file
	endif
  end subroutine makeSPIDERheader
end module fileio_sit_mod
