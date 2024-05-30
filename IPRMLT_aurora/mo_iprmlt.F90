module mo_iprmlt
!-------------------------------------------------------------------
!  module to read in ion pair produciton rates that are a function of
!  magnetic local time (MLT), L-shell and pressure level and transform them
!  onto the WACCM geographic grid 
! 
!  initial coding: DRM & PTV, October, 2017
!-------------------------------------------------------------------

  use shr_kind_mod,   only : r8 => shr_kind_r8
  use cam_abortutils, only : endrun
  use cam_logfile,    only : iulog
  use ppgrid,         only : pcols, pver
  use cam_history,    only : outfld, addfld
  use spmd_utils,     only : masterproc

  implicit none

  private

  public :: iprmlt_readnl
  public :: iprmlt_init
  public :: iprmlt_adv
  public :: iprmlt_ionization_noxhox

  save

  character(len=32)  :: specifier(1) = 'ipr'
  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0

  ! ... time related variables from file
  integer, allocatable :: datesec(:)
  integer, allocatable :: datearr(:) ! in case input includes array of dates, only first one will be used
  integer              :: filedate = 0

  logical :: has_iprmlt_ionization = .false.


  !-------------------------------------------------------------------
  !	... ipr data from netcdf file 
  !-------------------------------------------------------------------
  logical               :: has_dim_lshell = .false.        ! does the file include L-shell dimension (if not, geomagnetic latitude 'glat')
  character(len=256)    :: dim_y_name = ''                 ! name of the ipr y dimension, either 'lshell' or 'glat'
  character(len=256)    :: dim_y_bnds_name = ''            ! name of the ipr y dimension bounds, either 'lshell_width' or 'glat_bnds'
  real(r8), allocatable :: dim_y(:)                        ! ipr y dimension, either magnetic latitude in L-shell (dimensionless) or geomagnetic latitude
  real(r8), allocatable :: dim_y_bnds(:,:)                 ! ipr y dimension bounds, either magnetic latitude L-shell widths or geomagnetic latitude bounds
  real(r8), allocatable :: lshell(:)                       ! ipr magnetic latitude in L-shell (dimensionless)
  real(r8), allocatable :: lshell_width(:)                 ! ipr magnetic latitude L-shell widths 
  real(r8), allocatable :: mlt(:)                          ! ipr magnetic local time bins (sec)
  real(r8), allocatable :: mlt_width(:)                    ! ipr magnetic local time bin widths (sec)
  real(r8), allocatable :: lev(:)                          ! ipr pressure levels (hPa)
  real(r8), allocatable :: ipr(:,:,:,:)                    ! ipr (/cm^3/s), dimensions: MLT x L-shell x lev (x time)
  character(len=256)    :: ipr_units                       ! units of the ipr input: if /g/s, will be converted to /cm^3/s

  integer :: nlev, nydim, nshell, nmlt, ntime

contains

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine iprmlt_readnl(nlfile)
    !-------------------------------------------------------------------
    !	... routine to get IPR filename, etc. from namelist
    !-------------------------------------------------------------------

    use namelist_utils, only : find_group_name
    use units,          only : getunit, freeunit
    use mpishorthand

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    ! define a string with the name of this routine for error diagnostics
    character(len=*), parameter :: subname = 'iprmlt_readnl' 

    character(len=16)  ::  iprmlt_ionization_fldname
    character(len=256) ::  iprmlt_ionization_filename
    character(len=256) ::  iprmlt_ionization_datapath
    character(len=256) ::  iprmlt_ionization_filelist
    character(len=32)  ::  iprmlt_ionization_datatype
    integer            ::  iprmlt_ionization_cycle_yr
    integer            ::  iprmlt_ionization_fixed_ymd
    integer            ::  iprmlt_ionization_fixed_tod

    namelist /iprmlt_ionization_nl/ &
         iprmlt_ionization_fldname, &
         iprmlt_ionization_filename, &
         iprmlt_ionization_datapath, &
         iprmlt_ionization_filelist, &
         iprmlt_ionization_datatype, &
         iprmlt_ionization_cycle_yr, &
         iprmlt_ionization_fixed_ymd, &
         iprmlt_ionization_fixed_tod

    ! copy defaults into namelist

    iprmlt_ionization_fldname = specifier(1)
    iprmlt_ionization_filename = filename
    iprmlt_ionization_datapath = datapath
    iprmlt_ionization_filelist = filelist
    iprmlt_ionization_datatype = datatype
    iprmlt_ionization_cycle_yr = cycle_yr
    iprmlt_ionization_fixed_ymd = fixed_ymd
    iprmlt_ionization_fixed_tod = fixed_tod

    ! update namelist settings from namelist file
    if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'iprmlt_ionization_nl', status=ierr)
      if (ierr == 0) then
        read(unitn, iprmlt_ionization_nl, iostat=ierr)
        if (ierr /= 0) then
          call endrun(subname // ':: ERROR reading namelist')
        end if
      end if
      close(unitn)
      call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(iprmlt_ionization_fldname,  len(iprmlt_ionization_fldname),  mpichar, 0, mpicom)
    call mpibcast(iprmlt_ionization_filename, len(iprmlt_ionization_filename), mpichar, 0, mpicom)
    call mpibcast(iprmlt_ionization_filelist, len(iprmlt_ionization_filelist), mpichar, 0, mpicom)
    call mpibcast(iprmlt_ionization_datapath, len(iprmlt_ionization_datapath), mpichar, 0, mpicom)
    call mpibcast(iprmlt_ionization_datatype, len(iprmlt_ionization_datatype), mpichar, 0, mpicom)
    call mpibcast(iprmlt_ionization_cycle_yr, 1, mpiint,  0, mpicom)
    call mpibcast(iprmlt_ionization_fixed_ymd,1, mpiint,  0, mpicom)
    call mpibcast(iprmlt_ionization_fixed_tod,1, mpiint,  0, mpicom)
#endif

    ! Update module variables with user settings.
    specifier(1) = iprmlt_ionization_fldname
    filename  = iprmlt_ionization_filename
    filelist  = iprmlt_ionization_filelist
    datapath  = iprmlt_ionization_datapath
    datatype  = iprmlt_ionization_datatype
    cycle_yr  = iprmlt_ionization_cycle_yr
    fixed_ymd = iprmlt_ionization_fixed_ymd
    fixed_tod = iprmlt_ionization_fixed_tod

    ! Turn on MLT dependent IPR if user has specified an input dataset.
    if (len_trim(filename) > 0  .and. filename.ne.'NONE') has_iprmlt_ionization = .true.

  end subroutine iprmlt_readnl

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine iprmlt_init(first_call)
    !-------------------------------------------------------------------
    !	... read in MLT dependent ion production rates from netcdf file
    !-------------------------------------------------------------------

    use ioFileMod,     only : getfil
    use cam_pio_utils, only : cam_pio_openfile
    use pio,           only : file_desc_t, pio_get_var, pio_closefile, &
         pio_nowrite, pio_inq_varid, pio_inq_dimid, pio_inq_dimlen, pio_get_att, &
         pio_seterrorhandling, pio_bcast_error, pio_internal_error, pio_noerr
    use time_manager,  only : get_curr_date
    use cam_history,   only : horiz_only
    use shr_file_mod,  only : shr_file_getunit, shr_file_freeunit
    use mo_constants,  only : d2r

    implicit none

    ! date components             
    integer :: yr,    & ! year
               mon,   & ! month
               day,   & ! day of month
               tod      ! time of day (seconds past 0Z)

    logical, intent(in) :: first_call

    !-------------------------------------------------------------------
    !	... local variables
    !-------------------------------------------------------------------
    integer :: unitnumber, istat, ierr
    integer :: dim_id, var_id 
    type(file_desc_t) :: ncid
    character(len=256) :: filelistpath = ''
    character(len=256) :: filepath = ''
    character(len=256) :: locfn
    character(len=256) :: file_to_read
    character(len=256) :: line
    integer :: yrmonday = 0

    if (.not.has_iprmlt_ionization) return
    

    call get_curr_date(yr, mon, day, tod)
    yrmonday = 10000*yr + 100*mon + day
    if( masterproc ) then
      write(iulog,*) 'iprmlt_init: yr, mon, day, tod', yr, mon, day, tod
    end if


    !-------------------------------------------------------------------
    !     Read dimensions and allocate variables from the initial file
    !-------------------------------------------------------------------
    if ( first_call ) then

      !-------------------------------------------------------------------
      !	... open the initial netcdf file
      !-------------------------------------------------------------------
      if (len_trim(datapath) > 0 .and. datapath.ne.'NONE') then
        filepath = trim(datapath) // '/' // trim(filename)
      else
        filepath = trim(filename)
      end if
      call getfil(filepath, locfn, 0)
      call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE)

      !-------------------------------------------------------------------
      !	... check if file has L-shell dimension
      !-------------------------------------------------------------------
      call pio_seterrorhandling( ncid, PIO_BCAST_ERROR )
      ierr = pio_inq_dimid( ncid, 'lshell', dim_id )
      has_dim_lshell = (ierr==PIO_NOERR)
      call pio_seterrorhandling( ncid, PIO_INTERNAL_ERROR )

      if (masterproc) then
        write(iulog,*) 'iprmlt_init: has_dim_lshell', has_dim_lshell
      end if
      if ( has_dim_lshell ) then
        dim_y_name = 'lshell'
        dim_y_bnds_name = 'lshell_width'
      else
        dim_y_name = 'glat'
        dim_y_bnds_name = 'glat_bnds'
      end if

      !-------------------------------------------------------------------
      !	... read the dimensions
      !-------------------------------------------------------------------
      ierr = pio_inq_dimid( ncid, 'lev', dim_id )
      ierr = pio_inq_dimlen( ncid, dim_id, nlev )
      ierr = pio_inq_dimid( ncid, trim( dim_y_name ), dim_id )
      ierr = pio_inq_dimlen( ncid, dim_id, nydim )
      ierr = pio_inq_dimid( ncid, 'mlt', dim_id )
      ierr = pio_inq_dimlen( ncid, dim_id, nmlt )
      ierr = pio_inq_dimid( ncid, 'time', dim_id )
      ierr = pio_inq_dimlen( ncid, dim_id, ntime )

      !-------------------------------------------------------------------
      !	... allocate ipr variables
      !-------------------------------------------------------------------
      allocate( lev(nlev), dim_y(nydim), mlt(nmlt), mlt_width(nmlt), & 
            datearr(ntime), datesec(ntime), stat=istat )
      if( istat /= 0 ) then
        write(iulog,*) 'iprmlt_init: failed to allocate dimension arrays, error = ',istat
        call endrun
      end if

      if ( has_dim_lshell ) then
        allocate( dim_y_bnds(nydim,1), stat=istat )
        if( istat /= 0 ) then
          write(iulog,*) 'iprmlt_init: failed to allocate dimension array, error = ',istat
          call endrun
        end if
      else
        allocate( dim_y_bnds(nydim,2), stat=istat )
        if( istat /= 0 ) then
          write(iulog,*) 'iprmlt_init: failed to allocate dimension array, error = ',istat
          call endrun
        end if
      end if

      allocate( ipr(nlev,nydim,nmlt,ntime), stat=istat )
      if( istat /= 0 ) then
        write(iulog,*) 'iprmlt_init: failed to allocate ipr array, error = ',istat
        call endrun
      end if

      !-------------------------------------------------------------------
      !	... update file date
      !-------------------------------------------------------------------
      ierr = pio_inq_varid( ncid, 'date', var_id )
      ierr = pio_get_var( ncid, var_id, datearr )
      if( ierr /= 0 ) then
        write(iulog,*) 'iprmlt_init: failed to allocate filedate, error = ',ierr
        call endrun
      end if
      ! if there are multiple date values, only first will be used
      filedate = datearr(1)

      !	... close the netcdf file
      call pio_closefile( ncid )

    end if


    !-------------------------------------------------------------------
    !     Check that dates match!
    !-------------------------------------------------------------------
    if ( yrmonday .ne. filedate ) then
      write(iulog,*) ' '
      write(iulog,*) '***************'
      write(iulog,*) 'iprmlt_init: detected mismatching dates'
      write(iulog,*) 'iprmlt_init: model date: ', yrmonday
      write(iulog,*) 'iprmlt_init: file date: ', filedate
      write(iulog,*) 'iprmlt_init: attempting to find correct IPRMLT forcing file'
      write(iulog,*) '***************'
      write(iulog,*) ' '

      !-------------------------------------------------------------------
      !     Open the filelist of IPRMLT input.
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! ... read the input filelist (ASCII)
      !-------------------------------------------------------------------
      if ( len_trim(datapath) > 0 .and. filename.ne.'NONE' ) then
        filelistpath = trim(datapath) //'/'// trim(filelist)
      else
        filelistpath = trim(filelist)
      endif

      unitnumber = shr_file_getUnit()
      open(unit=unitnumber,FILE=filelistpath,iostat=ierr,STATUS='OLD',ACTION='READ')
      if (ierr /= 0) then
        call endrun('iprmlt_init: not able to open file: '//trim(filelistpath))
      endif
      if( masterproc ) then
        write(iulog,*) 'iprmlt_init: opened list file, unit', unitnumber
      end if

      !-------------------------------------------------------------------
      !	... read first line
      !-------------------------------------------------------------------
      read( unit=unitnumber, fmt='(A)' ,iostat=ierr ) line
      if (ierr /= 0) then
        call endrun('iprmlt_init: not able to read first file name from filelist: '//trim(filelist))
      end if

      !-------------------------------------------------------------------
      !	    If file date is (somehow) ahead of model date, start from the
      !     first file on the list
      !-------------------------------------------------------------------
      if ( yrmonday .lt. filedate ) then
        file_to_read = trim( line )
      else

        !-------------------------------------------------------------------
        !	... otherwise skip to current file
        !-------------------------------------------------------------------
        do while( trim( line ) /= trim( filename ) )
          read( unit=unitnumber, fmt='(A)', iostat=ierr ) line
          if (ierr /= 0) then
            call endrun('iprmlt_init: not able to increment file name from filelist: '//trim(filelist))
          end if
        end do

        !-------------------------------------------------------------------
        !	... and read next line and assign
        !-------------------------------------------------------------------
        read( unit=unitnumber, fmt='(A)') line
        file_to_read = trim( line )
      end if

      !-------------------------------------------------------------------
      !	    Loop until the correct date is found
      !-------------------------------------------------------------------
      do
        if (len_trim(datapath) > 0 .and. datapath.ne.'NONE') then
          filepath = trim(datapath) // '/' // trim(file_to_read)
        else
          filepath = trim(file_to_read)
        end if
        call getfil(filepath, locfn, 0)
        call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE)
        !-------------------------------------------------------------------
        !	... read file date
        !-------------------------------------------------------------------
        ierr = pio_inq_varid( ncid, 'date', var_id )
        ierr = pio_get_var( ncid, var_id, datearr )
        if( ierr /= 0 ) then
          write(iulog,*) 'iprmlt_init: failed to read filedate, error = ',ierr
          call endrun
        end if
        ! if there are multiple date values, only first will be used
        filedate = datearr(1)

        !	... close the netcdf file
        call pio_closefile( ncid )

        if  ( yrmonday == filedate ) then
          !-------------------------------------------------------------------
          !	... exit if dates match
          !-------------------------------------------------------------------
          write(iulog,*) 'iprmlt_init: found correct file for model date ', yrmonday
          write(iulog,*) 'iprmlt_init: new file to read: ', trim( file_to_read )
          write(iulog,*) 'iprmlt_init: file date: ', filedate
          write(iulog,*) '***************'
          write(iulog,*) ' '
          filename = trim( file_to_read )
          exit
        else
          !-------------------------------------------------------------------
          !	... else move on to next file on list
          !-------------------------------------------------------------------
          read( unit=unitnumber, fmt='(A)', iostat=ierr ) line
          if (ierr /= 0) then
            call endrun('iprmlt_init: not able to increment file name from filelist: '//trim(filelist))
          end if
          file_to_read = trim(line)
        end if

      end do
      !-------------------------------------------------------------------
      !	... close the filelist
      !-------------------------------------------------------------------
      close(unit=unitnumber)
      call shr_file_freeUnit(unitnumber)

    end if


    !-------------------------------------------------------------------
    !	... open the (correct) netcdf file
    !-------------------------------------------------------------------
    if (len_trim(datapath) > 0 .and. datapath.ne.'NONE') then
      filepath = trim(datapath) // '/' // trim(filename)
    else
      filepath = trim(filename)
    end if
    call getfil(filepath, locfn, 0)
    call cam_pio_openfile( ncid, trim(locfn), PIO_NOWRITE)


    !-------------------------------------------------------------------
    !	... update file date
    !-------------------------------------------------------------------
    ierr = pio_inq_varid( ncid, 'date', var_id )
    ierr = pio_get_var( ncid, var_id, datearr )
    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate filedate, error = ',ierr
      call endrun
    end if
    ! if there are multiple date values, only first will be used
    filedate = datearr(1)

    !-------------------------------------------------------------------
    !	... read the ipr variables
    !-------------------------------------------------------------------

    ierr = pio_inq_varid( ncid, trim( dim_y_name ), var_id )
    ierr = pio_get_var( ncid, var_id, dim_y )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate dimension y array, error = ',ierr
      call endrun
    end if

    ierr = pio_inq_varid( ncid, trim( dim_y_bnds_name ), var_id )
    ierr = pio_get_var( ncid, var_id, dim_y_bnds )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate dimension y bounds array, error = ',ierr
      call endrun
    end if

    ! convert glat from deg to rad
    if ( .not. has_dim_lshell ) then
      dim_y = dim_y * d2r
      dim_y_bnds = dim_y_bnds * d2r
    end if

    ierr = pio_inq_varid( ncid, 'mlt', var_id )
    ierr = pio_get_var( ncid, var_id, mlt )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate mlt array, error = ',ierr
      call endrun
    end if
  
    ierr = pio_inq_varid( ncid, 'mlt_width', var_id )
    ierr = pio_get_var( ncid, var_id, mlt_width )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate mlt_width array, error = ',ierr
      call endrun
    end if

    ierr = pio_inq_varid( ncid, 'lev', var_id )
    ierr = pio_get_var( ncid, var_id, lev )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate lev array, error = ',ierr
      call endrun
    end if

    ierr = pio_inq_varid( ncid, 'datesec', var_id )
    ierr = pio_get_var( ncid, var_id, datesec )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate datesec array, error = ',ierr
      call endrun
    end if

    ierr = pio_inq_varid( ncid, specifier(1), var_id ) 
    ierr = pio_get_var( ncid, var_id, (/1,1,1,1/), (/nlev, nydim, nmlt, ntime/), ipr )

    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to allocate ipr array, error = ',ierr
      call endrun
    end if

    !----------------------------------------------------------------------------
    ! CHECK: LET'S SEE THE UNITS, ALSO TEST
    !----------------------------------------------------------------------------
    ierr = pio_get_att( ncid, var_id, 'units', ipr_units)
    if( ierr /= 0 ) then
      write(iulog,*) 'iprmlt_init: failed to read ipr units = ',ierr
      call endrun
    end if
    if( masterproc ) then
	    write(iulog,*) ' '
	    write(iulog,*) 'iprmlt_init: ipr units:', trim( ipr_units )
      write(iulog,*) ' '
    end if


    !	... close the netcdf file
    call pio_closefile( ncid)

    if( masterproc ) then
	    write(iulog,*) ' '
	    write(iulog,*) 'iprmlt_init: filelist', filelist
      write(iulog,*) 'iprmlt_init: filename', filename
      write(iulog,*) 'iprmlt_init: filedate', filedate
	    write(iulog,*) 'iprmlt_init: dim_y', trim( dim_y_name ), dim_y
      ! write(iulog,*) 'iprmlt_init: mlt', mlt
      ! write(iulog,*) 'iprmlt_init: mlt_width', mlt_width
      ! write(iulog,*) 'iprmlt_init: dim_y_bnds', trim( dim_y_bnds_name ), dim_y_bnds
      ! write(iulog,*) 'iprmlt_init: lev', lev
      ! write(iulog,*) 'iprmlt_init: datesec', datesec
      ! write(iulog,*) 'iprmlt_init: ipr (lev,dim_y,mlt,datesec) (10,24,1:8,1)', ipr(10,24,1:8,1)
      ! write(iulog,*) 'iprmlt_init: ipr (lev,lshell,mlt)', ipr
	    write(iulog,*) ' '
    end if

    if (first_call) then
      if ( has_dim_lshell ) then
        call addfld('LSHELL', horiz_only, 'I', 'Earth radii', 'L shell' )
        call addfld('LSHELLI', horiz_only, 'I', ' ', 'L shell index' )
      end if
      call addfld('MLT', horiz_only, 'I', 'hour', 'Magnetic local time' )
      call addfld('MLTI', horiz_only, 'I', ' ', 'MLT index' )
      call addfld('IPRMLT', (/ 'lev' /), 'I', 'ions/cm^3/s', 'MLT dependent ionization' )
      ! for testing
      if ( index( trim(ipr_units), 'g^-1' ) .gt. 0 .OR. index( trim(ipr_units), '/g' ) .gt. 0) then
        call addfld('IPRgs', (/ 'lev' /), 'I', 'ions/g/s', 'MLT dependent ionization' )
      end if
    end if

  end subroutine iprmlt_init


  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine iprmlt_adv()
    use time_manager, only : get_curr_date

    ! date components
    integer :: yr,    & ! year
               mon,   & ! month
               day,   & ! day of month
               tod     ! time of day (seconds past 0Z)
    integer, parameter :: start_of_day = 0
    logical, parameter :: first_call = .false.

    if (.not.has_iprmlt_ionization) return

    ! Obtain date components valid at end of current timestep
    call get_curr_date(yr, mon, day, tod)

    if (tod.eq.start_of_day) then
      if( masterproc ) then
        write(iulog,*) 'iprmlt_adv: mon,day,tod', mon,day,tod
      end if

      call iprmlt_init(first_call)
    end if

  end subroutine iprmlt_adv

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine iprmlt_ionization_get( ncol, lchnk, pmid, temp, ionpairs )
    use physconst,    only : pi, rairv
    use mo_apex,      only : alatm, alonm   ! magnetic latitude grid (radians)
    use time_manager, only : get_curr_date


    integer, intent(in)   :: lchnk
    integer, intent(in)   :: ncol
    real(r8), intent(in)  :: pmid(:,:)  ! midpoint pressure (Pa) used for switching units if needed
    real(r8), intent(in)  :: temp(:,:)  ! midpoint temperature (K) used for switching units if needed

    real(r8), intent(out) :: ionpairs(:,:)

    real(r8) :: w_mlt(ncol)     ! magnetic local time at each WACCM gridpoint
    integer  :: w_mlti(ncol)    ! index of MLT array correstponding to WACCM mlt 
    !real(r8) :: w_lshell(ncol)  ! lshell value at each WACCM gridpoint
    !integer  :: w_lshelli(ncol) ! index of lshell array corresponding to WACCM mlat
    real(r8) :: w_ydim(ncol)    ! y dimension value at each WACCM gridpoint
    integer  :: w_ydimi(ncol)   ! index of y dimension array corresponding to WACCM mlat
    real(r8) :: mlat,mlon       ! local copy magnetic latitude and longitude in radians
    real(r8) :: mlond           ! local copy magnetic longitude in degrees
    ! for testing
    real(r8) :: ionpairs_gs(pcols,pver)


    ! Altitude at 90 km, in Earth Radii = (6371+90)/6371 = 1.0141.
    real(r8), parameter :: altre = 1.0141_r8
    real(r8), parameter :: s2h = 1._r8/3600._r8
    real(r8), parameter :: magplon = -72.8_r8 ! Longitude of the magnetic north pole in deg

    ! date components
    integer :: yr,    &! year
               mon,   &! month
    	         day,   &! day of month
    	         tod     ! time of day (seconds past 0Z)
    integer :: tod_ind = 1 ! index for datesec

    integer :: i, j, k

    ionpairs(:,:) = 0._r8
    ! for testing
    ionpairs_gs(:,:) = 0._r8

    if (.not.has_iprmlt_ionization) return

    ! Obtain date components valid at end of current timestep
    call get_curr_date(yr, mon, day, tod)

    ! Index for datesec dimension
    if ( ntime .ge. 2) then
      k = 1
      do while ( k .le. ntime )
        if ( tod .ge. datesec(k)) then
          tod_ind = k
        end if
        k = k + 1
      end do
    end if


    ! Calculate the magnetic local time at each location
    do i = 1,ncol
      mlon = alonm(i,lchnk)
      mlond = 180._r8*mlon/pi

      ! MLT calculated with longitude of the magnetic pole
      w_mlt(i) = (tod/3600._r8) + (mlond+magplon)/15._r8

      ! Making sure MLT is between 0 and 24
      w_mlt(i) = mod(w_mlt(i) + 24._r8, 24._r8)

      !w_mlt(i) = (tod/3600._r8) + (mlond+magplon)/15._r8
      !if ( w_mlt(i) .gt. 12._r8 ) then
        ! MLT going into the next day.
        !  w_mlt(i) = w_mlt(i) - 24._r8
      !end if
      ! Change into 0 - 24 range
      !w_mlt(i) = w_mlt(i) + 12._r8

      w_mlti(i) = -999

      ! index of MLT array corresponding to WACCM mlt
      ! Old style
      !do j = 1,nmlt
      !  if ( abs(w_mlt(i)-(mlt(j)*s2h)) .le. mlt_width(j)*0.5_r8*s2h ) then
      !    w_mlti(i) = j
      !    exit
      !  end if
      !enddo

      ! index of MLT array corresponding to WACCM mlt
      ! New style, without exit
      j = 1
      do while (j .le. nmlt )
        if ( abs(w_mlt(i)-(mlt(j)*s2h)) .le. mlt_width(j)*0.5_r8*s2h ) then
          w_mlti(i) = j
          j = nmlt
        end if
        j = j + 1
      end do
      ! If no index, check for cyclical (to account for the cyclic nature of MLT)
      if ( w_mlti(i) .lt. 0 ) then
        j = 1
        do while ( j .le. nmlt )
          if ( abs(w_mlt(i) - 24 -(mlt(j)*s2h)) .le. mlt_width(j)*0.5_r8*s2h ) then
            w_mlti(i) = j
            j = nmlt
          end if
          j = j + 1
        end do
      end if


      ! y dimension
      if ( has_dim_lshell ) then
        ! calculate L shell from magnetic latitude
        mlat = alatm(i,lchnk)
        w_ydim(i) = altre / (cos(mlat))**2
        if (w_ydim(i) .gt. 10.0_r8) then
          w_ydim(i) = 0.0_r8
        end if
        if (w_ydim(i) .lt. 2.0_r8) then
          w_ydim(i) = 0.0_r8
        end if

        w_ydimi(i) = -999

        ! index of L shell array corresponding to WACCM L shell
        ! Old style
        !do j = 1,nshell
        !  if ( abs(w_lshell(i)-lshell(j)) .le. lshell_width(j)*0.5_r8 ) then
        !    w_lshelli(i) = j
        !    exit
        !  end if
        !enddo

        ! index of L shell array corresponding to WACCM L shell
        ! New style, without exit
        j = 1
        do while ( j .le. nydim )
          if ( abs(w_ydim(i)-dim_y(j)) .le. dim_y_bnds(j,1)*0.5_r8 ) then
            w_ydimi(i) = j
            j = nydim
          end if
          j = j + 1
        end do
      else
        ! Geomagnetic latitude at grid points
        mlat = alatm(i,lchnk) ! radians

        w_ydimi(i) = -999
        j = 1
        do while ( j .le. nydim )
          if ( mlat .ge. dim_y_bnds(j,1) .and. mlat .lt. dim_y_bnds(j,2) ) then
            w_ydimi(i) = j
            j = nydim
          end if
          j = j + 1
        end do
      end if


      if (w_mlti(i) .gt. 0 .AND. w_ydimi(i) .gt. 0) then
        ionpairs(i,:) = ipr(:,w_ydimi(i),w_mlti(i),tod_ind)
        
        ! If ipr is in units /g/s, convert to /cm^3/s
        if ( index( trim(ipr_units), 'g^-1' ) .gt. 0 .OR. index( trim(ipr_units), '/g' ) .gt. 0) then
          ionpairs_gs(i,:) = ipr(:,w_ydimi(i),w_mlti(i),tod_ind)
          ! This conversion formula is directly copied from epp_ionization.F90
          ionpairs(i,:pver) = ionpairs(i,:pver) *(1.e-3_r8*pmid(i,:pver)/(rairv(i,:pver,lchnk)*temp(i,:pver)))
        end if
      end if
    enddo

    call outfld( 'MLT', w_mlt(:), ncol, lchnk )
    call outfld( 'MLTI', real(w_mlti(:),8) , ncol, lchnk )
    if ( has_dim_lshell ) then
      call outfld( 'LSHELL', w_ydim(:), ncol, lchnk )
      call outfld( 'LSHELLI', real(w_ydimi(:),8) , ncol, lchnk )
    end if
    call outfld( 'IPRMLT', ionpairs(:ncol,:pver), ncol, lchnk )

    ! for testing
    if ( index( trim(ipr_units), 'g^-1' ) .gt. 0 .OR. index( trim(ipr_units), '/g' ) .gt. 0) then
      call outfld( 'IPRgs', ionpairs_gs(:ncol,:pver), ncol, lchnk )
    end if

  end subroutine iprmlt_ionization_get

  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  subroutine iprmlt_ionization_noxhox( ncol, lchnk, zmid, pmid, temp, pr_nox, pr_hox )
    use spehox,  only : hox_prod_factor

    integer, intent(in)  :: ncol, lchnk
    real(r8), intent(in) :: zmid(:,:)	    ! altitude at mid-point used in HOX lookup table
    real(r8), intent(in) :: pmid(:,:)            ! midpoint pressure (Pa)
    real(r8), intent(in) :: temp(:,:)            ! midpoint temperature (K)

    real(r8), intent(out) :: pr_nox(:,:)  ! NOx production rate
    real(r8), intent(out) :: pr_hox(:,:)  ! HOx production rate

    real(r8) :: hoxprod_factor(pver)
    real(r8) :: ipr(pcols,pver)           ! ion-pair production rate

    integer :: i

    pr_nox(:,:) = 0._r8
    pr_hox(:,:) = 0._r8

    if (.not.has_iprmlt_ionization) return

    call  iprmlt_ionization_get( ncol, lchnk, pmid, temp, ipr )

    pr_nox(:ncol,:) = ipr(:ncol,:)

    do i = 1,ncol
      hoxprod_factor(:pver) = hox_prod_factor( ipr(i,:pver), zmid(i,:pver) )
      pr_hox(i,:pver) = hoxprod_factor(:pver) * ipr(i,:pver)
    end do

  end subroutine iprmlt_ionization_noxhox

end module mo_iprmlt
