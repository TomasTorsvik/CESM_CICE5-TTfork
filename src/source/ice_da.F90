!  SVN:$Id: ice_da.F90 2019-01-08 $
!=======================================================================
!
! Perform data assimilate for sea ice state, using
! 1) Combined Optimal Interpolation and Nudging (COIN) scheme
! authors: Keguang Wang, MET.no
!
! 2) Large revision to using the NorESM shared i/o
! infrastructure and implementation of nudging methods for
! pamip_short, and pamip_long
!
! Rev: 2020 Spring,  Jens Boldingh Debernard, met.no


module ice_da

! !USES:  
  use ice_kinds_mod
  use ice_constants, only: c0, c1, p1, p2, p01, p05, c10, secday, puny, &
       field_loc_center, field_type_scalar
  use ice_blocks, only: nx_block, ny_block
  use ice_domain_size, only: max_blocks, max_ntrcr
  use ice_domain_size, only: nilyr, nslyr, ncat 
  use ice_communicate, only: my_task, master_task, MPI_COMM_ICE
  use ice_broadcast
  use ice_fileunits  ! , only: nu_diag
  use ice_forcing, only: dbug
  use ice_timers, only: ice_timer_start, ice_timer_stop, timer_da
  use ice_calendar, only: istep, idate, new_day, yday, dt
  use ice_read_write, only: ice_open_nc, ice_read_nc, ice_close_nc

  !jd pio reading. What is nescessary?
  use shr_nl_mod, only : shr_nl_find_group_name
  use shr_strdata_mod
  use shr_dmodel_mod
  use shr_string_mod
  use shr_ncread_mod
  use shr_sys_mod
  use shr_mct_mod
  use mct_mod
  use pio
!jd 

      
  implicit none
  private
  public :: ice_da_init, ice_da_update, ice_da_prep, ice_da_init_streams
  save

  logical (kind=log_kind), public :: &
       da_ice ,       & ! perform data assimilation if true
       da_sic ,       & ! perform da of sic if true
       da_sit ,       & ! perform da of sea ice thickess if true
       da_sno ,       &  ! for snow depth if true
       da_update_ocn_heat, &
       da_update_ocn_mass 

  real (kind=dbl_kind), public :: da_timescale = c1   ! Nudging timescale in days
  
  character (char_len_long), public :: &
       da_method        ! data assimilation method
  
  character (char_len_long), public :: &
       da_data_dir       ! top directory for data to assimilate 

  !-----------------------------------------------------------------
  ! observed & model ice/snow variables & uncertainties
  !-----------------------------------------------------------------
  
  real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
       aice_obs     ,  & ! observed SIC 
       aice_obs_err ,  & ! observed SIC std
       aice_inc     ,  & ! Sic increment
       vice_obs     ,  & ! observed ice volume (m) 
       vice_obs_err ,  & ! observed ice volume std (m)
       vice_inc     ,  & ! Ice volume increment
       thice_obs     ,  & ! observed ice thickness (m) 
       thice_obs_err ,  & ! observed ice thickness std (m) 
       thice_inc     ,  & ! Ice thickness increment
       vsno_obs     ,  & ! observed snow volume (m) 
       vsno_obs_err ,  & ! observed snow volume std (m)
       vsno_inc       ,& ! Snow volume increment     
       thsnow_obs     ,  & ! observed snow thickness (m) 
       thsnow_obs_err ,  & ! observed snow thickness std (m)
       thsnow_inc          ! Snow thickness increment

!jd         trcr_obs, trcr_err, trcr_inc     ! tracers


  !jd Diagnostics of heat and freshwater fluxes due to 
  real (kind=dbl_kind), dimension (:,:,:), allocatable, public :: &
       da_fresh,   & ! Freshwater flux from nudging
       da_fsalt,   & ! Salt flux from nudging
       da_fheat      ! Heat flux loss from nudgning, more negative when ice
                     ! mass increase


  real (kind=dbl_kind), parameter :: min_aice_obs = p05
  real (kind=dbl_kind), parameter :: max_aice_obs = c1

  
      !jd For data stream
  integer(SHR_KIND_IN),parameter :: nFilesMaximum = 400 ! max number of files
  integer(kind=int_kind)         :: da_stream_year_first   ! first year in stream to use
  integer(kind=int_kind)         :: da_stream_year_last    ! last year in stream to use
  integer(kind=int_kind)         :: da_model_year_align    ! align stream_year_first
                                                        ! with this model year

  character(len=char_len_long)   :: da_stream_fldVarName
  character(len=char_len_long)   :: da_stream_fldFileName(nFilesMaximum)
  character(len=char_len_long)   :: da_stream_domTvarName
  character(len=char_len_long)   :: da_stream_domXvarName
  character(len=char_len_long)   :: da_stream_domYvarName
  character(len=char_len_long)   :: da_stream_domAreaName
  character(len=char_len_long)   :: da_stream_domMaskName
  character(len=char_len_long)   :: da_stream_domFileName
  character(len=char_len_long)   :: da_stream_mapread
  character(len=char_len_long)   :: da_stream_fillalgo
  logical(kind=log_kind)         :: da_ice_fill        ! true if data fill required

  type(shr_strdata_type)       :: sdat         ! prescribed data stream
  character(len=char_len_long) :: fldList      ! list of fields in data stream
      
      
!=======================================================================

contains

!=======================================================================

!  Allocates and initializes arrays needed for data assimilation, & 
!  read observation data.


  subroutine ice_da_init

    use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
    use ice_constants, only: c0,c1,c2
    use ice_grid, only: TLAT, tmask

    call ice_timer_start(timer_da)

!-----------------------------------------------------------------------
!     allocate variables for sea ice observations
!-----------------------------------------------------------------------

    if (trim(da_method)=='pamip_short') then
       da_sic=.true.
       da_sit=.true.
    else if (trim(da_method)=='pamip_long') then
       da_sic=.true.
       da_sit=.false.
    endif

    if (da_sit) da_sic=.true.
    
    if (da_sic) then
       allocate (aice_obs(nx_block,ny_block,max_blocks), &
            aice_obs_err(nx_block,ny_block,max_blocks), &
            aice_inc(nx_block,ny_block,max_blocks),  &
            ! Increment for volumes are always used
            vice_inc(nx_block,ny_block,max_blocks),  & 
            vsno_inc(nx_block,ny_block,max_blocks) )
        
       aice_obs = c0
       aice_obs_err = p1
       aice_inc = c0
       vice_inc = c0
       vsno_inc = c0
       ! we only assimilate at the ice edge?
       where( tmask .and. (TLAT >= 80. .or. TLAT <= -70 ))
          aice_obs = 0.99*c1
          aice_obs_err=p01*p01
       endwhere
       allocate(da_fresh(nx_block,ny_block,max_blocks), &
            da_fsalt(nx_block,ny_block,max_blocks),     &
            da_fheat(nx_block,ny_block,max_blocks) )
       da_fresh = c0
       da_fsalt = c0
       da_fheat = c0
    endif

    if (da_sit) then
       allocate (vice_obs(nx_block,ny_block,max_blocks), &
            vice_obs_err(nx_block,ny_block,max_blocks), &
            thice_obs(nx_block,ny_block,max_blocks), &
            thice_obs_err(nx_block,ny_block,max_blocks), &
            thice_inc(nx_block,ny_block,max_blocks) )
       thice_obs = c0
       thice_obs_err = p1
       thice_inc = c0
    end if
    if (da_sno) then
       allocate (vsno_obs(nx_block,ny_block,max_blocks), &
            vsno_obs_err(nx_block,ny_block,max_blocks))
       vsno_obs = c0
       vsno_obs_err = p1
       thsnow_obs = c0
       thsnow_obs_err = p1
       thsnow_inc = c0
    endif

    ! The rest of initialization is done in ice_da_init_streams, which is called
    ! from ice_comp_mct, after the mct initialization.
    ! call ice_da_init_stream

    call ice_timer_stop(timer_da)

           

  end subroutine ice_da_init

!=======================================================================

!  This initialize the pio data streams and finishes the initialization
      
  subroutine ice_da_init_streams(compid, gsmap, dom)
    use shr_pio_mod,     only : shr_pio_getiotype, shr_pio_getiosys, shr_pio_getioformat
    use ice_constants,   only: c0,c1,c2
    use ice_calendar,    only : calendar_type
    use ice_domain_size, only : nx_global,ny_global
    use ice_grid,        only:  TLAT
    use ice_itd,         only : hin_max

! !DESCRIPTION:
!    Prescribed ice initialization - needed to
!    work with new shr_strdata module derived type
!
! !Based on ice_prescribed_init, revised by 2009-Oct-12 - M. Vertenstein
!
! !Frist version:
!    2020-Apr-07 - J. B. Debernard   
!
! !INPUT/OUTPUT PARAMETERS:
!
   implicit none
   include 'mpif.h'
   integer(kind=int_kind), intent(in) :: compid
   type(mct_gsMap) :: gsmap
   type(mct_gGrid) :: dom

   !EOP

   !----- Local ------
   integer(kind=int_kind) :: nml_error ! namelist i/o error flag
   integer(kind=int_kind) :: n, nFile, ierr
   character(len=8)       :: fillalgo
   character(*),parameter :: subName = "('ice_da_init_streams')"


   namelist /da_nml/  &
        da_ice, da_sic, da_sit, da_sno, &
        da_timescale, da_method, da_data_dir,&
        da_update_ocn_heat,    &
        da_update_ocn_mass,    &
        da_stream_year_first  ,&
        da_stream_year_last   ,&
        da_model_year_align   ,&
        da_stream_fldVarName  ,&
        da_stream_fldFileName  ,&
        da_stream_domTvarName  ,&
        da_stream_domXvarName  ,&
        da_stream_domYvarName  ,&
        da_stream_domAreaName  ,&
        da_stream_domMaskName  ,&
        da_stream_domFileName  ,&
        da_stream_mapread  ,&
        da_ice_fill
        
! Default values 

   da_update_ocn_heat = .false.
   da_update_ocn_mass = .true.
   
   da_ice=.false.
   da_sic      = .true.
   da_sit      = .false.
   da_sno      = .false.
   da_timescale=5.   ! Timescale in days
   da_method   = 'pamip_long'
   da_data_dir = '.'

   da_stream_year_first      = 1                ! first year in  da_ice stream to use
   da_stream_year_last       = 1                ! last  year in  da_ice stream to use
   da_model_year_align       = 1                ! align stream_year_first with this model year
   da_stream_fldVarName      = 'ice_cov'
   da_stream_fldFileName(:)  = ' '
   da_stream_domTvarName     = 'time'
   da_stream_domXvarName     = 'lon'
   da_stream_domYvarName     = 'lat'
   da_stream_domAreaName     = 'area'
   da_stream_domMaskName     = 'mask'
   da_stream_domFileName     = ' '
   da_stream_mapread         = 'NOT_SET'
   da_ice_fill    = .false.   ! true if undefined assimilation data should be filled

   
   ! read from input file
   call get_fileunit(nu_nml)
   if (my_task == master_task) then
      open (nu_nml, file=nml_filename, status='old',iostat=nml_error)
      call shr_nl_find_group_name(nu_nml, 'da_nml', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, da_nml, iostat=nml_error)
         if (nml_error > 0) then
            call shr_sys_abort( 'problem on read of da_ice namelist in ice_da_init_streams' )
         endif
      endif
   end if
   call release_fileunit(nu_nml)
   call broadcast_scalar(da_ice,             master_task)

   if (my_task == master_task) then
      write(nu_diag,*) 'ice_da :: init, da_ice ::  ', da_ice
   endif
   ! *** If not da_ice then return ***
   if (.not. da_ice) RETURN
   call broadcast_scalar(da_sic,             master_task)
   call broadcast_scalar(da_sit,             master_task)
   call broadcast_scalar(da_sno,             master_task)
   call broadcast_scalar(da_timescale,       master_task)
   call broadcast_scalar(da_method,          master_task)
   call broadcast_scalar(da_data_dir,        master_task)
   
   call broadcast_scalar(da_model_year_align,master_task)
   call broadcast_scalar(da_stream_year_first,master_task)
   call broadcast_scalar(da_stream_year_last,master_task)
   call broadcast_scalar(da_stream_fldVarName,master_task)
   call broadcast_scalar(da_stream_domTvarName,master_task)
   call broadcast_scalar(da_stream_domXvarName,master_task)
   call broadcast_scalar(da_stream_domYvarName,master_task)
   call broadcast_scalar(da_stream_domAreaName,master_task)
   call broadcast_scalar(da_stream_domMaskName,master_task)
   call broadcast_scalar(da_stream_domFileName,master_task)
   call broadcast_scalar(da_stream_mapread,master_task)
   call broadcast_scalar(da_ice_fill,master_task)
   call mpi_bcast(da_stream_fldFileName, len(da_stream_fldFileName(1))*NFilesMaximum, &
        MPI_CHARACTER, 0, MPI_COMM_ICE, ierr)

   nFile = 0
   do n=1,nFilesMaximum
      if (da_stream_fldFileName(n) /= ' ') nFile = nFile + 1
   end do

   ! Read shr_strdata_nml namelist
   if (da_ice_fill) then
      fillalgo='nn'
   else
      fillalgo='none'
   endif

   if (my_task == master_task) then
      write(nu_diag,*) ' '
      write(nu_diag,*) '  This is the ice_da ice coverage option.'
      write(nu_diag,*) '  da_stream_year_first  = ',da_stream_year_first
      write(nu_diag,*) '  da_stream_year_last   = ',da_stream_year_last
      write(nu_diag,*) '  da_model_year_align   = ',da_model_year_align
      write(nu_diag,*) '  da_stream_fldVarName  = ',trim(da_stream_fldVarName)
      do n = 1,nFile
         write(nu_diag,*) '  da_stream_fldFileName = ',trim(da_stream_fldFileName(n)),n
      end do
      write(nu_diag,*) '  da_stream_domTvarName = ',trim(da_stream_domTvarName)
      write(nu_diag,*) '  da_stream_domXvarName = ',trim(da_stream_domXvarName)
      write(nu_diag,*) '  da_stream_domYvarName = ',trim(da_stream_domYvarName)
      write(nu_diag,*) '  da_stream_domFileName = ',trim(da_stream_domFileName)
      write(nu_diag,*) '  da_stream_mapread     = ',trim(da_stream_mapread)
      write(nu_diag,*) '  da_stream_fillalgo    = ',trim(fillalgo)
      write(nu_diag,*) ' '
   endif

   call shr_strdata_create(sdat,name="cice_da", &
        mpicom=MPI_COMM_ICE, compid=compid, &
        gsmap=gsmap, ggrid=dom,          &
        nxg=nx_global,nyg=ny_global,     &
        yearFirst=da_stream_year_first,     &
        yearLast=da_stream_year_last,       &
        yearAlign=da_model_year_align,      &
        offset=0,                        &
        domFilePath='',                  &
        domFileName=trim(da_stream_domFileName), &
        domTvarName=da_stream_domTvarName,  &
        domXvarName=da_stream_domXvarName,  &
        domYvarName=da_stream_domYvarName,  &
        domAreaName=da_stream_domAreaName,  &
        domMaskName=da_stream_domMaskName,  &
        filePath='',                     &
        filename=da_stream_fldFileName(1:nFile), &
        fldListFile=da_stream_fldVarName,   &
        fldListModel=da_stream_fldVarName,  &
        fillalgo=trim(fillalgo),       &
        calendar=trim(calendar_type),  &
        mapalgo='none',                &
        mapread=trim(da_stream_mapread))

   if (my_task == master_task) then
      call shr_strdata_print(sdat,'cice_da data')
   endif

   !-----------------------------------------------------------------
   ! For one ice category, set hin_max(1) to something big
   !-----------------------------------------------------------------
   if (ncat == 1) then
      hin_max(1) = 999._dbl_kind
   end if


   call ice_da_init
   
   if (trim(da_method)=='pamip_short') then
      if (my_task == master_task) then
         write(nu_diag,*) &
              subName, ':: DA PAMIP SHORT :: Ice thickness NH : 2.0 m , SH, 1.0 m '
!jd         write(nu_diag,*) &
!jd              'da_ice_init_streams :: DA PAMIP SHORT :: Sets aice = 0.99 north of 80N and South of 70 S'
      end if
      where ((TLAT >= c0))
         thice_obs = c2
         thice_obs_err = p01
      elsewhere
         thice_obs = c2
         thice_obs_err = p01
      endwhere

      vice_obs = thice_obs*aice_obs
      vice_obs_err = p1
      vice_inc = c0
   end if
   
    
 end subroutine ice_da_init_streams


      
!=======================================================================

!  This subroutine read observations, and call assimilation subroutines,

 subroutine ice_da_update

   use ice_blocks,   only: block, get_block, nblocks_x, nblocks_y
   use ice_domain,   only: ew_boundary_type, ns_boundary_type, &
        nblocks, blocks_ice
   use ice_state,    only: aicen, vicen, vsnon, trcrn, ntrcr, bound_state, &
        aice_init, aice0, aice, vice, vsno, trcr, trcr_depend
   use ice_itd,      only: aggregate
   use ice_grid,     only:  TLAT,tmask
   use ice_flux, only: Tf, sss, salinz, fresh, fsalt, fhocn

   
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
        i,j,iblk,nt,n,      & ! dummy loop indices
        ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
        iglob(nx_block),    & ! global indices
        jglob(ny_block),    & ! global indices
        iblock, jblock,     & ! block indices
        ibc,                & ! ghost cell column or row
        npad,               & ! padding column/row counter
        fid                   ! file id for netCDF routines

   character (char_len_long) :: &
        data_file             ! data file for observations

   character (char_len) :: &
        da_date,            & ! date for data assimilation
        fieldname             ! field name in netcdf file

   type (block) :: &
        this_block            ! block info for current block

   call ice_timer_start(timer_da)

!-----------------------------------------------------------------------
!  Read observations
!-----------------------------------------------------------------------

   if (trim(da_method)=='coin' .and. ( (istep == 1) .or. new_day ) ) then

      ! sea ice concentration & uncertainties
      write(da_date,'(i4)') idate/10000
      data_file = trim(da_data_dir)//'osisaf_'//trim(da_date)//'.nc'
      if (my_task == master_task) write(nu_diag,*) 'DA data file = ', data_file

      call ice_open_nc(data_file,fid)

      fieldname = 'obsAice'
      call ice_read_nc (fid, int(yday), fieldname, aice_obs, dbug, &
           field_loc_center, field_type_scalar)
        
      fieldname = 'obsAerr'
      call ice_read_nc (fid, int(yday), fieldname, aice_obs_err, dbug, &
           field_loc_center, field_type_scalar)
      
      where ((aice_obs >= c0) .and. (aice_obs) <= 100)
         aice_obs = aice_obs * p01
         aice_obs_err = aice_obs_err * p01
      elsewhere
         aice_obs = c0
         aice_obs_err = c0
      endwhere

      call ice_close_nc(fid)

   endif
   
   if (trim(da_method) == 'coin') then

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     iglob,jglob,iblock,jblock)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iglob = this_block%i_glob
         jglob = this_block%j_glob
         iblock = this_block%iblock
         jblock = this_block%jblock

         call da_coin  (nx_block,            ny_block,      &
                        ilo, ihi,            jlo, jhi,      &
                        iglob,               jglob,         &
                        iblock,              jblock,        &
                        tmask(:,:,    iblk),                &
                        aice(:,:,     iblk),                &
                        aice_obs(:,:, iblk),                &
                        aice_obs_err(:,:, iblk),            &
                        aicen(:,:,  :,iblk),                &
                        vicen(:,:,  :,iblk),                &
                        vsnon(:,:,  :,iblk),                &
                        trcrn(:,:,:,:,iblk), ntrcr)         
      enddo ! iblk
      !$OMP END PARALLEL DO
   else if (trim(da_method) == 'pamip_short' .or. trim(da_method) == 'pamip_long') then

      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     iglob,jglob,iblock,jblock)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call da_pamip_update (nx_block, ny_block, &
              ilo, ihi, jlo, jhi,                  &
              TLAT(:,:,  iblk),      &
              tmask(:,:,  iblk),     &
              sss(:,:,iblk),         &
              Tf(:,:,iblk),          &
              salinz(:,:,:,iblk),    &
              aice(:,:,  iblk),      &
              aice_inc(:,:,iblk),    &
              thice_inc(:,:,iblk),   &
              vice_inc(:,:,iblk),    &
              vsno_inc(:,:,iblk),           &
              da_fresh(:,:,iblk),           &
              da_fsalt(:,:,iblk),           &
              da_fheat(:,:,iblk),           &
              aicen(:,:,  :,iblk),          &
              vicen(:,:,  :,iblk),          &
              vsnon(:,:,  :,iblk),          &
              trcrn(:,:,:,:,iblk), ntrcr)         

         if (da_update_ocn_heat) then
            !jd not recommended in the present implementation, increasing ice
            ! would imply warming of the ocean. Could implement a redistribution
            ! before transfer to ocean?
            fhocn(:,:,iblk) = fhocn(:,:,iblk) &
                 + da_fheat(:,:,iblk)
         end if
         if (da_update_ocn_mass) then
            fresh(:,:,iblk) = fresh(:,:,iblk) &
                 + da_fresh(:,:,iblk)
            fsalt(:,:,iblk) = fsalt(:,:,iblk) &
                 + da_fsalt(:,:,iblk)
         end if
         

      enddo ! iblk
      !$OMP END PARALLEL DO



      
   endif

   
   !-----------------------------------------------------------------
   ! aggregate tracers
   !-----------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblk)
   do iblk = 1, nblocks

      call aggregate (nx_block, ny_block, &
                      aicen(:,:,:,iblk),  &
                      trcrn(:,:,:,:,iblk),&
                      vicen(:,:,:,iblk),  &
                      vsnon(:,:,:,iblk),  &
                      aice (:,:,  iblk),  &
                      trcr (:,:,:,iblk),  &
                      vice (:,:,  iblk),  &
                      vsno (:,:,  iblk),  &
                      aice0(:,:,  iblk),  &
                      tmask(:,:,  iblk),  &
                      max_ntrcr,          &
                      trcr_depend)

   enddo
   !$OMP END PARALLEL DO

   call ice_timer_stop(timer_da)

 end subroutine ice_da_update

 subroutine ice_da_prep(mDateIn,secIn)
   use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
   use ice_domain, only: ew_boundary_type, ns_boundary_type, &
        nblocks, blocks_ice
   use ice_state, only: aicen, vicen, vsnon, trcrn, ntrcr, bound_state, &
        aice_init, aice0, aice, vice, vsno, trcr, trcr_depend
   use ice_itd, only: aggregate
   use ice_grid, only:  tmask

   implicit none


! !INPUT/OUTPUT PARAMETERS:

   integer(kind=int_kind), intent(in) :: mDateIn  ! Current model date (yyyymmdd)
   integer(kind=int_kind), intent(in) :: secIn    ! Elapsed seconds on model date

!EOP
   
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,iblk,nt,n,      & ! dummy loop indices
     ilo,ihi,jlo,jhi,    & ! beginning and end of physical domain
     iglob(nx_block),    & ! global indices
     jglob(ny_block),    & ! global indices
     iblock, jblock,     & ! block indices
     ibc,                & ! ghost cell column or row
     npad,               & ! padding column/row counter
     fid                   ! file id for netCDF routines

   character (char_len_long) :: &
     data_file             ! data file for observations

   character (char_len) :: &
     da_date,            & ! date for data assimilation
     fieldname             ! field name in netcdf file

   type (block) :: &
     this_block            ! block info for current block

   call ice_timer_start(timer_da)

   if (trim(da_method) == 'pamip_short' .or. trim(da_method) == 'pamip_long') then
   call shr_strdata_advance(sdat,mDateIn,SecIn,MPI_COMM_ICE,'cice_da')
      
   n=0
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      do j = jlo, jhi
         do i = ilo, ihi
            n = n+1
            aice_obs(i,j,iblk) = sdat%avs(1)%rAttr(1,n)
         end do
      end do
   end do
   
!$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
!$OMP                     iglob,jglob,iblock,jblock)
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      call da_pamip_prep (nx_block,            ny_block,      &
                          ilo, ihi,            jlo, jhi,      &
                          tmask(:,:,    iblk),                &
                          aice(:,:,     iblk),                &
                          aice_obs(:,:, iblk),                &
                          thice_obs(:,:, iblk),               &
                          aice_inc(:,:, iblk),                &
                          thice_inc(:,:,iblk) )
   enddo ! iblk
   !$OMP END PARALLEL DO

   end if


   call ice_timer_stop(timer_da)
 end subroutine ice_da_prep



!=======================================================================

 subroutine da_pamip_prep (nx_block,            ny_block,      &
                          ilo, ihi,            jlo, jhi,      &
                          tmask,               aice,          &
                          o_aice,              o_thice,  &
                          inc_aice,            inc_thice )

!del   use ice_blocks, only: nblocks_x, nblocks_y

  integer (kind=int_kind), intent(in) :: &
       nx_block, ny_block, & ! block dimensions
       ilo, ihi          , & ! physical domain indices
       jlo, jhi              !

  logical (kind=log_kind), dimension (nx_block,ny_block), &
       intent(in) :: &
       tmask                 ! true for ice/ocean cells

  real (kind=dbl_kind), dimension (nx_block,ny_block), &
       intent(in) :: &
       aice,               & ! model aggregate sic
       o_aice,             & ! observed aggregate sic
       o_thice               ! observed ice thickness
         
  real (kind=dbl_kind), dimension (nx_block,ny_block), &
       intent(inout) :: &
       inc_aice , & ! increment in ice concentration
       inc_thice    ! increment in ice thickness

 
  ! local variables

  integer (kind=int_kind) :: &
       i, j        !, & ! horizontal indices
!del       ij          , & ! horizontal index, combines i and j loops
!del       icells          ! number of cells initialized with ice

!del  integer (kind=int_kind), dimension(nx_block*ny_block) :: &
!del       indxi, indxj    ! compressed indices for cells with restoring

  real (kind=dbl_kind) :: &
       obs_ai,       & ! local version of thickness used
       rda,          & ! dt/dT, where dT is observation time step
       radd            ! incremental ratio

  !-----------------------------------------------------------------
  ! Calculate increment in ice concentration and thickness
  !-----------------------------------------------------------------

  rda = dt / (da_timescale*real(secday,kind=dbl_kind))

  if (da_sic == .true.) then

     do j = jlo,jhi
        do i = ilo,ihi
           if (tmask(i,j) ) then
              obs_ai = min(o_aice(i,j), max_aice_obs)
              if (obs_ai < min_aice_obs) obs_ai = c0
              inc_aice(i,j)=rda*(obs_ai - aice(i,j))
           else
              inc_aice(i,j)=c0
              
           endif
           
        enddo
     enddo
  endif
  
end subroutine da_pamip_prep
!=======================================================================

subroutine da_pamip_update (nx_block,ny_block,   &
                          ilo, ihi,  jlo, jhi,   &
                          TLAT,                  &
                          tmask,     sss,        &
                          Tf,        salinz,     &
                          aice,                  &
                          inc_aice,  inc_thice,  &
                          inc_vice,  inc_vsno,   &
                          da_fresh,  da_fsalt,  da_fheat, &
                          aicen,     vicen,     vsnon,    &
                          trcrn,     ntrcr)         

  use ice_constants, only: rhoi, rhos, p001, ice_ref_salinity,&
                           c2,c5
  use ice_state, only: nt_qice, nt_qsno, nt_aero, tr_aero
  
  use ice_domain_size, only: n_aero
!  use ice_grid, only:  TLAT
    
  integer (kind=int_kind), intent(in) :: &
       nx_block, ny_block, & ! block dimensions
       ilo, ihi          , & ! physical domain indices
       jlo, jhi          , & !
       ntrcr                 ! number of tracers in use

  logical (kind=log_kind), dimension (nx_block,ny_block), &
       intent(in) :: &
       tmask                ! true for ice/ocean cells

  real (kind=dbl_kind), dimension (nx_block,ny_block), &
       intent(in) :: &
       Tf,           & ! freezing temperature (C)
       sss,          & ! sea surface salinity (ppt)
       aice,         & ! model aggregate sic
       inc_aice,     & ! increment in aggregated ice concentration
       inc_thice,    & ! incrmement in mean sea ice thickness
       TLAT

  real (kind=dbl_kind),dimension(nx_block,ny_block,nilyr+1),&
       intent(in) :: &
       salinz     ! initial salinity profile

  
  real (kind=dbl_kind), dimension (nx_block,ny_block), &
       intent(inout) :: &
       inc_vice,     & ! increment in ice volume
       inc_vsno,     & ! increment in snow volume
       da_fresh,     & ! freshwater change due to nudging
       da_fsalt,     & ! Salt flux due to nudgning
       da_fheat        ! heat change due to nudgning

  
  real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
       intent(inout) :: &
       aicen , & ! concentration of ice
       vicen , & ! volume per unit area of ice          (m)
       vsnon     ! volume per unit area of snow         (m)
  
  real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
       intent(inout) :: &
       trcrn     ! ice tracers
       ! 1: surface temperature of ice/snow (C)
  
  ! local variables

  integer (kind=int_kind) :: &
       i, j        , & ! horizontal indices
       ij          , & ! horizontal index, combines i and j loops
       ibc         , & ! ghost cell column or row
       npad        , & ! padding column/row counter
       k           , & ! ice layer index
       n           , & ! thickness category index
       it              ! tracer index

  real (kind=dbl_kind) :: &
       dvice,    &
       dvsno,    &
       dfresh,   &
       dfsalt,    &
       dfheat,    &
       vice0,    &
       vsno0,    &
       radd            ! incremental ratio

real (kind=dbl_kind), dimension(ncat) :: &
       radn_aero, &      ! Increment used for aerosols in each cathegory
       radn              ! Increment used in each cathegory

  
  !-----------------------------------------------------------------
  ! assimilate sic on grid
  !-----------------------------------------------------------------
  
  if (da_sic == .true.) then
     do j = jlo,jhi
        do i = ilo,ihi
!jd           if ((tmask(i,j)).and.(TLAT(i,j)>c0)) then
           if (tmask(i,j)) then
              da_fresh(i,j) = c0
              da_fsalt(i,j) = c0
              da_fheat(i,j) = c0
              
              !jd              if (aice(i,j) > puny) then
              if (aice(i,j) > c0) then

! Calculate increments for each category
! Limit ice growth factor to maximum 10.m.  ( 0 <= radd <= 10 )
                 radd = min(max(c0,c1 + inc_aice(i,j)/max(aice(i,j),puny)),c10)
                 radn(:)=radd
                 radn_aero(:)=c1
                 !remove ice in all categories
                 !but only increase it in relatively
                 !thin categories (below 10 m thickness)
                 if (radd>c1) then
                    do n=1, ncat
                       if ((vicen(i,j,n)/aicen(i,j,n))>c10) radn(n)=c1
                    enddo ! ncat
!jd !Reduce aerosol consentration when only small amount of ice exists before artificially increasing the ice concentration. 
                    if (aice(i,j) < p2) radn_aero(:)=c1/radn(:)
                 endif

! Update state
                 do n=1, ncat
                    vice0 = vicen(i,j,n)
                    vsno0 = vsnon(i,j,n)
                    aicen(i,j,n) = aicen(i,j,n) * radn(n)
                    vicen(i,j,n) = vicen(i,j,n) * radn(n)
                    vsnon(i,j,n) = vsnon(i,j,n) * radn(n)
                    if (tr_aero) then
                          trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,n)=&
                          trcrn(i,j,nt_aero:nt_aero+4*n_aero-1,n)*radn_aero(n)
                    end if
                    !
                    dvice = vicen(i,j,n) - vice0
                    dvsno = vsnon(i,j,n) - vsno0
                    ! From ice
                    dfresh = -dvice*rhoi/dt
                    dfsalt = ice_ref_salinity*p001*dfresh
                    ! From snow, no salt
                    dfresh = dfresh - dvsno*rhos/dt
                    !
                    da_fresh(i,j) = da_fresh(i,j) + dfresh
                    da_fsalt(i,j) = da_fsalt(i,j) + dfsalt
                    !
                    do k = 1, nilyr
                       ! enthalpy tracers do not change (e/v constant)
                       ! heat flux loss from ice increase
                       dfheat = trcrn(i,j,nt_qice+k-1,n)* dvice / (dt*real(nilyr,kind=dbl_kind))
                       da_fheat(i,j) = da_fheat(i,j) + dfheat
                    enddo   ! nilyr

                    do k = 1, nslyr
                       dfheat = trcrn(i,j,nt_qsno+k-1,n)*dvsno / (dt*real(nslyr,kind=dbl_kind))
                       da_fheat(i,j) = da_fheat(i,j) + dfheat
                    enddo                  ! nslyr

                 enddo  ! ncat
                 
              else if (inc_aice(i,j) > puny) then
                 !jd Add new ice based on similar approach  as in ice_therm_itd :: add_new_ice
                 call da_add_new_ice_pt(ntrcr, dt,   &
                           Tf        (i,j),          &
                           sss       (i,j),          &
                           salinz    (i,j,:),        &
                           inc_aice  (i,j),          &         
                           inc_vice  (i,j),          &
                           aicen     (i,j,:),        &
                           trcrn     (i,j,1:ntrcr,:),&
                           vicen     (i,j,:),        &
                           da_fresh  (i,j),          &
                           da_fsalt  (i,j),          &
                           da_fheat  (i,j) )
              endif
           endif  ! tmask
         enddo
      enddo

      
      endif  

    end subroutine da_pamip_update

!=======================================================================

    subroutine da_add_new_ice_pt(ntrcr,  dt,         &
                              Tf,        sss,        &
                              salinz,                &
                              inc_aice,  inc_vice,   &
                              aicen,     trcrn,      &
                              vicen,                 &
                              fresh,     fsalt,      &
                              fheat)
!jd                              fiso_ocn,              &
!jd                              HDO_ocn,             &
!jd                              H2_16O_ocn,          &
!jd                              H2_18O_ocn,          &
!jd                              nbtrcr,    flux_bio,   &
!jd                              ocean_bio, &


      use ice_constants, only: rhoi, Lfresh, ice_ref_salinity, p001,p5,c2,c4
      use ice_itd, only: hin_max, column_sum, &
                         column_conservation_check 
      use ice_state, only: nt_Tsfc, nt_iage, nt_FY, nt_alvl, nt_vlvl, nt_aero, &
                           nt_iso, nt_sice, nt_qice, &
                           nt_apnd, tr_pond_cesm, tr_pond_lvl, tr_pond_topo, &
                           tr_iage, tr_FY, tr_lvl, tr_aero, tr_iso, tr_brine
      use ice_therm_mushy, only: liquidus_temperature_mush, enthalpy_mush
      use ice_therm_shared, only: ktherm, hfrazilmin
      use ice_therm_vertical, only:  phi_init, dSin0_frazil
      
      integer (kind=int_kind), intent(in) :: &
         ntrcr              ! number of tracers in use

      real (kind=dbl_kind), intent(in) :: &
         dt        ! time step (s)

      real (kind=dbl_kind), intent(in) :: &
           inc_aice, & ! positive increment of sea ice consentration
           Tf    , & ! freezing temperature (C)
           sss       ! sea surface salinity (ppt)

      real (kind=dbl_kind), dimension(nilyr+1), intent(in) :: &
         salinz     ! initial salinity profile

      real (kind=dbl_kind), dimension (ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen     ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature

      real (kind=dbl_kind),  intent(inout) :: &
         inc_vice, & ! positive increment of sea ice consentration           
         fresh     , & ! fresh water flux to ocean (kg/m^2/s)
         fsalt     , & ! salt flux to ocean (kg/m^2/s)
         fheat         ! Heat flux change due to nudging


      ! local variables

      integer (kind=int_kind) :: &
         n            , & ! ice category index
         k            , & ! ice layer index
         it               ! aerosol tracer index

      real (kind=dbl_kind)  :: &
         thice_da,    &   ! thickness of ice added by nudgning in ice free grid points
         ai0new       , & ! area of new ice added to cat 1
         vi0new           ! volume of new ice added to cat 1

      real (kind=dbl_kind) :: &
         vice1        , & ! starting volume of existing ice
         vice_init, vice_final, & ! ice volume summed over categories
         eice_init, eice_final    ! ice energy summed over categories

      real (kind=dbl_kind) :: &
!jd         area1        , & ! starting fractional area of existing ice
!jd         alvl         , & ! starting level ice area
         dfresh       , & ! change in fresh
         dfsalt       , & ! change in fsalt
         dfheat        , & ! change in heat
         Ti               ! frazil temperature
      

      real (kind=dbl_kind) :: &
         qi0new       , & ! frazil ice enthalpy
         Si0new           ! frazil ice bulk salinity

      real (kind=dbl_kind), dimension (nilyr) :: &
         Sprofile         ! salinity profile used for new ice additions




      if (ktherm == 2) then  ! mushy
         if (sss  > c2 * dSin0_frazil) then
            Si0new = sss - dSin0_frazil
         else
            Si0new = sss**2 / (c4*dSin0_frazil)
         endif
         do k = 1, nilyr
            Sprofile(k) = Si0new
         enddo
         Ti = min(liquidus_temperature_mush(Si0new/phi_init), -p1)
         qi0new = enthalpy_mush(Ti, Si0new)
         
      else
         do k = 1, nilyr
            Sprofile(k) = salinz(k)
         enddo
         qi0new = -rhoi*Lfresh
      endif    ! ktherm


      !jd
      ! The main assumtion here is that inc_aice > 0,
      ! aice < puny, (open water)
      ! All ice goes to the thinnest cathegory
      ! We use the hfrazilmin value for the ice included.
      ! This is perhaps not optimal, and escpecially in
      ! with rapid melting, the new ice might melt before
      ! next time-step. Could increase this thickness.
      ! Also in the case with only ice cathegory (ncat=1),
      ! a larger value could be more approriate.

      thice_da = hfrazilmin
      if (ncat == 1) thice_da = p5
      
      inc_vice = thice_da*inc_aice

      dfresh = -rhoi*inc_vice / dt
      dfsalt = ice_ref_salinity*p001*dfresh
      
      fresh = fresh + dfresh
      fsalt = fsalt + dfsalt

      dfheat= qi0new* inc_vice / dt
      fheat = fheat + dfheat

      !jd Assume this are zero due to the nudge update logic
      !jd      vice1 = vicen(1)  ! Store previous value
      !jd      area1 = aicen(1)  !
      

      !jd      aicen(1) = aicen(1) + inc_aice
      !jd      vicen(1) = vicen(1) + inc_vice
      aicen(1) = inc_aice
      vicen(1) = inc_vice
      

      trcrn(nt_Tsfc,1) = Tf
      trcrn(nt_Tsfc,1) = min (trcrn(nt_Tsfc,1), c0)

      if (tr_FY) then
         trcrn(nt_FY,1) = c1
      endif

      if (tr_iage) &
           trcrn(nt_iage,1) = dt

      if (tr_aero) then
         do it = 1, nt_aero
            trcrn(nt_aero+2+4*(it-1),1) = c0
            trcrn(nt_aero+3+4*(it-1),1) = c0
         enddo
      endif

!jd      if (tr_iso) then  ! Not implemented
!jd      end if

      if (tr_lvl) then
         trcrn(nt_alvl,1) = c1
         trcrn(nt_vlvl,1) = c1
      end if

      !jd Ponds should be zero allready

      ! Update vertical profile
      do k = 1, nilyr
         if (vicen(1) > c0) then
            ! factor of nilyr cancels out
            ! enthalpy
            trcrn(nt_qice+k-1,1) = qi0new
            ! salinity
            trcrn(nt_sice+k-1,1) = Sprofile(k)
         endif
      enddo


      
    end subroutine da_add_new_ice_pt
      

!=======================================================================


subroutine da_coin    (nx_block,            ny_block,      &
                       ilo, ihi,            jlo, jhi,      &
                       iglob,               jglob,         &
                       iblock,              jblock,        &
                       tmask,               aice,          &
                       aice_obs,            aice_obs_err,  &
                       aicen,     vicen,    vsnon,         &
                       trcrn,     ntrcr)         

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)   , & !
         iblock            , & ! block indices
         jblock            , & !
         ntrcr                 ! number of tracers in use

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(inout) :: &
         tmask                 ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         aice,               & ! model aggregate sic
         aice_obs,           & ! observed aggregate sic
         aice_obs_err          ! observed aggregate sic_err
         
      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), &
         intent(inout) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), &
         intent(inout) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         ibc         , & ! ghost cell column or row
         npad        , & ! padding column/row counter
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with restoring

      real (kind=dbl_kind) :: &
         mod_err,      & ! model error
         mod_err2,     & ! model error square
         gain,         & ! Kalman gain, optimal estimated
         weight,       & ! nudging weight, increamental Kalman gain
         weightn,      & ! nudging weigth distributed among categories
         rda,          & ! dt/dT, where dT is observation time step
         radd            ! incremental ratio

      !-----------------------------------------------------------------
      ! assimilate sic on grid
      !-----------------------------------------------------------------

      rda = dt / real(secday,kind=dbl_kind)
!jd      if (my_task == master_task) write(nu_diag,*) 'rda = ', rda

      if (da_sic == .true.) then

         do j = 1, ny_block
         do i = 1, nx_block
            if (tmask(i,j)) then
               mod_err  = aice(i,j) - aice_obs(i,j)
               mod_err2 = mod_err**2 + aice_obs_err(i,j)**2
               gain = mod_err2/(mod_err2+puny+aice_obs_err(i,j)**2)
               weight = c1 - (c1 - gain)**rda
            else
               weight = c0
            endif

            if ((aice(i,j) > puny) .or. (aice_obs(i,j) > p1)) then
               radd = c1 + weight * (aice_obs(i,j)/max(aice(i,j),p1) - c1)
               do n=1, ncat
                  aicen(i,j,n) = aicen(i,j,n) * radd
                  vicen(i,j,n) = vicen(i,j,n) * radd
  
                  do it = 1, ntrcr
                     trcrn(i,j,it,n) = trcrn(i,j,it,n) * radd
                  enddo
               enddo
            endif
         enddo
         enddo
      endif

end subroutine da_coin

!=======================================================================

      end module ice_da

!=======================================================================
