module module_pet
 use netcdf
 use omp_lib
 implicit none
! integer intent(in):: ncid
contains
!----------------------------------------------------------------------------------
  subroutine rd_2d_var(vname, ncid, d_array, st_lat, en_lat, st_lon, en_lon, iserr)
    use netcdf    
    character(len=*), intent(in) :: vname
    integer, intent(in)  :: ncid
    integer, intent(in)  :: st_lat, en_lat, st_lon, en_lon
    real, dimension(st_lon:en_lon,st_lat:en_lat), intent(out) :: d_array
    integer              :: ierr, varid
    integer, intent(out) :: iserr

    ierr = nf90_inq_varid(ncid,  vname,  varid)
    !print*,varid
    if (ierr /= 0) then
          print*, 'ncid = ', ncid, "Error Reading Variable "//trim(vname)//" in input NetCDF file."
    else
          iserr = ierr
    endif

    ierr = nf90_get_var(ncid, varid, values=d_array, start=(/st_lon,st_lat,1/), count=(/en_lon-st_lon+1,en_lat-st_lat+1,1/))
    if (ierr /= 0) then
          print*, 'ncid = ', ncid, "Error Reading Variable "//trim(vname)//" in input NetCDF file."
    else
          iserr = ierr
          return
    endif

  end subroutine rd_2d_var
!----------------------------------------------------------------------------------
subroutine list_files(path,wrfFileNames,NFiles)
   character(len=100), dimension(:),allocatable, intent(inout) :: wrfFileNames
   character(*), intent(in) :: path
   integer :: cnt,iserr,cntwrf
   integer,intent(out) :: NFiles
   real    :: irer
   ! get the files
   call system("ls "//trim(path)//"OBJ1_wrfout_d01* >filelist.txt")
   open(41,FILE='filelist.txt',action="read")
   !counts of file
   cnt = 0
   do
     read(41,FMT='(a)',iostat=iserr) irer
     if (iserr/=0) EXIT
     cnt = cnt+1
   end do
   NFiles = cnt
   print*, "Number of NetCDF files: " , NFiles
   allocate(wrfFileNames(NFiles))
   rewind(41)
   do cntwrf = 1,NFiles
    read(41,'(a)') wrfFileNames(cntwrf)
   !print*,wrfFileNames(cntwrf)
   end do
   close(41)
   return
end subroutine list_files

!---------------------------------------------------------------------------------------

 subroutine wrt_pet_var(pet_out,outfname,nlons,nlats)
  use netcdf 
  character (len = *), intent(in) ::outfname   
  character (len = *), parameter :: FILE_NAME = "PET_"
  integer, intent(in):: nlons,nlats
  integer, parameter :: nrecs=1,ndims=3

  integer    :: lats(nlats), lons(nlons)
  integer :: lon_varid  , lat_varid

  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  character (len = *), parameter :: REC_NAME = "time"
  integer                        :: lon_dimid, lat_dimid, rec_dimid
  integer, dimension(3)          :: start, counts
  character (len = *), parameter :: PET_NAME="Penman_PET"
  integer :: pet_varid
  integer :: dimids(2)!dimids(3)

  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PET_UNITS = "mmh-1"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"
  real,dimension(1:nlons, 1:nlats), intent(in):: pet_out
  integer, parameter :: st_lat = 1, st_lon = 1
  integer :: lat, lon, rec, i
  integer :: ncid

  do lat = 1, nlats
     lats(lat) = lat
  end do
  do lon = 1, nlons
     lons(lon) = lon
  end do
  !print*,FILE_NAME//outfname(16:)

  call iserror( nf90_create(FILE_NAME//outfname(16:), nf90_clobber, ncid) )

  call iserror( nf90_def_dim(ncid, LAT_NAME, nlats, lat_dimid) )
  call iserror( nf90_def_dim(ncid, LON_NAME, nlons, lon_dimid) )
  !call iserror( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

  call iserror( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
  call iserror( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )

  call iserror( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
  call iserror( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

  dimids = (/ lon_dimid, lat_dimid /)

  call iserror( nf90_def_var(ncid, PET_NAME, NF90_REAL, dimids, pet_varid) )

  call iserror( nf90_put_att(ncid, pet_varid, UNITS, PET_UNITS) )

  call iserror( nf90_enddef(ncid) )

  call iserror( nf90_put_var(ncid, lat_varid, lats) )
  call iserror( nf90_put_var(ncid, lon_varid, lons) )
  counts = (/ nlons, nlats, 1 /)
  start = (/ 1, 1, 1 /)

  !do rec = 1, nrecs
     !start(3) = rec
     call iserror( nf90_put_var(ncid, pet_varid, pet_out, start = (/1,1/)) )
  !end do


  contains
    subroutine iserror(status)
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
    end subroutine iserror  
end subroutine wrt_pet_var

end module module_pet
