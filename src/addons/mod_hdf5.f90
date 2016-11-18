module mod_hdf5

interface hdf5_write
  module procedure hdf5_write_i4,hdf5_write_d,hdf5_write_z
end interface

interface hdf5_read
  module procedure hdf5_read_i4,hdf5_read_d,hdf5_read_z
end interface

public hdf5_initialize
public hdf5_finalize
public hdf5_create_file
public hdf5_create_group

contains

subroutine hdf5_initialize
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5open_f(ierr)
#endif
end subroutine

subroutine hdf5_finalize
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5close_f(ierr)
#endif
end subroutine

subroutine hdf5_create_file(fname)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
#ifdef _HDF5_
integer ierr
integer(hid_t) h5_root_id
call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_file): h5fcreate_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_file): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'("  fname: ",A)')trim(fname)
stop
#endif
end subroutine

subroutine hdf5_create_group(fname,path,gname)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: gname
#ifdef _HDF5_
integer(hid_t) h5_root_id,h5_group_id,h5_new_group_id
integer ierr

call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5gcreate_f(h5_group_id,trim(gname),h5_new_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gcreate_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(h5_new_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gclose_f for the new group returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(h5_group_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5gclose_f for the existing path returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,'("Error(hdf5_create_group): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path  : ",A)')trim(path)
write(*,'("  gname : ",A)')trim(gname)  
stop
#endif
end subroutine

!----------------------!
!     hdf5_write_i4    !
!----------------------!
subroutine hdf5_write_i4(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
integer, intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_write_array_i4(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_write_d    !
!---------------------!
subroutine hdf5_write_d(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
real(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_write_z    !
!---------------------!
subroutine hdf5_write_z(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)+1
  allocate(dims_(ndims))
  dims_(1)=2
  dims_(2:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=2
endif
call hdf5_write_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_read_i4    !
!---------------------!
subroutine hdf5_read_i4(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
integer, intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_read_array_i4(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_d    !
!--------------------!
subroutine hdf5_read_d(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
real(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_z    !
!--------------------!
subroutine hdf5_read_z(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)+1
  allocate(dims_(ndims))
  dims_(1)=2
  dims_(2:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=2
endif
call hdf5_read_array_d(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

end module

#ifdef _HDF5_


subroutine hdf5_write_array_i4(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer(4), intent(in) :: a(*)
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(hsize_t), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5screate_simple_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_INTEGER,dataspace_id,dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5dcreate_f returned ",I6)')ierr
  goto 10
endif 
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5dwrite_f returned ",I6)')ierr
  goto 10
endif 
call h5dclose_f(dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5dclose_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5gclose_f returned ",I6)')ierr
  goto 10
endif
call h5sclose_f(dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5sclose_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_i4): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims  : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path  : ",A)')trim(path)
write(*,'("  nm    : ",A)')trim(nm)
stop
end subroutine

subroutine hdf5_write_array_d(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
real(8), intent(in) :: a(*)
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(hsize_t), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5screate_simple_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_DOUBLE,dataspace_id,dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5dcreate_f returned ",I6)')ierr
  goto 10
endif 
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5dwrite_f returned ",I6)')ierr
  goto 10
endif 
call h5dclose_f(dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5dclose_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5gclose_f returned ",I6)')ierr
  goto 10
endif
call h5sclose_f(dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5sclose_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_write_array_d): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims  : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path  : ",A)')trim(path)
write(*,'("  nm    : ",A)')trim(nm)
stop
end subroutine


subroutine hdf5_read_array_i4(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer(4), intent(out) :: a(*)
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5dopen_f returned ",I6)')ierr
  goto 10
endif
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5dread_f returned ",I6)')ierr
  goto 10
endif
call h5dclose_f(dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5dclose_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5gclose_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_i4): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,*)
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path : ",A)')trim(path)
write(*,'("  nm : ",A)')trim(nm)
stop
end subroutine

subroutine hdf5_read_array_d(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
real(8), intent(out) :: a(*)
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5fopen_f returned ",I6)')ierr
  goto 10
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5gopen_f returned ",I6)')ierr
  goto 10
endif
call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5dopen_f returned ",I6)')ierr
  goto 10
endif
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5dread_f returned ",I6)')ierr
  goto 10
endif
call h5dclose_f(dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5dclose_f returned ",I6)')ierr
  goto 10
endif
call h5gclose_f(group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5gclose_f returned ",I6)')ierr
  goto 10
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(hdf5_read_array_d): h5fclose_f returned ",I6)')ierr
  goto 10
endif
return
10 continue
write(*,*)
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path : ",A)')trim(path)
write(*,'("  nm : ",A)')trim(nm)
stop
end subroutine


#endif

