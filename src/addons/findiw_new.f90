subroutine findiw_new(nw0,iwmap,dw)
use modmain
use mod_linresp
!
implicit none
integer,intent(in) :: nw0
integer,intent(out) :: iwmap(nw0,nw0)
real(8),intent(out) :: dw(nw0,nw0)
!
integer iw,iw1,iw2
real(8) :: dw1
!
! find the index t0 for t1-t2, i.e. A(t0)=A(t1-t2)
iwmap=-100
dw=0.d0

do iw2=1,nw0
 do iw1=1,nw0
  if (iw1.eq.iw2) then
   iwmap(iw1,iw2)=1
   dw(iw1,iw2)=0.d0
   cycle
  endif

  dw1=dreal(lr_w(iw1)-lr_w(iw2))
  if (iw1.lt.iw2) dw1=dw1+bhbar

  do iw=1,nw0-1
   if (dreal(lr_w(iw))-1d-8.le.dw1.and.dreal(lr_w(iw+1))+1d-8.ge.dw1) then
    iwmap(iw1,iw2)=iw
    dw(iw1,iw2)=(dw1-dreal(lr_w(iw)))/dreal(lr_w(iw+1)-lr_w(iw))
    exit
   endif
  enddo
 enddo !iw1
enddo !iw2

!if (mpi_grid_root()) then
! do iw2=1,nw0
!  do iw1=1,nw0
!   write(*,*) "iw1,iw2,iwmap:",iw1,iw2,iwmap(iw1,iw2)
!   write(*,*) "dw:",dw(iw1,iw2)
!  enddo
! enddo
!endif

return
end subroutine
