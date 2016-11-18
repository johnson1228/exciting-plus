integer function findiw(iw1,iw2)
use modmain
use mod_linresp
!
implicit none
integer,intent(in) :: iw1
integer,intent(in) :: iw2
!
integer iw,w0
real(8) dw
!
! find the index t0 for t1-t2, i.e. A(t0)=A(t1-t2)
w0=-100
dw=dreal(lr_w(iw1)-lr_w(iw2))
if (iw1.eq.iw2) then
 findiw=1
 return
endif
if (iw1.lt.iw2) dw=dw+bhbar

do iw=1,lr_nw-1
 if (dreal(lr_w(iw))-1d-8.le.dw.and.dreal(lr_w(iw+1))+1d-8.ge.dw) then
  if (dw.ge.dreal(lr_w(iw)+lr_w(iw+1))*0.5d0) then
   w0=iw+1
  else
   w0=iw
  endif
  exit
 endif
enddo

findiw=w0
return
end function findiw
