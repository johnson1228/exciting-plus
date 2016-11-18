subroutine cal_zaux(nw,gf0,sigc,z1,z2,zaux)
use modmain
use mod_linresp
!
implicit none
integer, intent(in) :: nw
real(8),intent(in) :: gf0(lr_nw)
real(8),intent(in) :: sigc(lr_nw)
real(8),intent(in) :: z1(lr_nw)
real(8),intent(in) :: z2(lr_nw)
real(8),intent(out) :: zaux(lr_nw)
!
integer iw,iw1,sgn
real(8) :: dw,dw0,gtmp,sctmp
real(8) :: w(nw)
real(8),external :: spline3_eval
!
dw0=dreal(lr_w(2)-lr_w(1))
w=0.d0

! define a dense uniform grid 
do iw=1,nw
 w(iw)=dw0*(iw-1)
enddo

do iw1=1,lr_nw
 zaux(iw1)=zaux(iw1)+gf0(iw1)*sigc(1)*dw0*0.5d0
 do iw=2,nw
  sgn=1
  dw=dreal(lr_w(iw1))-w(iw)

  if (dw.lt.-1.d-8) then 
   dw=dw+bhbar
   sgn=-1
  endif

  ! in case -1.d-8 < dw < 0
  if (abs(dw).lt.1.d-8) dw=1.d-10

  gtmp=spline3_eval(lr_nw-1,dreal(lr_w(:)),gf0(:),z1(:),dw)
  sctmp=spline3_eval(lr_nw-1,dreal(lr_w(:)),sigc(:),z2(:),w(iw))
  if (iw.lt.nw) then 
    zaux(iw1)=zaux(iw1)+sgn*gtmp*sctmp*dw0
  else !iw=nw
    zaux(iw1)=zaux(iw1)+sgn*gtmp*sctmp*dw0*0.5d0
  endif
 enddo !iw
enddo !iw1

return
end subroutine
