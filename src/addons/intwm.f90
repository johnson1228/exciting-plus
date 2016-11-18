complex(8) function intwm(w,wl,wh,dimwp)
! calculate \int^{wh}_{wl} ds 1/(w-s+i*dimwp)
implicit none
complex(8), intent(in) :: w
real(8), intent(in) :: wl
real(8), intent(in) :: wh
real(8), intent(in) :: dimwp
!local variable
real(8) :: re, im, arg

arg=((dreal(w)-wh)**2+dimwp**2)/((dreal(w)-wl)**2+dimwp**2)

re=-0.5d0*log(arg)

if (abs(dimwp).gt.1.d-10) then
 im=atan((dreal(w)-wh)/dimwp)-atan((dreal(w)-wl)/dimwp)
else
 im=0.d0
endif

intwm=dcmplx(re,im)

return
end function

