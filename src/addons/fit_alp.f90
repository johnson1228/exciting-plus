real(8) function fit_alp(x1,x2,x3,f1,f2,f3)
use modmain
use mod_linresp
!
implicit none
real(8),intent(in) :: x1
real(8),intent(in) :: x2
real(8),intent(in) :: x3
real(8),intent(in) :: f1
real(8),intent(in) :: f2
real(8),intent(in) :: f3
!
real(8) :: a1,a2,a3

a1=(f3-f2)/(x3-x2)-(f2-f1)/(x2-x1)
a1=a1/(x3-x1)

if (a1.lt.-1.d-10) then
 fit_alp=x3/2.d0
 return
endif

a2=(f3-f2)/(x3-x2)-a1*(x3+x2)
a3=f3-a1*x3**2-a2*x3

fit_alp=-a2/(2.d0*a1)
return
end function fit_alp
