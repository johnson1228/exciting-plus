!
! Numerical Mathematics and Computing, Fifth Edition
! Ward Cheney & David Kincaid
! Brooks/Cole Publ. Co.
! Copyright (c) 2003.  All rights reserved.
! For educational use with the Cheney-Kincaid textbook.
! Absolutely no warranty implied or expressed.
!
! Section 9.2
!
subroutine spline3_coef(n,t,y,z) 
implicit none
!
integer, intent(in) :: n
real(8), intent(in) :: t(n+1)
real(8), intent(in) :: y(n+1)
real(8), intent(out) :: z(n+1)
!
real(8) :: h(n),b(n)
real(8) :: u(n-1),v(n-1)
integer :: i
!
do i = 1,n
 h(i) = t(i+1)-t(i)
 b(i) = (y(i+1)-y(i))/h(i)    
enddo

u(1) = 2.0*(h(1) + h(2))
v(1) = 6.0*(b(2) - b(1))

do i = 2,n-1
 u(i) = 2.d0*(h(i+1) + h(i)) - h(i)**2/u(i-1)     
 v(i) = 6.d0*(b(i+1) - b(i)) - h(i)*v(i-1)/u(i-1) 
enddo

z(n+1)=0.0  

do i= n,2,-1     
 z(i) = (v(i-1)-h(i-1)*z(i+1))/u(i-1)
enddo

z(1)=0.0

return
end subroutine spline3_coef 
  
real(8) function spline3_eval(n,t,y,z,x)
!
integer, intent(in) :: n
real(8), intent(in) :: t(n+1)
real(8), intent(in) :: y(n+1)
real(8), intent(in) :: z(n+1)       
real(8), intent(in):: x
!
real(8) :: h, temp
integer :: i

do i=n,2,-1     
 if (x-t(i).ge.0.d0) exit    
end do
h = t(i+1) - t(i)     
temp = 0.5*z(i)+(x-t(i))*(z(i+1)-z(i))/(6.0*h) 
temp = (y(i+1)-y(i))/h-h*(z(i+1) + 2.0*z(i))/6.0 + (x- t(i))*temp     

spline3_eval = y(i) + (x-t(i))*temp  
return
end function spline3_eval 
