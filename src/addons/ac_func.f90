complex(8) function ac_func(n,coe,wn,z)
use modmain
use mod_linresp
!
implicit none
!
integer,intent(in) :: n
complex(8),intent(in) :: coe(n)
complex(8),intent(in) :: wn(n)
complex(8),intent(in) :: z
!
integer :: i
complex(8) :: tmp1,tmp2
!
!b=zzero
tmp1=zzero

do i=1,n-1
 tmp2=coe(n-i+1)*(z-wn(n-i))/(zone+tmp1)
 tmp1=tmp2
enddo

ac_func=coe(1)/(zone+tmp1)

return
end function

