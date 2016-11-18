subroutine computeNorm(n,m1,l,g,amat,sig1,alp,norm)
use modmain
!
implicit none
!
integer,intent(in) :: n
integer,intent(in) :: m1
real(8),intent(in) :: l
complex(8),intent(in) :: g(n)
complex(8),intent(in) :: amat(m1,n)
complex(8),intent(in) :: sig1(m1)
complex(8),intent(inout) :: alp(n)
real(8),intent(out) :: norm
!
complex(8),external :: zdotc
complex(8) :: tmp1,aux3(m1)
integer :: in
!
do in=1,n
  tmp1=alp(in)-l*g(in)
  if (dreal(tmp1).gt.0.d0) then
    alp(in)=dcmplx(dreal(tmp1),0.d0)
  else
    alp(in)=zzero
  endif
enddo

! compute the new norm
aux3(:)=sig1(:)
call zgemv('N',m1,n,zone,amat(1,1),m1,alp(1),1,zone,aux3(1),1)
tmp1=zdotc(m1,aux3(1),1,aux3(1),1)
norm=dreal(tmp1)

return
end subroutine
