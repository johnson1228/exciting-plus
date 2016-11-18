subroutine get_maxdel(nst,A,B,maxdel)
use modmain
! get the maximum difference (absolute value) between two arrays
implicit none
!
integer, intent(in) :: nst
real(kind=8),intent(in) :: A(nst,nspinor,nkptnr)
real(kind=8),intent(in) :: B(nst,nspinor,nkptnr)
real(kind=8),intent(out) :: maxdel
!
! local variables
integer :: i,isp1,ik
real(kind=8) :: tmp

maxdel=-1.0d6

do ik=1,nkptnr
 do isp1=1,nspinor
  do i=1,nst
   tmp=abs(A(i,isp1,ik)-B(i,isp1,ik))
   if (maxdel.le.tmp) then
     maxdel=tmp
   endif   
  enddo
 enddo
enddo

return
end subroutine
