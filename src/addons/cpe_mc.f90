subroutine cpe_mc(m, om, sigm, n, on, ns, os, sigman)
use modmain
!
implicit none
integer, intent(in) :: m ! number of input
complex(8), intent(in) :: om(m)  ! omega_m
complex(8), intent(in) :: sigm(m) !sig_m
integer, intent(in) :: n ! number of output
complex(8), intent(in) :: on(n)  ! omega_n
integer, intent(in) :: ns   ! # of sigma
complex(8), intent(in) :: os(ns) ! # of real frequency points for sigma
complex(8), intent(inout) :: sigman(ns) ! sig_n
!
! local var
!
complex*16, allocatable :: a(:,:)
complex*16 :: tmp1, tmp2, tmp3, tmp4
real*8 :: oldA, newA, stepmax, stepmin, step
real*8 :: sig0, oldsig0
real*8, allocatable :: alp(:), oldalp(:)
integer :: im, in
integer :: mc, imc
real*8 :: ran1, weight, ran2
complex*16 :: z
!
mc = 40000
stepmax = 0.2d0
stepmin = 0.01d0
allocate(a(m,n))
allocate(alp(n))
allocate(oldalp(n))
!
! prepare A(m,n)
!
do in=2,n-1
  do im=1,m
    tmp1 = ( on(in-1) - om(im) )/( on(in-1) - on(in) )
    tmp2 = ( on(in-1) - om(im) )/( on(in)   - om(im) )
    tmp3 = ( on(in+1) - om(im) )/( on(in)   - on(in+1) )
    tmp4 = ( om(im)   - on(in) )/( om(im)   - on(in+1) )
    a(im,in) = tmp1 * log(tmp2) - tmp3 * log(tmp4)
!    write(*,*) "im,in,a(im,in)", im,in,a(im,in)
  end do
end do
!
! init sig0 and alp
! 
!call init_random_seed()
!call random_number(sig0)
sig0 = -0.1d0
do in=1,n
  call random_number(alp(in))
  alp(in)=alp(in)*0.001d0
end do
!alp(n/4)=1.0d0
!
! mc loop
!
do imc = 1,mc
  step = stepmin + (stepmax-stepmin)/mc*(mc-imc)
  call random_number(ran1)
  oldsig0 = sig0
  sig0 = oldsig0+(ran1-0.5d0)*step*(oldsig0+0.0d0) !/ sqrt(float(imc))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do in=1,n
    call random_number(ran1)
    oldalp(in) = alp(in)
    alp(in) = oldalp(in)+(ran1-0.5d0)*step*(oldalp(in)+0.0d0) !/sqrt(float(imc))
    if (alp(in).lt.0.0d0) alp(in) = 0.0d0
  end do
  !
  ! calculate the target function
  !
  newA=0.0d0
  do im=1,m
    tmp1 = dcmplx(sig0,0.0d0)
    do in=2,n-1
      tmp1 = tmp1 + a(im,in)*alp(in)
    end do
    if (mpi_grid_root().and.(mc-imc.lt.10.or.imc.lt.10)) &
     write(*,'(A16,6F19.10)')"#om,sigma~,sigma",dreal(om(im)),&
                           &dimag(om(im)),dreal(tmp1),dimag(tmp1),&
                           &dreal(sigm(im)),dimag(sigm(im))
    tmp1 = tmp1-sigm(im)
    weight = 1.0d0
    if (dimag(om(im)).lt.0.1d0) weight = 10.0d0  ! need to check
    newA=newA+tmp1*dconjg(tmp1)*weight
  enddo
  if(mpi_grid_root().and.mod(imc,mc/100).eq.0) &
      write(*,'(A14,2I10,2F19.10)') "#imc,oldA,newA", imc,mc, oldA, newA
  if (imc==1) then
    oldA = newA
  else
    if (newA.lt.oldA) then
      oldA = newA
      oldsig0 = sig0
      oldalp(:)=alp(:)
    else
      sig0 = oldsig0
      alp(:) = oldalp(:)
    end if
  endif
enddo ! imc

do in=1,ns
  z=os(in)
  sigman(in)=dcmplx(oldsig0,0.0d0)
  do im=2,n-1
    tmp1 = (dreal(on(im-1))-z)/(dreal(on(im-1))-dreal(on(im)))
    tmp2 = (dreal(on(im-1))-z)/(dreal(on(im))-z)
    tmp3 = (dreal(on(im+1))-z)/(dreal(on(im))-dreal(on(im+1)))
    tmp4 = (z-dreal(on(im)))/(z-dreal(on(im+1)))
    sigman(in)=sigman(in)+oldalp(im)*(tmp1*log(tmp2)-tmp3*log(tmp4))
  end do
!  write(*,*) "#in,on(in),z,sigman(in)", in, on(in),z, sigman(in)
enddo

deallocate(a,alp,oldalp)
return
end subroutine cpe_mc
