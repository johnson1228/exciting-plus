subroutine cpe_init(ist,isp,ik,ikloc, m0, m, om, sigm, n, on, ns, os,&
                   &sigx_, vxcnk_,vclnk_, sigman)
use modmain
use mod_nrkp
use mod_linresp
!
implicit none
integer, intent(in) :: ist
integer, intent(in) :: isp
integer, intent(in) :: ik
integer, intent(in) :: ikloc
integer, intent(in) :: m0
integer, intent(in) :: m ! number of input
complex(8), intent(in) :: om(m)  ! omega_m
complex(8), intent(in) :: sigm(m) !sig_m
integer, intent(in) :: n ! number of output
complex(8), intent(in) :: on(n)  ! omega_n
integer, intent(in) :: ns   ! # of sigma
complex(8), intent(in) :: os(ns) ! # of real frequency points for sigma
real(8),intent(in) :: sigx_
real(8),intent(in) :: vxcnk_
real(8),intent(in) :: vclnk_
complex(8), intent(inout) :: sigman(ns) ! sig_n
!
! local var
!
complex(8), allocatable :: amat(:,:),sig1(:)
complex(8) :: tmp1, tmp2, tmp3, tmp4
real(8) :: norm_new,norm0,l_min,x(20),y(20)
complex(8) :: sig0
complex(8), allocatable :: alp(:),alp0(:),grad(:)
complex(8), allocatable :: a_aux1(:,:),a_aux2(:),aux3(:),aux4(:)
integer :: im, in, iter, is, i,m1,j, nls
logical :: exst
integer,parameter :: maxstep=1000
real(8),parameter :: m0wt=4.d0
complex(8) :: z,sig_tmp,g_tmp,g2
complex(8),external :: zdotc
real(8),external :: fit_alp
character*100 :: fname
!
m1=m-m0+1
allocate(amat(m1,n))
allocate(sig1(m1))
allocate(alp(n))
allocate(alp0(n))
allocate(a_aux1(n,n))
allocate(a_aux2(n))
allocate(aux3(m1))
allocate(aux4(n))
allocate(grad(n))
amat=zzero
! init sig0 and alp
alp0=zzero
norm_new=0.d0
x=0.d0
y=0.d0
aux3=zzero
aux4=zzero
a_aux1=zzero
a_aux2=zzero
!
!initial guess for alp and sig0, either read in or using default values
write(fname,'("Sig_c_files/alp_tmp_k",I4.4,"_ist",I4.4)') ik,ist
inquire(file=fname, exist=exst)
if (exst) then
 open(172,file=trim(adjustl(fname)),action='read',form="unformatted",&
      &status="old")
 read(172) alp
 read(172) sig0
 close(172)
else
 write(*,'("alp_tmp file not found!",2I5)') ik,ist
 alp=zzero
 in=int(n/4)
 alp(in:3*in)=dcmplx(cpe_alp0,0.d0)
 sig0=cpe_sig0
endif
!
! construct A(m,n)
do i=1,m1
 im=m0+i-1
 do in=2,n-1
  tmp1=(on(in-1)-om(im))/(on(in-1)-on(in))
  tmp2=(on(in-1)-om(im))/(on(in)-om(im))
  tmp3=(on(in+1)-om(im))/(on(in)-on(in+1))
  tmp4=(om(im)-on(in))/(om(im)-on(in+1))
  amat(i,in)=tmp1*log(tmp2)-tmp3*log(tmp4)
 enddo
enddo

! define sig1=sig0-sigm
do i=1,m1
 im=m0+i-1
 sig1(i)=-sigm(im)+sig0
enddo

! The norm is summed over square deviations / # of iwn
amat(:,:)=amat(:,:)/sqrt(dble(m1))
sig1(:)=sig1(:)/sqrt(dble(m1))

! increase weight for the first four data points
if (m1.gt.4) then
  amat(1:4,:)=amat(1:4,:)*sqrt(m0wt)
  sig1(1:4)=sig1(1:4)*sqrt(m0wt)
endif

! test
if (mpi_grid_root()) then
 write(fname,'("Sig_c_files/sig_test_ib",I4.4)') ist
 open(169,file=fname,status='replace')
 write(fname,'("Sig_c_files/sigr_test_ib",I4.4)') ist
 open(170,file=fname,status='replace') 
endif

write(fname,'("Sig_c_files/greenfunc_test_ik",I4.4,"_ib",I4.4)') ik,ist
open(173,file=fname,status='replace')

! (A^H)*A
call zgemm('C','N',n,n,m1,zone,amat(1,1),m1,amat(1,1),m1,zzero,a_aux1(1,1),n)

! (A^H)*sig1
call zgemv('C',m1,n,zone,amat(1,1),m1,sig1(1),1,zzero,a_aux2(1),1)

! main loop for searching alp
do iter=1,maxstep
  if (iter.eq.1) then
   ! calculate the initial norm
   aux3(:)=sig1(:)
   call zgemv('N',m1,n,zone,amat(1,1),m1,alp(1),1,zone,aux3(1),1)
   tmp1=zdotc(m1,aux3(1),1,aux3(1),1)
   norm_new=dreal(tmp1)
   if (mpi_grid_root()) write(*,'("initial norm:",G19.10)') norm_new
  endif
  ! compute gradient d(norm)/d(alp)
  aux4(:)=a_aux2(:)
  call zgemv('N',n,n,zone,a_aux1(1,1),n,alp(1),1,zone,aux4(1),1)
  do in=1,n
   grad(in)=(2.d0,0.d0)*dreal(aux4(in))
  enddo

 ! prepare for line search
  g2=zdotc(n,grad(1),1,grad(1),1)
  alp0(:)=alp(:)
    
  l_min=10.d0
  if (dreal(g2).lt.1.d-5.and.dreal(g2).gt.1.d-7) then
    l_min=l_min*4.d0
    nls=4
  elseif (dreal(g2).lt.1.d-7) then
    l_min=l_min*8.d0
    nls=4
  else
    nls=6
  endif

  i=0
  norm0=norm_new

! line search for lambda_min
  do while (.true.)
! initialization
   alp(:)=alp0(:)
   i=i+1
   call computeNorm(n,m1,l_min,grad,amat,sig1,alp,norm_new)
   x(i)=l_min
   y(i)=norm_new

!   if (mpi_grid_root()) write(*,*) "i,x,y:",i,x(i),y(i)

   if (i.lt.nls) then
     if (y(i).lt.norm0) then
       if (iter.gt.10) then
         norm0=y(i)
         alp0(:)=alp(:)
         exit
       else
         l_min=l_min/2.d0
       endif
     else
       l_min=l_min/2.d0
     endif
   else
     ! inter-/extra- polation for lambda_min
     l_min=fit_alp(x(i-2),x(i-1),x(i),y(i-2),y(i-1),y(i))
     if (abs(l_min).lt.1.d-09) l_min=0.01d0

     ! compute the new norm
     alp(:)=alp0(:)
     call computeNorm(n,m1,l_min,grad,amat,sig1,alp,norm0)

     if (norm0.gt.y(i)) then
       alp(:)=alp0(:)
       l_min=x(i)
       call computeNorm(n,m1,l_min,grad,amat,sig1,alp,norm0)
     endif

     alp0(:)=alp(:)
     if (mpi_grid_root()) write(*,*) "lambda_min:",l_min
     exit
   endif 
  enddo ! line search

  alp(:)=alp0(:)
  norm_new=norm0
  ! store the array alp
  if (mod(iter,100).eq.0) then
   write(fname,'("Sig_c_files/alp_tmp_k",I4.4,"_ist",I4.4)') ik, ist
   open(172,file=trim(adjustl(fname)),form="unformatted",status="replace")
   write(172) alp
   write(172) sig0
   close(172)
  endif

  if (mpi_grid_root()) &
     &write(*,'(I6,1X,2(G18.12,1X))') iter,dreal(g2),norm_new

  if (iter.gt.1.and.dreal(g2).gt.1.d+4) then
    write(*,*) "gradient is too large!", dreal(g2),ist,ik
    call pstop
  endif

  if ((dreal(g2).lt.cpe_gtor.and.norm_new.lt.cpe_gtor*10.d0)&
      .or.iter.eq.maxstep) then
   if (mpi_grid_root()) write(*,*) "CPE process is done!", iter
   ! warning for not converged
   if (iter.eq.maxstep) write(*,'("iter.eq.maxstep! ",2(I5,1X))') ist,ik
   exit !exit the loop
  endif
enddo !iter

! calculate sigma on real-w axis

if (mpi_grid_root()) then
 write(*,*) "g2,sig0,norm:",dreal(g2),dreal(sig0),norm_new
 do in=1,n
  write(*,*) "in,alp:",in,dreal(alp(in)),dimag(alp(in))
 enddo
endif

! self-energy at E_{ks}
do is=1,ns
 z=os(is)
 sigman(is)=sig0
 do in=2,n-1
  tmp1 = (on(in-1)-z)/(on(in-1)-on(in))
  tmp2 = (on(in-1)-z)/(on(in)-z)
  tmp3 = (on(in+1)-z)/(on(in)-on(in+1))
  tmp4 = (z-on(in))/(z-on(in+1))
  sigman(is)=sigman(is)+alp(in)*(tmp1*log(tmp2)-tmp3*log(tmp4))
 enddo
enddo

! self-energy at iw_m
do i=1,m1
 im=m0+i-1
 z=om(im)
 sig_tmp=sig0
 do in=2,n-1
  tmp1 = (on(in-1)-z)/(on(in-1)-on(in))
  tmp2 = (on(in-1)-z)/(on(in)-z)
  tmp3 = (on(in+1)-z)/(on(in)-on(in+1))
  tmp4 = (z-on(in))/(z-on(in+1))
  sig_tmp=sig_tmp+alp(in)*(tmp1*log(tmp2)-tmp3*log(tmp4))
 enddo
 if (mpi_grid_root()) write(169,'(I4,1X,2G18.6)') im,dreal(sig_tmp)*ha2ev&
                                          &,dimag(sig_tmp)*ha2ev
enddo

! along real-w axis
do im=1,2*n
 !z=(on(im)+on(im+1))*0.5d0+dcmplx(0.d0,lr_eta)
 z=dcmplx(-5.d0,lr_eta)+(im-1)*10.d0/(2*n-1)
 sig_tmp=sig0 
 do in=2,n-1 
   tmp1 = (on(in-1)-z)/(on(in-1)-on(in))
   tmp2 = (on(in-1)-z)/(on(in)-z) 
   tmp3 = (on(in+1)-z)/(on(in)-on(in+1))
   tmp4 = (z-on(in))/(z-on(in+1))
   sig_tmp=sig_tmp+alp(in)*(tmp1*log(tmp2)-tmp3*log(tmp4))
 enddo
 if (mpi_grid_root()) write(170,'(G12.6,1X,2G18.6)') dreal(z)*ha2ev,&
                         & dreal(sig_tmp)*ha2ev,dimag(sig_tmp)*ha2ev

 ! green function along real-w axis!
 i=ist-qpnb(1)+1
 g_tmp=zone/(z-evalsvnr(ist,ik)+vxcnk_-sigx_-vclnk_-sig_tmp)

 write(173,'(G12.6,1X,2G18.6)') dreal(z)*ha2ev,dreal(g_tmp)/ha2ev,&
                                &dimag(g_tmp)/ha2ev
enddo

if (mpi_grid_root()) then
 close(169)
 close(170)
endif

close(173)

deallocate(amat,sig1,alp,a_aux1,a_aux2,aux3,aux4,grad)
deallocate(alp0)
return
end subroutine cpe_init
