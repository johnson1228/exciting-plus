subroutine gen_gw_wmesh(gw_type,emax,nw_se)
use modmain
use mod_linresp
use mod_addons
use mod_addons_q
!
implicit none
integer,intent(in) :: gw_type
real*8,intent(in) :: emax
integer,intent(out) :: nw_se
integer iw,i,nw1
!
! gw_type = 0  <--> task= 802
!         = 1  <--> task= 8022
!
nw_se=2 !default

if (gw_type.eq.0) then
 if (caltype.eq.0) then  !ppa
  lr_nw=2
  if (allocated(lr_w)) deallocate(lr_w)
  allocate(lr_w(lr_nw))
  lr_w=0.d0
 elseif (caltype.eq.1) then !real-axis integration
  if (allocated(lr_w)) deallocate(lr_w)
  allocate(lr_w(lr_nw))
  lr_dw=(emax+2.d0)/(lr_nw-1)
  do iw=1,lr_nw
    lr_w(iw)=dcmplx((iw-1)*lr_dw,lr_eta)
  enddo
 elseif (caltype.eq.2) then ! full self-energy
  if (allocated(lr_w)) deallocate(lr_w)
  allocate(lr_w(lr_nw))
  nw_se=51
  lr_dw=(emax+2.d0)/(lr_nw-1)
  do iw=1,lr_nw
   lr_w(iw)=dcmplx((iw-1)*lr_dw,lr_eta)
  enddo
 elseif (caltype.eq.3) then ! new version of real-axis integration
  if (allocated(lr_w)) deallocate(lr_w)
  allocate(lr_w(lr_nw))
  ! low and high frequency cutoffs
  if (raicut(1).lt.1.d-8.or.raicut(2).lt.1.d-8) then
    raicut(1)=min(emax,2.d0)
    raicut(2)=max(emax+2.d0,gqmax**2)
  endif
  ! in case h cutoff is lower than l cutoff
  if (raicut(2).lt.raicut(1)+1.d-8) raicut(2)=raicut(1)+2.d0

  nw1=lr_nw/4*3
  lr_dw=raicut(1)/(nw1-1)  
  do iw=1,nw1
    lr_w(iw)=dcmplx((iw-1)*lr_dw,lr_eta)
  enddo
  ! high frequency region
  lr_dw=(raicut(2)-raicut(1))/(lr_nw-nw1)
  do iw=1,lr_nw-nw1
    lr_w(iw+nw1)=dcmplx(iw*lr_dw+raicut(1),lr_eta)
  enddo 

 endif !caltype 0:ppa 1:rai 2:full SE 3: new rai

elseif (gw_type.eq.1) then    ! \tau-GW method 
! using uniform power mesh
 lr_nw=2*upm(1)*upm(2)+1
 if (allocated(lr_w)) deallocate(lr_w)
 allocate(lr_w(lr_nw))
 call upmesh(0.d0,bhbar,upm(1),upm(2),lr_w)

 if (mpi_grid_root()) then
  do i=1,2*upm(1)*upm(2)+1
   write(*,*) "i,w(i):",i,lr_w(i)
  enddo
 endif

endif

return
end subroutine

subroutine upmesh(t1,t2,p,u,w)
use modmain
!
implicit none
!
real(8),intent(in) :: t1
real(8),intent(in) :: t2
integer,intent(in) :: p
integer,intent(in) :: u
complex(8),intent(out) :: w(2*u*p+1)
!
integer :: i,j,k
real(8) :: b,b1,wj

w(:)=zzero
b=t2-t1 !original interval
w(1)=t1
k=1
b1=b

do i=1,p-1
 b1=b1*0.5d0
enddo

do i=1,p-1
 do j=1,2*u
  wj=(j*b1)/(2*u)
  if ((i.eq.1).or.(wj.gt.b1*0.5d0+1.d-8)) then
   k=k+1
   w(k)=dcmplx(wj)
  endif
 enddo
 b1=b1*2.d0
enddo

do i=1,u*p
 j=u*p+1+i
 w(j)=b-w(u*p-i+1)
enddo

return
end subroutine 
