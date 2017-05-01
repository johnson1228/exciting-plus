subroutine dyson_gf(iter,bndrg,gf0,sigc,exxnk,exxvc,vxcnk,vclnk,gf)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp

implicit none
integer,intent(in) :: iter
integer,intent(in) :: bndrg(2,nspinor)
real(8),intent(in) :: gf0(nstsv,lr_nw,nkptnrloc)
real(8),intent(in) :: sigc(nbnd,nspinor,lr_nw,nkptnrloc)
real(8),intent(in) :: exxnk(nbnd,nspinor,nkptnrloc)
real(8),intent(in) :: exxvc(nbnd,nspinor,nkptnrloc)
real(8),intent(in) :: vxcnk(nbnd,nspinor,nkptnrloc)
real(8),intent(in) :: vclnk(nbnd,nspinor,nkptnrloc)
real(8), intent(out) :: gf(nstsv,lr_nw,nkptnrloc)
! local
integer :: iw1,iw2,iw3,iwp,iwr,ist,ikloc,info,i,nw,sgn
complex(8) :: gtmp,ctmp
real(8) :: ztmp,dw1,dele
! begin test
real(8) :: gf0_aux(nbnd,lr_nw,nkptnrloc),gf_aux(nbnd,lr_nw,nkptnrloc)
real(8) :: nelec0,nelec1,dw0,d
complex(8) :: gf0_auxc(nbnd,lr_nw,nkptnrloc),gf_auxc(nbnd,lr_nw,nkptnrloc)
complex(8) :: gf0_w(nbnd,niw,nkptnrloc),gf_w(nbnd,niw,nkptnrloc)
complex(8) :: sigc_w(nbnd,niw,nkptnrloc)
complex(8) :: sigc1(nbnd,lr_nw,nkptnrloc)
real(8) :: z1(lr_nw),z2(lr_nw)
real(8),external :: spline3_eval
real(8),allocatable :: amat(:,:),bmat(:),z_aux(:)
real(8),allocatable :: wd(:),dw(:),zint(:)
integer,allocatable :: ipiv(:)
integer :: iwmap(lr_nw)
integer :: iw,ik,isp1,ibnd,fbnd,ie,signdf,signdn
character*100 :: fname,fname_tot,fspn
!end test
!
! assume the green function is diagonal at this moment
! beyond qpnb(1) and qpnb(2), Green's function is not updated
!
! finer grid points for z_aux
dw0=dreal(lr_w(2)-lr_w(1))
nw=int(dreal(lr_w(lr_nw))/dw0+0.1d0)+1
allocate(wd(nw))
allocate(zint(nw))
allocate(dw(nw))
allocate(amat(nw,nw))
allocate(bmat(nw))
allocate(ipiv(nw))
allocate(z_aux(lr_nw))
!
amat=0.d0
bmat=0.d0
ipiv=0
z_aux=0.d0
iwmap=0
!
! define a dense uniform grid
do iw=1,nw
 wd(iw)=dw0*(iw-1)
enddo

! mapping from coarse grid to dense grid
do iw=1,lr_nw
 do iw1=1,nw
   if (abs(dreal(lr_w(iw))-wd(iw1)).lt.1.d-8) then
     iwmap(iw)=iw1
     exit
   endif
 enddo
 if (iwmap(iw).eq.0) then
  write(*,*) "iwmap .eq. 0, Error!"
  call pstop
 endif
enddo

! calculate dw^{r}
dw(:)=dw0
dw(1)=dw0*0.5d0
dw(nw)=dw0*0.5d0
!
gf(:,:,:)=gf0(:,:,:)
!
! calculate the total number of electrons, should be conserved
nelec0=0.d0

do ikloc=1,nkptnrloc
 do ist=1,nstsv
   if (abs(gf0(ist,lr_nw,ikloc)).gt.occ_tor) &
     nelec0=nelec0-gf0(ist,lr_nw,ikloc)/nkptnr
 enddo
enddo

nelec0=nelec0*occmax
call mpi_grid_reduce(nelec0,1,dims=(/dim_k/),all=.true.)
if (mpi_grid_root()) write(*,'("initial # of electrons:",f8.4)') nelec0

ie=0
signdf=0
signdn=0

do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   ibnd=bndrg(1,isp1)
   fbnd=bndrg(2,isp1)
   do ist=ibnd,fbnd
    i=ist-ibnd+1
    ! evaluate Z(\tau1,\tau2)=Z(\tau1-\tau2)
    ! exchange part
    do iw1=1,lr_nw
      z_aux(iw1)=gf0(ist,iw1,ikloc)*(exxnk(i,isp1,ikloc)+exxvc(i,isp1,ikloc)&
               & -vxcnk(i,isp1,ikloc)+vclnk(i,isp1,ikloc)-dyson_mu)
    enddo !iw1

    ! cubic spline interpolation for gf0 and sigc
    call spline3_coef(lr_nw-1,dreal(lr_w(:)),gf0(ist,:,ikloc),z1(:))
    call spline3_coef(lr_nw-1,dreal(lr_w(:)),sigc(i,isp1,:,ikloc),z2(:))

    ! correlation part
    call cal_zaux(nw,gf0(ist,:,ikloc),sigc(i,isp1,:,ikloc),z1(:),&
                & z2(:),z_aux)
   
    ! cubic spline for z_aux as well
    call spline3_coef(lr_nw-1,dreal(lr_w(:)),z_aux(:),z2(:))
 
    do iw1=1,nw
     bmat(iw1)=spline3_eval(lr_nw-1,dreal(lr_w(:)),gf0(ist,:,ikloc),z1(:),&
                & wd(iw1))
     zint(iw1)=spline3_eval(lr_nw-1,dreal(lr_w(:)),z_aux(:),z2(:),wd(iw1))
    enddo

    ! amat(\tau1,\tau2)
    do iw2=1,nw
      do iw1=1,nw
       sgn=1
       dw1=wd(iw1)-wd(iw2)
       if (dw1.lt.-1.d-8) then
        dw1=dw1+bhbar
        sgn=-1
       endif
       iw=int(dw1/dw0+0.1d0)+1
       amat(iw1,iw2)=-sgn*zint(iw)*dw(iw2)
       if (iw1.eq.iw2) amat(iw1,iw2)=amat(iw1,iw2)+1.d0
      enddo ! iw1
    enddo !iw2

    ! solve for dressed Green function
    call dgesv(nw,1,amat(1,1),nw,ipiv,bmat(1),nw,info)
    if (info.ne.0) write(*,*) "info.ne.0!",info

    do iw=1,lr_nw
     gf(ist,iw,ikloc)=bmat(iwmap(iw))
    enddo

    if (mpi_grid_root()) write(*,'("ist,ikloc:",2I4)') ist,ikloc
   enddo !ist
  enddo !isp1
enddo !ikloc

! calculate the number of electrons from gf
nelec1=0.d0
do ikloc=1,nkptnrloc
  do ist=1,nstsv
    if (abs(gf(ist,lr_nw,ikloc)).gt.occ_tor) &
      & nelec1=nelec1-gf(ist,lr_nw,ikloc)/nkptnr
  enddo
enddo

nelec1=nelec1*occmax
call mpi_grid_reduce(nelec1,1,dims=(/dim_k/),all=.true.)
if (mpi_grid_root()) write(*,'("# of electrons:",f8.4,1X,I4)') nelec1,iter

! Change of total number of electrons
dele=nelec1-nelec0
if (mpi_grid_root()) write(*,'("signdf:",I4,1X,I4)') signdf,ie
if (mpi_grid_root()) write(*,'("signdn:",I4,1X,I4)') signdn,ie
if (mpi_grid_root()) write(*,'("dyson_mu:",f8.4,1X,I4)') dyson_mu,ie
if (mpi_grid_root()) write(*,'("nelec1-nelec0:",f8.4,1X,I4)') dele,ie

! Updated Green's functin gf
if (mpi_grid_root((/dim_q/))) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  write(fname,'("gf_tau_I",I1.1,"_k",I4.4)') iter,ik
  do ist=qpnb(1),qpnb(2)
   write(fspn,'("_ist",I3.3)') ist
   fname_tot="Green_func/"//trim(adjustl(fname))//trim(adjustl(fspn))
   open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
   do iw=1,lr_nw
     write(168,'(2G18.8)') dreal(lr_w(iw))/ha2ev,gf(ist,iw,ikloc)
   enddo
   close(168)
  enddo !ist
 enddo !ikloc
endif

! Get Green's function in Matsubara frequency domain, gf_w
do ist=qpnb(1),qpnb(2)
 i=ist-qpnb(1)+1
 gf_aux(i,:,:)=gf(ist,:,:)
enddo

!G(iw)
gf_auxc(:,:,:)=dcmplx(gf_aux(:,:,:))
call ft_tw3(1,-1,nbnd,lr_nw,nkptnrloc,niw,gf_auxc(:,:,:),gf_w)

! G(iw): gf_w
if (mpi_grid_root((/dim_q/))) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  write(fname,'("gf_iw1_I",I1.1,"_k",I4.4)') iter,ik
  do ist=qpnb(1),qpnb(2)
   i=ist-qpnb(1)+1
   write(fspn,'("_ist",I3.3)') ist
   fname_tot="Green_func/"//trim(adjustl(fname))//trim(adjustl(fspn))
   open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
   do iw=1,niw
     write(168,'(3G18.8)') iw,dreal(gf_w(i,iw,ikloc)),&
                         &dimag(gf_w(i,iw,ikloc))
   enddo
   close(168)
  enddo !ist
 enddo !ikloc
endif

deallocate(amat,z_aux,bmat,ipiv)
deallocate(wd,zint,dw)
return
end subroutine
