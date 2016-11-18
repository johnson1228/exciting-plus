subroutine sigma_ac(fnum,n1,n2,n3,nk,bndrg,evalmap,neval,sigx,vxcnk,&
                     &vclnk,sigc,sigcr)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer,intent(in) :: fnum
integer,intent(in) :: n1
integer,intent(in) :: n2
integer,intent(in) :: n3
integer,intent(in) :: nk
integer,intent(in) :: bndrg(2,nspinor)
integer,intent(in) :: evalmap(n1,n2,nk)
integer,intent(in) :: neval(n2,nk)
real(8),intent(in) :: sigx(n1,n2,nk)
real(8),intent(in) :: vxcnk(n1,n2,nk)
real(8),intent(in) :: vclnk(n1,n2,nk)
complex(8),intent(in) :: sigc(n1,n2,n3,nk)
complex(8),intent(out) :: sigcr(2,n1,n2,nk)
!
integer :: iw1,iw2,ist1,isp1,ikloc,n,ik,ibnd,fbnd,i,j
integer :: iloc,nevalloc(n2,nk)
complex(8),allocatable :: gn(:,:),ctmp(:,:,:,:)
complex(8),allocatable :: sigcoe(:,:,:,:)
complex(8) :: enk,wn(n3),zn(2)
!test
integer :: n0
complex(8) :: rwn(2*cpe_N+1)  ! change it later
!complex(8),allocatable :: tmp1(:,:,:,:)
!character*100 :: fname,fspn,fname_tot
complex(8),external :: ac_func
!
allocate(gn(n3,n3))
allocate(sigcoe(n3,n1,n2,nk))
allocate(ctmp(n3,n1,n2,nk))
!
!test
!allocate(tmp1(n4,n1,n2,n3))

sigcoe=zzero
sigcr=zzero
gn=zzero
ctmp=zzero
wn=zzero
!tmp1=zzero

!change order of indices
do ikloc=1,nkptnrloc
 do iw1=1,n3
  ctmp(iw1,:,:,ikloc)=sigc(:,:,iw1,ikloc)
 enddo
enddo

! iw mesh
n=0
do iw1=-int((n3-1)/2),int(n3/2)
 n=n+1
 wn(n)=zi*(2.d0*iw1+1)*pi/bhbar
 if (iw1.eq.0) n0=n
enddo

if (ac_sigma.eq.1) then
! pade approximation
  if (mpi_grid_root()) write(fnum,'("Perform Pade approximation for &
                        &the self-energy")')
  do ikloc=1,nk
   do isp1=1,n2
    do ist1=1,n1
     do iw1=1,n3
      do iw2=iw1,n3
       if (iw1.eq.1) then
        gn(1,iw2)=ctmp(iw2,ist1,isp1,ikloc)
       else
        gn(iw1,iw2)=(gn(iw1-1,iw1-1)-gn(iw1-1,iw2))/&
                 &(wn(iw2)-wn(iw1-1))/gn(iw1-1,iw2)
       endif
      enddo !iw2
! gn(i,i) are the coefficients of the fitted self-energy
      sigcoe(iw1,ist1,isp1,ikloc)=gn(iw1,iw1)
     enddo !iw1
    enddo !ist1
   enddo !ist2
  enddo !ikloc

  do ikloc=1,nk
   ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
   do isp1=1,n2
    ibnd=bndrg(1,isp1)
    do ist1=1,n1
     do iw1=1,2
      enk=dcmplx(evalsvnr(ibnd+ist1-1,ik)+(iw1-1)*del_e,lr_eta)
      sigcr(iw1,ist1,isp1,ikloc)=ac_func(n3,sigcoe(:,ist1,isp1,ikloc),wn,enk)
     enddo
    enddo
   enddo
  enddo
elseif (ac_sigma.eq.2.or.ac_sigma.eq.3) then
! continuous pole expansion
 if (mpi_grid_root()) write(fnum,'("Perform continuous-pole expansion for &
                        &the self-energy")')

! bands are parallelized along dim_q
 do ikloc=1,nk
  do isp1=1,n2
    nevalloc(isp1,ikloc)=mpi_grid_map(neval(isp1,ikloc),dim_q)
    if (mpi_grid_root()) &
    write(*,*) "isp1,ikloc,nevalloc:",isp1,ikloc,nevalloc(isp1,ikloc)
  enddo
 enddo
! define \omega_n along the real-w axis
 n=0
 do iw1=-cpe_N,cpe_N
  n=n+1
  rwn(n)=dcmplx(iw1)*cpe_delta/cpe_N
 enddo

 do ikloc=1,nk
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,n2
   do iloc=1,nevalloc(isp1,ikloc)
     ! parallel over j along q_dim
     i=mpi_grid_map(neval(isp1,ikloc),dim_q,loc=iloc)
     ibnd=evalmap(i,isp1,ikloc)

     if (ibnd.gt.bndrg(2,isp1)) cycle
     if (mpi_grid_root()) write(*,*) "ibnd:",ibnd

     j=ibnd-bndrg(1,isp1)+1
     if (mpi_grid_root()) write(*,*) "j:",j

     zn(1)=dcmplx(evalsvnr(ibnd,ik),lr_eta)
     zn(2)=dcmplx(evalsvnr(ibnd,ik)+del_e,lr_eta)

     if (ac_sigma.eq.2) then
      ! Monte-Carlo like algorithm of CPE
      if (mpi_grid_root()) write(fnum,'("Monte-Carlo like algorithm!")')

      call cpe_mc(n3,wn,ctmp(:,j,isp1,ikloc),n,rwn,2,zn,&
                 &sigcr(:,j,isp1,ikloc))
     elseif (ac_sigma.eq.3) then
      ! original algorithm of CPE
      if (mpi_grid_root()) write(fnum,'("Frank-Wolfe like algorithm!")')
      !
      call cpe_init(ibnd,isp1,ik,ikloc,n0,n3,wn,ctmp(:,j,isp1,ikloc),n,&
                 &rwn,2,zn,sigx(j,isp1,ikloc),vxcnk(j,isp1,ikloc),&
                 &vclnk(j,isp1,ikloc),sigcr(:,j,isp1,ikloc))
     endif

     ! perform CPE once for degenerate states
     if (i.lt.neval(isp1,ikloc)) then
      fbnd=evalmap(i+1,isp1,ikloc)-1
      if (mpi_grid_root()) write(*,*) "fbnd:",fbnd
      if (fbnd.gt.ibnd) then
       do ist1=j+1,fbnd-bndrg(1,isp1)+1
        sigcr(:,ist1,isp1,ikloc)=sigcr(:,j,isp1,ikloc)
       enddo
      endif
     endif

   enddo !iloc
  enddo !isp1
 enddo !ikloc

 call mpi_grid_reduce(sigcr(1,1,1,1),2*n1*n2*nk,dims=(/dim_q/),all=.true.)
endif

!test
!do ikloc=1,n3
! do isp1=1,n2
!  do ist1=1,n1
!   do iw1=1,n4
!    enk=wn(iw1)
!    tmp1(iw1,ist1,isp1,ikloc)=ac_func(n4,sigcoe(:,ist1,isp1,ikloc),wn,enk)
!   enddo
!  enddo
! enddo
!enddo

! do ikloc=1,n3
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  write(fname,'("sigc_iw2_k",I4.4)') ik
!  do ist1=qpnb(1),qpnb(2)
!   write(fspn,'("_ist",I3.3)') ist1
!   fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
!   open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
!   i=ist1-qpnb(1)+1
!   do iw1=1,n4
!    write(168,'(3G18.8)') iw1,dreal(tmp1(iw1,i,1,ikloc))*ha2ev,&
!          & dimag(tmp1(iw1,i,1,ikloc))*ha2ev
!   enddo
!   close(168)
!  enddo !ist1
! enddo !ikloc

deallocate(gn,ctmp,sigcoe)
!deallocate(tmp1)
return
end subroutine

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
