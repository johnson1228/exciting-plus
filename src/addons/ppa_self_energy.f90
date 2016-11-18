subroutine ppa_self_energy(iq,iqloc,istep,nw_se,corr_se_aux)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer, intent(in) :: iq
integer, intent(in) :: iqloc
integer, intent(in) :: istep
integer, intent(in) :: nw_se
complex(8), intent(inout) :: corr_se_aux(nw_se,nbnd,nspinor,nkptnr,nkpt+nvq0-1)
!
integer ikloc,iw,ig,ig1,ig2,n,jk,i,j,ik
integer ist1,ist2
real(8) ppa_e0
complex(8) zw,zt1
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: ppa_w(:,:)
complex(8), allocatable :: ppa_r(:,:)
complex(8), allocatable :: zm(:,:)
complex(8), allocatable :: ame(:,:)
integer :: isp1, ibnd, fbnd
integer :: nst
integer :: bndrg(2,nspinor)
character*100 fname
!
allocate(chi0(ngq(iq),ngq(iq)))
allocate(epsinv(ngq(iq),ngq(iq),2))
allocate(ppa_w(ngq(iq),ngq(iq)))
allocate(ppa_r(ngq(iq),ngq(iq)))
allocate(krnl(ngq(iq),ngq(iq)))
allocate(zm(ngq(iq),ngq(iq)))
allocate(ame(ngq(iq),nstsv*nstsv))
!
ppa_e0=1.d0

nst=(int(chgval/2.0)+nebd_se)*nspinor
nst=min(nst,nstsv)

do isp1=1,nspinor
 if (isp1.eq.1) then
  bndrg(1,isp1)=qpnb(1)
  bndrg(2,isp1)=qpnb(2)
 elseif (isp1.eq.2) then   !spin-polarized
  bndrg(1,isp1)=qpnb(1)+int(nstsv/2)
  bndrg(2,isp1)=qpnb(2)+int(nstsv/2)
 endif
enddo

! for the use of test, change them later
!
krnl=zzero
do ig=1,ngq(iq)
  krnl(ig,ig)=vhgq(ig,iq)
enddo
epsinv=zzero

!calculate the dielectric matrix only when istep=1
if (istep.eq.1) then 
 do iw=1,2
   if (iw.eq.1) then
     zw=dcmplx(0.d0,lr_eta)
   else
     zw=dcmplx(0.d0,ppa_e0)
   endif
   chi0=zzero
   do ikloc=1,nkptnrloc
     if (nmegqblh(ikloc).gt.0) then
! for each k-point : sum over interband transitions
       call genchi0blh(ikloc,ngq(iq),zw,chi0)
     endif
   enddo
   call mpi_grid_reduce(chi0(1,1),ngq(iq)*ngq(iq),dims=(/dim_k/),all=.true.)
   chi0(:,:)=chi0(:,:)/nkptnr/omega
! compute epsilon=1-chi0*v
   do ig=1,ngq(iq)
     epsinv(ig,ig,iw)=zone
   enddo
   call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),-zone,chi0,ngq(iq),&
     &krnl,ngq(iq),zone,epsinv(1,1,iw),ngq(iq))
! inverse epsilon: eps^{-1}=(1-chi0*v)^{-1} = 1+chi*v    !check
   call invzge(epsinv(1,1,iw),ngq(iq))
! find frequency dependent part of inverse epsilon
   do ig=1,ngq(iq)
     epsinv(ig,ig,iw)=epsinv(ig,ig,iw)-zone
   enddo
 enddo

! compute coefficients of the plasmon-pole approximation
 do ig1=1,ngq(iq)
   do ig2=1,ngq(iq)
     ppa_w(ig1,ig2)=ppa_e0/sqrt(epsinv(ig1,ig2,1)/epsinv(ig1,ig2,2)-1.d0)
     ppa_w(ig1,ig2)=dreal(ppa_w(ig1,ig2))
     ppa_r(ig1,ig2)=-epsinv(ig1,ig2,1)*ppa_w(ig1,ig2)/2.d0
     if (real(epsinv(ig1,ig2,1)/epsinv(ig1,ig2,2)).le.1.d0) &
       ppa_w(ig1,ig2)=(1.d0,0.d0) !yambo
   enddo
 enddo
endif !istep=1

if (scgwni.gt.1) then
 write(fname,'("Temp_files/ppa_mat_q",I4.4)') iq
 if (mpi_grid_root().and.(istep.eq.1)) then
   open(171,file=trim(adjustl(fname)),action='write',form='unformatted', &
&       status='replace')
   write(171) ppa_w,ppa_r
   close(171)
 elseif (istep.gt.1) then
   open(171,file=trim(adjustl(fname)),action='read',form='unformatted', &
&       status='old')
   read(171) ppa_w,ppa_r
   close(171)
 endif
endif

if (iq.eq.1.and.mpi_grid_root()) then
 do ig=1,ngq(iq)
   write(*,*) "ppa_w(ig,ig),ig:",ppa_w(ig,ig),ig
 enddo
endif

! wait for all the processors
call mpi_grid_barrier()

! restore q->0 matrix elements
 if (vq_gamma(iq)) then
   do ikloc=1,nkptnrloc
     ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
     do ig=1,ngq(iq)
       if (igqig(ig,iq).eq.1) then
         amegqblh(:,ig,ikloc)=zzero
         do i=1,nmegqblh(ikloc)
           ist1=bmegqblh(1,i,ikloc)
           ist2=bmegqblh(2,i,ikloc)
           if (ist1.eq.ist2) then
             amegqblh(i,ig,ikloc)=zone
           endif
         enddo
       endif
     enddo
   enddo
 endif

do ikloc=1,nkptnrloc
   ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
   jk=idxkq(3,ik)  !k'=k-q
   ! change order of indices
   do ig=1,ngq(iq)
     ame(ig,:)=amegqblh(:,ig,ikloc)
   enddo

   do isp1=1,nspinor
     ibnd=bndrg(1,isp1)
     fbnd=bndrg(2,isp1)

     do i=1,namegqblh(ikloc)
       n=bamegqblh(2,i,ikloc)
       j=bamegqblh(1,i,ikloc)
       if (n.gt.fbnd.or.n.lt.ibnd) cycle    !
       if (j.gt.nst) cycle
 ! G0W0
       if (caltype.eq.0) then
         lr_w(1)=dcmplx(evalsvnr(n,ik),lr_eta)
         lr_w(2)=dcmplx(evalsvnr(n,ik)+del_e,lr_eta)
       endif

       do iw=1,nw_se
         j=-1
         if (bamegqblh(1,i,ikloc).ne.j) then
           j=bamegqblh(1,i,ikloc)

           do ig1=1,ngq(iq)
             do ig2=1,ngq(iq)
               zm(ig1,ig2)=-wtvhgq(ig1,iq)*ppa_r(ig1,ig2)*(occsvnr(j,jk)/(evalsvnr(j,jk)-dconjg(lr_w(iw))-ppa_w(ig1,ig2))+&
                 &(occmax-occsvnr(j,jk))/(evalsvnr(j,jk)-lr_w(iw)+ppa_w(ig1,ig2)))/occmax
             enddo
           enddo
         endif
         do ig2=1,ngq(iq)
           zt1=zzero
           do ig1=1,ngq(iq)
             zt1=zt1+dconjg(ame(ig1,i))*zm(ig1,ig2)
           enddo
           corr_se_aux(iw,n-ibnd+1,isp1,ik,iqloc)= &
                 corr_se_aux(iw,n-ibnd+1,isp1,ik,iqloc)+zt1*ame(ig2,i)
         enddo
       enddo  !iw
     enddo  !i
   enddo !isp1
enddo !ikloc

deallocate(chi0)
deallocate(ppa_w)
deallocate(ppa_r)
deallocate(krnl)
deallocate(epsinv)
deallocate(zm)
deallocate(ame)
return
end subroutine
