subroutine full_self_energy(iq,iqloc,istep,nw_se,corr_se_aux)
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
integer ikloc,iw,ig,ig1,ig2,n,jk,i,j,ik,iz
integer ist1,ist2
complex(8) zw(nw_se),zt1,intw_aux1,intw_aux2
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsinv(:,:,:)
complex(8), allocatable :: ame(:,:)
complex(8), allocatable :: amat(:,:,:)
complex(8), allocatable :: zm(:,:)
complex(8), external :: intwp,intwm
complex(8) :: pp,pm,mp,mm
complex(8) :: w
real(8) :: wh, wl
integer :: isp1, ibnd, fbnd
integer :: bndrg(2,nspinor)
integer :: nst, st
!

allocate(chi0(ngq(iq),ngq(iq)))
allocate(epsinv(ngq(iq),ngq(iq),2))  ! remove frequency dependency
allocate(krnl(ngq(iq),ngq(iq)))
allocate(ame(ngq(iq),nstsv*nstsv))
allocate(amat(ngq(iq),ngq(iq),2))
allocate(zm(ngq(iq),ngq(iq)))
!
krnl=zzero
epsinv=zzero
amat=zzero
zm=zzero
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

do ig=1,ngq(iq)
  krnl(ig,ig)=dcmplx(vhgq(ig,iq))
enddo

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

do iz=1,nw_se
  zw(iz)=dcmplx((-(nw_se-1)/2.d0+iz-1)/ha2ev,lr_eta)
enddo

do iw=1,lr_nw   ! calculate eps^{-1} at lr_nw points
 chi0=zzero
 do ikloc=1,nkptnrloc
  if (nmegqblh(ikloc).gt.0) then
! for each k-point : sum over interband transitions
   call genchi0blh(ikloc,ngq(iq),lr_w(iw),chi0)
  endif
 enddo
 call mpi_grid_reduce(chi0(1,1),ngq(iq)*ngq(iq),dims=(/dim_k/),all=.true.)
 chi0(:,:)=chi0(:,:)/nkptnr/omega
! compute epsilon=1-chi0*v
 do ig=1,ngq(iq)
  epsinv(ig,ig,2)=zone
 enddo
 call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),-zone,chi0,ngq(iq),&
           &krnl,ngq(iq),zone,epsinv(1,1,2),ngq(iq))

! inverse epsilon: eps^{-1}=(1-chi0*v)^{-1} = 1+chi*v    !check
 call invzge(epsinv(1,1,2),ngq(iq))

 do ig=1,ngq(iq)
! find frequency dependent part of inverse epsilon
  epsinv(ig,ig,2)=epsinv(ig,ig,2)-zone
 enddo
  
 if (iw.gt.1) then
  !slope
   amat(:,:,1)=(epsinv(:,:,2)-epsinv(:,:,1))/&
                    (lr_w(iw)-lr_w(iw-1))
  !y-intercept
   amat(:,:,2)=-amat(:,:,1)*lr_w(iw)+epsinv(:,:,2)
 endif
 !epsinv(iw-1)
 epsinv(:,:,1)=epsinv(:,:,2)

! A,B are obtained now!
 if (iw.eq.1) cycle ! start with iw > 1

! range for the convolution 
 wh=dreal(lr_w(iw))
 wl=dreal(lr_w(iw-1))

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
    if ((n.gt.fbnd).or.(n.lt.ibnd)) cycle    ! nk

    j=bamegqblh(1,i,ikloc)  !m
    if ((j.gt.nst/nspinor).and.(j.le.nstsv/nspinor)) cycle
    if (j.gt.(nst/nspinor+nstsv/nspinor)) cycle

 ! G0W0 and real-axis integreal
    do iz=1,nw_se  !w  
     w=zw(iz)-evalsvnr(j,jk)
     pp=intwp(w,wl,wh,lr_eta)
     mp=intwm(w,wl,wh,lr_eta)
     intw_aux1=occsvnr(j,jk)*dconjg(w)*dconjg(mp-pp)+&
              &(occmax-occsvnr(j,jk))*w*(mp-pp)
     intw_aux2=occsvnr(j,jk)*dconjg(pp+mp)+(occmax-occsvnr(j,jk))*(pp+mp)

     do ig2=1,ngq(iq)
       do ig1=1,ngq(iq)
        zm(ig1,ig2)=(intw_aux1*amat(ig1,ig2,1)+intw_aux2*amat(ig1,ig2,2))&
                  &*wtvhgq(ig1,iq)
       enddo
     enddo

     do ig2=1,ngq(iq)
       zt1=zzero
       do ig1=1,ngq(iq)
        zt1=zt1+dconjg(ame(ig1,i))*zm(ig1,ig2)
      enddo
      corr_se_aux(iz,n-ibnd+1,isp1,ik,iqloc)=&
          corr_se_aux(iz,n-ibnd+1,isp1,ik,iqloc)+zt1*ame(ig2,i)
     enddo
    enddo !iz  
   enddo  !i
  enddo ! isp1
 enddo !ikloc

enddo !iw

deallocate(chi0)
deallocate(krnl)
deallocate(epsinv)
deallocate(ame)
deallocate(amat)
deallocate(zm)
return
end subroutine

