subroutine rai_self_energy_new(iq,iqloc,istep,nw_se,corr_se_aux)
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
integer ikloc,iw,ig,ig1,ig2,n,jk,i,j,ik,iz,iwloc
integer ist1,ist2,nwloc,nvqr
complex(8) zw(2),zt1,gf0_aux
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsinv(:,:)
complex(8), allocatable :: ame(:,:)
complex(8), allocatable :: zm(:,:)
!complex(8), allocatable :: corr_aux2(:,:,:,:)
complex(8) :: w
integer :: isp1, ibnd, fbnd
integer :: nst
integer :: bndrg(2,nspinor)
real(8) :: dw(lr_nw)
character*100 fname
!
allocate(chi0(ngq(iq),ngq(iq)))
allocate(epsinv(ngq(iq),ngq(iq))) 
allocate(krnl(ngq(iq),ngq(iq)))
allocate(ame(ngq(iq),nstsv*nstsv))
allocate(zm(ngq(iq),ngq(iq)))
!allocate(corr_aux2(nw_se,nbnd,nspinor,nkptnr))
!
krnl=zzero
epsinv=zzero
zm=zzero
!corr_aux2=zzero
nst=(int(chgval/2.0)+nebd_se)*nspinor
nst=min(nst,nstsv)
nwloc=mpi_grid_map(lr_nw,dim_q)
nvqr=nkpt+nvq0-1

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

! GW0,need to modify for better performance
if (scgwni.gt.1) then
 write(fname,'("Temp_files/rai_eps_q",I4.4,"_ipq",I4.4)') iq, &
      & mpi_grid_dim_pos(dim_q)
 if (istep.eq.1) then
  open(171,file=trim(adjustl(fname)),action='write',form='unformatted', &
&       status='replace')
 elseif (istep.gt.1) then
  open(171,file=trim(adjustl(fname)),action='read',form='unformatted', &
&       status='old')
 endif
endif

! interval used in the integration, applicable for non-uniform mesh
dw(1)=dreal(lr_w(2)-lr_w(1))/2.d0
dw(lr_nw)=dreal(lr_w(lr_nw)-lr_w(lr_nw-1))/2.d0
do iw=2,lr_nw-1
  dw(iw)=dreal(lr_w(iw+1)-lr_w(iw-1))/2.d0
enddo

if (mpi_grid_root()) then
 do iw=1,lr_nw
  write(*,*) iw,dreal(lr_w(iw))*ha2ev,dw(iw)*ha2ev
 enddo
endif

del_e=dw(1)*2.d0
if (mpi_grid_root()) write(*,*) "del_e:",del_e*ha2ev,iqloc

! parallelize along dim_q
do iwloc=1,nwloc   ! calculate eps^{-1} at lr_nw points
 iw=mpi_grid_map(lr_nw,dim_q,loc=iwloc)
 if (mpi_grid_root()) write(*,*) "iw:",iw, nwloc
 chi0=zzero
 !calculate the dielectric matrix only when istep=1
 if (istep.eq.1) then
  do ikloc=1,nkptnrloc
   if (nmegqblh(ikloc).gt.0) then
! for each k-point : sum over interband transitions
    call genchi0blh(ikloc,ngq(iq),lr_w(iw),chi0)
   endif
  enddo
  call mpi_grid_reduce(chi0(1,1),ngq(iq)*ngq(iq),dims=(/dim_k/),all=.true.)
  chi0(:,:)=chi0(:,:)/nkptnr/omega
! compute epsilon=1-chi0*v
  epsinv(:,:)=zzero
  do ig=1,ngq(iq)
   epsinv(ig,ig)=zone
  enddo
  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),-zone,chi0,ngq(iq),&
             &krnl,ngq(iq),zone,epsinv(1,1),ngq(iq))

! inverse epsilon: eps^{-1}=(1-chi0*v)^{-1} = 1+chi*v    !check
  call invzge(epsinv(1,1),ngq(iq))

  do ig=1,ngq(iq)
! find frequency dependent part of inverse epsilon
   epsinv(ig,ig)=epsinv(ig,ig)-zone
  enddo
 endif ! istep=1

! for GW0 calculation 
 if (scgwni.gt.1) then
  if (mpi_grid_root((/dim_k/)).and.(istep.eq.1)) then
   write(171) epsinv
  elseif (istep.gt.1) then
   read(171) epsinv
  endif
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
    if ((n.gt.fbnd).or.(n.lt.ibnd)) cycle    ! nk

    zw(1)=dcmplx(evalsvnr(n,ik),0.d0)   !w
    zw(2)=dcmplx(evalsvnr(n,ik)+del_e,0.d0)

    j=bamegqblh(1,i,ikloc)  !m
    if ((j.gt.nst/nspinor).and.(j.le.nstsv/nspinor)) cycle
    if (j.gt.(nst/nspinor+nstsv/nspinor)) cycle

   ! G0W0 and real-axis integreal
    do iz=1,nw_se  !w 
      w=zw(iz)+lr_w(iw)-evalsvnr(j,jk)
      gf0_aux=occsvnr(j,jk)/dconjg(w)+(occmax-occsvnr(j,jk))/w
      w=zw(iz)-lr_w(iw)-evalsvnr(j,jk)
      gf0_aux=gf0_aux+occsvnr(j,jk)/w+(occmax-occsvnr(j,jk))/dconjg(w)

      do ig2=1,ngq(iq)
       do ig1=1,ngq(iq)
         zm(ig1,ig2)=gf0_aux*epsinv(ig1,ig2)*dw(iw)*wtvhgq(ig1,iq)
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

if (scgwni.gt.1) close(171)

deallocate(chi0)
deallocate(krnl)
deallocate(epsinv)
deallocate(ame)
deallocate(zm)
!deallocate(corr_aux2)
return
end subroutine
