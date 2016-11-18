subroutine genvclnk(nst,bndrg,vclnk)
use modmain
use mod_nrkp
use mod_addons_q
use mod_linresp,  only : qpnb,occ_tor
implicit none
!
integer, intent(in) :: nst
integer, intent(in) :: bndrg(2,nspinor)
real(8), intent(out) :: vclnk(nst,nspinor,nkptnrloc)
!
!local variables
integer jst,ik,ikloc,ias,is,ic,l1,l2,l3,io1,io2,ispn,lm1,lm2,lm3,ig,ir
integer m1, m2, m3, nr
integer ibnd, fbnd, isp1, lmax, indx
integer :: inbndrg(nstsv)
real(8) :: t1,t2,rt,cfq,sum1,sum2
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
real(8), allocatable :: jlgqr(:,:,:)
real(8), allocatable :: gqc(:)
real(8), allocatable :: tpgqc(:,:)
real(8), external :: gaunt
real(8), external :: rfinp
real(8) zn(nspecies)
!
real(8), allocatable :: rhomt_(:,:,:)
real(8), allocatable :: rhoir_(:)
real(8), allocatable :: rhonkmt(:,:,:,:,:,:)
real(8), allocatable :: rhonkir(:,:,:,:)
real(8), allocatable :: vclmt_(:,:,:)
real(8), allocatable :: vclir_(:)
real(8), allocatable :: rf(:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
complex(8), allocatable :: gntmp(:,:)
complex(8), allocatable :: zfft(:)
complex(8) zt, zrho0
!
! calculate the Hartree potential due to change in \rho
!
allocate(gntmp(lmmaxapw,lmmaxapw))
allocate(jlgqr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(gqc(ngvec))
allocate(tpgqc(2,ngvec))
allocate(zfft(ngrtot))
allocate(rf(nrmtmax))
if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxvr,ngvec))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngvec,natmtot))
!
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot))
allocate(rhomt_(lmmaxvr,nrmtmax,natmtot))
allocate(vclmt_(lmmaxvr,nrmtmax,natmtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(rhonkmt(lmmaxvr,nrmtmax,natmtot,nst,nspinor,nkptnrloc))
allocate(rhoir_(ngrtot))
allocate(zrhoir(ngrtot))
allocate(zvclir(ngrtot))
allocate(vclir_(ngrtot))
allocate(rhonkir(ngrtot,nst,nspinor,nkptnrloc))
!
vclnk=0.d0
zvclmt=zzero
zvclir=zzero
zrhomt=zzero
zrhoir=zzero
rhomt_=0.d0
rhoir_=0.d0
vclmt_=0.d0
vclir_=0.d0
rhonkmt=0.d0
rhonkir=0.d0
zfft=zzero
inbndrg(:)=0

! check if band index falls in bndrg
do jst=1,nstsv
 do isp1=1,nspinor
  if (jst.ge.bndrg(1,isp1).and.jst.le.bndrg(2,isp1)) then
   inbndrg(jst)=isp1
   exit
  endif
 enddo
enddo

if (mpi_grid_root()) then
 do jst=1,nstsv
  write(*,*) "jst,inbndrg:",jst,inbndrg(jst)
 enddo
endif

! construct the new charge density
do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)

 do jst=1,nstsv
   if (occsvnr(jst,ik).ge.occ_tor.or.inbndrg(jst).gt.0) then
    if (inbndrg(jst).gt.0) indx=jst-bndrg(1,inbndrg(jst))+1
    do lm3=1,lmmaxvr
      gntmp(:,:)=gntyry(lm3,:,:) !Gaunt coefficient
      l3=lm2l(lm3)
      do ias=1,natmtot
        is=ias2is(ias)
        ic=ias2ic(ias)
        rf=0.d0
        do l1=0,lmaxapw; do io1=1,nufr(l1,is)
          do l2=0,lmaxapw; do io2=1,nufr(l2,is)
            if (mod(l1+l2+l3,2).eq.0) then
              zt=zzero
              do ispn=1,nspinor
                do lm2=l2**2+1,(l2+1)**2
                  do lm1=l1**2+1,(l1+1)**2
                    zt=zt+dconjg(wfsvmtnrloc(lm1,io1,ias,ispn,jst,ikloc))*&
                      &wfsvmtnrloc(lm2,io2,ias,ispn,jst,ikloc)*gntmp(lm1,lm2)
                  enddo !lm1
                enddo !lm2
              enddo !ispn
              rf(:)=rf(:)+dreal(zt)*ufr(:,l1,io1,ic)*ufr(:,l2,io2,ic)
            endif
          enddo; enddo !l2, io2
        enddo; enddo !l1, io1
        ! store the orbital square if in bndrg
        if (inbndrg(jst).gt.0) then
         rhonkmt(lm3,:,ias,indx,inbndrg(jst),ikloc)=rf(:)
        endif
        if (occsvnr(jst,ik).ge.occ_tor) &
          &rhomt_(lm3,:,ias)=rhomt_(lm3,:,ias)+rf(:)*occsvnr(jst,ik)
      enddo !ias
    enddo !lm3

    do ispn=1,nspinor
      zfft=zzero
      do ig=1,ngknr(ikloc)
        zfft(igfft(igkignr(ig,ikloc)))=wfsvitnrloc(ig,ispn,jst,ikloc)
      enddo
      call zfftifc(3,ngrid,1,zfft)
      do ir=1,ngrtot
        rt=(abs(zfft(ir))**2)/omega
        if (inbndrg(jst).gt.0) rhonkir(ir,indx,inbndrg(jst),ikloc)=&
              & rhonkir(ir,indx,inbndrg(jst),ikloc)+rt
        if (occsvnr(jst,ik).ge.occ_tor) &
              & rhoir_(ir)=rhoir_(ir)+rt*occsvnr(jst,ik)
      enddo
    enddo !ispn

   endif
 enddo !jst
enddo !ikloc

call mpi_grid_reduce(rhomt_(1,1,1),lmmaxvr*nrmtmax*natmtot,dims=(/dim_k/),&
                      &all=.true.)
call mpi_grid_reduce(rhoir_(1),ngrtot,dims=(/dim_k/),all=.true.)

rhomt_(:,:,:)=rhomt_(:,:,:)/nkptnr
rhoir_(:)=rhoir_(:)/nkptnr

!add core charge density to total \rho
sum1=0.d0
sum2=0.d0
! MT charge density
do ias=1,natmtot
 is=ias2is(ias)
 nr=nrmt(is)
 do ispn=1,nspncr
   do ir=1,nr
     rhomt_(1,ir,ias)=rhomt_(1,ir,ias)+rhocr(ir,ias,ispn)/y00
     fr(ir)=rhocr(ir,ias,ispn)*spr(ir,is)**2
   enddo
   call fderiv(-1,nr,spr(:,is),fr,gr,cf)
   sum1=sum1+fourpi*gr(nr)
 enddo
enddo
! interstitial charge density
do is=1,nspecies
 sum2=sum2+dble(natoms(is))*(4.d0*pi/3.d0)*rmt(is)**3
enddo

rhoir_(:)=rhoir_(:)+(chgcr-sum1)/(omega-sum2)

!renormalize the charge density
!sum1=0.d0
!do ias=1,natmtot
!  is=ias2is(ias)
!  do ir=1,nrmt(is)
!   fr(ir)=rhomt_(1,ir,ias)*spr(ir,is)**2
!  enddo
!  call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
!  sum1=sum1+fourpi*y00*gr(nrmt(is))
!enddo
! find the interstitial charge
!sum2=0.d0
!do ir=1,ngrtot
!  sum2=sum2+rhoir_(ir)*cfunir(ir)
!enddo
! total charge
!sum1=sum1+sum2*omega/dble(ngrtot)

! error in average density
!t1=(chgtot-sum1)/omega
! add the constant difference to the density
!t2=t1/y00
!do ias=1,natmtot
! is=ias2is(ias)
! do ir=1,nrmt(is)
!   rhomt_(1,ir,ias)=rhomt_(1,ir,ias)+t2
! enddo
!enddo
! interstitial charge density
!rhoir_(:)=rhoir_(:)+t1

!convert real MT charge density to complex spherical harmonic expansion
do ias=1,natmtot
 is=ias2is(ias)
 do ir=1,nrmt(is)
  call rtozflm(lmaxvr,rhomt_(:,ir,ias),zrhomt(:,ir,ias))
 enddo
enddo
! convert to complex array
zrhoir(:)=dcmplx(rhoir_(:))

do ig=1,ngvec
! G+q-vector length and (theta, phi) coordinates
  call sphcrd(vgc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vectors
  call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
end do

!compute the required spherical Bessel functions
lmax=lmaxvr+npsden+1
call genjlgpr(lmax,gqc,jlgqr)
! structure factor
call gensfacgp(ngvec,vgc,ngvec,sfacgq)

!solve poisson equation for Hartree potential
call zpotcoul(nrmt,nrmtmax,spnrmax,spr,1,gqc,jlgqr,ylmgq,sfacgq,spzn,&
     &zrhomt,zrhoir,zvclmt,zvclir,zrho0)

! convert complex MT potential to real spherical harmonic expansion
do ias=1,natmtot
  is=ias2is(ias)
  do ir=1,nrmt(is)
   call ztorflm(lmaxvr,zvclmt(:,ir,ias),vclmt_(:,ir,ias))
  end do
end do
! store complex interstitial potential in real array
vclir_(:)=dble(zvclir(:))

! compute vclnk
do ikloc=1,nkptnrloc
 do isp1=1,nspinor
   do jst=1,nst
    vclnk(jst,isp1,ikloc)=rfinp(1,rhonkmt(:,:,:,jst,isp1,ikloc),&
                  vclmt_,rhonkir(:,jst,isp1,ikloc),vclir_)
    t1=rfinp(1,rhonkmt(:,:,:,jst,isp1,ikloc),&
                  vclmt,rhonkir(:,jst,isp1,ikloc),vclir)
    if (mpi_grid_root()) write(*,*) "t1:",t1,jst,ikloc
    vclnk(jst,isp1,ikloc)=vclnk(jst,isp1,ikloc)-t1
   enddo !jst
 enddo !isp1
enddo !ikloc

! wait for all processors
call mpi_grid_barrier()

deallocate(gntmp,zfft,zrhomt,zrhoir,zvclmt,zvclir,rhonkmt,rhonkir)
deallocate(rhomt_,rhoir_,vclmt_,vclir_)
deallocate(gqc,jlgqr,ylmgq,sfacgq,tpgqc,rf)
return
end subroutine
