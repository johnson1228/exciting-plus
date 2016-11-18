subroutine genexxnk(iter,nst,ex_nk,ex_vc)
use modmain
use mod_nrkp
use mod_linresp
use mod_addons_q
implicit none
integer, intent(in) :: iter
integer, intent(in) :: nst
real(8), intent(out) :: ex_nk(nst,nspinor,nkptnrloc) 
real(8), intent(out) :: ex_vc(nst,nspinor,nkptnrloc)
!
integer ist,jst,ik,ikloc,ias,is,ic,l1,l2,l3,io1,io2,ispn,lm1,lm2,lm3,ig,ir
integer lmax, igq0, ik1, ik2, iq, ikloc1, m1, m2, m3,l,m0,lm,ia,m,isym,kmap
complex(8) :: zt,zt2
complex(8), allocatable :: zt1(:)
complex(8), allocatable :: zt3(:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zfftn(:)
complex(8), allocatable :: zfftm(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
complex(8), allocatable :: wfcr(:,:,:)
complex(8), allocatable :: wfcrsp(:,:,:)
complex(8) :: zrho0

integer, allocatable :: igkignr_ik1(:)
real(8), allocatable :: rrmt(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: vgqc1(:,:)
real(8), allocatable :: tpgqc(:,:)
real(8), allocatable :: gqc(:)
real(8), allocatable :: jlgqr(:,:,:)
real(8), allocatable :: extmp(:,:,:)
!
complex(8), external :: zfinp2,zfmtinp
real(8), external :: gaunt
real(8) :: cfq, t1, zn(nspecies),v(3),sum1
complex(8), allocatable :: wfmt(:,:,:,:,:), wfit(:,:,:)
integer :: wfmttd, wfittd, ngknr_ik1, iv(3)
character(len=256) :: fname
integer :: isp1, ibnd, fbnd, exst, ik0, ios, ist0, ist1
integer :: bndrg(2,nspinor)
integer,allocatable :: neval(:,:),evalmap(:,:,:)
!

if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxvr,ngvec))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngvec,natmtot))
!
allocate(igkignr_ik1(ngkmax))
allocate(gkcnr(ngkmax))
allocate(vgqc1(3,ngvec))
allocate(tpgqc(2,ngvec))
allocate(gqc(ngvec))
allocate(jlgqr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(wfmt(lmmaxapw,nufrmax,natmtot,nspinor,nstsv))
allocate(wfit(ngkmax,nspinor,nstsv))
allocate(zfftn(ngrtot))
allocate(zfftm(ngrtot))
allocate(zt1(nrmtmax))
allocate(zt3(nrmtmax))
allocate(zrhomt(lmmaxvr,nrmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrmtmax,natmtot))
allocate(zvclir(ngrtot))
allocate(rrmt(nrmtmax,nspecies))
allocate(wfcr(lmmaxvr,nrmtmax,2))
allocate(wfcrsp(lmmaxvr,nrmtmax,2))
allocate(extmp(nst,nspinor,nkpt))
! find out degenencies between bndrg(1,isp1) and bndrg(2,isp1)
allocate(evalmap(nbnd,nspinor,nkptnrloc))
! number of subintervals between bndrg(1,isp1) and bndrg(2,isp1)
allocate(neval(nspinor,nkptnrloc))
!
rrmt=0.d0
wfcr=zzero
wfcrsp=zzero
lmax=lmaxvr+npsden+1
wfmttd=lmmaxapw*nufrmax*natmtot*nspinor*nstsv
wfittd=ngkmax*nspinor*nstsv
ik0=0
extmp=0.d0
ex_vc=0.d0
ist0=0

call genwiq2_new

! coefficient for long-range term
cfq=0.5d0*(omega/pi)**2
! set the nuclear charges to zero
zn(:)=0.d0

do is=1,nspecies
 do ir=1,nrmt(is)
   rrmt(ir,is)=spr(ir,is)
 enddo
enddo

do isp1=1,nspinor
 if (isp1.eq.1) then
  bndrg(1,isp1)=qpnb(1)
  bndrg(2,isp1)=qpnb(2)
 elseif (isp1.eq.2) then   !spin-polarized
  bndrg(1,isp1)=qpnb(1)+int(nstsv/2)
  bndrg(2,isp1)=qpnb(2)+int(nstsv/2)
 endif
enddo

if (gw_restart) then
  ! exchange part from valence electrons
  if (exxtype.eq.0) then
    write(fname,'("Temp_files/tmp_Sx_ip",I4.4)') mpi_grid_dim_pos(dim_k)
    inquire(file=trim(adjustl(fname)),exist=exst)

    if (exst) then
     open(166,file=trim(adjustl(fname)),action='read',form='unformatted',&
             &status='old')
     read(166) ik0
     read(166) extmp
     close(166)
    else
     write(*,'("Error reading tmp_Sx",I4)') mpi_grid_dim_pos(dim_k)
     ik0=0
    endif
  endif

! read-in the ex_vc array
 if (.not.rho_val) then
   write(fname,'("Temp_files/tmp_Sx_core_ip",I4.4)') mpi_grid_dim_pos(dim_k)
   inquire(file=trim(adjustl(fname)),exist=exst)

   if (exst) then
    open(166,file=trim(adjustl(fname)),form='unformatted',status='old')
    read(166) ex_vc
    read(166) ist0
    close(166)
   else
    write(*,'("Error reading tmp_Sc_core!!",I4)') mpi_grid_dim_pos(dim_k)
    ist0=0
   endif
 endif

endif

!----------------------------
! first we calculate the valence-valence contribution to \Sigma^x_{jk}
!----------------------------

if (exxtype.eq.0) then

do ik1=1,nkpt  ! kIBZ
 ! in case when gw_restart=.true.
 if (ik1.le.ik0) cycle
 if (mpi_grid_root()) write(*,*) "exx ik1:",ik1

! find out index of kIBZ in the non-reduced k-mesh
 do ik2=1,nkptnr
  if (abs(vklnr(1,ik2)-vkl(1,ik1))+abs(vklnr(2,ik2)-vkl(2,ik1))+&
      &abs(vklnr(3,ik2)-vkl(3,ik1)).lt.epslat) then
    kmap=ik2
    exit
  endif
 enddo

 wfmt=zzero
 wfit=zzero
 ngknr_ik1=0
 igkignr_ik1=0
 do ikloc1=1,nkptnrloc
   ik2=mpi_grid_map(nkptnr,dim_k,loc=ikloc1)
   if (kmap.eq.ik2) then
     wfmt(:,:,:,:,:)=wfsvmtnrloc(:,:,:,:,:,ikloc1)
     wfit(:,:,:)=wfsvitnrloc(:,:,:,ikloc1)
     ngknr_ik1=ngknr(ikloc1)
     igkignr_ik1(:)=igkignr(:,ikloc1)
     exit
   endif
 enddo
 call mpi_grid_reduce(wfmt(1,1,1,1,1),wfmttd,dims=(/dim_k/),all=.true.)
 call mpi_grid_reduce(wfit(1,1,1),wfittd,dims=(/dim_k/),all=.true.)
 call mpi_grid_reduce(ngknr_ik1,1,dims=(/dim_k/),all=.true.)
 call mpi_grid_reduce(igkignr_ik1(1),ngkmax,dims=(/dim_k/),all=.true.)

 do ikloc=1,nkptnrloc   ! k'
   ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)

! determine q-vector
   iv(:)=ivknr(:,kmap)-ivknr(:,ik)    ! k-k'
   iv(:)=modulo(iv(:),ngridk(:))
   iq=ikmap(iv(1),iv(2),iv(3))   
   v(:)=vkcnr(:,kmap)-vkcnr(:,ik)     ! +q

   do ig=1,ngvec
! determine G+q vectors
     vgqc1(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
     call sphcrd(vgqc1(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vectors
     call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
   end do

! structure factor for G+q
   call gensfacgp(ngvec,vgqc1,ngvec,sfacgq)
! find the shortest G+q-vector
   call findigp0(ngvec,gqc,igq0)
! compute the required spherical Bessel functions
   call genjlgpr(lmax,gqc,jlgqr)

   do isp1=1,nspinor
    ibnd=bndrg(1,isp1)
    fbnd=bndrg(2,isp1)

    do jst=ibnd,fbnd
     do ist=1,nstsv
       if (occsvnr(ist,ik).ge.occ_tor) then
         zrhomt=zzero
         zrhoir=zzero
         do lm3=1,lmmaxvr
           l3=lm2l(lm3)
           m3=lm2m(lm3)
           do ias=1,natmtot
             is=ias2is(ias)
             ic=ias2ic(ias)
             zt1=zzero
             do l1=0,lmaxapw; do io1=1,nufr(l1,is)
               do l2=0,lmaxapw; do io2=1,nufr(l2,is)
                 if (mod(l1+l2+l3,2).eq.0) then
                   zt=zzero
                   do ispn=1,nspinor
                     do m2= -l2, l2
                       lm2=idxlm(l2, m2)
                       do m1=-l1, l1    
                        lm1=idxlm(l1, m1)                      
                        zt=zt+dconjg(wfsvmtnrloc(lm1,io1,ias,ispn,ist,ikloc))*&
                          &wfmt(lm2,io2,ias,ispn,jst)*gaunt(l2,l3,l1,m2,m3,m1)
                       enddo !m1
                     enddo !m2
                   enddo !ispn
                   zt1(:)=zt1(:)+zt*ufr(:,l1,io1,ic)*ufr(:,l2,io2,ic)
                 endif
               enddo; enddo !l2, io2
             enddo; enddo !l1, io1
             zrhomt(lm3,:,ias)=zt1(:)*occsvnr(ist,ik)/occmax
           enddo !ias
         enddo !lm3

         do ispn=1,nspinor
           zfftn=zzero
           zfftm=zzero
           do ig=1,ngknr(ikloc) 
             zfftn(igfft(igkignr(ig,ikloc)))=wfsvitnrloc(ig,ispn,ist,ikloc)
           enddo
           do ig=1,ngknr_ik1
             zfftm(igfft(igkignr_ik1(ig)))=wfit(ig,ispn,jst) 
           enddo
           call zfftifc(3,ngrid,1,zfftn)
           call zfftifc(3,ngrid,1,zfftm)

           do ir=1,ngrtot
             zrhoir(ir)=zrhoir(ir)+dconjg(zfftn(ir))*zfftm(ir)/omega*occsvnr(ist,ik)/occmax
           enddo
         enddo !ispn

! solve the complex Poisson's equation

         call zpotcoul(nrmt,nrmtmax,nrmtmax,rrmt,igq0,gqc,jlgqr,&
             ylmgq,sfacgq,zn,zrhomt,zrhoir,zvclmt,zvclir,zrho0)

         zt2=zfinp2(.true.,zrhomt,zvclmt,zrhoir,zvclir,nrmt,rrmt)
         t1=cfq*wiq2(iq)*(dble(zrho0)**2+aimag(zrho0)**2)
         extmp(jst-ibnd+1,isp1,ik1)=extmp(jst-ibnd+1,isp1,ik1)&
                                   & -dble(zt2/nkptnr)-t1

       endif
     enddo !ist
    enddo !jst
   enddo !isp1
 enddo !ikloc

! save the intermediate Exx terms
 if (mpi_grid_root((/dim_q/))) then
  write(fname,'("Temp_files/tmp_Sx_ip",I4.4)') mpi_grid_dim_pos(dim_k)
  open(166,file=trim(adjustl(fname)),form='unformatted',status='replace')
  write(166) ik1
  write(166) extmp
  close(166)
 endif

enddo !ik1

call mpi_grid_reduce(extmp(1,1,1),nst*nspinor*nkpt,dims=(/dim_k/),all=.true.)

! valence-valence part of \Sigma^x_{jk}
do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 call findkpt(vklnr(:,ik),isym,ik1)
 ex_nk(:,:,ikloc)=extmp(:,:,ik1)
enddo

endif !exxtype.eq.0

!----------------------------
! valence-core contribution to \Sigma^x_{jk}
!----------------------------
if (.not.rho_val.and.iter.eq.1) then
 ist1=0

 do is=1,nspecies
  do ia=1,natoms(is)
   ias=idxas(ia,is)
   ic=ias2ic(ias)
   do ist=1,spnst(is)
     if (spcore(ist,is)) then
      if (pt_core) then
       if (ist.lt.crst_in(is)) cycle
      endif
      m0=spk(ist,is)
      do m=-m0,m0-1
       ist1=ist1+1
       if (ist1.le.ist0) cycle
       if (mpi_grid_root()) write(*,'("ist1: ",I4)') ist1

! pass m-1/2 to wavefcr
       call wavefcr(1,is,ia,ist,m,nrmtmax,wfcr)

! convert from spherical coordinates to spherical harmonic
       do ispn=1,nspinor
         call zgemm('N','N',lmmaxvr,nrmtmax,lmmaxvr, zone,zfshtvr, &
              lmmaxvr,wfcr(1,1,ispn),lmmaxvr,zzero,wfcrsp(1,1,ispn),lmmaxvr)
       enddo ! ispn

       do ikloc=1,nkptnrloc
        do isp1=1,nspinor
         ibnd=bndrg(1,isp1)
         fbnd=bndrg(2,isp1)

         do jst=ibnd,fbnd
          zrhomt=zzero
          zvclmt=zzero

          do lm3=1,lmmaxvr
           l3=lm2l(lm3)
           m3=lm2m(lm3)
           zt3=zzero
           do l1=0,lmaxapw; do io1=1,nufr(l1,is)
            zt1=zzero
            do l2=0,lmaxapw
              if (mod(l1+l2+l3,2).eq.0) then
                do ispn=1,nspinor
                 do m2= -l2, l2
                  lm2=idxlm(l2,m2)
                  do m1=-l1, l1
                   lm1=idxlm(l1,m1)
                   zt1(:)=zt1(:)+dconjg(wfsvmtnrloc(lm1,io1,ias,ispn,&
                         jst,ikloc))*wfcrsp(lm2,:,ispn)* &
                         gaunt(l2,l3,l1,m2,m3,m1)
                  enddo !m1
                 enddo !m2
                enddo !ispn
              endif
            enddo !l2  
            zt3(:)=zt3(:)+zt1(:)*ufr(:,l1,io1,ic)            
           enddo; enddo !l1, io1
           zrhomt(lm3,:,ias)=zt3(:)
          enddo !lm3

! calculate the Coulomb potential in MT

          call zpotclmt(.false.,lmaxvr,nrmt(is),rrmt(:,is),0.d0,lmmaxvr, &
                  zrhomt(:,:,ias),zvclmt(:,:,ias))
          zt2=zfmtinp(.true.,lmaxvr,nrmt(is),rrmt(:,is),lmmaxvr, &
                  zrhomt(:,:,ias),zvclmt(:,:,ias))
          ex_vc(jst-ibnd+1,isp1,ikloc)=ex_vc(jst-ibnd+1,isp1,ikloc)-dble(zt2)
         enddo !jst
        enddo !isp1
       enddo !ikloc  

       ! store the ex_vc array
       if (mpi_grid_root((/dim_q/))) then
         write(fname,'("Temp_files/tmp_Sx_core_ip",I4.4)') mpi_grid_dim_pos(dim_k)
         open(166,file=trim(adjustl(fname)),form='unformatted',status='replace')
         write(166) ex_vc
         write(166) ist1
         close(166)
       endif
    
      enddo !m
     endif
   enddo !ist
  enddo !ia
 enddo !is

endif

! wait for all processors
call mpi_grid_barrier()

! find out degenerate states
call find_degenency(bndrg,evalmap,neval)

! symmetrization of ex_nk
do ikloc=1,nkptnrloc
 do isp1=1,nspinor
  do m=1,neval(isp1,ikloc)-1
   ibnd=evalmap(m,isp1,ikloc)
   fbnd=evalmap(m+1,isp1,ikloc)-1
   if (fbnd.eq.ibnd) cycle

   sum1=0.d0
   do ist=ibnd,fbnd
    sum1=sum1+ex_nk(ist-bndrg(1,isp1)+1,isp1,ikloc)
   enddo
   sum1=sum1/(fbnd-ibnd+1)
   do ist=ibnd,fbnd
    ex_nk(ist-bndrg(1,isp1)+1,isp1,ikloc)=sum1
   enddo

   if (.not.rho_val.and.iter.eq.1) then
    sum1=0.d0
    do ist=ibnd,fbnd
     sum1=sum1+ex_vc(ist-bndrg(1,isp1)+1,isp1,ikloc)
    enddo
    sum1=sum1/(fbnd-ibnd+1)
    do ist=ibnd,fbnd
     ex_vc(ist-bndrg(1,isp1)+1,isp1,ikloc)=sum1
    enddo
   endif

  enddo !m
 enddo !isp1
enddo !ikloc

deallocate(zfftn,zfftm,zrhomt,zrhoir,zt1,zt3,zvclmt,zvclir)
deallocate(igkignr_ik1,gkcnr)
deallocate(vgqc1,tpgqc,gqc,jlgqr,ylmgq,wfmt,wfit)
deallocate(wfcrsp,wfcr)
deallocate(extmp,sfacgq,rrmt)
deallocate(evalmap,neval)
return
end subroutine
