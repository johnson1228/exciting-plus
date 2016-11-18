subroutine gwmain
use modmain
use mod_addons_q
use mod_nrkp
use mod_hdf5
use mod_wannier
use mod_linresp
use mod_expigqr
implicit none
integer iq,i,tmp,ndk,ib,exind
integer nvqrloc,nvqr,iqloc,it,ikloc,ik,ist,is
integer nw_se
character*100 qnm,fname,fspn,fgf,fname_tot
integer iw
logical exst
integer,allocatable :: bndrg(:,:)
real(8), allocatable :: vxcnk(:,:,:)
real(8), allocatable :: vclnk(:,:,:)
real(8), allocatable :: exxnk(:,:,:)
real(8), allocatable :: exxvc(:,:,:)
real(8), allocatable :: eqp_aux(:,:,:)
!gw band
real(8), allocatable :: Zf(:,:,:),Eqp(:,:,:),prev_Eqp(:,:,:)
complex(8) :: dSe, factor
real(8) :: emin, emax
! test on plotting green function including QP correction
logical :: gfdos
real(8) :: ww, dw_gfdos, step, maxdel
complex(8) :: gf, delta, spec_f ! spectral function of dressed G
integer :: nw_gfdos,istep,step0
integer :: isp1, ibnd, fbnd, iz, q0, ios, ierr1,ierr2
integer,allocatable :: iqrmap(:,:)
integer,allocatable :: qqnrmap(:,:,:),rkmap(:,:),evalmap(:,:,:),neval(:,:)
integer,allocatable :: kknrmap(:,:,:),kmap(:,:)
integer :: nqkmax,irk,ik1,isym,iqbz
complex(8), allocatable :: corr_se_aux(:,:,:,:,:),sum1(:),se_aux2(:,:,:,:)
!
call init0
call init1
if (.not.mpi_grid_in()) return
if (mpi_grid_root()) call timestamp(6,"[gwmain] done init")
!call init_q_mesh(8)
! if vq0c not used, set nvq0=8
if (nvq0.ne.1) nvq0=8
call init_q_mesh(nvq0)   ! test
call genvq
call genvgq
! read the density and potentials from file
call readstate
! generate radial functions
call genradf
! generate core states
if (.not. rho_val) call gencore

if (mpi_grid_root()) then
  open(151,file="GW.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",8I6)')&
    &(mpi_grid_dim_size(i),i=1,mpi_grid_nd)
  write(151,*)
  do i=1,nvq0
    write(151,'(" vqc : ",3G18.10)')vqc(:,i)
  enddo
  call flushifc(151)

  !create folders for temp files
  call system("mkdir -p Temp_files")

endif
! generate wave-functions for entire BZ
call genwfnr(151,tq0bz)

! allocate array for band range
allocate(bndrg(2,nspinor))
! find out mapping from q_{IBZ} to q_{BZ}, and k1 to k2=R^{-1}*k1
nvqr=nkpt+nvq0-1
allocate(iqrmap(2,nvqr))
allocate(qqnrmap(nvq,2,nvqr))
allocate(rkmap(48,nkptnr))
allocate(kmap(2,nkpt))
allocate(kknrmap(nkptnr,2,nkpt))
call kq_map(iqrmap,qqnrmap,rkmap,kmap,kknrmap)

if (mpi_grid_root()) write(*,*) "after kq_map"

! set Efermi=0
evalsvnr(:,:)=evalsvnr(:,:)-efermi
! find out the energy range
emin=minval(evalsvnr)
emax=maxval(evalsvnr)

!checking qpnb(1) and qpnb(2)
if ((qpnb(1).gt.qpnb(2)).or.(qpnb(1).le.0)) qpnb(1)=1
if ((qpnb(2).le.0).or.(qpnb(2).gt.nstsv)) qpnb(2)=nstsv
nbnd=qpnb(2)-qpnb(1)+1
!band range
do isp1=1,nspinor
 if (isp1.eq.1) then
  bndrg(1,isp1)=qpnb(1)
  bndrg(2,isp1)=qpnb(2)
 elseif (isp1.eq.2) then   !spin-polarized
  bndrg(1,isp1)=qpnb(1)+int(nstsv/2)
  bndrg(2,isp1)=qpnb(2)+int(nstsv/2)
 endif
enddo
! find out degenencies between bndrg(1,isp1) and bndrg(2,isp1)
allocate(evalmap(nbnd,nspinor,nkptnrloc))
! number of subintervals between bndrg(1,isp1) and bndrg(2,isp1)
allocate(neval(nspinor,nkptnrloc))
call find_degenency(bndrg,evalmap,neval)
! energy mesh for the GW calculation
call gen_gw_wmesh(0,emax,nw_se)
! allocate arrays for self-energies
allocate(corr_se_aux(nw_se,nbnd,nspinor,nkptnr,nvqr))
allocate(gw_self_energy(nw_se,nbnd,nspinor,nkptnrloc))
allocate(se_aux2(nw_se,nbnd,nspinor,nkpt))
allocate(vxcnk(nbnd,nspinor,nkptnrloc))
allocate(vclnk(nbnd,nspinor,nkptnrloc))
allocate(exxnk(nbnd,nspinor,nkptnrloc))
allocate(exxvc(nbnd,nspinor,nkptnrloc))
allocate(eqp_aux(nbnd,nspinor,nkptnrloc))
allocate(sum1(nw_se))
! array for quasiparticle energies
allocate(Eqp(nbnd,nspinor,nkptnr))
allocate(prev_Eqp(nbnd,nspinor,nkptnr))
allocate(Zf(nbnd,nspinor,nkptnrloc))
!
!initial values
exxnk=0.d0  
exxvc=0.d0
vxcnk=0.d0
vclnk=0.d0
gw_self_energy=zzero
corr_se_aux=zzero
se_aux2=zzero
prev_Eqp=0.d0
Zf=0.d0
q0=0
step0=0
ierr1=1
ierr2=1
istep=1
exind=0

! gw band plot
if (gwband) then
  call gwband_dat(151)
  if (mpi_grid_root()) call flushifc(151)
! go to the last line
  goto 102 
endif ! gwband

!reset nebd_chi and nebd_se
if (nebd_chi.eq.1) nebd_chi=nempty
if (nebd_se.eq.1) nebd_se=nempty

! check exxtype
if (exxtype.ne.0) then
 write(*,'("exxtype needs to be zero!")')
 call pstop
endif

! printing out variables
if (mpi_grid_root()) then
 do isp1=1,nspinor
   write(151,'("spin component:",I1)') isp1
   write(151,'("band range:",2(I4))') bndrg(:,isp1)
 enddo
 write(151,'("emin,emax:",2f16.8)') emin,emax
 write(151,'("lr_eta:",f16.8)') lr_eta
 write(151,'("lr_nw:",I5)') lr_nw
 write(151,'("niw:",I5)') niw
 write(151,'("gf_niw:",I6)') gf_niw
 write(151,'("nw_se:",I5)') nw_se
 write(151,'("nstsv:",I5)') nstsv
 write(151,'("nempty:",I5)') nempty
 write(151,'("nebd_chi:",I5)') nebd_chi
 write(151,'("nebd_se:",I5)') nebd_se
 write(151,'("ngq(1):",I5)') ngq(1)
 write(151,*) "rho_val:",rho_val
 write(151,*) "pt_core:",pt_core
 write(151,'("raicut:",2(f16.8,1X))') raicut(1),raicut(2)
 write(151,'("scgwni:",I5)') scgwni
 call timestamp(151)
 call flushifc(151)
endif

! generate matrix elements <nk|Vxc|nk>
if (gw_restart) then
! read data from old files
 call read_exxvxc(nbnd,exxnk,exxvc,vxcnk,vclnk,exind,ierr1)
 if (ierr1.ne.0) goto 202
else
 202 continue
 call genvxcnk(nbnd,vxcnk)
 call genexxnk(istep,nbnd,exxnk,exxvc)
! outputing exxnk and vxcnk
 if (mpi_grid_root((/dim_q/))) &
    &call write_exxvxc(1,nbnd,bndrg,exxnk,exxvc,vxcnk,vclnk)
endif !gw_restart

! get quantity E1-Vxc
do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isp1=1,nspinor
  ibnd=bndrg(1,isp1)
  fbnd=bndrg(2,isp1)
  do ist=ibnd,fbnd
   i=ist-ibnd+1
   eqp_aux(i,isp1,ikloc)=evalsvnr(ist,ik)-vxcnk(i,isp1,ikloc)
  enddo
 enddo
enddo

megq_include_bands=chi0_include_bands
! check for tmp_Sc files if ierr1=0
if ((gw_restart).and.(ierr1.eq.0)) then
! read in the stored correlation self-energy elements
  call read_gwse(nw_se,nbnd,corr_se_aux,prev_Eqp,ierr2,step0,q0)
endif !gw_restart and ierr1=0

! printing variables
if (mpi_grid_root()) then
 write(151,'("q0: ",I4)') q0
 write(151,'("step0: ",I4)') step0
 write(151,'("ierr1,ierr2:",2(I4,1X))') ierr1,ierr2
 write(151,'(" ")')
 write(151,'("start to calculate Sc:")')
 call timestamp(151)
 call flushifc(151)
endif

! for GW0 calculation, evalsvnr needs to be updated when restarting!
do ik=1,nkptnr
 do isp1=1,nspinor
  ibnd=bndrg(1,isp1)
  fbnd=bndrg(2,isp1)
  do i=ibnd,fbnd
   if ((step0.eq.0).or.(step0.eq.1)) then
   ! E_KS
    prev_Eqp(i-ibnd+1,isp1,ik)=evalsvnr(i,ik)  
   elseif (step0.ge.2) then
   ! prev_Eqp is read from the file
    evalsvnr(i,ik)=prev_Eqp(i-ibnd+1,isp1,ik)
   endif
  enddo !i
 enddo !isp1
enddo !ik

!total exchange self-energy
exxnk(:,:,:)=exxnk(:,:,:)+exxvc(:,:,:)

! outmost loop for GW0 calculations, exx and vxc don't need to be updated
! when the system has a nonzero band gap
do istep=1,scgwni 
 if (istep.le.(step0-1)) cycle
 ! need to reset zero every time
 Eqp=0.d0
 ! nonzero only when restarting properly at istep=step0!
 if ((ierr2.ne.0).or.(istep.gt.step0)) then
  corr_se_aux=zzero
  se_aux2=zzero
  gw_self_energy=zzero
 endif

 if (mpi_grid_root()) then
  write(151,'("G ",I2,"W 0 calculation starts:")') istep-1
  write(151,'(" ")')
  call flushifc(151)
 endif
! main loop over q-points
 do iqloc=1,nvqr
  iq=iqrmap(1,iqloc)
  if ((iqloc.le.q0).and.(istep.le.step0)) cycle
 
  if (istep.eq.1) then
   call genmegq(iq,.false.,.true.,.true.)
   call get_adjoint_megqblh(iq)
   if (scgwni.gt.1.and.mpi_grid_root((/dim_q/))) call write_amegqblh(iq)
  elseif (istep.gt.1) then
! need to initialize idxkq(3,jk)
   call init_kq(iq)
! initialize nmegqblh
   call init_band_trans(.true.)
! read amegqblh
   call read_amegqblh(iq)
  endif
 
! type of GW calculations
  if (caltype.eq.0) then !ppa
   call ppa_self_energy(iq,iqloc,istep,nw_se,corr_se_aux)
  elseif (caltype.eq.1) then !rai
   call rai_self_energy(iq,iqloc,istep,nw_se,corr_se_aux)
  elseif (caltype.eq.2) then !full se, need check!
   call full_self_energy(iq,iqloc,istep,nw_se,corr_se_aux)
  elseif (caltype.eq.3) then ! new version of rai
   call rai_self_energy_new(iq,iqloc,istep,nw_se,corr_se_aux)
  endif

! creating tmp files for Sc    
  write(fname,'("Temp_files/tmp_Sc_ik",I4.4,"_iq",I4.4)') &
       & mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
  open(166,file=trim(adjustl(fname)),form='unformatted',status='replace')
  ! create the other tmp files, just in case!
  if (mod(iqloc,2).eq.0) then !even
   write(fname,'("Temp_files/e_tmp_Sc_ik",I4.4,"_iq",I4.4)') &
       & mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
  else  !odd
   write(fname,'("Temp_files/o_tmp_Sc_ik",I4.4,"_iq",I4.4)') &
       & mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
  endif

  open(1660,file=trim(adjustl(fname)),form='unformatted',status='replace')
  write(166) istep,iqloc
  write(166) corr_se_aux  ! in Ha
  write(166) prev_Eqp
  write(1660) istep,iqloc
  write(1660) corr_se_aux  ! in Ha
  write(1660) prev_Eqp

  close(166)
  close(1660)

  if (mpi_grid_root()) then
   write(151,'("iq : ",I4," out of ",I4)') iqloc,nvqr
   call timestamp(151)
   call flushifc(151)
  endif

 enddo !iqloc

 if (caltype.eq.3) then
  call mpi_grid_reduce(corr_se_aux(1,1,1,1,1),nw_se*nbnd*nspinor*nkptnr*nvqr,&
                     dims=(/dim_k,dim_q/),all=.true.)
 else
  call mpi_grid_reduce(corr_se_aux(1,1,1,1,1),nw_se*nbnd*nspinor*nkptnr*nvqr,&
                     dims=(/dim_k/),all=.true.)
 endif

! impose symmetry between k_{BZ} <-> k_{IBZ}
 do ik1=1,nkpt
  ik=kmap(1,ik1)
  do iq=1,nvqr
   do iqbz=1,iqrmap(2,iq)
    irk=rkmap(qqnrmap(iqbz,2,iq),ik)
    se_aux2(:,:,:,ik1)=se_aux2(:,:,:,ik1)+&
                         &corr_se_aux(:,:,:,irk,iq)
   enddo
  enddo

  do ikloc=1,nkptnrloc
   ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
   do i=1,kmap(2,ik1)
    if (ik.eq.kknrmap(i,1,ik1)) gw_self_energy(:,:,:,ikloc)=se_aux2(:,:,:,ik1)
   enddo
  enddo
 enddo !ik1

! symmetrization of self-energy matrix elements for degenerate states
 if (nbnd.gt.1) then
  do ikloc=1,nkptnrloc
   do isp1=1,nspinor
    do i=1,neval(isp1,ikloc)-1
     ibnd=evalmap(i,isp1,ikloc)
     fbnd=evalmap(i+1,isp1,ikloc)-1
     if (fbnd.eq.ibnd) cycle
    
     sum1(:)=zzero
     do ist=ibnd,fbnd
      sum1(:)=sum1(:)+gw_self_energy(:,ist-bndrg(1,isp1)+1,isp1,ikloc)
     enddo
     sum1(:)=sum1(:)/(fbnd-ibnd+1)
     do ist=ibnd,fbnd
      gw_self_energy(:,ist-bndrg(1,isp1)+1,isp1,ikloc)=sum1(:)
     enddo

    enddo !i
   enddo !isp1
  enddo !ikloc
 endif

 if (caltype.gt.0) then !rai or full_se
  factor=zi/twopi/occmax/nkptnr/omega
 else ! ppa
  factor=zone/nkptnr/omega
 endif
 gw_self_energy(:,:,:,:)=gw_self_energy(:,:,:,:)*factor

 ! total self-energy
 do iw=1,nw_se
  gw_self_energy(iw,:,:,:)=gw_self_energy(iw,:,:,:)+exxnk(:,:,:)
 enddo

 ! update the quasiparticle energies
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   ibnd=bndrg(1,isp1)
   fbnd=bndrg(2,isp1)
   do ist=ibnd,fbnd
    i=ist-ibnd+1
    dSe=gw_self_energy(2,i,isp1,ikloc)-gw_self_energy(1,i,isp1,ikloc)
    dSe=dSe/del_e
    Zf(i,isp1,ikloc)=1.d0/(1.d0-dreal(dSe))
    Eqp(i,isp1,ik)=prev_Eqp(i,isp1,ik)+Zf(i,isp1,ikloc)*&
                 &dreal(gw_self_energy(1,i,isp1,ikloc)+&
                 &eqp_aux(i,isp1,ikloc)-prev_Eqp(i,isp1,ik))
   enddo !i
  enddo !isp1
 enddo !ikloc
 
 call mpi_grid_reduce(Eqp(1,1,1),nbnd*nspinor*nkptnr,dims=(/dim_k/), &
                    & all=.true.)

 if (scgwni.gt.1) then
  call get_maxdel(nbnd,Eqp,prev_Eqp,maxdel)
  if (mpi_grid_root()) then
   write(151,'("maximum delta is, ",f8.4," eV.")') maxdel*ha2ev
   call flushifc(151)
  endif
! check for convergence
  if (maxdel*ha2ev.le.0.01d0) then
   if (mpi_grid_root()) then
    write(151,'("maximum delta is less than 0.01 eV, ",f6.3)') maxdel*ha2ev
    write(151,'("GW0 calculation is converged!")') 
    write(151,'(" ")')
    call flushifc(151)
   endif
   exit ! exit the GW0 loop
  else
! store the quasiparticle energies, do not update in the last step
   if (istep.le.scgwni-1) then
    do ik=1,nkptnr
     do isp1=1,nspinor
      ibnd=bndrg(1,isp1)
      fbnd=bndrg(2,isp1)
      do i=ibnd,fbnd
       evalsvnr(i,ik)=Eqp(i-ibnd+1,isp1,ik)
       if (mpi_grid_root()) then
         write(151,*) i,isp1,ik,Eqp(i-ibnd+1,isp1,ik)*ha2ev
         call flushifc(151)
       endif
      enddo !i
     enddo !isp1
    enddo !ik
    prev_Eqp(:,:,:)=Eqp(:,:,:)
   endif
  endif
 endif !scgwni.gt.1

enddo !istep

! G0W0 part, modified by I.H. Chu, H.P. Cheng, Jun,2012
 if ((caltype.eq.0).or.(caltype.eq.1).or.(caltype.eq.3)) then !ppa or rai
  if (mpi_grid_root((/dim_q/))) then
   do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do isp1=1,nspinor
     ibnd=bndrg(1,isp1)
     fbnd=bndrg(2,isp1)

     if (nspinor.eq.1) then
      fspn=''
     else
      write(fspn,'("_spn",I1.1)') isp1
     endif

     write(fname,'("Eqp_k",I4.4)')ik
     fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
     open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
     do i=ibnd,fbnd
      write(168,'(I5,1X,6(G16.6,1X))') i,Eqp(i-ibnd+1,isp1,ik)*ha2ev, &
             &                prev_Eqp(i-ibnd+1,isp1,ik)*ha2ev,&
             &                exxnk(i-ibnd+1,isp1,ikloc)*ha2ev,&
             &                vxcnk(i-ibnd+1,isp1,ikloc)*ha2ev, &
             &                dreal(gw_self_energy(1,i-ibnd+1,isp1,ikloc) &
             &                -exxnk(i-ibnd+1,isp1,ikloc))*ha2ev, &
             &                Zf(i-ibnd+1,isp1,ikloc)
     enddo !i
     write(168,'("Done!")')
     close(168)
    enddo !isp1
   enddo !ikloc
  endif
!!!!!!!!!!!!!!!!!!!!!
 elseif (caltype.eq.2) then ! full self-energy
  if (mpi_grid_root((/dim_q/))) then
   do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do isp1=1,nspinor
      ibnd=bndrg(1,isp1)
      fbnd=bndrg(2,isp1)

      if (nspinor.eq.1) then
       fspn=''
      else
       write(fspn,'("_spn",I1.1)') isp1
      endif   

      do i=ibnd,fbnd
       write(fname,'("self_energy_k",I4.4,"_b",I4.4)')ik,i
       fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
       open(170,file=fname_tot,form="FORMATTED",status="REPLACE")
       do iw=1,nw_se
         write(170,'(3G16.6)')dble(iw),dreal(gw_self_energy(iw,i-ibnd+1,isp1,ikloc)),dimag(gw_self_energy(iw,i-ibnd+1,isp1,ikloc))
       enddo
       write(170,*) "Done!"
       close(170)
      enddo  !i
    enddo !isp1
   enddo !ikloc
  endif
 endif !caltype=2

! wait for all the processors
call mpi_grid_barrier()

if (mpi_grid_root()) then
 call timestamp(151)
 write(151,*) 
 write(151,'("Done.")')
 close(151)
endif
deallocate(lr_w)
deallocate(gw_self_energy)
deallocate(vxcnk,vclnk,exxnk,exxvc,Eqp,Zf)
deallocate(prev_Eqp,bndrg,eqp_aux)
deallocate(iqrmap,qqnrmap,rkmap,evalmap,neval,sum1)

! delete all the tmp files using the root processor
if (mpi_grid_root()) then
 call system('rm -f Temp_files/*tmp_* amegq_*')
 if (caltype.eq.0) call system('rm -f Temp_files/ppa_mat*')
 if (caltype.eq.1) call system('rm -f Temp_files/rai_amat_*')
 if (caltype.eq.3) call system('rm -f Temp_files/rai_eps*')
endif

102 continue
return
end subroutine
