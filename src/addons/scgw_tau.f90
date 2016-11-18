subroutine scgw_tau
use modmain
use mod_addons_q
use mod_nrkp
use mod_hdf5
use mod_wannier
use mod_linresp
use mod_expigqr
!
implicit none
integer :: ibnd,fbnd,iter,nvqr,iqloc,iq,exind,ierr1,ierr2,n
integer :: ist1,ist,ik1,iqbz,irk
integer :: i,ik,ikloc,iw,isp1
integer :: nw_se,nwloc
integer :: iter0, q0, iter1
real(8) :: diff,maxdel,maxdel_prev
real(8) :: Zf,emin,emax
real(8) :: factor,c,tmem
complex(8) :: dSe
character*100 fname,fspn,fname_tot
!
integer,allocatable :: bndrg(:,:)
! symmetry mappings for k- and q-points
integer,allocatable :: iqrmap(:,:)
integer,allocatable :: qqnrmap(:,:,:),rkmap(:,:),evalmap(:,:,:),neval(:,:)
integer,allocatable :: kknrmap(:,:,:),kmap(:,:)
! exchange self-energy matrices
real(8),allocatable :: exxvc(:,:,:),exxnk(:,:,:)
! K-S Exc potential, Hartree potential difference
real(8),allocatable :: vxcnk(:,:,:),vclnk(:,:,:)
! quasiparticle energies
real(8),allocatable :: Eqp(:,:,:)
! dressed Green's function
real(8),allocatable :: gf_tau(:,:,:)
! K-S Green's function
real(8),allocatable :: gf0_tau(:,:,:)
! correlation self-energy
real(8),allocatable :: sig_c(:,:,:,:)
! auxiliary sig_x and sig_c matrices
real(8),allocatable :: sig_c_aux(:,:,:,:,:)
real(8),allocatable :: sig_x_aux(:,:,:,:)
real(8),allocatable :: se_aux(:,:,:,:)
real(8),allocatable :: sx_aux(:,:,:)
real(8),allocatable :: sum1(:)
! correlation self-energy in frequency domain
complex(8),allocatable :: sig_cw(:,:,:,:)
! correlation self-energy in real-frequency domain
complex(8),allocatable :: sig_cr(:,:,:,:)
complex(8),allocatable :: sig_cc(:,:,:,:)
!
call init0
call init1
if (.not.mpi_grid_in()) return
if (mpi_grid_root()) call timestamp(6,"[gwmain] done init")
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
 
  !create folders for output files
  call system("mkdir -p Green_func")
  call system("mkdir -p Sig_c_files")
  call system("mkdir -p Temp_files")
endif
!
! generate wave-functions for entire BZ
call genwfnr(151,tq0bz)
!
! find out mapping from q_{IBZ} to q_{BZ}, and k1 to k2=R^{-1}*k1
nvqr=nkpt+nvq0-1
allocate(iqrmap(2,nvqr))
allocate(qqnrmap(nvq,2,nvqr))
allocate(rkmap(48,nkptnr))
allocate(kmap(2,nkpt))
allocate(kknrmap(nkptnr,2,nkpt))
call kq_map(iqrmap,qqnrmap,rkmap,kmap,kknrmap)
!
! set Efermi=0
evalsvnr(:,:)=evalsvnr(:,:)-efermi
!
! find out the energy range
allocate(bndrg(2,nspinor))
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
!
! energy mesh for the GW calculation
call gen_gw_wmesh(1,emax,nw_se)
! distribute frequency along dim_q
nwloc=mpi_grid_map(lr_nw,dim_q)
!
! allocate arrays for self-energies
allocate(sig_c(nbnd,nspinor,lr_nw,nkptnrloc))
allocate(sig_cc(nbnd,nspinor,lr_nw,nkptnrloc))
allocate(sig_c_aux(nbnd,nspinor,lr_nw,nkptnr,nvqr))
allocate(sig_x_aux(nbnd,nspinor,nkptnr,nvqr))
allocate(se_aux(nbnd,nspinor,lr_nw,nkpt))
allocate(sx_aux(nbnd,nspinor,nkpt))
allocate(sig_cw(nbnd,nspinor,niw,nkptnrloc))
allocate(sig_cr(nw_se,nbnd,nspinor,nkptnrloc))
allocate(exxnk(nbnd,nspinor,nkptnrloc))
allocate(exxvc(nbnd,nspinor,nkptnrloc))
allocate(vxcnk(nbnd,nspinor,nkptnrloc))
allocate(vclnk(nbnd,nspinor,nkptnrloc))
allocate(gf_tau(nstsv,lr_nw,nkptnrloc))
allocate(gf0_tau(nstsv,lr_nw,nkptnrloc))
allocate(Eqp(nbnd,nspinor,nkptnrloc))
allocate(sum1(lr_nw))
! find out degenencies between bndrg(1,isp1) and bndrg(2,isp1)
allocate(evalmap(nbnd,nspinor,nkptnrloc))
! number of subintervals between bndrg(1,isp1) and bndrg(2,isp1)
allocate(neval(nspinor,nkptnrloc))
!
! initial values
gf_tau=0.d0
gf0_tau=0.d0
sig_c=0.d0
sig_c_aux=0.d0
sig_x_aux=0.d0
se_aux=0.d0
sx_aux=0.d0
exxnk=0.d0
exxvc=0.d0
vxcnk=0.d0
vclnk=0.d0
Eqp=0.d0
sig_cc=zzero
sig_cw=zzero
q0=0
iter0=0
iter1=0
ierr1=1
ierr2=1
exind=0
evalmap=0
neval=0
maxdel=0.d0
maxdel_prev=0.d0
! find out degenencies between qpnb(1) and qpnb(2)
call find_degenency(bndrg,evalmap,neval)

! printing out variables
if (mpi_grid_root()) then
 do isp1=1,nspinor
   write(151,'("spin component:",I1)') isp1
   write(151,'("band range:",2(I4))') bndrg(:,isp1)
 enddo
 write(151,'("emin,emax:",2f16.8)') emin,emax
 write(151,'("gw_mode:",I3)') gw_mode
 write(151,'("lr_eta:",f16.8)') lr_eta
 write(151,'("lr_nw:",I5)') lr_nw
 write(151,'("niw:",I5)') niw
 write(151,'("nw_se:",I5)') nw_se
 write(151,'("nstsv:",I5)') nstsv
 write(151,'("nempty:",I5)') nempty
 write(151,'("nebd_chi:",I5)') nebd_chi
 write(151,'("nebd_se:",I5)') nebd_se
 write(151,'("ngq(1):",I5)') ngq(1)
 write(151,'("rho_val:",I1)') rho_val
 write(151,'("pt_core:",I1)') pt_core
 write(151,'("exxtype:",I1)') exxtype
 write(151,'("ac_sigma:",I5)') ac_sigma
 write(151,'("cpe_alp0:",G18.10)') cpe_alp0
 write(151,'("cpe_sig0:",G18.10)') cpe_sig0
 write(151,'("cpe_gtor:",G18.10)') cpe_gtor
 write(151,'("cpe_N:",I6)') cpe_N
 write(151,'("scgwni:",I5)') scgwni
 write(151,'("upm:",2I5)') upm(:)
 write(151,'("occ_tor:",f8.4)') occ_tor
 write(151,'(" ")')
 ! memory estimate
 ! two green functions, real matrices
 write(151,'("Memory allocations for the following matrices:")')
 tmem=lr_nw*nstsv*nkptnrloc*16.d0
 tmem=tmem/1024.d0/1024.d0
 write(151,'("Green functions:",f10.3, " MB")') tmem
 ! sig_c_aux
 tmem=nbnd*nspinor*lr_nw*nkptnr*nvqr*8.d0
 tmem=tmem/1024.d0/1024.d0
 write(151,'("Auxiliary sig_c:",f10.3, " MB")') tmem
 ! svq and P
 tmem=ngq(1)*ngq(1)*(lr_nw+nwloc)*16.d0
 tmem=tmem/1024.d0/1024.d0
 write(151,'("W and P:",f10.3, " MB")') tmem
 ! megqblh and amegqblh
 tmem=ngq(1)*nkptnrloc*nstsv*nstsv*32.d0
 tmem=tmem/1024.d0/1024.d0
 write(151,'("megq and amegq:",f10.3, " MB")') tmem
 ! refresh the output file
 call timestamp(151)
 call flushifc(151)
endif

! scGW
do iter=1,scgwni

 if (mpi_grid_root()) then
   write(151,'(" ")')
   write(151,'("Iteration ",I3,":")') iter
 endif

 if (iter.eq.1) then
   !initial Green's function
   call init_gf_tau(bndrg,gf0_tau)
   if (gw_restart) then
    ! read exxnk, vxcnk and vclnk from old files
    call read_exxvxc(nbnd,exxnk,exxvc,vxcnk,vclnk,exind,ierr1)
    ! read sig_c,gf_tau from file
    if (ierr1.eq.0) &
    call read_gwsetau(151,nbnd,sig_c_aux,sig_x_aux,gf_tau,sig_c,exxnk,&
                       &ierr2,iter0,q0,iter1,maxdel_prev)
   endif

   if (mpi_grid_root()) then
    write(151,'("ierr1,ierr2:",2I4)') ierr1,ierr2
    write(151,'("exind:",I4)') exind
    write(151,'("iter0,q0,iter1:",3I4)') iter0,q0,iter1
    write(151,'("bhbar:",G12.6)') bhbar
    write(151,'("dyson_mu:",G12.6)') dyson_mu
    write(151,'("maxdel_prev:",G12.6)') maxdel_prev
    call timestamp(151)
    call flushifc(151)
   endif
 endif !iter.eq.1

 if (iter.le.iter0-1) cycle

 ! if gf_tau and sig_c are already updated in tmp_Sc
 if (iter.eq.iter1-1) then
  sig_c_aux=0.d0
  se_aux=0.d0
  if (exxtype.eq.1) then
    sig_x_aux=0.d0
    sx_aux=0.d0
  endif
  q0=0
  cycle
 endif

 ! update occupation number of each state
 if (iter.eq.1) then
  call update_occ(gf0_tau)
 else 
  call update_occ(gf_tau)
 endif

 if (mpi_grid_root()) then
  write(151,'(" ")')
  write(151,'("Self-Consistent GW calculation starts: ",I2)') iter
  write(151,'(" ")')
  call flushifc(151)
 endif

 ! if no read-in data, compute exxnk, vxcnk and vclnk
 if ((ierr1.ne.0).and.iter.eq.1) then
   call genvxcnk(nbnd,vxcnk)   ! vxc
 endif

 if ((ierr1.ne.0).or.(ierr2.eq.0.and.iter.gt.iter0).or.&
    &(ierr2.ne.0.and.iter.gt.1)) then
   if (mpi_grid_root()) then
    write(151,'("Calculating exchange self-energy:",I4)') iter
    write(151,'(" ")')
    call flushifc(151)
   endif

   if ((iter.gt.exind.and.exxtype.eq.0).or.exxtype.eq.1) then
     !exx 
     ! exxtype=0: exxnk is the total exchange self-energy
     ! exxtype=1: exxnk is the core part of exchange self-energy
     if (exxtype.eq.0.or.iter.eq.1) call genexxnk(iter,nbnd,exxnk,exxvc)

     !vcl only when iter > 1
     if (iter.gt.1) then
       call genvclnk(nbnd,bndrg,vclnk)
     endif

     ! output exxnk and vxcnk
     if (mpi_grid_root((/dim_q/))) &
       &call write_exxvxc(iter,nbnd,bndrg,exxnk,exxvc,vxcnk,vclnk)
     if (mpi_grid_root()) then
       write(151,'("Outputing exchange self-energy:",I4)') iter
       call timestamp(151)
       call flushifc(151)
     endif
   endif

 endif

 ! calculate the correlation self-energy
 do iqloc=1,nvqr
  if ((iqloc.le.q0).and.(iter.le.iter0)) cycle
  iq=iqrmap(1,iqloc)

  if (mpi_grid_root()) then
   write(151,'(" ")')
   write(151,'("iq : ",I4," out of ",I4)') iqloc,nvqr
   call timestamp(151)
   call flushifc(151)
  endif
  
  if (iter.eq.1) then
   call genmegq(iq,.false.,.true.,.true.)
   if (scgwni.gt.1.and.mpi_grid_root((/dim_q/))) call write_megqblh(iq)
  elseif (iter.gt.1) then
   ! need to initialize idxkq(3,jk)
   call init_kq(iq)
   ! initialize nmegqblh
   call init_band_trans(.true.)
   call read_megqblh(iq)
  endif
 
  call get_adjoint_megqblh(iq)

  if (iter.gt.1) then
    ! GW
    call cal_sigma(iq,iqloc,nwloc,iter,bndrg,gf_tau,sig_c_aux,sig_x_aux)
  elseif (iter.eq.1) then
    ! G0W0
    call cal_sigma(iq,iqloc,nwloc,iter,bndrg,gf0_tau,sig_c_aux,sig_x_aux)
  endif

! create tmp files for sig_c_aux  
  write(fname,'("Temp_files/tmp_Sc_ik",I4.4,"_iq",I4.4)') &
         &mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
  open(166,file=trim(adjustl(fname)),form='unformatted',status='replace')
  write(166) iter,iqloc,iter
  write(166) sig_c_aux  ! in Ha
  write(166) gf_tau
  write(166) sig_c  ! in Ha
  write(166) dyson_mu,maxdel_prev
  if (exxtype.eq.1) write(166) sig_x_aux,exxnk
  close(166)

 enddo !iqloc
 
 call mpi_grid_reduce(sig_c_aux(1,1,1,1,1),nbnd*nspinor*lr_nw*nkptnr*nvqr,&
                   &dims=(/dim_k,dim_q/),all=.true.)
 
 if (exxtype.eq.1) then
   call mpi_grid_reduce(sig_x_aux(1,1,1,1),nbnd*nspinor*nkptnr*nvqr,&
                     &dims=(/dim_k/),all=.true.)
   sig_x_aux(:,:,:,:)=sig_x_aux(:,:,:,:)/nkptnr/omega
 endif

 do ik1=1,nkpt
  ik=kmap(1,ik1)
  do iq=1,nvqr
   do iqbz=1,iqrmap(2,iq)
    irk=rkmap(qqnrmap(iqbz,2,iq),ik)
    se_aux(:,:,:,ik1)=se_aux(:,:,:,ik1)+sig_c_aux(:,:,:,irk,iq)
    if (exxtype.eq.1) &
      sx_aux(:,:,ik1)=sx_aux(:,:,ik1)+sig_x_aux(:,:,irk,iq)
   enddo
  enddo
  do ikloc=1,nkptnrloc
   ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
   do i=1,kmap(2,ik1)
    if (ik.eq.kknrmap(i,1,ik1)) sig_c(:,:,:,ikloc)=se_aux(:,:,:,ik1)
    if (ik.eq.kknrmap(i,1,ik1).and.exxtype.eq.1) &
      & exxnk(:,:,ikloc)=sx_aux(:,:,ik1)
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

     sum1(:)=0.d0
     do ist=ibnd,fbnd
      sum1(:)=sum1(:)+sig_c(ist-bndrg(1,isp1)+1,isp1,:,ikloc)
     enddo
     sum1(:)=sum1(:)/(fbnd-ibnd+1)
     do ist=ibnd,fbnd
      sig_c(ist-bndrg(1,isp1)+1,isp1,:,ikloc)=sum1(:)
     enddo

     if (exxtype.eq.1) then
       sum1(1)=0.d0
       do ist=ibnd,fbnd
        sum1(1)=sum1(1)+exxnk(ist-bndrg(1,isp1)+1,isp1,ikloc)
       enddo
       sum1(1)=sum1(1)/(fbnd-ibnd+1)
       do ist=ibnd,fbnd
        exxnk(ist-bndrg(1,isp1)+1,isp1,ikloc)=sum1(1)
       enddo
     endif

    enddo !i
   enddo !isp1
  enddo !ikloc
 endif

 factor=-1.d0/nkptnr/omega
 sig_c(:,:,:,:)=sig_c(:,:,:,:)*factor

 if (mpi_grid_root((/dim_q/))) then
 ! sig_c isp1=1
  do ikloc=1,nkptnrloc
   ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
   write(fname,'("sigc_tau_I",I1.1,"_k",I4.4)') iter,ik
   do ist1=qpnb(1),qpnb(2)
    write(fspn,'("_ist",I3.3)') ist1
    fname_tot="Sig_c_files/"//trim(adjustl(fname))//trim(adjustl(fspn))
    open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
    i=ist1-qpnb(1)+1
    do iw=1,lr_nw
     write(168,'(2G18.8)') dreal(lr_w(iw))/ha2ev,&
            & sig_c(i,1,iw,ikloc)*ha2ev
    enddo
    close(168)
   enddo !ist1
  enddo !ikloc
 endif
 
! dyson equation
 call dyson_gf(iter,bndrg,gf0_tau,sig_c,exxnk,exxvc,vxcnk,vclnk,gf_tau)

! check for convergence
 maxdel=-100.d0
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   ibnd=bndrg(1,isp1)
   fbnd=bndrg(2,isp1)
   do ist1=ibnd,fbnd
     do iw=1,lr_nw
      c=gf_tau(ist1,iw,ikloc)-gf0_tau(ist1,iw,ikloc)
      diff=c**2
      if (diff.ge.maxdel) then
       maxdel=diff
      endif
     enddo
   enddo
  enddo
 enddo

 call mpi_grid_reduce(maxdel,1,dims=(/dim_k,dim_q/),all=.true.,op=op_max)
! in 1/eV
 maxdel=maxdel/ha2ev

 if (mpi_grid_root()) then
   write(151,'("maximum delta is, ",G12.6," 1/eV.")') maxdel
   write(151,'("maxde_prev is, ",G12.6," 1/eV.")') maxdel_prev
   call flushifc(151)
 endif

! store maxdel
 maxdel_prev=maxdel

! create tmp files for sig_c_aux  
 write(fname,'("Temp_files/tmp_Sc_ik",I4.4,"_iq",I4.4)') &
         &mpi_grid_dim_pos(dim_k),mpi_grid_dim_pos(dim_q)
 open(166,file=trim(adjustl(fname)),form='unformatted',status='replace')
 write(166) iter,iqloc,iter+1
 write(166) sig_c_aux  ! in Ha
 write(166) gf_tau
 write(166) sig_c  ! in Ha
 write(166) dyson_mu,maxdel_prev
 if (exxtype.eq.1) write(166) sig_x_aux,exxnk
 close(166)

 !reset zero
 sig_c_aux=0.d0
 se_aux=0.d0
 if (exxtype.eq.1) then
   sig_x_aux=0.d0
   sx_aux=0.d0
 endif
enddo !scgwni

! analytic continuation for Green's function using Pade approximation
call pade_ac_gf(bndrg,2*niw+1,gf_tau)
!sigma_c(iw)
sig_cc(:,:,:,:)=dcmplx(sig_c(:,:,:,:))
call ft_tw4(1,-1,nbnd,nspinor,lr_nw,nkptnrloc,niw,sig_cc(:,:,:,:),&
          &sig_cw)

! sig_cw
if (mpi_grid_root((/dim_q/))) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do isp1=1,nspinor
   write(fname,'("sigc_iw_k",I4.4,"_isp",I1.1)') ik,isp1
   do ist1=bndrg(1,isp1),bndrg(2,isp1)
    write(fspn,'("_ist",I3.3)') ist1
    fname_tot="Sig_c_files/"//trim(adjustl(fname))//trim(adjustl(fspn))
    open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")
    i=ist1-bndrg(1,isp1)+1
    do iw=1,niw
     write(168,'(3G18.8)') iw,dreal(sig_cw(i,1,iw,ikloc))*ha2ev,&
          & dimag(sig_cw(i,1,iw,ikloc))*ha2ev
    enddo
    close(168)
   enddo !ist1
  enddo !isp1
 enddo !ikloc
endif

if (mpi_grid_root()) then
 write(151,'("Now perform analytic continuation for self-energy!")')
 write(151,'(" ")')
 call flushifc(151)
endif

! total exchange self-energy
exxnk(:,:,:)=exxnk(:,:,:)+exxvc(:,:,:)

!analytic continuation for self-energy
call sigma_ac(151,nbnd,nspinor,niw,nkptnrloc,bndrg,evalmap,neval,&
             &exxnk,vxcnk,vclnk,sig_cw,sig_cr)

! shift according to the new Fermi energy
sig_cr(:,:,:,:)=sig_cr(:,:,:,:)-dyson_mu

! output quasiparticle energies
if (mpi_grid_root((/dim_q/))) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  write(fname,'("Eqp_k",I4.4)')ik
  do isp1=1,nspinor
   ibnd=bndrg(1,isp1)
   fbnd=bndrg(2,isp1)

   if (nspinor.eq.1) then
    fspn=''
   else
    write(fspn,'("_spn",I1.1)') isp1
   endif

   fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
   open(168,file=fname_tot,action='write',form="FORMATTED",status="REPLACE")

   do ist1=ibnd,fbnd
    i=ist1-ibnd+1
    dSe=sig_cr(2,i,isp1,ikloc)-sig_cr(1,i,isp1,ikloc)
    Zf=1.d0/(1.d0-dreal(dSe)/del_e)
    Eqp(i,isp1,ikloc)=evalsvnr(ist1,ik)+&
         & Zf*dreal(sig_cr(1,i,isp1,ikloc)+exxnk(i,isp1,ikloc)- &
         &  vxcnk(i,isp1,ikloc)+vclnk(i,isp1,ikloc))
    write(168,'(I5,1X,6(G16.6,1X))') ist1,Eqp(i,isp1,ikloc)*ha2ev,&
          & evalsvnr(ist1,ik)*ha2ev,exxnk(i,isp1,ikloc)*ha2ev,&
          & vxcnk(i,isp1,ikloc)*ha2ev,&
          & dreal(sig_cr(1,i,isp1,ikloc))*ha2ev,Zf
   enddo !ist1
   close(168)
  enddo !isp1
 enddo !ikloc
endif

! wait for all the processors
call mpi_grid_barrier()

if (mpi_grid_root()) then
 call timestamp(151)
 write(151,*)
 write(151,'("Done.")')
 close(151)
endif

deallocate(lr_w,bndrg)
deallocate(gf_tau)
deallocate(sig_c,exxnk,sig_cw,sig_cr)
deallocate(sig_c_aux,se_aux)
deallocate(sig_x_aux,sx_aux)
deallocate(vxcnk,vclnk,sum1)
deallocate(Eqp)
deallocate(kknrmap,kmap)
deallocate(iqrmap,qqnrmap,rkmap,evalmap,neval)
deallocate(sig_cc)

!if (mpi_grid_root().and.scgwni.gt.1) then
! call system('rm -f Temp_files/megq_*')
!endif

return
end subroutine
