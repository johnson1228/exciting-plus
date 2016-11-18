subroutine writewann
use modmain
use modldapu
use mod_nrkp
use mod_linresp
implicit none
integer i,j,ik,ikloc,n1,n2,j1,j2,nwan,ias,n,is,ia
integer t1,t2,t3,ntv,ibnd,fbnd,npp1dloc,bndtmp,it,ik1,isp1
integer t1tot,t2tot,t3tot
real(8) eng,mag
complex(8) arg
logical lpmat
complex(8), allocatable :: zm(:,:)
complex(8), allocatable :: tb_ht(:,:,:),tb_hk(:,:)
real(8), allocatable :: dm(:,:),eval(:)
real(8), allocatable :: gw_bs(:,:)
character*100 fname, fname_tot, fspn
call init0
call init1

wproc=mpi_grid_root()
if (.not.wannier) then
  write(*,*)
  write(*,'("Error(writewann_h) : WF generation is switched off")')
  write(*,*)
  call pstop
endif
! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp
lpmat=.false.
if (task.eq.808) lpmat=.true.
call genwfnr(6,.false.)
if (allocated(wann_h)) deallocate(wann_h)
allocate(wann_h(nwantot,nwantot,nkptnr))
wann_h=zzero
if (allocated(wann_e)) deallocate(wann_e)
allocate(wann_e(nwantot,nkptnr))
wann_e=0.d0

call occupy
!print out the old Fermi energy 
if (mpi_grid_root()) write(*,'("Old Fermi Energy:",G18.10)') efermi*ha2ev
! old magnetization 
if (nspinor.eq.2) then 
 mag=0.d0 
 do ik=1,nkptnr 
  do i=1,int(nstsv/2) 
   mag=mag+(occsvnr(i,ik)-occsvnr(i+int(nstsv/2),ik))/nkptnr
  enddo 
 enddo 
 if (mpi_grid_root()) write(*,'("Old magnetization:",G18.10)') mag
endif 

! update evalsvnr using GW quasiparticle energies
if (gwband.and.eval_udt) then
 do isp1=1,nspinor
   if (isp1.eq.1) then
     ibnd=qpnb(1)
     fbnd=qpnb(2)
     if (nspinor.eq.1) then !not spin-polarized
       fspn=''
     elseif (nspinor.eq.2) then
       write(fspn,'("_spn",I1.1)') isp1
     endif
   elseif (isp1.eq.2) then  !spin down
     ibnd=qpnb(1)+int(nstsv/2)
     fbnd=qpnb(2)+int(nstsv/2)
     write(fspn,'("_spn",I1.1)') isp1
   endif

   do ik=1,nkptnr
    write(fname,'("Eqp_k",I4.4)')ik
    fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
    open(unit=50,file=fname_tot,action='read',form="formatted",status='old')
    do i=ibnd,fbnd
     read(50,*) bndtmp,eng
     evalsvnr(i,ik)=eng/ha2ev+efermi !efermi is the DFT Fermi energy
    enddo
    close(50)
   enddo !ik
 enddo !isp1

! make sure evalsvnr is identical among the processors
 call mpi_grid_bcast(evalsvnr(1,1),nstsv*nkptnr)
! make sure svalsv are updated at the same time, for occupy.f90
 do ik=1,nkpt
  do ik1=1,nkptnr
   if (abs(vkl(1,ik)-vklnr(1,ik1))+abs(vkl(2,ik)-vklnr(2,ik1))+&
     & abs(vkl(3,ik)-vklnr(3,ik1)).le.1d-10) then
     evalsv(:,ik)=evalsvnr(:,ik1)
     if (mpi_grid_root()) write(*,*) "ik,ik1:",ik,ik1
     exit
   endif
  enddo
 enddo

! find out the GW Fermi energy
 call occupy
 if (mpi_grid_root()) write(*,'("New Fermi Energy:",G18.10)') efermi*ha2ev
! new magnetization
 if (nspinor.eq.2) then
  mag=0.d0
  do ik=1,nkpt
   do i=1,int(nstsv/2)
    mag=mag+(occsv(i,ik)-occsv(i+int(nstsv/2),ik))*wkpt(ik)
   enddo
  enddo
  if (mpi_grid_root()) write(*,'("New magnetic moment:",G18.10)') mag
 endif

endif !gwband

! for each k, obtain H(k), i.e. wann_h(:,:,k)
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call genwann_h(.true.,evalsvnr(1,ik),wanncnrloc(1,1,ikloc),&
    & wann_h(1,1,ik),wann_e(1,ik))
enddo

!allocate(wann_ene_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
!allocate(wann_occ_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
!call wann_ene_occ_(wann_ene_m,wann_occ_m)

! get h(n,n',k) for all the processors
call mpi_grid_reduce(wann_h(1,1,1),nwantot*nwantot*nkptnr,dims=(/dim_k/),side=.true.)
call mpi_grid_reduce(wann_e(1,1),nwantot*nkptnr,dims=(/dim_k/),side=.true.)
!call mpi_grid_reduce(wann_p(1,1,1,1),3*nwantot*nwantot*nkpt,dims=(/dim_k/),side=.true.)

! GW band structure
if (gwband) then
 t1tot=ngridk(1)
 t2tot=ngridk(2)
 t3tot=ngridk(3)
 ntv=t1tot*t2tot*t3tot
! get the k-vector along the k path
 call connect(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
! solve the band structure in parallel
 npp1dloc=mpi_grid_map(npp1d,dim_k)
 allocate(tb_ht(nwantot,nwantot,ntv))
 allocate(tb_hk(nwantot,nwantot))
 allocate(gw_bs(nwantot,npp1d))
 tb_ht=zzero
 gw_bs=0.d0
 it=0
 !first, obtain h(T) using h(k)
 do t3=-int((t3tot-1)/2),int(t3tot/2)
  do t2=-int((t2tot-1)/2),int(t2tot/2)
   do t1=-int((t1tot-1)/2),int(t1tot/2)
    it=it+1
    do ik=1,nkptnr
     arg=-zi*twopi*(vklnr(1,ik)*t1+vklnr(2,ik)*t2+vklnr(3,ik)*t3)
     tb_ht(:,:,it)=tb_ht(:,:,it)+exp(arg)*wann_h(:,:,ik)
    enddo
   enddo
  enddo
 enddo

! tb_ht(T) is the tight-binding Hamiltonian H(T)
 call mpi_grid_bcast(tb_ht(1,1,1),nwantot*nwantot*ntv)
 tb_ht(:,:,:)=tb_ht(:,:,:)/nkptnr 

! solve the band structure in parallel
 if (mpi_grid_root()) write(*,'("npp1d:",I4,1X,I4)') npp1dloc, npp1d
 do ikloc=1,npp1dloc
  ik=mpi_grid_map(npp1d,dim_k,loc=ikloc)
  write(*,'("Info(GW_bandstr): ",I6," of ",I6," k-points, ")') ik,npp1d
  it=0
  tb_hk=zzero
  do t3=-int((t3tot-1)/2),int(t3tot/2)
   do t2=-int((t2tot-1)/2),int(t2tot/2)
    do t1=-int((t1tot-1)/2),int(t1tot/2)
     it=it+1
     arg=zi*twopi*(vplp1d(1,ik)*t1+vplp1d(2,ik)*t2+vplp1d(3,ik)*t3)
     tb_hk(:,:)=tb_hk(:,:)+exp(arg)*tb_ht(:,:,it)
    enddo
   enddo
  enddo
! tb_hk will chanage after calling diagzhe
  call diagzhe(nwantot,tb_hk(1,1),gw_bs(1,ik))
 enddo !ikloc

! collect gw_bs from all the processors
 call mpi_grid_reduce(gw_bs(1,1),nwantot*npp1d,dims=(/dim_k/),&
                    & side=.true.)

! output GW-corrected band structure
 if (mpi_grid_root()) then
  write(fname,'("gw_bs.dat")')
  if (nspinor.eq.1) then
   fspn=''
  else
   write(fspn,'("_spn",I1.1)') gwband_sp
  endif

  fname_tot=trim(adjustl(fname))//trim(adjustl(fspn))
  open(50,file=fname_tot,form='formatted',status='replace')

  do i=1,nwantot
   do ik=1,npp1d
    write(50,'(2G18.10)') dpp1d(ik),(gw_bs(i,ik)-efermi)*ha2ev
   enddo
   write(50,'(   )')
  enddo
  close(50)

! output the GW-corrected tight-binding Hamiltonian
  open(50,file="GW_TB_HT.OUT",form="FORMATTED",status="REPLACE")
  it=0
  do t3=-int((t3tot-1)/2),int(t3tot/2)
   do t2=-int((t2tot-1)/2),int(t2tot/2)
    do t1=-int((t1tot-1)/2),int(t1tot/2)
     it=it+1
     write(50,'("# T-vector: ",I8)')it
     write(50,'("# lattice coordinates")')
     write(50,'(3G18.10)') t1,t2,t3
     write(50,'("# real part of tb_ht")')
     do i=1,nwantot
       write(50,'(255G18.10)')(dreal(tb_ht(i,j,it)),j=1,nwantot)
     enddo
     write(50,'("# imaginary part of H")')
     do i=1,nwantot
       write(50,'(255G18.10)')(dimag(tb_ht(i,j,it)),j=1,nwantot)
     enddo
    enddo
   enddo
  enddo
  close(50)
 endif

 deallocate(tb_ht,tb_hk,gw_bs)

endif !gwband

if (mpi_grid_root().and.task.eq.807) then
  call readfermi
  open(200,file="WANN_H.OUT",form="FORMATTED",status="REPLACE")
  write(200,'("# units of energy are Hartree, 1 Ha=",F18.10," eV")')ha2ev
  write(200,'("# fermi energy")')
  write(200,'(G18.10)')efermi !new fermi energy if gwband=.true.
  write(200,'("# lattice vectors (3 rows)")')
  do i=1,3
    write(200,'(3G18.10)')avec(:,i)
  enddo
  write(200,'("# reciprocal lattice vectors (3 rows)")')
  do i=1,3
    write(200,'(3G18.10)')bvec(:,i)
  enddo
  write(200,'("# k-grid size")')
  write(200,'(3I6)')ngridk
  write(200,'("# number of k-points")')
  write(200,'(I8)')nkptnr
  write(200,'("# number of Wannier functions")')
  write(200,'(I8)')nwantot
  write(200,'("# number of atoms")')
  write(200,'(I8)')natmtot 
  write(200,'("# number of species")')
  write(200,'(I8)')nspecies  
  write(200,'("# wf -> atom mapping")')
  do n=1,nwantot
    write(200,'(2I8)')n,wan_info(wi_atom,n)
  enddo
  write(200,'("# atom -> species mapping")')
  do ias=1,natmtot
    is=ias2is(ias)
    write(200,'(2I8)')ias,is
  enddo
  do ias=1,natmtot
    ia=ias2ia(ias)
    is=ias2is(ias)
    write(200,'("# atom : ",I8)')ias
    write(200,'("# Cartesian coordinates")')
    write(200,'(3G18.10)')atposc(:,ia,is)
    write(200,'("# lattice coordinates")')
    write(200,'(3G18.10)')atposl(:,ia,is)
  enddo

  do ik=1,nkptnr
    write(200,'("# k-point : ",I8)')ik
    write(200,'("# weight")')
    write(200,'(G18.10)')wkptnr(ik)
    write(200,'("# lattice coordinates")')
    write(200,'(3G18.10)')vklnr(:,ik)
    write(200,'("# Cartesian coordinates")')
    write(200,'(3G18.10)')vkcnr(:,ik)
    write(200,'("# real part of H")')
    do i=1,nwantot
      write(200,'(255G18.10)')(dreal(wann_h(i,j,ik)),j=1,nwantot)
    enddo
    write(200,'("# imaginary part of H")')
    do i=1,nwantot
      write(200,'(255G18.10)')(dimag(wann_h(i,j,ik)),j=1,nwantot)
    enddo
    write(200,'("# eigen-values of H")')
    write(200,'(255G18.10)')(wann_e(j,ik),j=1,nwantot)
  enddo
  close(200)
endif

if (mpi_grid_root()) then
  do ias=1,natmtot
    nwan=nwannias(ias)
    if (nwan.ne.0) then
      write(*,*)"ias : ",ias,"  nwan : ",nwan
      allocate(zm(nwan,nwan))
      allocate(dm(nwan,nwan))
      allocate(eval(nwan))
      zm=zzero
      dm=zzero
      j1=0
      do n1=1,nwantot
        if (wan_info(wi_atom,n1).eq.ias) then
          j1=j1+1
          j2=0
          do n2=1,nwantot
            if (wan_info(wi_atom,n2).eq.ias) then
              j2=j2+1
              zm(j1,j2)=wann_h(n1,n2,1)
              dm(j1,j2)=dreal(zm(j1,j2))
            endif
          enddo
        endif
      enddo
      write(*,*)"Hamiltonian at gamma point :"
      do j1=1,nwan
        write(*,'(255F12.6)')(dreal(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)
      do j1=1,nwan
        write(*,'(255F12.6)')(dimag(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)"eigen-vectors : "
      call diagdsy(nwan,dm,eval)
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j1,j2),j2=1,nwan)
      enddo
      write(*,*)
      write(*,'(2X,7G18.10)')(eval(j1),j1=1,nwan)
      write(*,*)
      write(*,*)"transpose of eigen-vectors : "
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j2,j1),j2=1,nwan)
      enddo
      write(*,*)
      deallocate(zm,dm,eval)
    endif
  enddo
  write(*,*)
  do ias=1,natmtot
    nwan=nwannias(ias)
    if (nwan.ne.0) then
      write(*,*)"ias : ",ias,"  nwan : ",nwan
      allocate(zm(nwan,nwan))
      zm=zzero
      allocate(dm(nwan,nwan))
      allocate(eval(nwan))
      j1=0
      do n1=1,nwantot
        if (wan_info(wi_atom,n1).eq.ias) then
          j1=j1+1
          j2=0
          do n2=1,nwantot
            if (wan_info(wi_atom,n2).eq.ias) then
              j2=j2+1
              do ik=1,nkptnr
                zm(j1,j2)=zm(j1,j2) + wann_h(n1,n2,ik)/nkptnr
              enddo
              dm(j1,j2)=dreal(zm(j1,j2))
            endif
          enddo
        endif
      enddo
      write(*,*)"Average Hamiltonian :"
      do j1=1,nwan
        write(*,'(255F12.6)')(dreal(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)
      do j1=1,nwan
        write(*,'(255F12.6)')(dimag(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)"eigen-vectors : "
      call diagdsy(nwan,dm,eval)
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j1,j2),j2=1,nwan)
      enddo
      write(*,*)
      write(*,'(2X,7G18.10)')(eval(j1),j1=1,nwan)
      write(*,*)
      write(*,*)"transpose of eigen-vectors : "
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j2,j1),j2=1,nwan)
      enddo
      write(*,*)
      deallocate(zm,dm,eval)
    endif
  enddo
endif

!deallocate(wann_ene_m,wann_occ_m)
return
end
