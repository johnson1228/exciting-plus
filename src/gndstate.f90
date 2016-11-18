
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gndstate
! !INTERFACE:
subroutine gndstate
! !USES:
use modmain
use modldapu
use mod_wannier
use mod_sic
use mod_seceqn
use mod_libapw
! !DESCRIPTION:
!   Computes the self-consistent Kohn-Sham ground-state. General information is
!   written to the file {\tt INFO.OUT}. First- and second-variational
!   eigenvalues, eigenvectors and occupancies are written to the unformatted
!   files {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT}
!   and {\tt OCCSV.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
logical exist
integer ik,is,ia,ias,ist,i,ikloc
integer n,nwork
real(8) dv,etp,de,timetot
! allocatable arrays
real(8), allocatable :: v(:)
real(8), allocatable :: work(:)
real(8), allocatable :: evalfv(:,:,:)
! initialise global variables
call timer_start(t_init,reset=.true.)
call init0
call init1
call timer_stop(t_init)
if (.not.mpi_grid_in()) return
! only root processor writes
wproc=mpi_grid_root()
! allocate arrays for eigen-values/-vectors
if (tsveqn) then
  allocate(evalfv(nstfv,nspnfv,nkptloc))
  allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc))
  allocate(evecsvloc(nstsv,nstsv,nkptloc))
else
  allocate(evecfdloc(nspinor*nmatmax,nstsv,nkptloc))
endif
! initialise OEP variables if required
if (xctype(1).lt.0) call init2
if (wproc) then
! write the real and reciprocal lattice vectors to file
  call writelat
! write interatomic distances to file
  call writeiad(.false.)
! write symmetry matrices to file
  call writesym
! output the k-point set to file
  call writekpts
! write lattice vectors and atomic positions to file
  call writegeom(.false.)
! open INFO.OUT file
  open(ioinfo,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
! open TOTENERGY.OUT
  open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open FERMIDOS.OUT
  open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
   &form='FORMATTED')
! open FORCEMAX.OUT if required
!  if (tforce) open(64,file='FORCEMAX'//trim(filext),action='WRITE', &
!   &form='FORMATTED')
! open RMSDVEFF.OUT
  open(65,file='RMSDVEFF'//trim(filext),action='WRITE',form='FORMATTED')
! open DTOTENERGY.OUT
  open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open TENSMOM.OUT
  if (tmomlu) open(67,file='TENSMOM'//trim(filext),action='WRITE', &
   &form='FORMATTED')
! write out general information to INFO.OUT
  call writeinfo(ioinfo)
  write(ioinfo,*)
  write(ioinfo,'("initialization time : ",F10.2," sec.")')timer_get_value(t_init)
  call writenn
  call writegshells
endif
! size of mixing vector
n=lmmaxvr*nrmtmax*natmtot+ngrtot
if (spinpol) n=n*(1+ndmag)
if (ldapu.ne.0) n=n+2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
! allocate mixing arrays
allocate(v(n))
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,n,v,dv,nwork,v)
allocate(work(nwork))
! initialise or read the charge density and potentials from file
iscl=0
if (wproc) write(ioinfo,*)
if ((task.eq.1).or.(sic.and.isclsic.gt.1)) then
  call readstate
  if (wproc) write(ioinfo,'("Potential read in from STATE.OUT")')
  if (autolinengy) call readfermi
else if (task.eq.200) then
  call phveff
  if (wproc) write(ioinfo,'("Supercell potential constructed from STATE.OUT")')
else
  call timer_start(t_rhoinit,reset=.true.)
  call allatoms
  do ias=1,natmtot
    is=ias2is(ias)
    do ist=1,spnst(is)
      evalcr(ist,ias)=speval(ist,is)
    enddo
  enddo
  call rhoinit
  call poteff
  call timer_stop(t_rhoinit)
  if (wproc) then
    write(ioinfo,'("Density and potential initialised from atomic data")')
    write(ioinfo,'(" done in : ",F10.2," sec.")')timer_get_value(t_rhoinit)
    write(ioinfo,*)
    do is=1,nspecies
      write(ioinfo,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
      do ist=1,spnst(is)
        write(ioinfo,'("  n : ",I1,"   l : ",I1,"   k : ",I1,"   energy : ",G18.10)')&
          &spn(ist,is),spl(ist,is),spk(ist,is),speval(ist,is)
      enddo
      write(ioinfo,*)
    enddo
  endif
end if
if (wproc) call flushifc(ioinfo)
if (sic) then
  call sic_read_data(.true.)
  call sic_genblochsum(.true.,.true.)
endif
! set stop flag
tstop=.false.
!10 continue
! set last iteration flag
tlast=.false.
etp=0.d0
! delete any existing eigenvector files
!if (wproc.and.(task.eq.0).or.(task.eq.2)) call delevec
! begin the self-consistent loop
if (wproc) then
  write(ioinfo,*)
  write(ioinfo,'("+------------------------------+")')
  write(ioinfo,'("| Self-consistent loop started |")')
  write(ioinfo,'("+------------------------------+")')
endif
do iscl=1,maxscl
! reset all timers
  call timer_reset(0)
  call timer_start(t_iter_tot)
  if (wproc) then
    write(ioinfo,*)
    write(ioinfo,'("+-------------------------+")')
    write(ioinfo,'("| Iteration number : ",I4," |")') iscl
    write(ioinfo,'("+-------------------------+")')
    call flushifc(ioinfo)
  endif
  if (iscl.ge.maxscl) then
    if (wproc) then
      write(ioinfo,*)
      write(ioinfo,'("Reached self-consistent loops maximum")')
    endif
    tlast=.true.
  end if
  if (wproc) call flushifc(ioinfo)
! Fourier transform effective potential to G-space
  call genveffig  
! generate the core wavefunctions and densities
  call gencore
! generate radial functions
  call genradf
! write out the linearisation energies
  if (wproc) call writelinen
! compute the overlap radial integrals
  call timer_start(t_hbo_rad)
  call olprad
#ifdef _LIBAPW_
  call libapw_seceqn_init
#else
! compute the Hamiltonian radial integrals
  call hmlrad
! generate effective magntic field integrals for full diagonalization
  call genbeff
! generate muffin-tin effective magnetic fields and s.o. coupling functions
  if (texactrho) then
    call seceqnsv_init
  else
    call genbeffmt
  endif 
#endif
  call timer_stop(t_hbo_rad)
  evalsv(:,:)=0.d0
  if (sic) sic_evalsum=0.d0
  if (iscl.eq.1.and.ldapu.ne.0) then
    call init_vmatlu
  endif
! begin parallel loop over k-points
  do ikloc=1,nkptloc
! solve the secular equation
    if (tsveqn) then
      call seceqn(ikloc,evalfv(1,1,ikloc),evecfvloc(1,1,1,ikloc),&
        &evecsvloc(1,1,ikloc))
    else
      call seceqnfd(ikloc,evecfdloc(1,1,ikloc))
    endif
  end do  
  call mpi_grid_reduce(evalsv(1,1),nstsv*nkpt,dims=(/dim_k/),all=.true.)
  if (sic) call sic_write_eval
  if (wproc) then
! find the occupation numbers and Fermi energy
    call occupy
    if (autoswidth) then
      write(ioinfo,*)
      write(ioinfo,'("New smearing width : ",G18.10)') swidth
    end if    
! write out the eigenvalues and occupation numbers
    call writeeval
! write the Fermi energy to file
    call writefermi
    call flushifc(ioinfo)
  endif
  call mpi_grid_bcast(swidth)
  call mpi_grid_bcast(occsv(1,1),nstsv*nkpt)
  if (wannier) call wann_ene_occ
#ifdef _LIBAPW_
  call libapw_rhomag
#else  
  if (texactrho.or..not.tsveqn) then
    call rhomag_exact
  else
    call rhomag
  endif
#endif
  if (sic) call mpi_grid_reduce(sic_evalsum,dims=(/dim_k/))
! LDA+U
  if (iscl.gt.1.and.ldapu.ne.0) then
! generate the LDA+U density matrix
    call gendmatrsh
! generate the LDA+U potential matrix
    call genvmatlu
! write the LDA+U matrices to file
    if (wproc) call writeldapu
! calculate and write tensor moments to file
    if (tmomlu.and.wproc) then
      call tensmom(67)
      call flushifc(67)
    end if
  end if
! compute the effective potential
  call poteff
! pack interstitial and muffin-tin effective potential and field into one array
  call mixpack(.true.,n,v)
! mix in the old potential and field with the new
  call mixerifc(mixtype,n,v,dv,nwork,work)
! unpack potential and field
  call mixpack(.false.,n,v)
! add the fixed spin moment effect field
  if (fixspin.ne.0) call fsmfield
! reduce the external magnetic fields if required
  if (reducebf.lt.1.d0) then
    bfieldc(:)=bfieldc(:)*reducebf
    bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
  end if
! compute the energy components
  call energy
  call timer_stop(t_iter_tot)
! output energy components
  if (wproc) then
    call writeengy(ioinfo)
    write(ioinfo,*)
    write(ioinfo,'("Density of states at Fermi energy : ",G18.10)') fermidos
    write(ioinfo,'(" (states/Hartree/unit cell)")')
    write(ioinfo,'("Estimated band gap (eV) : ",G18.10)') bandgap(1)*ha2ev
    write(ioinfo,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
    write(ioinfo,'("Estimated direct band gap (eV) : ",G18.10)') bandgap(2)*ha2ev
    write(ioinfo,'(" at k-point ",I6)') ikgap(3)
    call flushifc(ioinfo)
! write total energy to TOTENERGY.OUT and flush
    write(61,'(G22.12)') engytot
    call flushifc(61)
! write DOS at Fermi energy to FERMIDOS.OUT and flush
    write(62,'(G18.10)') fermidos
    call flushifc(62)
! output charges and moments
    call writechg(60)
! write total moment to MOMENT.OUT and flush
    if (spinpol) then
      write(63,'(3G18.10)') momtot(1:ndmag)
      call flushifc(63)
    end if
! output effective fields for fixed spin moment calculations
    if (fixspin.ne.0) call writefsm(60)
! check for WRITE file
    inquire(file='WRITE',exist=exist)
    if (exist) then
      write(ioinfo,*)
      write(ioinfo,'("WRITE file exists - writing STATE.OUT")')
      call writestate
      open(50,file='WRITE')
      close(50,status='DELETE')
    end if
! write STATE.OUT file if required
    if (nwrite.ge.1) then
      if (mod(iscl,nwrite).eq.0) then
        call writestate
        write(ioinfo,*)
        write(ioinfo,'("Wrote STATE.OUT")')
      end if
    end if
  endif ! wproc
  call mpi_grid_bcast(tlast)
! exit self-consistent loop if last iteration is complete
  if (tlast) goto 20
! check for convergence
  if (wproc) then
    if (iscl.ge.2) then
      write(ioinfo,*)
      write(ioinfo,'("RMS change in effective potential (target) : ",G18.10," (",&
       &G18.10,")")') dv,epspot
      write(65,'(G18.10)') dv
      call flushifc(65)
      de=abs(engytot-etp)
      write(60,'("Absolute change in total energy (target)   : ",G18.10," (",&
       &G18.10,")")') de,epsengy
      write(66,'(G18.10)') de
      call flushifc(66)
      if ((dv.lt.epspot).and.(de.lt.epsengy)) then
        write(ioinfo,*)
        write(ioinfo,'("Convergence targets achieved")')
        tlast=.true.
      end if
    end if
    etp=engytot
    if (xctype(1).lt.0) then
      write(ioinfo,*)
      write(ioinfo,'("Magnitude of OEP residual : ",G18.10)') resoep
    end if
! check for STOP file
    inquire(file='STOP',exist=exist)
    if (exist) then
      write(ioinfo,*)
      write(ioinfo,'("STOP file exists - stopping self-consistent loop")')
      tstop=.true.
      tlast=.true.
      open(50,file='STOP')
      close(50,status='DELETE')
    end if
! output the current total CPU time
    timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
    write(ioinfo,*)
    write(ioinfo,'("Time (CPU seconds) : ",F12.2)') timetot
    write(ioinfo,*)
    write(ioinfo,'("iteration time (seconds)                    : ",F12.2)')&
      &timer_get_value(t_iter_tot)
    write(ioinfo,'("  radial setup (linen, APW, radint)         : ",3F12.2)')&
      &timer_get_value(t_lin_en),&
      &timer_get_value(t_apw_rad),&
      &timer_get_value(t_hbo_rad)
    write(ioinfo,'("  total for secular equation                : ",F12.2)')&
      &timer_get_value(t_seceqn)
    write(ioinfo,'("    first-variational                       : ",F12.2)')&
      &timer_get_value(t_seceqnfv)
    write(ioinfo,'("      setup                                 : ",F12.2)')&
      &timer_get_value(t_seceqnfv_setup)
    write(ioinfo,'("        setup H (total, MT, IT)             : ",3F12.2)')&
      &timer_get_value(t_seceqnfv_setup_h),&
      &timer_get_value(t_seceqnfv_setup_h_mt),&
      &timer_get_value(t_seceqnfv_setup_h_it)
    write(ioinfo,'("        setup O (total, MT, IT)             : ",3F12.2)')&
      &timer_get_value(t_seceqnfv_setup_o),&
      &timer_get_value(t_seceqnfv_setup_o_mt),&
      &timer_get_value(t_seceqnfv_setup_o_it)
    write(ioinfo,'("      diagonalization                       : ",F12.2)')&
      &timer_get_value(t_seceqnfv_diag)
    write(ioinfo,'("    second-variational                      : ",F12.2)')&
      &timer_get_value(t_seceqnsv)
    write(ioinfo,'("      setup (total, MT, IT)                 : ",3F12.2)')&
      &timer_get_value(t_seceqnsv_setup),&
      &timer_get_value(t_seceqnsv_setup_mt),&
      &timer_get_value(t_seceqnsv_setup_it)
    write(ioinfo,'("      diagonalization                       : ",F12.2)')&
      &timer_get_value(t_seceqnsv_diag)
    write(ioinfo,'("  total for charge and magnetization        : ",F12.2)')&
      &timer_get_value(t_rho_mag_tot)
    write(ioinfo,'("    k-point summation (total, MT, IT)       : ",3F12.2)')&
      &timer_get_value(t_rho_mag_sum),&
      &timer_get_value(t_rho_mag_mt),&
      &timer_get_value(t_rho_mag_it)
    if (texactrho) then
      write(ioinfo,'("      wave-function setup                   : ",F12.2)')&
        &timer_get_value(t_rho_wf)
      write(ioinfo,'("      convert to r-mesh                     : ",F12.2)')&
        &timer_get_value(t_rho_mag_conv)
    endif
    write(ioinfo,'("    symmetrization                          : ",F12.2)')&
      &timer_get_value(t_rho_mag_sym)
    write(ioinfo,'("  potential (total, XC, Ha)                 : ",3F12.2)')&
      &timer_get_value(t_pot),&
      &timer_get_value(t_pot_xc),&
      &timer_get_value(t_pot_ha)
    write(ioinfo,'("  density matrix setup                      : ",F12.2)')&
      &timer_get_value(t_dmat)
    if (sic) then
      write(ioinfo,'("  sic_genfvprj                              : ",F12.2)')&
        &timer_get_value(t_sic_genfvprj)
      write(ioinfo,'("  sic_hunif                                 : ",F12.2)')&
        &timer_get_value(t_sic_hunif)
    endif
  endif !wproc
! end the self-consistent loop
end do !iscl
20 continue
if (wproc) then
  write(ioinfo,*)
  write(ioinfo,'("+------------------------------+")')
  write(ioinfo,'("| Self-consistent loop stopped |")')
  write(ioinfo,'("+------------------------------+")')
  call flushifc(ioinfo)
! write density and potentials to file only if maxscl > 1
  if (maxscl.gt.1) then
    call writestate
    write(ioinfo,*)
    write(ioinfo,'("Wrote STATE.OUT")')
    call flushifc(ioinfo)
  end if
endif !wproc
! write eigenvalues/vectors and occupancies to file
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
        if (tsveqn) then
          call putevalfv(ik,evalfv(1,1,ikloc))
          call putevecfv(ik,evecfvloc(1,1,1,ikloc))
          call putevecsv(ik,evecsvloc(1,1,ikloc))
        else
          call putevecfd(ikloc,evecfdloc(1,1,ikloc))
        endif
        call putevalsv(ik,evalsv(1,ik))
        call putoccsv(ik,occsv(1,ik))
      end do
    end if
    call mpi_grid_barrier(dims=(/dim_k/))
  end do
endif
! write SIC related functions
if (sic) then
  call sic_write_blochsum
endif
call mpi_grid_bcast(tstop)
!-----------------------!
!     compute forces    !
!-----------------------!
if (tforce) then
  call force
! output forces to INFO.OUT
  if (wproc) then
    call writeforce(ioinfo)
    call flushifc(ioinfo)
! write maximum force magnitude to FORCEMAX.OUT
!    write(64,'(G18.10)') forcemax
!    call flushifc(64)
  endif
end if
! output timing information
if (wproc) then
  write(ioinfo,*)
  write(ioinfo,'("Timings (CPU seconds) :")')
  write(ioinfo,'(" initialisation",T40,": ",F12.2)') timeinit
  write(ioinfo,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
  write(ioinfo,'(" first-variational secular equation",T40,": ",F12.2)') timefv
  if (spinpol) then
    write(ioinfo,'(" second-variational calculation",T40,": ",F12.2)') timesv
  end if
  write(ioinfo,'(" charge density calculation",T40,": ",F12.2)') timerho
  write(ioinfo,'(" potential calculation",T40,": ",F12.2)') timepot
  if (tforce) then
    write(ioinfo,'(" force calculation",T40,": ",F12.2)') timefor
  end if
  timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
  write(ioinfo,'(" total",T40,": ",F12.2)') timetot
  write(ioinfo,*)
  write(ioinfo,'("+-----------------------------+")')
  write(ioinfo,'("| Elk version ",I1.1,".",I1.1,".",I3.3," stopped |")') version
  write(ioinfo,'("+-----------------------------+")')
! close the INFO.OUT file
  close(ioinfo)
! close the TOTENERGY.OUT file
  close(61)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT file
  if (spinpol) close(63)
! close the FORCEMAX.OUT file
!  if (tforce) close(64)
! close the RMSDVEFF.OUT file
  close(65)
! close the DTOTENERGY.OUT file
  close(66)
! close TENSMOM.OUT file
  if (tmomlu) close(67)
endif
deallocate(v,work)
if (tsveqn) then
  deallocate(evalfv)
  deallocate(evecfvloc)
  deallocate(evecsvloc)
else
  deallocate(evecfdloc)
endif
#ifdef _LIBAPW_
call lapw_timers
#endif
return
end subroutine
!EOC
