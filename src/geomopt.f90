
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine geomopt
use modmain
use mod_mpi_grid
use mod_seceqn

implicit none
! local variables
integer istp,jstp,i,j,fnum1
real(8) t1
logical :: wproc1
! initialise global variables
call init0
call init1
call initmpigrid 
! store orginal volume
if (.not.mpi_grid_in()) return
wproc1=mpi_grid_root()
omega0=omega
! atomic forces are required
tforce=.true.

! file unit number for 'RELAX.OUT'
fnum1=1111

! initial atomic step sizes
if (allocated(tauatm)) deallocate(tauatm)
allocate(tauatm(natmtot))
tauatm(:)=tau0atm
! initialise the previous total force on each atom
if (allocated(forcetp)) deallocate(forcetp)
allocate(forcetp(3,natmtot))
forcetp(:,:)=0.d0
! initial lattice vector step size
!taulatv(:)=tau0latv
! initialise previous stress matrix
!stressp(:,:)=0.d0
! open INFO.OUT file

if (wproc1) then
 open(fnum1,file='RELAX.OUT',form='FORMATTED',status='replace')
 ! open TOTENERGY.OUT
 open(71,file='TOTENERGY_OPT.OUT',form='FORMATTED',status='replace')
 ! open FORCEMAX.OUT
 open(72,file='FORCEMAX.OUT',form='FORMATTED',status='replace')
 ! open GEOMETRY_OPT.OUT
 open(73,file='GEOMETRY_OPT.OUT',form='FORMATTED',status='replace')
 ! open IADIST_OPT.OUT
 open(74,file='IADIST_OPT.OUT',form='FORMATTED',status='replace')
 ! open FORCES_OPT.OUT
 open(75,file='FORCES_OPT.OUT',form='FORMATTED',status='replace')
! open STRESSMAX.OUT and STRESS_OPT.OUT if required
!if (latvopt.ne.0) then
!  open(76,file='STRESSMAX.OUT',action='WRITE',form='FORMATTED')
!  open(77,file='STRESS_OPT.OUT',action='WRITE',form='FORMATTED')
!end if
 write(fnum1,*)
 call flushifc(fnum1)
endif

! change it later
!if (task.eq.3) then
!  call readstate
!  if (wproc) write(fnum1,'("Potential read in from STATE.OUT")')
!  if (autolinengy) call readfermi
!end if

do istp=1,maxlatvstp
  do jstp=1,maxatmstp
    if (wproc1) then
      write(fnum1,'(" ")')
      write(fnum1,'("Info(geomopt): atomic position optimisation step : ",I6)') jstp
      call timestamp(fnum1)
      call flushifc(fnum1)
    endif
! ground-state and forces calculation
    call gndstate
! subsequent calculations will read in the potential from STATE.OUT
!    trdstate=.true.
! update the atomic positions
!    if (wproc) then
!     write(fnum1,*) "before updatpos"
!     call flushifc(fnum1)
!    endif
    call updatpos
! write total energy, forces, atomic positions, interatomic distances to file
    if (wproc1) then
      call flushifc(fnum1)
      write(71,'(G22.12)') engytot
      call flushifc(71)
      write(72,'(G18.10)') forcemax
      call flushifc(72)
      write(73,*)
      write(73,*)
      write(73,'("! Lattice and atomic position optimisation steps : ",2I6)') istp,jstp
      call writegeom(73)
      call flushifc(73)
      write(74,*)
      write(74,*)
      write(74,'("Lattice and atomic position optimisation steps : ",2I6)') istp,jstp
      call writeiad(73)
      call flushifc(74)
      write(75,*)
      write(75,*)
      write(75,'("Lattice and atomic position optimisation steps : ",2I6)') istp,jstp
      call writeforce(75)
      write(75,*)
      write(75,'("Maximum force magnitude over all atoms (target) : ",G18.10," (",G18.10,")")') forcemax,epsforce
      call flushifc(75)
    end if
! check force convergence
    if (forcemax.le.epsforce) then
      if (wproc1) then
        write(75,*)
        write(75,'("Force convergence target achieved")')
      end if
      exit
    end if
    if (wproc1.and.(jstp.eq.maxatmstp)) then
      write(fnum1,*)
      write(fnum1,'("Warning(geomopt): atomic position optimisation failed to &
       &converge in ",I6," steps")') maxatmstp
    end if

! store the current forces array
    forcetp(:,:)=forcetot(:,:)
! end loop over atomic position optimisation
  end do
! exit lattice optimisation loop if required
  if (latvopt.eq.0) exit
  if (wproc1) then
    write(fnum1,'("Info(geomopt): lattice vector optimisation step : ",I6)') istp
  end if
! generate the stress matrix
!  call genstress
! update the lattice vectors
!  call latvstep
! write stress magnitude and matrix to file
!  if (mp_mpi) then
!    write(76,'(G18.10)') stressmax
!    call flushifc(76)
!    write(77,*)
!    write(77,'("Lattice vector optimisation step : ",I6)') istp
!    do j=1,3
!      write(77,'(3G18.10)') (stress(i,j),i=1,3)
!    end do
!    call flushifc(77)
!  end if
! check for stress convergence; stress may be non-zero because of volume
! constraint; checking change in stress matrix instead
!  t1=sum(abs(stress(:,:)-stressp(:,:)))
!  if (t1.le.epsstress*tau0latv) then
!    if (mp_mpi) then
!      write(77,*)
!      write(77,'("Stress convergence target achieved")')
!    end if
!    exit
!  end if
!  if (mp_mpi.and.(istp.eq.maxlatvstp)) then
!    write(*,*)
!    write(*,'("Warning(geomopt): lattice vector optimisation failed to &
!     &converge in ",I6," steps")') maxlatvstp
!  end if
! store the current stress matrix
!  stressp(:,:)=stress(:,:)

  call mpi_grid_barrier(dims=(/dim_k/))
! end loop over lattice optimisation
end do

if (wproc1) then
 write(fnum1,'("Structual optimization is done!")')
 call timestamp(fnum1)
 call flushifc(fnum1)
 close(fnum1)
 close(71); close(72); close(73); close(74); close(75)
endif
!if (latvopt.ne.0) then
!  close(76); close(77)
!end if

! ground-state should be run again after lattice vector optimisation
!if (latvopt.ne.0) call gndstate
return
end subroutine

