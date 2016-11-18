
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: updatpos
! !INTERFACE:
subroutine updatpos
! !USES:
use modmain
! !DESCRIPTION:
!   Updates the current atomic positions according to the force on each atom. If
!   ${\bf r}_{ij}^m$ is the position and ${\bf F}_{ij}^m$ is the force acting on
!   it for atom $j$ of species $i$ and after time step $m$, then the new
!   position is calculated by
!   $$ {\bf r}_{ij}^{m+1}={\bf r}_{ij}^m+\tau_{ij}^m\left({\bf F}_{ij}^m
!    +{\bf F}_{ij}^{m-1}\right), $$
!   where $\tau_{ij}^m$ is a parameter governing the size of the displacement.
!   If ${\bf F}_{ij}^m\cdot{\bf F}_{ij}^{m-1}>0$ then $\tau_{ij}^m$ is
!   increased, otherwise it is decreased.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,ispn,is,ia,ias
integer n
real(8) t1
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the dot-product between the current and previous total force
    t1=dot_product(forcetot(:,ias),forcetp(:,ias))
! if the force is in the same direction then increase step size parameter
    if (t1.gt.0.d0) then
      tauatm(ias)=tauatm(ias)+tau0atm
    else
      tauatm(ias)=tau0atm
    end if
! check for negative mass
!    if (spmass(is).gt.0.d0) then
      atposc(:,ia,is)=atposc(:,ia,is)+tauatm(ias)*(forcetot(:,ias) &
       +forcetp(:,ias))
!    end if
  end do
end do

! make sure the atomic positions are the same in all the cpu
n=3*maxatoms*maxspecies
call mpi_grid_bcast(atposc(1,1,1),n)

do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the lattice coordinates of the atomic positions
    call r3mv(ainv,atposc(:,ia,is),atposl(:,ia,is))
  end do
end do
! check for overlapping muffin-tins
!call checkmt
return
end subroutine
!EOC
