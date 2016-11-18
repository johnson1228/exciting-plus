
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: charge
! !INTERFACE:
subroutine charge_val
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Computes the muffin-tin, interstitial and total charges by integrating the
!   density.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,ir
real(8) sum,t1
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
! find the muffin-tin charges
chgmttot_val=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      fr(ir)=rhomt_val(1,ir,ias)*spr(ir,is)**2
    end do
    call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
    chgmt_val(ias)=fourpi*y00*gr(nrmt(is))
    chgmttot_val=chgmttot_val+chgmt_val(ias)
  end do
end do
! find the interstitial charge
sum=0.d0
do ir=1,ngrtot
  sum=sum+rhoir_val(ir)*cfunir(ir)
end do
chgir_val=sum*omega/dble(ngrtot)
! total calculated charge
chgcalc_val=chgmttot_val+chgir_val
if (.not. pt_core) then
 t1=chgval/chgcalc_val
else
 t1=(chgval+ptchgcr)/chgcalc_val
endif
if (abs(t1-1.d0).gt.epschg) then
  write(*,*)
  write(*,'("Warning(charge): total valence charge density incorrect for s.c. &
   &loop ",I5)') iscl
  write(*,'(" Calculated : ",G18.10)') chgcalc_val
  write(*,'(" Required   : ",G18.10)') chgval
end if
! write calculated total charge to test file
!call writetest(400,'calculated total valence charge',tol=1.d-6,rv=chgval)
return
end subroutine
!EOC
