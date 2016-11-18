
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: addrhocr
! !INTERFACE:
subroutine addptrhocr
! !USES:
use modmain
! !DESCRIPTION:
!   Adds the core density to the muffin-tin and interstitial densities. A
!   uniform background density is added in the interstitial region to take into
!   account leakage of core charge from the muffin-tin spheres.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
integer nr,ir,ispn
real(8) sum1,sum2
real(8) v(3),t1,t2
! automatic arrays
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
real(8) ptchgcrlk
sum1=0.d0
sum2=0.d0
do is=1,nspecies
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! loop over spin channels
    do ispn=1,nspncr
      do ir=1,nr
! add the core density to the muffin-tin density
        rhomt_val(1,ir,ias)=rhomt_val(1,ir,ias)+ptrhocr(ir,ias,ispn)/y00
        fr(ir)=ptrhocr(ir,ias,ispn)*spr(ir,is)**2
      end do
! compute the core charge inside the muffin-tins
      call fderiv(-1,nr,spr(:,is),fr,gr,cf)
      sum1=sum1+fourpi*gr(nr)
    end do
  end do
  sum2=sum2+dble(natoms(is))*(4.d0*pi/3.d0)*rmt(is)**3
end do
! add remaining core charge to interstitial density
ptchgcrlk=ptchgcr-sum1
t1=ptchgcrlk/(omega-sum2)
rhoir_val(:)=rhoir_val(:)+t1
return
end subroutine
!EOC

