subroutine update_sigxc(gf0_tau,iq,iqloc,iter_,bndrg,sigc_aux)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
real(8),intent(in) :: gf0_tau(nstsv,lr_nw,nkptnrloc)
integer,intent(in) :: iq
integer,intent(in) :: iqloc
integer,intent(in) :: iter_
integer,intent(in) :: bndrg(2,nspinor)
real(8),intent(inout) :: sigc_aux(nbnd,nspinor,lr_nw,nkptnr,nkpt+nvq0-1)
!
complex(8),allocatable :: pol(:,:,:)
integer :: nwloc
!
!parallelize along the frequency mesh along dim_q
nwloc=mpi_grid_map(lr_nw,dim_q)
allocate(pol(ngq(iq),ngq(iq),nwloc))
pol=zzero
! gw_mode == 0    GW0
!         == 1    full GW

if (gw_mode.eq.1.or.(gw_mode.eq.0.and.iter_.eq.1)) then
 call cal_pol(iq,gf0_tau,pol)
endif

! calculate self-energy at iq
call cal_sigma(iq,iqloc,nwloc,iter_,bndrg,gf0_tau,pol,sigc_aux)

deallocate(pol)
return
end subroutine
