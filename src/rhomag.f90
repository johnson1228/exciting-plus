subroutine rhomag
use modmain
implicit none
integer ikloc,idm
call timer_start(t_rho_mag_tot)
! set the charge density and magnetisation to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
do ikloc=1,nkptloc
  call rhomagk(ikloc,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
end do
call mpi_grid_reduce(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,&
  dims=(/dim_k,dim2/))
call mpi_grid_reduce(rhoir(1),ngrtot,dims=(/dim_k,dim2/))
if (spinpol) then
  call mpi_grid_reduce(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag,&
    dims=(/dim_k,dim2/))
  call mpi_grid_reduce(magir(1,1),ngrtot*ndmag,dims=(/dim_k,dim2/))
endif
if (mpi_grid_root(dims=(/dim_k,dim2/))) then
! convert muffin-tin density/magnetisation to spherical harmonics
  call rhomagsh
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(lradstp,magmt,magir)
! convert the density from a coarse to a fine radial mesh
  call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
  do idm=1,ndmag
    call rfmtctof(magmt(:,:,:,idm))
  end do
! save valence charge density----new
  if (rho_val.or.pt_core) then
   rhomt_val(:,:,:)=rhomt(:,:,:)
   rhoir_val(:)=rhoir(:)
  endif
! add the core density to the total density
  call addrhocr
  if (pt_core) call addptrhocr
endif
call mpi_grid_bcast(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,&
  dims=(/dim_k,dim2/))
call mpi_grid_bcast(rhoir(1),ngrtot,dims=(/dim_k,dim2/))
! bcast for valence charge density as well-new
if (rho_val.or.pt_core) then
 call mpi_grid_bcast(rhomt_val(1,1,1),lmmaxvr*nrmtmax*natmtot,&
   dims=(/dim_k,dim2/))
 call mpi_grid_bcast(rhoir_val(1),ngrtot,dims=(/dim_k,dim2/))
endif

if (spinpol) then
  call mpi_grid_bcast(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag,&
    dims=(/dim_k,dim2/))
  call mpi_grid_bcast(magir(1,1),ngrtot*ndmag,dims=(/dim_k,dim2/))
endif
! calculate the charges
call charge
! calculate the valence charges
if (rho_val.or.pt_core) then
  call charge_val
endif
! calculate the moments
if (spinpol) call moment
! normalise the density
call rhonorm
!normalize the valence density
if (rho_val.or.pt_core) then
  call rhovalnorm
endif

call timer_stop(t_rho_mag_tot)
return
end
