subroutine cal_sigma(iq,iqloc,nwloc,iter_,bndrg,gf_tau,sigc_aux,sigx_aux)
use modmain
use mod_linresp
use mod_addons_q
use mod_nrkp
use mod_expigqr
implicit none
!
integer,intent(in) :: iq
integer,intent(in) :: iqloc
integer,intent(in) :: nwloc
integer,intent(in) :: iter_
integer,intent(in) :: bndrg(2,nspinor)
real(8),intent(in) :: gf_tau(nstsv,lr_nw,nkptnrloc)
real(8),intent(inout) :: sigc_aux(nbnd,nspinor,lr_nw,nkptnr,nkpt+nvq0-1)
real(8),intent(inout) :: sigx_aux(nbnd,nspinor,nkptnr,nkpt+nvq0-1)
!
integer :: ig1,ig2,iw,i,ikloc,ist1,ist2,ik,j,n
integer :: iwloc, ibnd, fbnd, isp1
real(8),allocatable :: gfkq_tau(:,:,:)
complex(8),allocatable :: svq(:,:,:)
complex(8),allocatable :: ame(:,:)
complex(8) :: zt1,zt2
character*100 :: fname
logical :: exst
!
allocate(gfkq_tau(nstsv,lr_nw,nkptnrloc))
allocate(ame(ngq(iq),nstsv*nstsv))
allocate(svq(ngq(iq),ngq(iq),nwloc))
!
gfkq_tau=0.d0
svq=zzero

! restore q->0 matrix elements
if (vq_gamma(iq)) then
 do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do ig1=1,ngq(iq)
   if (igqig(ig1,iq).eq.1) then
     amegqblh(:,ig1,ikloc)=zzero
     do i=1,nmegqblh(ikloc)
       ist1=bmegqblh(1,i,ikloc)
       ist2=bmegqblh(2,i,ikloc)
       if (ist1.eq.ist2) then
         amegqblh(i,ig1,ikloc)=zone
       endif
     enddo
   endif
  enddo
 enddo
endif

! get G(k-q) from other processors
call get_gf0kq(3,gf_tau,gfkq_tau)

!Screened Coulomb potential W(tau)
if (gw_mode.eq.1.or.(gw_mode.eq.0.and.iter_.eq.1)) then
 ! solve W
 call psolve_wtau(iq,nwloc,gf_tau,svq)
 !
 if (gw_mode.eq.0.and.scgwni.gt.1.and.mpi_grid_root((/dim_k/))) then
  if (mpi_grid_root()) then
    write(151,'("writing svq_mat_q",I4.4)') iq
    write(151,'(" ")')
    call flushifc(151)
  endif

  write(fname,'("Temp_files/svq_mat_q",I4.4,"_ipq",I4.4)') iq, &
       & mpi_grid_dim_pos(dim_q)
  open(171,file=trim(adjustl(fname)),action='write',form='unformatted', &
       & status='replace')
  write(171) svq
  close(171)
 endif
elseif (gw_mode.eq.0.and.iter_.gt.1) then
 write(fname,'("Temp_files/svq_mat_q",I4.4,"_ipq",I4.4)') iq, mpi_grid_dim_pos(dim_q)
 inquire(file=trim(adjustl(fname)),exist=exst)
 
 if (exst) then
   open(171,file=trim(adjustl(fname)),action='read',form='unformatted', &
        & status='old')
   read(171) svq
   close(171)
 else
   if (mpi_grid_root()) then
    write(151,'("cannot open file svq_mat_q",I4.4,"_ipq",I4.4)') iq, &
           & mpi_grid_dim_pos(dim_q)
    write(151,'("Error!")')
    call flushifc(151)
    call pstop
   endif
 endif
endif

! wait for all processors
call mpi_grid_barrier()

if (mpi_grid_root()) then
 write(151,*) "Now calculate the self-energy elements..."
 call timestamp(151)
 call flushifc(151)
endif

!calculate the self-energy
do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 !note that amegqblh M(k-q), without flipping jj'
 do ig1=1,ngq(iq)
   ame(ig1,:)=amegqblh(:,ig1,ikloc)
 enddo

 do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_q,loc=iwloc)
  do isp1=1,nspinor
   ibnd=bndrg(1,isp1)
   fbnd=bndrg(2,isp1)
   do i=1,namegqblh(ikloc)
     n=bamegqblh(2,i,ikloc)
     if ((n.gt.fbnd).or.(n.lt.ibnd)) cycle  !nk  
     j=bamegqblh(1,i,ikloc)

     zt2=zzero
     do ig2=1,ngq(iq)
      zt1=zzero
      do ig1=1,ngq(iq)   
        zt1=zt1+dconjg(ame(ig1,i))*svq(ig1,ig2,iwloc)
      enddo
      zt2=zt2+zt1*ame(ig2,i)
     enddo
  
     sigc_aux(n-ibnd+1,isp1,iw,ik,iqloc)=sigc_aux(n-ibnd+1,isp1,iw,ik,iqloc)+&
                            & gfkq_tau(j,iw,ikloc)*dreal(zt2)

   enddo !i
  enddo !isp1
 enddo !iwloc

 !evaluate valence contribution of exxnk in k-space
 if (exxtype.eq.1) then
  do isp1=1,nspinor
   ibnd=bndrg(1,isp1)
   fbnd=bndrg(2,isp1)
   do i=1,namegqblh(ikloc)
     n=bamegqblh(2,i,ikloc)
     if ((n.gt.fbnd).or.(n.lt.ibnd)) cycle  !nk  
     j=bamegqblh(1,i,ikloc)

     do ig1=1,ngq(iq)
       sigx_aux(n-ibnd+1,isp1,ik,iqloc)=sigx_aux(n-ibnd+1,isp1,ik,iqloc)+&
                  & ame(ig1,i)*dconjg(ame(ig1,i))*gfkq_tau(j,lr_nw,ikloc)*&
                  & wtvhgq(ig1,iq)
     enddo !ig1
   enddo !i
  enddo !isp1
 endif

enddo !ikloc

deallocate(ame)
deallocate(gfkq_tau,svq)
return
end subroutine
