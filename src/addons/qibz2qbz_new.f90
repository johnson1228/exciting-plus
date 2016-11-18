subroutine qibz2qbz_new(nqk,iqkrmap,nqkmax)
use modmain
use mod_addons_q
use mod_linresp
!
implicit none
!
integer,intent(out) :: nqk(nkptnrloc)
integer,intent(out) :: iqkrmap(2,nvq,nkptnrloc)
integer,intent(out) :: nqkmax
integer :: symk(3,3,nsymcrys,nkptnrloc)
integer :: nsymk(nkptnrloc)
integer :: ik,iq,iq1,isym,ik1,lspl,iv(3),ikloc
real(8) :: t1,v1(3),s(3,3)
integer :: iqkmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
integer :: ivqk(3,ngridk(1)*ngridk(2)*ngridk(3))
integer :: iktmp(nvq,nkptnrloc)
real(8) :: vqkc(3,ngridk(1)*ngridk(2)*ngridk(3))
real(8) :: wppt(ngridk(1)*ngridk(2)*ngridk(3),nkptnrloc)
real(8) :: vqkl(3,ngridk(1)*ngridk(2)*ngridk(3),nkptnrloc)
real(8) :: boxl(3,4)
character*100 :: fname
!
! generate mapping from q_IBZ to q_BZ, for a given k-point (small group of k)
! iqkrmap(1,iqibz,ik): include in the calculation or not (1/0)
! iqkrmap(2,iqibz,ik): # of qbz related to qibz
nqk=0
iqkrmap=0
symk=0
nsymk=0
vqkl=0.d0
wppt=0.d0
nqkmax=0

do iq=1,nvq0
 iqkrmap(:,iq,:)=1
enddo

! setup the default k-point box
boxl(:,1)=vkloff(:)/dble(ngridk(:))
boxl(:,2)=boxl(:,1); boxl(:,3)=boxl(:,1); boxl(:,4)=boxl(:,1)
boxl(1,2)=boxl(1,2)+1.d0
boxl(2,3)=boxl(2,3)+1.d0
boxl(3,4)=boxl(3,4)+1.d0

! find out the small group of a given k
do ikloc=1,nkptnrloc
 ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
 do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,vklnr(:,ik),v1)
  call r3frac(epslat,v1,iv)
  t1=abs(v1(1)-vklnr(1,ik))+abs(v1(2)-vklnr(2,ik))+abs(v1(3)-vklnr(3,ik))
  if (t1.lt.epslat) then
   nsymk(ikloc)=nsymk(ikloc)+1
   symk(:,:,nsymk(ikloc),ikloc)=s(:,:)
  endif
 enddo

 call genppts(.false.,nsymk(ikloc),symk(:,:,:,ikloc),ngridk,epslat,&
            &bvec,boxl,nqk(ikloc),iqkmap,ivqk,vqkl(:,:,ikloc),vqkc,&
            &wppt(:,ikloc))

 do iq1=1,nvq
  do iq=1,nqk(ikloc)
   if ((abs(vqkl(1,iq,ikloc)).lt.epslat).and.&
     &(abs(vqkl(2,iq,ikloc)).lt.epslat).and.&
     &(abs(vqkl(3,iq,ikloc)).lt.epslat)) cycle
   t1=abs(vqkl(1,iq,ikloc)-vql(1,iq1))+abs(vqkl(2,iq,ikloc)-vql(2,iq1))&
      &+abs(vqkl(3,iq,ikloc)-vql(3,iq1))
   if (t1.lt.epslat) then
     iqkrmap(1,iq1,ikloc)=1
     iqkrmap(2,iq1,ikloc)=int(wppt(iq,ikloc)*ngridk(1)*ngridk(2)*ngridk(3))
     exit
   endif
  enddo
 enddo

 do iq1=1,nvq
  iktmp(iq1,ikloc)=iqkrmap(1,iq1,ikloc)
 enddo

 call mpi_grid_reduce(iktmp(1,1),nvq*nkptnrloc,dims=(/dim_k/),all=.true.,&
                      &op=op_max)

 do iq1=nvq0+1,nvq
  iqkrmap(1,iq1,ikloc)=iktmp(iq1,ikloc)
 enddo

! test
 if (mpi_grid_root((/dim_q/))) then
   write(fname,'("iqkrmap_k",I4.4)') ik
   open(164,file=fname,action='write',status='replace')
   write(164,'("ik,nqk:",2(I4,1X))') ik,nqk(ikloc)+nvq0-1
   do iq=1,nvq
    write(164,*) "iq:",iq
    write(164,*) iqkrmap(1,iq,ikloc),iqkrmap(2,iq,ikloc)
    write(164,*) "vql:",vql(:,iq)
   enddo
   write(164,'(" ")')

   do iq=1,nqk(ikloc)
    write(164,*) "iq,vqkl:",iq,vqkl(:,iq,ikloc)
   enddo

   close(164)
 endif

enddo !ikloc

return
end subroutine
