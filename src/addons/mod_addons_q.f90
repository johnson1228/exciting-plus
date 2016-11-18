module mod_addons_q
!-------------------------!
!     q and G+q vectors   !
!-------------------------!
! number of q-vectors
integer nvq
! q-vectors in k-mesh coordinates
integer, allocatable :: vqm(:,:)
! non-reduced (to first BZ) q-vectors in lattice coordinates
real(8), allocatable :: vqlnr(:,:)
! non-reduced (to first BZ) q-vectors in Cartesian coordinates
real(8), allocatable :: vqcnr(:,:)
! reduced to first BZ q-vectors in lattice coordinates
real(8), allocatable :: vql(:,:)
! reduced to first BZ q-vectors in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! weights associated with the integral of 1/q^2
real(8), allocatable :: wiq2(:)
! cutoff for |G+q|
real(8) gqmax
data gqmax/2.d0/
integer gqsh
data gqsh/2/
logical tgqsh
data tgqsh/.false./
! maximum number of G+q vectors
integer ngqmax
! global index of Gq-vector, which brigs q-vector to first BZ
integer, allocatable :: ig0q(:)
! index of Gq vector in the range[1,ngq(iq)]
integer iig0q
! number of G+q vectors
integer, allocatable :: ngq(:)
! G+q vectors in Cartesian coordinates
real(8), allocatable ::  vgqc(:,:,:)
! length of |G+q| vectors
real(8), allocatable :: gq(:,:)
! global index of G of G+q vector  
integer, allocatable :: igqig(:,:)
! 4*Pi/|G+q|^2 (Fourier transform of Hartree potential)
real(8), allocatable :: vhgq(:,:)
! weighted G+q components of Hartee potential
real(8), allocatable :: wtvhgq(:,:)

integer, allocatable :: ngqsh(:)
integer, allocatable :: gqshidx(:,:)
real(8), allocatable :: gqshlen(:,:)

real(8), allocatable :: tpgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: ylmgq(:,:)

! true if Gamma point included in the integration over Brillouin zone
logical tq0bz
data tq0bz/.true./
integer nvq0
data nvq0/0/
real(8) vq0c(3)
real(8) q0wt

contains 

logical function vq_gamma(iq)
implicit none
integer, intent(in) :: iq
!
vq_gamma=.false.
if (sum(vql(:,iq)**2).lt.1d-10) vq_gamma=.true.
return
end function

integer function getngvecme()
use modmain
implicit none
integer ngsh
integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)
allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))
call getgshells(ngsh,igishell,ishellng)
getngvecme=ishellng(gqsh,2)
deallocate(igishell)
deallocate(ishellng)
return
end function

subroutine init_q_mesh(nvq0_)
use modmain
implicit none
integer, intent(in) :: nvq0_
integer j,i1,i2,i3
!
nvq0=nvq0_
if (allocated(vqm)) deallocate(vqm)
nvq=nkptnr-1+nvq0
allocate(vqm(3,nvq))
vqm=0
j=nvq0
do i1=0,ngridk(1)-1
  do i2=0,ngridk(2)-1
    do i3=0,ngridk(3)-1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        j=j+1
        vqm(:,j)=(/i1,i2,i3/)
      endif
    enddo
  enddo
enddo
return
end subroutine

! generate q-vector related data
subroutine genvq
use modmain
implicit none
integer iq,ig,v1(3),i,i1,i2,i3
logical f
!
if (allocated(vqlnr)) deallocate(vqlnr)
allocate(vqlnr(3,nvq))
if (allocated(vqcnr)) deallocate(vqcnr)
allocate(vqcnr(3,nvq))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nvq))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nvq))
if (allocated(ig0q)) deallocate(ig0q)
allocate(ig0q(nvq))

do iq=1,nvq
! non-reduced q-vector in lattice coordinates
  vqlnr(:,iq)=dble(vqm(:,iq))/dble(ngridk(:))
! non-reduced q-vector in Cartesian coordinates 
  call r3mv(bvec,vqlnr(:,iq),vqcnr(:,iq))
! find Gq vector which reduces q to first BZ
  f=.false.
  do ig=1,ngvec
! q-G vector in k-mesh coordinates 
    v1(:)=vqm(:,iq)-ngridk(:)*ivg(:,ig)
    if (v1(1).ge.0.and.v1(1).lt.ngridk(1).and.&
        v1(2).ge.0.and.v1(2).lt.ngridk(2).and.&
        v1(3).ge.0.and.v1(3).lt.ngridk(3).and..not.f) then
      ig0q(iq)=ig
      f=.true.
    endif
  enddo !ig
  if (.not.f) then
    write(*,*)
    write(*,'("Error(genvq): Gq-vector is not found because &
      &q-vector is too large")')
    call pstop
  endif
! reduced q-vector in lattice coordinates
  vql(:,iq)=vqlnr(:,iq)-dble(ivg(:,ig0q(iq)))
! reduced q-vector in Cartesian coordinates
  call r3mv(bvec,vql(:,iq),vqc(:,iq))
enddo
! init q=0 points
if (nvq0.eq.1) then
  if (.not.vq_gamma(1)) then
    write(*,'("[genvq] : firt q-pont in the list is not a Gamma")')
    call pstop
  endif
  vqc(:,1)=vq0c(:)
else if (nvq0.eq.8) then
  i=0
  do i1=0,1
    do i2=0,1
      do i3=0,1
        i=i+1
! take small vector on the diagonal of the eight nearest corner-sharing micro-cells
        vqc(:,i)=0.01d0*((i1-0.5d0)*bvec(:,1)/ngridk(1)+(i2-0.5d0)*bvec(:,2)/ngridk(2)+&
          &(i3-0.5d0)*bvec(:,3)/ngridk(3))
      enddo
    enddo
  enddo
else if (nvq0.ne.0) then
  write(*,'("[genvq] : nvq0= ",I4," is not implemented")')nvq0
  call pstop
endif
return
end subroutine

! generate G+q-vector related data
subroutine genvgq
use modmain
!use mod_expigqr
implicit none
real(8) t2,gqmax2
integer iq,ig,i
real(8) v2(3)
logical tautogqmax,tfound
!
! find |G+q| cutoff
tautogqmax=.true.
if (tautogqmax) then
  do iq=1,nvq
    t2=sqrt(sum(vqcnr(:,iq)**2))
    gqmax=max(gqmax,t2*1.01)
  enddo
endif
! find maximum number of G+q vectors
ngqmax=0
gqmax2=gqmax**2
do iq=1,nvq
  i=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if (t2.le.gqmax2) i=i+1
  enddo
  ngqmax=max(ngqmax,i)
enddo
if (tgqsh) then
  ngqmax=getngvecme()
endif
! generate G+q vectors
if (allocated(ngq)) deallocate(ngq)
allocate(ngq(nvq))
if (allocated(gq)) deallocate(gq)
allocate(gq(ngqmax,nvq))
if (allocated(vgqc)) deallocate(vgqc)
allocate(vgqc(3,ngqmax,nvq))
if (allocated(igqig)) deallocate(igqig)
allocate(igqig(ngqmax,nvq))
if (allocated(vhgq)) deallocate(vhgq)
allocate(vhgq(ngqmax,nvq))
vhgq(:,:)=0.d0
do iq=1,nvq
  ngq(iq)=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if ((.not.tgqsh.and.t2.le.gqmax2).or.(tgqsh.and.ig.le.ngqmax)) then
      ngq(iq)=ngq(iq)+1
      gq(ngq(iq),iq)=sqrt(t2)
      vgqc(:,ngq(iq),iq)=v2(:)
      igqig(ngq(iq),iq)=ig
      vhgq(ngq(iq),iq)=fourpi/t2
    endif
  enddo !ig
enddo
! compute weights for 1/|q|^2 integral
if (nvq0.ne.0) call genwtvhgq
! find |G+q| shells
if (allocated(ngqsh)) deallocate(ngqsh)
allocate(ngqsh(nvq))
ngqsh=0
if (allocated(gqshidx)) deallocate(gqshidx)
allocate(gqshidx(ngqmax,nvq))
gqshidx=0
if (allocated(gqshlen)) deallocate(gqshlen)
allocate(gqshlen(ngqmax,nvq))
gqshlen=0.d0
do iq=1,nvq
  do ig=1,ngq(iq)
    tfound=.false.
    do i=1,ngqsh(iq)
      if (abs(gqshlen(i,iq)-gq(ig,iq)).lt.1d-10) then
        tfound=.true.
        gqshidx(ig,iq)=i
        exit
      endif
    enddo
    if (.not.tfound) then
      ngqsh(iq)=ngqsh(iq)+1
      gqshlen(ngqsh(iq),iq)=gq(ig,iq)
      gqshidx(ig,iq)=ngqsh(iq)
    endif
  enddo
enddo
return
end subroutine

subroutine genwtvhgq
use modmain
implicit none
integer i,n,i1,i2,i3,ig,iq
real(8) x(2),y(2)
real(8) vq(3),p1,p2,p3,vol
!
if (allocated(wtvhgq)) deallocate(wtvhgq)
allocate(wtvhgq(ngqmax,nvq))
!
n=100
y(:)=0.d0
do i=1,2
  do i1=-n,n-1
    p1=(0.5d0+i1)/dble(2*n)
    do i2=-n,n-1
      p2=(0.5d0+i2)/dble(2*n)
      do i3=-n,n-1
        p3=(0.5d0+i3)/dble(2*n)
        vq(:)=(p1*bvec(:,1)/ngridk(1)+p2*bvec(:,2)/ngridk(2)+p3*bvec(:,3)/ngridk(3))
        y(i)=y(i)+1.d0/dot_product(vq,vq)
      enddo
    enddo
  enddo
  vol=(twopi**3)/omega/nkptnr/((2*n)**3)
  y(i)=y(i)*vol
  x(i)=vol**(1/3.d0)
  n=n+40
enddo
! numerical derivative to get y(x=0), where x is the effective size of the micro-cell
q0wt=y(2)-x(2)*(y(1)-y(2))/(x(1)-x(2))

wtvhgq(:,:)=vhgq(:,:)
do iq=1,nvq
  if (vq_gamma(iq)) then
    do ig=1,ngq(iq)
      if (igqig(ig,iq).eq.1) then
        wtvhgq(ig,iq)=fourpi*q0wt*nkptnr*omega/(twopi**3)
      endif
    enddo
    wtvhgq(:,iq)=wtvhgq(:,iq)/dble(nvq0)
  endif
enddo

return
end subroutine

subroutine genwiq2_new
use modmain
implicit none
! local variables
integer, parameter :: np=5
integer, parameter :: ns0=10, nss=20
integer nsymqpt,isym,nqpt,reduceq
integer symqpt(3,3,48)
logical lsym(48)
integer, allocatable :: iqmaptmp(:,:,:),ivqtmp(:,:)
real(8), allocatable :: vqltmp(:,:),vqctmp(:,:),wqpttmp(:)
!real(8) :: vtmp(3)
!real(8) :: wtq
integer ns,iq,i1,i2,i3,i,ip
real(8) d(3),dv,sum,t1,t2
real(8) v1(3),v2(3),v3(3)
real(8) xa(np),ya(np),c(np)
real(8) boxl(3,4)
! external functions
real(8) polynom
external polynom

! allocate global wiq2 array
if (allocated(wiq2)) deallocate(wiq2)
allocate(wiq2(nkpt))

! first construct the q-point set according to init2
reduceq=1   !default value in readinput.f90

if (reduceq.eq.0) then
  nsymqpt=1
  symqpt(:,:,1)=symlat(:,:,1)
else
  lsym(:)=.false.
  do isym=1,nsymcrys
    lsym(lsplsymc(isym))=.true.
  end do
  nsymqpt=0
  do isym=1,nsymlat
    if (lsym(isym)) then
      nsymqpt=nsymqpt+1
      symqpt(:,:,nsymqpt)=symlat(:,:,isym)
    end if
  end do
end if

reduceq=0 ! similar to Hartree-Fock

allocate(ivqtmp(3,ngridk(1)*ngridk(2)*ngridk(3)))
allocate(vqltmp(3,ngridk(1)*ngridk(2)*ngridk(3)))
allocate(vqctmp(3,ngridk(1)*ngridk(2)*ngridk(3)))
allocate(wqpttmp(ngridk(1)*ngridk(2)*ngridk(3)))
allocate(iqmaptmp(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1))

! setup the q-point box (offset should always be zero)
boxl(:,:)=0.d0
boxl(1,2)=1.d0; boxl(2,3)=1.d0; boxl(3,4)=1.d0;

call genppts(.true.,nsymqpt,symqpt,ngridk,epslat,bvec,boxl,nqpt,&
         iqmaptmp,ivqtmp,vqltmp,vqctmp,wqpttmp)

if (mpi_grid_root()) then

  do iq=1,ngridk(1)*ngridk(2)*ngridk(3)
    write(*,*) vqltmp(:,iq),iq
    write(*,*) vqctmp(:,iq),iq
  enddo

endif
!----end of q-point setting

! if system is a molecule wiq2 should be zero
if (molecule) then
  wiq2(:)=0
  return
end if

! begin loop over q-points, note that the vectors vqc are assumed to be in the
! first Brillouin zone
do iq=1,nqpt
! loop over different subdivisions
  ns=ns0
  do ip=1,np
! subdivision vectors in lattice coordinates
    do i=1,3
      d(i)=1.d0/(dble(ngridk(i)*2*ns))
    end do
! smallest volume element
    dv=((twopi**3)/omega)*d(1)*d(2)*d(3)
! compute the integral of 1/q^2
    sum=0.d0
    do i1=-ns,ns-1
      t1=dble(i1)*d(1)
!      v1(:)=vqc(:,iq)+t1*bvec(:,1)
       v1(:)=vqctmp(:,iq)+t1*bvec(:,1)
      do i2=-ns,ns-1
        t1=dble(i2)*d(2)
        v2(:)=v1(:)+t1*bvec(:,2)
        do i3=-ns,ns-1
          t1=dble(i3)*d(3)
          v3(:)=v2(:)+t1*bvec(:,3)
          t2=v3(1)**2+v3(2)**2+v3(3)**2
          if (t2.gt.1.d-14) then
            sum=sum+1.d0/t2
          end if
        end do
      end do
    end do
    sum=sum*dv
    xa(ip)=dv**(1.d0/3.d0)
    ya(ip)=sum
! increment number of subdivisions
    ns=ns+nss
  end do
! extrapolate the volume element to zero with a polynomial
   wiq2(iq)=polynom(0,np,xa,ya,c,0.d0)

  if (mpi_grid_root()) then
    write(*,*) "wiq2,iq:",wiq2(iq)
!    write(*,*) "vqc(:),iq:",vqc(:,iq),iq
    write(*,*) "vqltmp(:),iq:",vqltmp(:,iq),iq
  endif
end do
 if (mpi_grid_root()) then
  write(*,*) "nvq,molecule,omega:",nvq,molecule,omega
 endif

deallocate(ivqtmp,vqltmp,vqctmp,wqpttmp,iqmaptmp)

return
end subroutine

subroutine kq_map(iqrmap,qqnrmap,akmap,kmap,kknrmap)
use modmain
!
implicit none
!
integer,intent(out) :: iqrmap(2,nkpt+nvq0-1)
integer,intent(out) :: qqnrmap(nvq,2,nkpt+nvq0-1)
integer,intent(out) :: akmap(48,nkptnr)
integer,intent(out) :: kmap(2,nkpt)
integer,intent(out) :: kknrmap(nkptnr,2,nkpt)
! local variable
integer :: ik,iq,isym,ik1,iv(3),lspl
real(8) :: v1(3),vak(3),s(3,3),t1

! generate mapping from q_IBZ to q_BZ
! iqrmap(1,iqibz): index of qibz in the non-reduced BZ mesh
! iqrmap(2,iqibz): # of qbz related to qibz
iqrmap=0
! qqnrmap: set of q_{BZ} points generated by q_{IBZ} through 
!         R*q_{IBZ}+GR = q_{BZ} where {R} are the symmetry operations. 
!         Note that q_{BZ}s are in non-reduced q-mesh. 
! qqnrmap(iqbz,1,iqibz): index of each q_{BZ} point in the
!                        non-reduced q-mesh, related to q_{IBZ};
! qqnrmap(iqbz,2,iqibz): index of the symmetry operation "R" that yields 
!                        q_{BZ} = R*q_{IBZ}+GR
qqnrmap(:,:,:)=0
! akmap(isym,ik): mapping non-reduced k1 to k2=R^{-1}*k1 in non-reduced k-mesh
akmap(:,:)=0
! similar iqrmap and qqnrmap for k_{IBZ} <-> k_{BZ}
kmap(:,:)=0
kknrmap(:,:,:)=0
!
do iq=1,nvq
 if (vq_gamma(iq)) then
  iqrmap(1,iq)=iq
  iqrmap(2,iq)=1
  qqnrmap(iqrmap(2,iq),1,iq)=iq
  qqnrmap(iqrmap(2,iq),2,iq)=1
 else
  call findkpt(vql(:,iq),isym,ik1)
  iqrmap(2,ik1+nvq0-1)=iqrmap(2,ik1+nvq0-1)+1
  qqnrmap(iqrmap(2,ik1+nvq0-1),1,ik1+nvq0-1)=iq
  qqnrmap(iqrmap(2,ik1+nvq0-1),2,ik1+nvq0-1)=lsplsymc(isym)
! check
  if (lsplsymc(isym).eq.1) iqrmap(1,ik1+nvq0-1)=iq
 endif
enddo

if (mpi_grid_root()) then
 do iq=1,nkpt+nvq0-1
  write(*,*) "iq,iqrmap,nq_iq:",iq,iqrmap(1,iq),iqrmap(2,iq)
  write(*,*) "vql:",vql(:,iqrmap(1,iq))
 enddo
endif

do ik=1,nkptnr
 v1(:)=vklnr(:,ik)    ! k2
 do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,v1,vak)   ! R*k2
  call r3frac(epslat,vak,iv)   ! (R*k2)_{BZ}
  do ik1=1,nkptnr
   v1(:)=vklnr(:,ik1)  ! k1
   t1=abs(vak(1)-v1(1))+abs(vak(2)-v1(2))+abs(vak(3)-v1(3))
   if (t1.lt.epslat) then
    akmap(lspl,ik1)=ik   ! index of k2=R^{-1}k1
    exit
   endif
  enddo
 enddo
enddo

! find out kIBZ <-> kBZ
do ik=1,nkptnr
 call findkpt(vklnr(:,ik),isym,ik1)
 kmap(2,ik1)=kmap(2,ik1)+1
 kknrmap(kmap(2,ik1),1,ik1)=ik
 kknrmap(kmap(2,ik1),2,ik1)=lsplsymc(isym)
 if (lsplsymc(isym).eq.1) kmap(1,ik1)=ik
enddo

!if (mpi_grid_root()) then
! do ik=1,nkptnr
!  write(*,*) "ik:",ik
!  do isym=1,nsymcrys
!    write(*,*) "isym:",isym,ik,akmap(isym,ik)
!  enddo
! enddo
!endif

!if (mpi_grid_root()) then
! do isym=1,nsymcrys
!  write(*,*) "isym,lspl(isym):",isym,lsplsymc(isym)
! enddo
!endif

return
end subroutine

end module
