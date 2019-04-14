program hexaTB
! This is hexaTB version 1.0 (April 2019). Written by Viktor Zolyomi.
!
! This program computes numerically the band structure of 2D hexagonal crystals.
! It solves as an example the bilayer graphene tight-binding model published in
! E. McCann and M. Koshino, Rep. Prog. Phys. 76, 056503 (2013)
! To use this program with any other hexagonal system, simply replace the functions
! provided in the fun_mod module. The code will work for any hexagonal 2D
! material and can easily be generalized to non-hexagonal or 3D crystals.
! Support for nontrivial overlap matrices not yet implemented in this release.
  use :: dtype_mod, only: DP ! use double precision
  use :: fun_mod ! import functions and parameters from separate module
  use :: omp_lib ! openmp parallelization is used, see below
  use :: settings_mod ! import parameter settings from separate module
  implicit none
  real (kind=dp), parameter :: PI = 4.0_dp*atan(1.0_dp) ! double precision pi
  real (kind=dp),allocatable :: bb(:,:) ! array of reciprocal lattice vectors
  real (kind=dp),allocatable :: kkline(:,:) ! array of quasimomentum vectors for GKMG path
  real (kind=dp),allocatable :: kkgrid(:,:,:) ! array of quasimomentum vectors for 2D grid; you may want to redefine this if you have a 3D crystal
  real (kind=dp) :: ktrack ! total length of distance covered by array of quasimomentum vectors along GKMG path
  real (kind=dp),allocatable :: energline(:,:) ! array of band energies for GKMG path
  real (kind=dp),allocatable :: energgrid(:,:,:) ! array of band energies for 2D grid; you may want to redefine this if you have a 3D crystal
  complex (kind=dp),allocatable :: work(:) ! complex work array for lapack
  real (kind=dp),allocatable :: rwork(:) ! real work array for lapack
  integer :: ierr ! integer for error tracking
  integer :: istart ! integer storing start time
  integer :: iend ! integer storing end time
  real (kind=dp),allocatable :: params(:) ! tight-binding parameters
  character(len=120) :: msg ! error message variable
  integer :: num_threads ! number of threads in openmp
  integer :: thread_counter ! index of thread in openmp
  complex (kind=dp),allocatable :: hh(:,:) ! Hamiltonian matrix
!  complex (kind=dp),allocatable :: ss(:,:) ! overlap matrix; neglected in this version
  integer :: ii, jj ! integers used in loops throughout
  character(len=20) :: cfmt1, cfmt2 ! export format specification strings
  istart = time()
  write(*,*) 'This is hexaTB.x version 1.0 (April 2019). Written by Viktor Zolyomi.'
  write(*,*)
  write(*,*) 'This program computes the eigenvalues of the TB Hamiltonian'
  write(*,*) 'for bilayer graphene, following the model published in'
  write(*,*) 'E. McCann and M. Koshino, Rep. Prog. Phys. 76, 056503 (2013)'
  write(*,*)
  write(*,*) 'Allocating memory'
  allocate(params(paramsize), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  call tb_params(params)
  allocate(hh(n_size,n_size), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
!  allocate(ss(n_size,n_size), stat = ierr, errmsg = msg)
!  if(ierr/=0)then
!    write(*,*) 'Allocation error:',msg
!  endif
  allocate(bb(n_dim,n_dim), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  allocate(kkline(numk,n_dim), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  allocate(kkgrid(numk,numk,n_dim), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  allocate(energline(numk,n_size), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  allocate(energgrid(numk,numk,n_size), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  allocate(work(lwork), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  allocate(rwork(lwork), stat = ierr, errmsg = msg)
  if(ierr/=0)then
    write(*,*) 'Allocation error:',msg
  endif
  write(*,*) 'Memory allocation complete'

! openmp testing:
  if(lomptest)then
!$omp parallel default(shared), private(thread_counter)
!$omp single
    num_threads = omp_get_num_threads()
!$omp end single
    thread_counter = omp_get_thread_num()
    print *, 'This is thread', thread_counter, 'of', num_threads
!$omp end parallel
  endif

! Here we define reciprocal lattice vectors. Make sure your model uses the
! correct choice of lattice vectors, or update these to match your model:
  bb(:,1)=[2.0_dp*pi/aa,2.0_dp*pi/(aa*sqrt(3.0_dp))]
  bb(:,2)=[2.0_dp*pi/aa,-2.0_dp*pi/(aa*sqrt(3.0_dp))]
! Two sets of k-points for which to calculate energies. The first is the
! Gamma-K-M-Gamma path, the second is a grid spanning the entire Brillouin zone.
  do ii=1,2*numk/5
    kkline(ii,:)=((ii-1)*(bb(:,1)+bb(:,2))/3.0_dp)/(2*numk/5) ! Gamma to K
  enddo
  do ii=2*numk/5+1,3*numk/5
    kkline(ii,:)=((bb(:,1)+bb(:,2))/3.0_dp)+&
&((ii-2*numk/5-1)*( bb(:,1)*0.5_dp - (bb(:,1)+bb(:,2))/3.0_dp ))/(3*numk/5-2*numk/5) ! K to M
  enddo
  do ii=3*numk/5+1,numk
    kkline(ii,:)=((bb(:,1))*0.5_dp)-((ii-3*numk/5-1)*(bb(:,1))*0.5_dp)/(numk-3*numk/5) ! M to Gamma
  enddo
  do ii=1,numk
    do jj=1,numk
      kkgrid(ii,jj,:)=((ii-numk/2)*bb(:,1))/numk+((jj-numk/2)*bb(:,2))/numk ! 2D grid spanning full BZ
    enddo
  enddo

  write(*,*) 'GKMG path k loop commencing'
  energline=0.0_dp
!$omp parallel default(none), shared(kkline,energline,params), private(ii,hh,work,rwork,ierr)
!$omp do
  do ii = 1, numk
    hh = Ham(aa,kkline(ii,:),params,n_size)
    work=0.0_dp
    rwork=0.0_dp
    call ZHEEV('N','U',n_size,hh,n_size,energline(ii,:),work,lwork,rwork,ierr)
! the following line can be used to determine the optimal dimension of the work array:
!    write(*,*) 'optimal lwork: ',work(1)
    if(ierr/=0)then
      write(*,*) 'k-point number',ii,'encountered lapack error,  ierr=',ierr
    endif
  enddo
!$omp end do
!$omp end parallel
  write(*,*) 'GKMG path k loop complete'

  write(*,*) 'Full grid k loop commencing'
  energgrid=0.0_dp
!$omp parallel default(none), shared(kkgrid,energgrid,params), private(ii,hh,work,rwork,ierr,jj)
!$omp do
  do ii = 1, numk
    do jj = 1, numk
      hh = Ham(aa,kkgrid(ii,jj,:),params,n_size)
      work=0.0_dp
      rwork=0.0_dp
      call ZHEEV('N','U',n_size,hh,n_size,energgrid(ii,jj,:),work,lwork,rwork,ierr)
! the following line can be used to determine the optimal dimension of the work array:
!    write(*,*) 'optimal lwork: ',work(1)
      if(ierr/=0)then
        write(*,*) 'k-point number',ii,'encountered lapack error, ierr=',ierr
      endif
    enddo
  enddo
!$omp end do
!$omp end parallel
  write(*,*) 'Full grid k loop complete'

! Here we set export precision:
  if(lhiprec)then
    cfmt1='(i5,4f12.8,3f7.3)'
    cfmt2='(2f7.3,4f12.8)'
  else
    cfmt1='(i5,4f7.3,3f7.3)'
    cfmt2='(2f7.3,4f7.3)'
  endif
! Exporting eigenvalues to ascii file compatible with gnuplot:
  open(unit=11,file='bandstruct_GKMG.dat',status='unknown')
  open(unit=12,file='bandstruct_fullgrid.dat',status='unknown')
  ktrack=0.0_dp ! variable tracking the total distance traveled along GKM triangle
  write(11,*) '# k-point   energies(eV)                          k-distance(1/A)  k-vectors(1/A)'
  do ii=1,numk
    write(11,cfmt1) ii,energline(ii,:),ktrack,kkline(ii,:)
    if(ii==1)then
      ktrack=ktrack+norm2(kkline(ii,:))
    else
      ktrack=ktrack+norm2(kkline(ii,:)-kkline(ii-1,:))
    endif
  enddo
  write(12,*) '# k-vectors(1/A)      energies(eV)'
  do ii=1,numk
    do jj=1,numk
      write(12,cfmt2) kkgrid(ii,jj,:),energgrid(ii,jj,:)
    enddo
  enddo
  write(*,*) 'Band structure calculation complete. Bands saved in file: bandstruct.dat'
  close(11)
  close(12)
  deallocate(hh,kkline,energline,kkgrid,energgrid,work,rwork)
  iend = time()
  write(*,*) 'Elapsed time: ',iend-istart,' seconds'
  write(*,*) 'Program complete.'
1001 format(i5,4f12.8,3f7.3) ! used for high precision output of GKMG path energies
1002 format(i5,4f7.3,3f7.3) ! used for low precision output of GKMG path energies
1003 format(2f7.3,4f12.8) ! used for high precision output of full grid energies
1004 format(2f7.3,4f7.3) ! used for low precision output of full grid energies
end program
