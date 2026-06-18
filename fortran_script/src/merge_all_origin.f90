    Program Color_PDOS
    implicit none
    integer :: Natom
    real,dimension(3,3) :: AL
    integer, allocatable, dimension(:) :: Type_atom
    real, allocatable, dimension(:,:) :: Pos_frac, Pos_cart
    real, allocatable, dimension(:,:) :: ATOM 
    real, allocatable, dimension(:) :: ATOM_tot 
    real, allocatable, dimension(:) :: ATOM_ave
    real, allocatable, dimension(:) :: ATOM_z
    integer, allocatable, dimension(:,:) :: Move_atom
    integer :: status
    integer :: i,j,k
    integer :: ii,jj,kk
    integer :: ji,jf
    integer :: i2,j2
    character(10) :: chr_ave_num
    integer :: N_int, N_tail
    character(3) :: layer,ave
    integer :: layer_num,ave_num
    real :: ave_coor
    real, dimension(4000,50) :: DOS 
    logical  alive
    real, dimension(1000) :: merge_coor
    real, dimension(1000,2) :: merge_DOS
    integer  :: ios,ios2

    open(unit=30,file='merge_coor.dat',status='old',action='read',iostat=status)
    if (status/=0) then
         write(*,*) 'coor file open failed!'
    end if

    open(unit=31,file='merge_DOS.dat',status='old',action='read',iostat=status)
    if (status/=0) then
         write(*,*) 'DOS file open failed!'
    end if

    open(unit=32,file='merge_all_origin.dat',status='replace',action='readwrite',iostat=status)

    DO 
       read(30,*,iostat=ios) merge_coor(i)
       read(31,*,iostat=ios2) merge_DOS(i,1:2)
       write(32,"(3F16.9)"), merge_coor(i), merge_DOS(i,1:2)
          if (ios2<0) exit
    end DO

    write(*,*) 'Done!'
    end program
    
    
!2023/11/28, By Yueyang.
    
    
