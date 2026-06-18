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

    call getarg(1,layer)    !第几层
    read(layer,'(I3)') layer_num

    open(unit=20,file='DOS.totalspin',status='old',action='read',iostat=status)
    if (status/=0) then
         write(*,*) 'DOS file open failed!'
    end if
    read(20,*)

    DO i=1,2000 
       read(20,*) 
    end DO

    DO i=1,2000 
       read(20,*) DOS(i,1:2)
    end DO
   
    inquire(file='merge_DOS.dat',exist=alive)
    if (.NOT. alive) then
       open(unit=21,file='merge_DOS.dat',status='replace',action='readwrite',iostat=status)
    else
       open(unit=21,file='merge_DOS.dat',status='old',action='readwrite',iostat=status)
    end if

    DO i2=1,(layer_num-1)*2000
       read(21,*)
    end DO

    DO j2=1,2000
       write(21,"(2F16.9)") DOS(j2,1:2)
    end DO

101  format(I3,3F16.9,2X,3I2,I5)
    write(*,*) 'Done!'
    end program
    
    
!2023/11/28, By Yueyang.
    
    
