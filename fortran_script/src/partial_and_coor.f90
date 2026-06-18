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
    call getarg(2,ave)    !每层多少原子
    read(layer,'(I3)') layer_num
    read(ave,'(I3)') ave_num
    open(unit=9,file='atom.config',status='old',action='read',iostat=status)
    if (status/=0) then
         write(*,*) 'atom file open failed!'
    end if
    read(9,*) Natom
    allocate( ATOM(Natom,7))
    allocate( ATOM_tot(Natom))
    allocate( ATOM_ave(Natom))
    allocate( ATOM_z(Natom))
    read(9,*)
    read(9,*) AL(1,1), AL(1,2), AL(1,3)
    read(9,*) AL(2,1), AL(2,2), AL(2,3)
    read(9,*) AL(3,1), AL(3,2), AL(3,3) 
    read(9,*)

    DO i=1,Natom
        read(9,*) ATOM(i,1:7)
    end DO
    write(*,*) "Natom is: ",NATOM

    open(unit=10,file='atom.config_1',status='replace',action='readwrite',iostat=status)
    write(10,*) Natom
    write(10,*) 'Lattice'
    write(10,"(3F16.9)") AL(1,1), AL(1,2), AL(1,3)
    write(10,"(3F16.9)") AL(2,1), AL(2,2), AL(2,3)
    write(10,"(3F16.9)") AL(3,1), AL(3,2), AL(3,3) 
    write(10,*) 'Positions' 

    N_int=int(NATOM/ave_num)
    N_tail=mod(NATOM,ave_num)
    write(*,*) "N_int and N_tail is: ",N_int, N_tail
    ave_coor=0.0
    
    DO kk=1,(layer_num-1)*ave_num
       write(10,101) int(ATOM(kk,1)),ATOM(kk,2:4),int(ATOM(kk,5:7)),0
    end DO

    DO kk=(layer_num-1)*ave_num+1,layer_num*ave_num
       write(10,101) int(ATOM(kk,1)),ATOM(kk,2:4),int(ATOM(kk,5:7)),1
       ave_coor=ave_coor+ATOM(kk,4) 
    end DO
 !      ave_coor=ave_coor/ave_num
       ave_coor=ave_coor/ave_num*AL(3,3)
    write(*,*)'The ave_coor is: ', ave_coor
    DO kk=(layer_num*ave_num)+1,NATOM
       write(10,101) int(ATOM(kk,1)),ATOM(kk,2:4),int(ATOM(kk,5:7)),0
    end DO

 !  open(unit=11,file='DOS.totalspin',status='old',action='read',iostat=status)
 !  if (status/=0) then
 !       write(*,*) 'DOS file open failed!'
 !  end if
 !      
 !  read(11,*)
 !  DO i2=1,4000 
 !  read(11,*) DOS(i2,1:2)
 !  end DO
 !  write(*,*)'DOS(1,1) and DOS(1,2) is: ', DOS(1,1),DOS(1,2) 
   
    inquire(file='merge_coor.dat',exist=alive)
    if (.NOT. alive) then
       open(unit=11,file='merge_coor.dat',status='replace',action='readwrite',iostat=status)
    else 
       open(unit=11,file='merge_coor.dat',status='old',action='readwrite',iostat=status)
    end if

    DO i2=1,(layer_num-1)*2000
       read(11,*)
    end DO

    DO j2=1,2000
       write(11,"(F16.9)") ave_coor
    end DO

101  format(I3,3F16.9,2X,3I2,I5)
!    open(unit=11,file='layer_coordinate.dat',status='replace',action='readwrite',iostat=status)
    write(*,*) 'Done!'
    end program
    
    
!2023/11/28, By Yueyang.
    
    
