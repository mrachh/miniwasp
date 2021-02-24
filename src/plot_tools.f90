subroutine plot1D(x,y,N,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: x(N),y(N)

!List of local variables
!character (len=100) nombre
integer ( kind = 8 ) umio,count,flag
integer :: ierror
    open(8, FILE=nombre,STATUS='REPLACE')

    flag=1
    write(8,*) N
    do count=1,N
        write(8,*) x(count)
    enddo
    do count=1,N
        write(8,*) y(count)
    enddo
    close (8)
return
end
!
!
!
!
!
subroutine plot_curve_3D(x,y,z,N,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: x(N),y(N),z(N)

!List of local variables
!character (len=100) nombre
integer ( kind = 8 ) umio,count,flag
integer :: ierror
    open(8, FILE=nombre,STATUS='REPLACE')
    flag=3
    write(8,*) flag
    write(8,*) N
    do count=1,N
        write(8,*) x(count)
    enddo
    do count=1,N
        write(8,*) y(count)
    enddo
    do count=1,N
        write(8,*) z(count)
    enddo
    close (8)
return
end
!
!
!
!
!
subroutine plot2D(x,y,F,N,M,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer ( kind = 8 ), intent(in) :: N,M
real ( kind = 8 ), intent(in) :: x(M,N),y(M,N),F(M,N)

!List of local variables
!character (len=100) nombre
integer ( kind = 8 ) umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=nombre,STATUS='REPLACE')
    flag=2
    write(8,*) flag
    write(8,*) M
    write(8,*) N
    do count1=1,N
        do count2=1,M
            write(8,*) x(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) y(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) F(count2,count1)
        enddo
    enddo
    close (8)
return
end
!
!
!
!
!
subroutine plot2D_v2(F,N,M,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer ( kind = 8 ), intent(in) :: N,M
real ( kind = 8 ), intent(in) :: F(M,N)

!List of local variables
!character (len=100) nombre
integer ( kind = 8 ) umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=nombre,STATUS='REPLACE')
    flag=4
    write(8,*) flag
    write(8,*) M
    write(8,*) N
    do count1=1,N
        do count2=1,M
            write(8,*) F(count2,count1)
        enddo
    enddo
    close (8)
return
end
!
!
!
!
!
subroutine plot2D_v3(x,y,z,F,N,M,nombre)
implicit none

!List of calling arguments
character (len=100), intent(in) :: nombre
integer ( kind = 8 ), intent(in) :: N,M
real ( kind = 8 ), intent(in) :: x(M,N),y(M,N),z(M,N),F(M,N)

!List of local variables
integer ( kind = 8 ) umio,count1,count2,flag
integer :: ierror
!    nombre='fortran_plot'
    open(8, FILE=nombre,STATUS='REPLACE')
    flag=5
    write(8,*) flag
    write(8,*) M
    write(8,*) N
    do count1=1,N
        do count2=1,M
            write(8,*) x(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) y(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) z(count2,count1)
        enddo
    enddo
    do count1=1,N
        do count2=1,M
            write(8,*) F(count2,count1)
        enddo
    enddo
    close (8)
return
end
