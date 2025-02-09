module interpolation
    use, intrinsic :: iso_c_binding
    implicit none

    contains


        function INTERP_NEAREST(grid_values,xs,ys)result(interp_values)
            real(kind=c_double),dimension(:,:),allocatable,intent(in):: grid_values
            real(kind=c_double),dimension(:),allocatable,intent(in):: xs,ys
            real(kind=c_double),dimension(:),allocatable:: interp_values
            integer:: i,N,Npart,xindex,yindex

            Npart=SIZE(xs)
            allocate(interp_values(Npart))
            N=SIZE(grid_values,dim=1)
            do i=1,Npart
                xindex=CEILING(N*MODULO(xs(i),1.0))
                yindex=CEILING(N*MODULO(ys(i),1.0))
                if (xindex==0) then
                    xindex=N
                end if
                if (yindex==0) then
                    yindex=N
                end if
                interp_values(i)=grid_values(xindex,yindex)
            end do

        end function INTERP_NEAREST

        function INTERP_LINEAR(grid_values,xs,ys)result(interp_values)
            real(kind=c_double),dimension(:,:),allocatable,intent(in):: grid_values
            real(kind=c_double),dimension(:),allocatable:: xs,ys
            real(kind=c_double),dimension(:),allocatable:: interp_values
            real(kind=c_double):: x1,y1,x2,y2,x,y,f11,f12,f21,f22
            integer:: i,N,Npart,xindex,yindex,xindex2,yindex2

            Npart=SIZE(xs)
            allocate(interp_values(Npart))
            N=SIZE(grid_values,dim=1)
            do i=1,Npart
                xindex=CEILING(N*MODULO(xs(i),1.0))
                yindex=CEILING(N*MODULO(ys(i),1.0))
                if (xindex==0) then
                    xindex=N
                end if
                if (yindex==0) then
                    yindex=N
                end if
                xindex2=MODULO(xindex,N)+1
                yindex2=MODULO(yindex,N)+1
                x1=real(xindex-1,kind=c_double)/N
                y1=real(yindex-1,kind=c_double)/N
                x2=real(xindex2-1,kind=c_double)/N
                y2=real(yindex2-1,kind=c_double)/N
                x=N*MIN(ABS(xs(i)-x1),ABS(1-xs(i)-x1))
                y=N*MIN(ABS(ys(i)-y1),ABS(1-ys(i)-y1))
                f11=grid_values(xindex,yindex)
                f12=grid_values(xindex,yindex2)
                f21=grid_values(xindex2,yindex)
                f22=grid_values(xindex2,yindex2)
                interp_values(i)=f11*(1-x)*(1-y)+f21*x*(1-y)+f12*(1-x)*y+f22*x*y

            end do

        end function INTERP_LINEAR

        function INTERP_CUBIC_FD(grid_values,xs,ys)result(interp_values)
            real(kind=c_double),dimension(:,:),allocatable,intent(in):: grid_values
            real(kind=c_double),dimension(:,:),allocatable:: f,f_x,f_y,f_xy
            real(kind=c_double),dimension(:),allocatable,intent(in):: xs,ys
            real(kind=c_double),dimension(:),allocatable:: interp_values
            real(kind=c_double):: x1,x2,y1,y2,xbar,ybar
            integer:: i,j,k,N,Npart,xindex,yindex,xindex2,yindex2
            real(kind=c_double),dimension(4,4):: coeff,coeff2
            real(kind=c_double),dimension(4):: valuesx,valuesy
            integer,dimension(4,4):: mat1,mat2

            Npart=SIZE(xs)
            allocate(interp_values(Npart))
            N=SIZE(grid_values,dim=1)
            f=grid_values
            allocate(f_x(N,N))
            allocate(f_y(N,N))
            allocate(f_xy(N,N))
            
            mat1(:,1)=(/1,0,0,0/)
            mat1(:,2)=(/0,0,1,0/)
            mat1(:,3)=(/-3,3,-2,-1/)
            mat1(:,4)=(/2,-2,1,1/)
            
            mat2(:,1)=(/1,0,-3,2/)
            mat2(:,2)=(/0,0,3,-2/)
            mat2(:,3)=(/0,1,-2,1/)
            mat2(:,4)=(/0,0,-1,1/)

            do i=1,N
                do j=1,N
                    f_x(i,j)=0.5*N*(f(MODULO(i,N)+1,j)-f(N-MODULO(1-i,N),j))
                    f_y(i,j)=0.5*N*(f(i,MODULO(j,N)+1)-f(i,N-MODULO(1-j,N)))
                    f_xy(i,j)=0.25*N*N*(f(MODULO(i,N)+1,MODULO(j,N)+1)&
                    &-f(MODULO(i,N)+1,N-MODULO(1-j,N))-f(N-MODULO(1-i,N),MODULO(j,N)+1)&
                    &+f(N-MODULO(1-i,N),N-MODULO(1-j,N)))
                end do 
            end do

            do i=1,Npart
                xindex=CEILING(N*MODULO(xs(i),1.0))
                yindex=CEILING(N*MODULO(ys(i),1.0))
                if (xindex==0) then
                    xindex=N
                end if
                if (yindex==0) then
                    yindex=N
                end if
                xindex2=MODULO(xindex,N)+1
                yindex2=MODULO(yindex,N)+1

                x1=real(xindex-1,kind=c_double)/N
                y1=real(yindex-1,kind=c_double)/N
                x2=real(xindex2-1,kind=c_double)/N
                y2=real(yindex2-1,kind=c_double)/N

                coeff(:,1)=(/f(xindex,yindex),f(xindex2,yindex),f_x(xindex,yindex)/N,f_x(xindex2,yindex)/N/)
                coeff(:,2)=(/f(xindex,yindex2),f(xindex2,yindex2),f_x(xindex,yindex2)/N,f_x(xindex2,yindex2)/N/)
                coeff(:,3)=(/f_y(xindex,yindex)/N,f_y(xindex2,yindex)/N,f_xy(xindex,yindex)/(N*N),f_xy(xindex2,yindex)/(N*N)/)
                coeff(:,4)=(/f_y(xindex,yindex2)/N,f_y(xindex2,yindex2)/N,f_xy(xindex,yindex2)/(N*N),f_xy(xindex2,yindex2)/(N*N)/)

                coeff2=real(0.0,kind=c_double)
                do j=1,4
                    do k=1,4
                        coeff2(:,j)=coeff2(:,j)+coeff(:,k)*mat1(k,j)
                    end do
                end do
                coeff=real(0.0,kind=c_double)
                do j=1,4
                    do k=1,4
                        coeff(:,j)=coeff(:,j)+mat2(:,k)*coeff2(k,j)
                    end do
                end do
                xbar=N*MIN(ABS(xs(i)-x1),ABS(1-xs(i)-x1))
                ybar=N*MIN(ABS(ys(i)-y1),ABS(1-ys(i)-y1))
                valuesx=(/real(1.0,kind=c_double),xbar,xbar*xbar,xbar*xbar*xbar/)
                valuesy=(/real(1.0,kind=c_double),ybar,ybar*ybar,ybar*ybar*ybar/)
                interp_values(i)=real(0.0,kind=c_double)
                do j=1,4
                    do k=1,4
                        interp_values(i)=interp_values(i)+valuesx(j)*coeff(j,k)*valuesy(k)
                    end do
                end do
            end do

        end function INTERP_CUBIC_FD

end module interpolation