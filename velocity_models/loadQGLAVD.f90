module loadQGLAVD
    
    use, intrinsic :: iso_c_binding
    use velocitymodule
    implicit none
    include 'fftw3.f03'

    real(kind=c_double),parameter:: c2=0.13502027922908531468,&
    c3=0.58712036810181981280, c4=0.57332577194527991038, c5=1.0  

    TYPE QGfield
        integer:: N,timestep,timewindow
        TYPE(c_ptr):: planf,planb
        real(kind=c_double):: timescale,basinscale,dt
        real(kind=c_double),allocatable,dimension(:,:):: vorticity1,vorticity2,u1,u2,v1,v2,psi1,psi2
        TYPE(Velocity):: vel
        

        contains
            procedure:: init => QG_INIT
            procedure:: update => QG_UPDATE

    END TYPE QGfield

    contains

        subroutine QG_INIT(this,dt,timescale,N,basinscale,R_d,H1,H2)
            Class(QGfield):: this
            TYPE(c_ptr):: planf,planb
            integer:: i,j,j1,j2,k1,k2,k,N,timewindow
            real(kind=c_double):: basinscale,R_d,H1,H2,timescale,dt,dt_day,pi
            real(kind=c_double),allocatable,dimension(:):: wave_numbers
            real(kind=c_double),allocatable,dimension(:,:):: in,wave_x,wave_y,zeta1_d,zeta2_d
            complex(kind=c_double_complex),allocatable,dimension(:,:):: out
            open(21,file='/path/to/data/streamfunction.d',form='unformatted',access='sequential')
            this%N=N
            this%dt=dt
            timewindow=CEILING(0.5/dt)

            this%timewindow=timewindow
            this%timescale=timescale
            this%basinscale=basinscale


            call this%vel%create_vector('velocity')
            call this%vel%create_scalar('vorticity')
            call this%vel%create_scalar('mean_vort')

            allocate(this%psi1(N,N))
            allocate(this%psi2(N,N))
            allocate(this%u1(N,N))
            allocate(this%u2(N,N))
            allocate(this%v1(N,N))
            allocate(this%v2(N,N))
            allocate(this%vorticity1(N,N))
            allocate(this%vorticity2(N,N))

            read(21) this%psi1
            read(21) this%psi2

            this%psi1=(this%basinscale/this%N)*this%psi1
            this%psi2=(this%basinscale/this%N)*this%psi2

            do j=1,N
                do k=1,N
                    k1=MODULO(k,N)+1
		    j1=MODULO(j,N)+1
		    k2=k-1
                    j2=j-1
	            if (j==1) then
			j2=N
		    end if
                    if (k==1) then
                        k2=N
                    end if
                    this%u1(j,k)=0.5*N*(this%psi1(j,k2)-this%psi1(j,k1))/this%basinscale
                    this%v1(j,k)=0.5*N*(this%psi1(j1,k)-this%psi1(j2,k))/this%basinscale
                    this%u2(j,k)=0.5*N*(this%psi2(j,k2)-this%psi2(j,k1))/this%basinscale
                    this%v2(j,k)=0.5*N*(this%psi2(j1,k)-this%psi2(j2,k))/this%basinscale
                    this%vorticity1(j,k) = N*N*(this%psi1(j2,k)+this%psi1(j1,k)+this%psi1(j,k2)+this%psi1(j,k1)-4*this%psi1(j,k))/(this%basinscale*this%basinscale)
                    this%vorticity2(j,K) = N*N*(this%psi2(j2,k)+this%psi2(j1,k)+this%psi2(j,k2)+this%psi2(j,k1)-4*this%psi2(j,k))/(this%basinscale*this%basinscale)
                end do
            end do

            this%u1=(86400/this%basinscale)*this%u1
            this%v1=(86400/this%basinscale)*this%v1
            this%vorticity1=86400*this%vorticity1
            this%u2=(86400/this%basinscale)*this%u2
            this%v2=(86400/this%basinscale)*this%v2
            this%vorticity2=86400*this%vorticity2
        end subroutine QG_INIT

        subroutine QG_UPDATE(this)
            Class(QGfield):: this
            integer:: N,j,k,j1,j2,k1,k2
            real(kind=c_double):: T
            real(kind=c_double),allocatable,dimension(:,:):: zeta1_d,zeta2_d,grid_data1, &
            & grid_data2,grid_data3,grid_data4

            N=this%N
            allocate(grid_data4(N,N))
            T=REAL(MODULO(this%timestep,this%timewindow),kind=c_double)/this%timewindow
            grid_data1=T*this%u2+(real(1.0,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%vorticity2+(real(1.0,kind=c_double)-T)*this%vorticity1
            grid_data4=SUM(grid_data3)/(this%N*this%N)
            call this%vel%set_vector('velocity',1,grid_data1,grid_data2)
            call this%vel%set_scalar('vorticity',1,grid_data3)
            call this%vel%set_scalar('mean_vort',1,grid_data4)


            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c2)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%vorticity2+(real(1.0,kind=c_double)-T)*this%vorticity1
            grid_data4=SUM(grid_data3)/(this%N*this%N)
            call this%vel%set_vector('velocity',2,grid_data1,grid_data2)
            call this%vel%set_scalar('vorticity',2,grid_data3)
            call this%vel%set_scalar('mean_vort',2,grid_data4)
            


            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c3)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%vorticity2+(real(1.0,kind=c_double)-T)*this%vorticity1
            grid_data4=SUM(grid_data3)/(this%N*this%N)
            call this%vel%set_vector('velocity',3,grid_data1,grid_data2)
            call this%vel%set_scalar('vorticity',3,grid_data3)
            call this%vel%set_scalar('mean_vort',3,grid_data4)

            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c4)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%vorticity2+(real(1.0,kind=c_double)-T)*this%vorticity1
            grid_data4=SUM(grid_data3)/(this%N*this%N)
            call this%vel%set_vector('velocity',4,grid_data1,grid_data2)
            call this%vel%set_scalar('vorticity',4,grid_data3)
            call this%vel%set_scalar('mean_vort',4,grid_data4)

            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c5)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%vorticity2+(real(1.0,kind=c_double)-T)*this%vorticity1
            grid_data4=SUM(grid_data3)/(this%N*this%N)
            call this%vel%set_vector('velocity',5,grid_data1,grid_data2)
            call this%vel%set_scalar('vorticity',5,grid_data3)
            call this%vel%set_scalar('mean_vort',5,grid_data4)
            

            if (MODULO(this%timestep+1,this%timewindow)==0) then
                print*, "update"
                this%u1=this%u2
                this%v1=this%v2
                this%vorticity1=this%vorticity2
                this%psi1=this%psi2
                read(21) this%psi2

                this%psi2=(this%basinscale/this%N)*this%psi2

                do j=1,N
                    do k=1,N
			j1=MODULO(j,N)+1
			k1=MODULO(k,N)+1
			j2=j-1
			k2=k-1
			if (j == 1) then
			    j2=N	
			end if
			if (k == 1) then
			    k2=N
			end if 
                        this%u2(j,k)=0.5*N*(this%psi2(j,k2)-this%psi2(j,k1))/this%basinscale
                        this%v2(j,k)=0.5*N*(this%psi2(j1,k)-this%psi2(j2,k))/this%basinscale
                        this%vorticity2(j,K) = N*N*(this%psi2(j2,k)+this%psi2(j1,k)+this%psi2(j,k2)+this%psi2(j,k1)-4*this%psi2(j,k))/(this%basinscale*this%basinscale)
                    end do
                end do

                this%u2=(86400/this%basinscale)*this%u2
                this%v2=(86400/this%basinscale)*this%v2
                this%vorticity2=86400*this%vorticity2

            
            end if

            this%timestep = this%timestep+1


        end subroutine QG_UPDATE

end module loadQGLAVD
