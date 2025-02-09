module loadQG
    
    use, intrinsic :: iso_c_binding
    use velocitymodule
    implicit none
    include 'fftw3.f03'
    real(kind=c_double),parameter:: c2=0.13502027922908531468,&
    c3=0.58712036810181981280, c4=0.57332577194527991038, c5=1.0

    TYPE QGfield
        integer:: N,timestep,timewindow
        TYPE(c_ptr):: planf,planb
        real(kind=c_double):: timescale,basinscale
        real(kind=c_double),allocatable,dimension(:,:):: ageo_div1,ageo_div2,u1,u2,v1,v2,psi1,psi2
        TYPE(Velocity):: vel
        

        contains
            procedure:: init => QG_INIT
            procedure:: update => QG_UPDATE

    END TYPE QGfield

    contains

        subroutine QG_INIT(this,dt,timescale,N,basinscale,R_d,H1,H2)
            Class(QGfield):: this
            TYPE(c_ptr):: planf,planb
            integer:: i,j,k,k1,k2,j1,j2,N,timewindow
            real(4),allocatable,dimension(:,:):: zeta1,zeta2
            real(kind=c_double):: basinscale,R_d,H1,H2,timescale,dt,dt_day,pi
            real(kind=c_double),allocatable,dimension(:):: wave_numbers
            real(kind=c_double),allocatable,dimension(:,:):: in,wave_x,wave_y,zeta1_d,zeta2_d,psi1,psi2
            complex(kind=c_double_complex),allocatable,dimension(:,:):: out
            open(21,file='/path/to/data/x_velocity.unf',form='unformatted',access='sequential')
            open(22,file='/path/to/data/y_velocity.unf',form='unformatted',access='sequential')
            open(23,file='/path/to/data/divergence.unf',form='unformatted',access='sequential')

            timewindow=CEILING(0.5/dt)

            this%timewindow=timewindow
            this%timescale=timescale
            this%basinscale=basinscale
            this%N=N

            call this%vel%create_vector('velocity')
            call this%vel%create_scalar('divergence')

            allocate(this%u1(N,N))
            allocate(this%u2(N,N))
            allocate(this%v1(N,N))
            allocate(this%v2(N,N))
            allocate(this%ageo_div1(N,N))
            allocate(this%ageo_div2(N,N))

            read(21) this%u1
            read(21) this%u2
            read(22) this%v1
            read(22) this%v2
            read(23) this%ageo_div1
            read(23) this%ageo_div2


            this%u1=(86400/this%basinscale)*this%u1
            this%v1=(86400/this%basinscale)*this%v1
            this%ageo_div1=86400*this%ageo_div1
            this%u2=(86400/this%basinscale)*this%u2
            this%v2=(86400/this%basinscale)*this%v2
            this%ageo_div2=86400*this%ageo_div2
        end subroutine QG_INIT

        subroutine QG_UPDATE(this)
            Class(QGfield):: this
            integer:: k,j,N,k1,k2,j1,j2
            real(4),allocatable,dimension(:,:):: zeta1,zeta2
            real(kind=c_double):: T
            real(kind=c_double),allocatable,dimension(:,:):: zeta1_d,zeta2_d,grid_data1, &
            & grid_data2,grid_data3

            T=REAL(MODULO(this%timestep,this%timewindow),kind=c_double)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%ageo_div2+(real(1.0,kind=c_double)-T)*this%ageo_div1
            call this%vel%set_vector('velocity',1,grid_data1,grid_data2)
            call this%vel%set_scalar('divergence',1,grid_data3)



            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c2)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%ageo_div2+(real(1.0,kind=c_double)-T)*this%ageo_div1
            call this%vel%set_vector('velocity',2,grid_data1,grid_data2)
            call this%vel%set_scalar('divergence',2,grid_data3)
            


            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c3)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%ageo_div2+(real(1.0,kind=c_double)-T)*this%ageo_div1
            call this%vel%set_vector('velocity',3,grid_data1,grid_data2)
            call this%vel%set_scalar('divergence',3,grid_data3)

            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c4)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%ageo_div2+(real(1.0,kind=c_double)-T)*this%ageo_div1
            call this%vel%set_vector('velocity',4,grid_data1,grid_data2)
            call this%vel%set_scalar('divergence',4,grid_data3)

            T=(REAL(MODULO(this%timestep,this%timewindow),kind=c_double)+c5)/this%timewindow
            grid_data1=T*this%u2+(real(1,kind=c_double)-T)*this%u1
            grid_data2=T*this%v2+(real(1,kind=c_double)-T)*this%v1
            grid_data3=T*this%ageo_div2+(real(1.0,kind=c_double)-T)*this%ageo_div1
            call this%vel%set_vector('velocity',5,grid_data1,grid_data2)
            call this%vel%set_scalar('divergence',5,grid_data3)
            
            if (MODULO(this%timestep+1,this%timewindow)==0) then
                print*, "update"
                this%u1=this%u2
                this%v1=this%v2
                this%ageo_div1=this%ageo_div2

                read(21) this%u2
                read(22) this%v2
                read(23) this%ageo_div2

                this%u2=(86400/this%basinscale)*this%u2
                this%v2=(86400/this%basinscale)*this%v2
                this%ageo_div2=86400*this%ageo_div2



            
            end if

            this%timestep = this%timestep+1


        end subroutine QG_UPDATE

end module loadQG
