module deltacorrelated
    use, intrinsic :: iso_c_binding
    use velocitymodule
    implicit none
    include 'fftw3.f03'
    real(kind=c_double),parameter :: pi=4*ATAN(REAL(1.0,kind=c_double))

    TYPE DeltaCorrelatedField
        integer:: n_modes,n_steps,i
        TYPE(Velocity):: vel
        TYPE(c_ptr):: planb
        real(kind=c_double):: gamma,sigma_u,l_corr,t_mem,dt
        real(kind=c_double),dimension(:,:),allocatable:: a_s,a_p,b_s,b_p,en_spec,&
        &kernel_x,kernel_y,kernel_div

        contains
            procedure:: init => DELTACORR_INIT
            procedure:: update => DELTACORR_UPDATE
            procedure:: update_spectrum => DELTACORR_UPDATE_SPECTRUM
    END TYPE

    contains
        subroutine DELTACORR_INIT(this,dt,path)
            integer:: i
            character(len=*),intent(in):: path
            Class(DeltaCorrelatedField):: this
            TYPE(c_ptr):: planb
            real(kind=c_double):: dt
            real(kind=c_double),dimension(:),allocatable:: wave_numbers
            real(kind=c_double),dimension(:,:),allocatable:: initial,wave_sq,wave_x,wave_y,out
            complex(kind=c_double_complex),dimension(:,:),allocatable:: in

            call DELTACORR_LOAD_PARAMETERS(this,path)
            this%dt=0.5*dt
            this%i=0
            this%n_steps=NINT(this%t_mem/this%dt)     
            print*, this%n_steps     
            allocate(in(this%n_modes/2+1,this%n_modes))
            allocate(out(this%n_modes,this%n_modes))
            in=DCMPLX(0.0,0.0)
            out=real(0.0,kind=c_double)
            print*, "preparing FFT plan"
            call dfftw_plan_dft_c2r_2d(planb,this%n_modes,this%n_modes,in,out,FFTW_MEASURE)
            print*, "done"
            this%planb=planb
            allocate(initial(this%n_modes,this%n_modes))
            print*, "creating vector/scalar fields"
            call this%vel%create_vector('velocity')
            call this%vel%create_scalar('divergence')
            initial=REAL(0.0,kind=c_double)
            print*, "initialising grid values for all fields"
            do i=1,3
                call this%vel%set_vector('velocity',i,initial,initial)
                call this%vel%set_scalar('divergence',i,initial)
            end do

            allocate(this%a_s(this%n_modes/2+1,this%n_modes))
            allocate(this%a_p(this%n_modes/2+1,this%n_modes))
            allocate(this%b_s(this%n_modes/2+1,this%n_modes))
            allocate(this%b_p(this%n_modes/2+1,this%n_modes))
            this%a_s=0.0
            this%a_p=0.0
            this%b_s=0.0
            this%b_p=0.0

            allocate(wave_numbers(this%n_modes))
            allocate(wave_x(this%n_modes/2+1,this%n_modes))
            allocate(wave_y(this%n_modes/2+1,this%n_modes))
            allocate(wave_sq(this%n_modes/2+1,this%n_modes))

            do i=0,this%n_modes/2
                wave_numbers(i+1)=2*pi*i
            end do
            do i=this%n_modes/2+1,this%n_modes-1
                wave_numbers(i+1)=2*pi*(i-this%n_modes)
            end do

            do i=1,this%n_modes
                wave_y(:,i)=wave_numbers(i)
            end do

            do i=1,this%n_modes/2+1
                wave_x(i,:)=wave_numbers(i)
            end do

            wave_sq=wave_x*wave_x+wave_y*wave_y
            allocate(this%en_spec(this%n_modes/2+1,this%n_modes))

            this%en_spec=DSQRT((1/(8*pi))*this%l_corr*this%l_corr*this%l_corr*this%l_corr*&
            &wave_sq*DEXP(-0.5*wave_sq*this%l_corr*this%l_corr))

            this%kernel_div=DSQRT(wave_sq)

            wave_sq(1,1)=real(1.0,kind=c_double)
            this%kernel_x=wave_x/DSQRT(wave_sq)
            this%kernel_y=wave_y/DSQRT(wave_sq)
            
        end subroutine DELTACORR_INIT

        subroutine DELTACORR_UPDATE(this)
            Class(DeltaCorrelatedField):: this
            integer:: i,M

            
            M=SIZE(this%vel%scalar_fields)
            do i=1,M
                this%vel%scalar_fields(i)%grid_values1=this%vel%scalar_fields(i)%grid_values3
            end do
            if (MODULO(this%i,this%n_steps)==0) then
                call DELTACORR_UPDATE_SPECTRUM(this)   

                call DELTACORR_VELOCITY_OUTPUT(this,2)

                call DELTACORR_UPDATE_SPECTRUM(this)

                call DELTACORR_VELOCITY_OUTPUT(this,3)
                
                this%i=0
            end if
            this%i=this%i+1
        end subroutine DELTACORR_UPDATE

        subroutine DELTACORR_UPDATE_SPECTRUM(this)
            Class(DeltaCorrelatedField):: this
            
            this%a_p=this%sigma_u*this%en_spec*GAUSSIAN(this%n_modes)
            this%a_s=this%sigma_u*this%en_spec*GAUSSIAN(this%n_modes)
            this%b_p=this%sigma_u*this%en_spec*GAUSSIAN(this%n_modes)
            this%b_s=this%sigma_u*this%en_spec*GAUSSIAN(this%n_modes)

        end subroutine

        subroutine DELTACORR_VELOCITY_OUTPUT(this,timestep)
            Class(DeltaCorrelatedField):: this
            integer:: timestep
            real(kind=c_double),dimension(:,:),allocatable:: grid_data1,grid_data2,grid_data3 
            complex(kind=c_double_complex),dimension(:,:),allocatable:: spectrum

            spectrum=2*pi*(this%gamma*DCMPLX(this%a_p,this%b_p)*this%kernel_x+&
            & (1-this%gamma)*DCMPLX(this%a_s,this%b_s)*this%kernel_y)
            
            call DELTACORR_FFT(this%planb,spectrum,grid_data1)

            spectrum=2*pi*(this%gamma*DCMPLX(this%a_p,this%b_p)*this%kernel_y-&
            & (1-this%gamma)*DCMPLX(this%a_s,this%b_s)*this%kernel_x)

            call DELTACORR_FFT(this%planb,spectrum,grid_data2)

            spectrum=2*pi*this%gamma*DCMPLX(-1*this%b_p,this%a_p)*this%kernel_div

            call DELTACORR_FFT(this%planb,spectrum,grid_data3)

            call this%vel%set_vector('velocity',timestep,grid_data1,grid_data2)
            call this%vel%set_scalar('divergence',timestep,grid_data3)
            
        end subroutine DELTACORR_VELOCITY_OUTPUT

        subroutine DELTACORR_LOAD_PARAMETERS(this,path)
            Class(DeltaCorrelatedField):: this
            character(len=*),intent(in):: path
            character(len=100) :: buffer, label
            real(kind=c_double) :: sigma,gamma,lamb,lcorr
            integer :: pos,N
            integer, parameter :: fh = 15
            integer :: ios = 0
            integer :: line = 0

            open(fh, file=path)
      
            ! ios is negative if an end of record condition is encountered or if
            ! an endfile condition was detected.  It is positive if an error was
            ! detected.  ios is zero otherwise.
          
            do while (ios == 0)
               read(fh, '(A)', iostat=ios) buffer
               if (ios == 0) then
                  line = line + 1
          
                  ! Find the first instance of whitespace.  Split label and data.
                  pos = scan(buffer, '    ')
                  label = buffer(1:pos)
                  buffer = buffer(pos+1:)
          
                  select case (label)
                  case ('Nmodes')
                     read(buffer, *, iostat=ios) N
                     print *, 'Read Nmodes: ', N
                  case ('sigma')
                     read(buffer, *, iostat=ios) sigma
                     print *, 'Read sigma: ', sigma
        
                  case ('lcorr')
                        read(buffer, *, iostat=ios) lcorr
                        print *, 'Read lcorr: ', lcorr
                  case ('gamma')
                        read(buffer, *, iostat=ios) gamma
                        print *, 'Read gamma: ', gamma
                  case ('lamb')
                        read(buffer, *, iostat=ios) lamb
                        print *, 'Read lamb: ', lamb
                  case default
                     print *, 'Skipping invalid label at line', line
                  end select
               end if
            end do

            this%n_modes=N
            this%gamma=gamma
            this%sigma_u=sigma
            this%t_mem=1/lamb
            this%l_corr=lcorr
        end subroutine DELTACORR_LOAD_PARAMETERS

        function GAUSSIAN(N)result(normal)
            !function to generate an (N/2+1)xN array of independent standard normal random variable
            integer :: N,i,j
            real(c_double), dimension(:,:), allocatable :: normal
            real(c_double) :: U1,U2
            
            allocate(normal(N/2+1,N))
            
            do i=1,N/2+1
                do j=1,N/2
                    CALL RANDOM_NUMBER(U1)
                    CALL RANDOM_NUMBER(U2)
                    normal(i,j)=DSQRT(-2*LOG(U1))*COS(2*pi*U2)
                    normal(i,j+N/2)=DSQRT(-2*LOG(U1))*SIN(2*pi*U2)
                end do
            end do
        end function GAUSSIAN

        subroutine DELTACORR_FFT(plan,in,out)
            TYPE(c_ptr):: plan
            integer:: N
            real(kind=c_double),dimension(:,:),allocatable:: out
            complex(kind=c_double_complex),dimension(:,:),allocatable:: in
            N=SIZE(in,dim=2)
            if (ALLOCATED(out)) then
                deallocate(out)
            end if
            allocate(out(N,N))
            call dfftw_execute_dft_c2r(plan,in,out)


        end subroutine DELTACORR_FFT

end module deltacorrelated
