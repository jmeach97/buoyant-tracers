module particlemodule

    use, intrinsic :: iso_c_binding
    use velocitymodule
    implicit none

    TYPE ParticleAttribute
        character(len=:),allocatable:: name
        integer:: position

    END TYPE ParticleAttribute

    TYPE ParticleParameter
        character(len=:),allocatable:: name
        integer:: position
        real(kind=c_double):: param_value

    END TYPE ParticleParameter

    TYPE ParticleManager
        TYPE(ParticleAttribute),dimension(:),allocatable::  attributes
        TYPE(ParticleParameter),dimension(:),allocatable:: particle_parameters
        real(kind=c_double):: dt
        character(len=:),allocatable:: evaluate_type
        integer,dimension(:),allocatable:: ids
        integer:: num_attributes,n_part,timestep
        logical,dimension(:),allocatable:: alive
        real(kind=c_double),dimension(:,:),allocatable:: attribute_values,arg

        contains
            procedure:: init => PARTICLE_INIT
            procedure:: create_attribute => PARTICLE_CREATE_ATTRIBUTE
            procedure:: create_parameter => PARTICLE_CREATE_PARAMETER
            procedure:: get => PARTICLE_GET_ATTRIBUTE_VAlUE
            procedure:: get_param => PARTICLE_GET_PARAMETER
            procedure:: get_scalar => PARTICLE_GET_VELOCITY_SCALAR
            procedure:: get_vector => PARTICLE_GET_VELOCITY_VECTOR
            procedure:: get_pos => PARTICLE_GET_ATTRIBUTE
            procedure:: update => PARTICLE_RUNGE_KUTTA
            procedure:: save => PARTICLE_SAVE
    END TYPE ParticleManager

    contains
        subroutine PARTICLE_INIT(this,n_part,dt,evaluate_type)
            Class(ParticleManager):: this
            real(kind=c_double):: dt
            character(len=*),intent(in):: evaluate_type
            integer:: n_part
            this%dt=dt
            this%n_part=n_part
            this%evaluate_type=evaluate_type
            allocate(this%ids(this%n_part))
            allocate(this%alive(this%n_part))
            this%alive=.True.
            allocate(this%attribute_values(this%n_part,this%num_attributes))
        end subroutine

        subroutine PARTICLE_CREATE_PARAMETER(this,name,param_value)
            Class(ParticleManager):: this
            character(len=*),intent(in):: name
            TYPE(ParticleParameter),dimension(:),allocatable:: parameters_new
            integer:: M,i
            real(kind=c_double):: param_value
            if (allocated(this%particle_parameters)) then 
                M=SIZE(this%particle_parameters)
                allocate(parameters_new(M+1))
                do i=1,M
                    parameters_new(i)=this%particle_parameters(i)
                end do
                deallocate(this%particle_parameters)
            else
                M=0
                allocate(parameters_new(1))
            end if

            parameters_new(M+1)%name=name
            parameters_new(M+1)%position=M+1
            parameters_new(M+1)%param_value=param_value
            
            this%particle_parameters=parameters_new
        end subroutine PARTICLE_CREATE_PARAMETER

        function PARTICLE_GET_PARAMETER(this,name)result(param_value)
            Class(ParticleManager):: this
            character(len=*),intent(in):: name
            real(kind=c_double):: param_value
            integer:: i,M  
            M=SIZE(this%particle_parameters)

            do i=1,M
                if (this%particle_parameters(i)%name==name) then
                    param_value=this%particle_parameters(i)%param_value
                end if
            end do

        end function PARTICLE_GET_PARAMETER

        subroutine PARTICLE_CREATE_ATTRIBUTE(this,name)
            Class(ParticleManager):: this
            TYPE(ParticleAttribute),dimension(:),allocatable:: new_attributes
            integer:: i,M
            character(len=*),intent(in):: name
            if (allocated(this%attributes)) then
                M=SIZE(this%attributes)
                allocate(new_attributes(M+1))
                do i=1,M
                    new_attributes(i)=this%attributes(i)
                end do
                deallocate(this%attributes)
            else
                M=0
                allocate(new_attributes(1))
            end if

            
            new_attributes(M+1)%name=name
            new_attributes(M+1)%position=M+1


            this%attributes=new_attributes
            this%num_attributes=M+1

        end subroutine PARTICLE_CREATE_ATTRIBUTE

        function PARTICLE_GET_ATTRIBUTE(this,name)result(position)
            Class(ParticleManager):: this
            integer:: position,M,i
            character(len=*),intent(in):: name

            M=SIZE(this%attributes)
            do i=1,M 
                if (this%attributes(i)%name==name) then
                    position=i
                end if
            end do

        end function PARTICLE_GET_ATTRIBUTE

        function PARTICLE_GET_ATTRIBUTE_VALUE(this,name)result(att_values)
            Class(ParticleManager):: this
            character(len=*),intent(in):: name
            integer:: position
            real(kind=c_double),dimension(:),allocatable:: att_values

            position=PARTICLE_GET_ATTRIBUTE(this,name)
            att_values=this%arg(:,position)
        end function PARTICLE_GET_ATTRIBUTE_VALUE

        function PARTICLE_GET_VELOCITY_SCALAR(this,name,velocity_instance)result(vel_values)
            Class(ParticleManager):: this
            character(len=*),intent(in):: name
            TYPE(Velocity):: velocity_instance
            real(kind=c_double),dimension(:),allocatable:: vel_values,xs,ys


            xs=PARTICLE_GET_ATTRIBUTE_VALUE(this,'x')
            ys=PARTICLE_GET_ATTRIBUTE_VALUE(this,'y')

            vel_values=velocity_instance%scalar(name,this%timestep,this%evaluate_type,xs,ys)

        end function PARTICLE_GET_VELOCITY_SCALAR

        function PARTICLE_GET_VELOCITY_VECTOR(this,name,dimension,velocity_instance)result(vel_values)
            Class(ParticleManager):: this
            character(len=*),intent(in):: name,dimension
            TYPE(Velocity):: velocity_instance
            real(kind=c_double),dimension(:),allocatable:: vel_values,xs,ys

            xs=PARTICLE_GET_ATTRIBUTE_VALUE(this,'x')
            ys=PARTICLE_GET_ATTRIBUTE_VALUE(this,'y')

            vel_values=velocity_instance%vector(name,this%timestep,dimension,this%evaluate_type,xs,ys)

        end function PARTICLE_GET_VELOCITY_VECTOR

        subroutine PARTICLE_RUNGE_KUTTA(this,update_function,velocity_instance)
            Class(ParticleManager):: this
            TYPE(Velocity):: velocity_instance
            real(kind=c_double),dimension(:,:),allocatable:: k1,k2,k3,k4
            interface
                function update_function(this,velocity_instance)result(step_val)
                    import ParticleManager,Velocity,c_double
                    TYPE(ParticleManager),intent(in):: this
                    TYPE(Velocity),intent(in):: velocity_instance
                    real(kind=c_double),dimension(:,:),allocatable:: step_val
                end function update_function
            end interface


            this%arg=this%attribute_values
            this%timestep=1
            k1=update_function(this,velocity_instance)

            this%arg=this%attribute_values+0.5*this%dt*k1
            
            this%timestep=2
            k2=update_function(this,velocity_instance)
            
            this%arg=this%attribute_values+0.5*this%dt*k2

            k3=update_function(this,velocity_instance)

            this%arg=this%attribute_values+this%dt*k3
            
            this%timestep=3
            k4=update_function(this,velocity_instance)

            this%attribute_values=this%attribute_values+this%dt*(real(1.0,kind=c_double)/6)*(k1+2*k2+2*k3+k4)
        end subroutine PARTICLE_RUNGE_KUTTA

        subroutine PARTICLE_SAVE(this,path,fname)
            Class(ParticleManager):: this
            character(len=*),intent(in):: path,fname
            integer:: i,M

            M=this%n_part
            OPEN(UNIT=15, FILE=TRIM(path)//TRIM(ADJUSTL(fname))//".dat", ACTION="write", STATUS="replace")            
            do i=1,M
                WRITE(15,"(I0,',',*(G0,:,','))") this%ids(i),this%attribute_values(i,:)
            end do
            CLOSE(UNIT=15)
        end subroutine PARTICLE_SAVE
end module particlemodule