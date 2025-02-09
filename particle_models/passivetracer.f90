module passivetracer

    use, intrinsic :: iso_c_binding
    use velocitymodule
    use particlemodule
    implicit none

    TYPE PassiveTracers
        TYPE(ParticleManager):: particle_manager

        contains

            procedure:: init=> PASSIVE_TRACER_INIT
            procedure:: update => PASSIVE_UPDATE
    END TYPE PassiveTracers

    contains

        subroutine PASSIVE_TRACER_INIT(this,dt,sqrt_n_part)
            Class(PassiveTracers):: this
            real(kind=c_double):: dt
            integer:: sqrt_n_part,i,j
            integer,dimension(:),allocatable:: positions

            this%particle_manager%dt=dt
            this%particle_manager%n_part=sqrt_n_part*sqrt_n_part
            call this%particle_manager%create_attribute('x')
            call this%particle_manager%create_attribute('y')
            call this%particle_manager%create_attribute('C')

            call this%particle_manager%init(sqrt_n_part*sqrt_n_part,dt,'linear')

            allocate(positions(this%particle_manager%num_attributes))
            positions(1)=this%particle_manager%get_pos('x')
            positions(2)=this%particle_manager%get_pos('y')
            positions(3)=this%particle_manager%get_pos('C')
            do i=1,sqrt_n_part
                do j=1,sqrt_n_part
                    this%particle_manager%attribute_values(i+sqrt_n_part*(j-1),positions(1))=real(i-1,kind=c_double)/sqrt_n_part
                    this%particle_manager%attribute_values(i+sqrt_n_part*(j-1),positions(2))=real(j-1,kind=c_double)/sqrt_n_part
                    this%particle_manager%attribute_values(i+sqrt_n_part*(j-1),positions(3))=real(1.0,kind=c_double)
                    this%particle_manager%ids(i+sqrt_n_part*(j-1))=i+sqrt_n_part*(j-1)
                end do
            end do
        end subroutine PASSIVE_TRACER_INIT

        subroutine PASSIVE_UPDATE(this,update_function,velocity_instance)
            Class(PassiveTracers):: this
            TYPE(Velocity):: velocity_instance
            interface
                function update_function(this,velocity_instance)result(step_val)
                    import ParticleManager,Velocity,c_double
                    TYPE(ParticleManager),intent(in):: this
                    TYPE(Velocity),intent(in):: velocity_instance
                    real(kind=c_double),dimension(:,:),allocatable:: step_val
                end function update_function
            end interface

            call this%particle_manager%update(update_function,velocity_instance)

        end subroutine PASSIVE_UPDATE

        function PASSIVE_UPDATE_FUNCTION(this,velocity_instance)result(step_val)
            TYPE(ParticleManager),intent(in):: this
            TYPE(Velocity),intent(in):: velocity_instance
            integer,dimension(:),allocatable:: positions
            real(kind=c_double),dimension(:,:),allocatable:: step_val

            allocate(step_val(this%n_part,this%num_attributes))
            allocate(positions(this%num_attributes))

            step_val(:,this%get_pos('x'))=this%get_vector('velocity','x',velocity_instance)
            step_val(:,this%get_pos('y'))=this%get_vector('velocity','y',velocity_instance)
            step_val(:,this%get_pos('C'))=-1*this%get('C')*this%get_scalar('divergence',velocity_instance)
            
        end function PASSIVE_UPDATE_FUNCTION
end module passivetracer