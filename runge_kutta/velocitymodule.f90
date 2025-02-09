module velocitymodule

    use, intrinsic :: iso_c_binding
    use interpolation
    implicit none

    TYPE ScalarField
        character(len=:),allocatable:: name
        integer:: position
        real(kind=c_double),dimension(:,:),allocatable:: grid_values1,grid_values2,grid_values3

    END TYPE ScalarField

    TYPE VectorField
        character(len=:),allocatable:: name
        integer:: position
        real(kind=c_double),dimension(:,:),allocatable:: grid_x1,grid_x2,grid_x3,grid_y1,grid_y2,grid_y3
    END TYPE VectorField

    TYPE Velocity
        TYPE(ScalarField),dimension(:),allocatable:: scalar_fields
        TYPE(VectorField),dimension(:),allocatable:: vector_fields

        contains
            procedure:: scalar => VELOCITY_GET_SCALAR_VALUES
            procedure:: vector => VELOCITY_GET_VECTOR_VALUES
            procedure:: get_vector_grid => VELOCITY_GET_VECTOR
            procedure:: get_scalar_grid => VELOCITY_GET_SCALAR
            procedure:: create_scalar => VELOCITY_CREATE_SCALAR
            procedure:: create_vector => VELOCITY_CREATE_VECTOR
            procedure:: find_scalar => VELOCITY_FIND_SCALAR
            procedure:: find_vector => VELOCITY_FIND_VECTOR
            procedure:: set_scalar => VELOCITY_SET_SCALAR_VALUES
            procedure:: set_vector => VELOCITY_SET_VECTOR_VALUES
            procedure:: save_scalar => VELOCITY_SAVE_SCALAR
            procedure:: save_vector => VELOCITY_SAVE_VECTOR
    END TYPE Velocity


    contains

        subroutine VELOCITY_CREATE_SCALAR(this,name)
            Class(Velocity):: this
            integer:: i,M
            TYPE(ScalarField),dimension(:),allocatable:: new_scalars
            character(len=*),intent(in):: name

            if (ALLOCATED(this%scalar_fields)) then
                M=SIZE(this%scalar_fields)
                allocate(new_scalars(M+1))
                do i=1,M
                    new_scalars(i)=this%scalar_fields(i)
                end do
                deallocate(this%scalar_fields)
            else 
                M=0
                allocate(new_scalars(1))
            end if

            new_scalars(M+1)%name=name
            new_scalars(M+1)%position=M+1
            this%scalar_fields=new_scalars
        end subroutine VELOCITY_CREATE_SCALAR

        subroutine VELOCITY_CREATE_VECTOR(this,name)
            Class(Velocity):: this
            integer:: i,M
            TYPE(VectorField),dimension(:),allocatable:: new_vectors
            character(len=*),intent(in):: name

            if (ALLOCATED(this%vector_fields)) then
                M=SIZE(this%vector_fields)
                allocate(new_vectors(M+1))
                do i=1,M
                    new_vectors(i)=this%vector_fields(i)
                end do
                deallocate(this%vector_fields)
            else
                M=0
                allocate(new_vectors(1))
            end if

            new_vectors(M+1)%name=name
            new_vectors(M+1)%position=M+1
            this%vector_fields=new_vectors
        end subroutine VELOCITY_CREATE_VECTOR



        function VELOCITY_GET_SCALAR(this,position,timestep)result(grid)
            Class(Velocity):: this
            integer:: position,timestep
            real(kind=c_double),dimension(:,:),allocatable:: grid
            
            SELECT CASE(timestep)
                CASE(1)
                    grid=this%scalar_fields(position)%grid_values1
                CASE(2)
                    grid=this%scalar_fields(position)%grid_values2
                CASE(3)
                    grid=this%scalar_fields(position)%grid_values3            
            END SELECT
        end function VELOCITY_GET_SCALAR

        function VELOCITY_GET_VECTOR(this,position,timestep,dimension)result(grid)
            Class(Velocity):: this
            integer:: position,timestep,dimension
            real(kind=c_double),dimension(:,:),allocatable:: grid
            SELECT CASE(timestep)
                CASE(1)
                    if (dimension==1) then
                        grid=this%vector_fields(position)%grid_x1
                    else 
                        grid=this%vector_fields(position)%grid_y1
                    end if
                CASE(2)
                    if (dimension==1) then
                        grid=this%vector_fields(position)%grid_x2
                    else 
                        grid=this%vector_fields(position)%grid_y2
                    end if
                CASE(3)
                    if (dimension==1) then
                        grid=this%vector_fields(position)%grid_x3
                    else 
                        grid=this%vector_fields(position)%grid_y3
                    end if          
            END SELECT            
        end function VELOCITY_GET_VECTOR

        function VELOCITY_EVALUATE(grid_values,xs,ys,evaluate_type)result(interp_values)
            character(len=*),intent(in):: evaluate_type
            real(kind=c_double),dimension(:),allocatable:: xs,ys,interp_values
            real(kind=c_double),dimension(:,:),allocatable,intent(in):: grid_values

            SELECT CASE(evaluate_type)
            CASE('nearest')
                interp_values=INTERP_NEAREST(grid_values,xs,ys)                                      
            CASE('linear')
                interp_values=INTERP_LINEAR(grid_values,xs,ys)
            CASE('cubic')
                interp_values=INTERP_CUBIC_FD(grid_values,xs,ys)

            END SELECT
        end function VELOCITY_EVALUATE

        function VELOCITY_FIND_SCALAR(this,name)result(position)
            Class(Velocity):: this
            character(len=*),intent(in):: name
            integer:: i,M,position
            M=SIZE(this%scalar_fields)
            do i=1,M
                if (this%scalar_fields(i)%name==name) then 
                    position=i
                end if
            end do
        end function VELOCITY_FIND_SCALAR

        function VELOCITY_FIND_VECTOR(this,name)result(position)
            Class(Velocity):: this
            character(len=*),intent(in):: name
            integer:: i,M,position
            M=SIZE(this%vector_fields)
            do i=1,M
                if (this%vector_fields(i)%name==name) then 
                    position=i
                end if
            end do
        end function VELOCITY_FIND_VECTOR

        subroutine VELOCITY_SET_SCALAR_VALUES(this,name,timestep,grid_values)
            Class(Velocity):: this
            integer:: timestep,position
            character(len=*),intent(in):: name
            real(kind=c_double),dimension(:,:),allocatable,intent(in):: grid_values

            position=this%find_scalar(name)
            SELECT CASE(timestep)
                CASE(1)
                    this%scalar_fields(position)%grid_values1=grid_values
                CASE(2)
                    this%scalar_fields(position)%grid_values2=grid_values
                CASE(3)
                    this%scalar_fields(position)%grid_values3=grid_values
            END SELECT
        end subroutine VELOCITY_SET_SCALAR_VALUES

        subroutine VELOCITY_SET_VECTOR_VALUES(this,name,timestep,gridx_values,gridy_values)
            Class(Velocity):: this
            integer:: timestep,position
            character(len=*),intent(in):: name
            real(kind=c_double),dimension(:,:),allocatable,intent(in):: gridx_values,gridy_values

            position=this%find_vector(name)
            SELECT CASE(timestep)
                CASE(1)
                    this%vector_fields(position)%grid_x1=gridx_values
                    this%vector_fields(position)%grid_y1=gridy_values
                CASE(2)
                    this%vector_fields(position)%grid_x2=gridx_values
                    this%vector_fields(position)%grid_y2=gridy_values
                CASE(3)
                    this%vector_fields(position)%grid_x3=gridx_values
                    this%vector_fields(position)%grid_y3=gridy_values
                END SELECT
        end subroutine VELOCITY_SET_VECTOR_VALUES

        function VELOCITY_GET_SCALAR_VALUES(this,name,timestep,evaluate_type,xs,ys)result(interp_values)
            Class(Velocity):: this
            character(len=*),intent(in):: name,evaluate_type
            integer:: position,timestep
            real(kind=c_double),dimension(:),allocatable,intent(in):: xs,ys
            real(kind=c_double),dimension(:),allocatable:: interp_values
            real(kind=c_double),dimension(:,:),allocatable:: grid_values
            
            position=VELOCITY_FIND_SCALAR(this,name)
            grid_values=VELOCITY_GET_SCALAR(this,position,timestep)
            interp_values=VELOCITY_EVALUATE(grid_values,xs,ys,evaluate_type)
        end function VELOCITY_GET_SCALAR_VALUES

        function VELOCITY_GET_VECTOR_VALUES(this,name,timestep,dimension,evaluate_type,xs,ys)result(interp_values)
            Class(Velocity):: this
            character(len=*),intent(in):: name,evaluate_type,dimension
            integer:: position,timestep
            real(kind=c_double),dimension(:),allocatable,intent(in):: xs,ys
            real(kind=c_double),dimension(:),allocatable:: interp_values
            real(kind=c_double),dimension(:,:),allocatable:: grid_values
            
            position=VELOCITY_FIND_VECTOR(this,name)
            SELECT CASE(dimension)
                CASE('x')
                    grid_values=VELOCITY_GET_VECTOR(this,position,timestep,1)
                CASE('y')
                    grid_values=VELOCITY_GET_VECTOR(this,position,timestep,2)
            END SELECT
            interp_values=VELOCITY_EVALUATE(grid_values,xs,ys,evaluate_type)
        end function VELOCITY_GET_VECTOR_VALUES

        subroutine VELOCITY_SAVE_SCALAR(this,name,path,fname)
            Class(Velocity):: this
            integer:: i,N
            character(len=*),intent(in):: path,name,fname
            real(kind=c_double),dimension(:,:),allocatable:: out

            out=VELOCITY_GET_SCALAR(this,VELOCITY_FIND_SCALAR(this,name),2)
            N=SIZE(out,dim=1)
            OPEN(UNIT=12, FILE=TRIM(ADJUSTL(path))//TRIM(ADJUSTL(fname))//".dat", ACTION="write", STATUS="replace")
            DO i=1,N
                WRITE(12,"(*(G0,:,','))") out(i,:)
            END DO
            CLOSE(UNIT=12)


        end subroutine VELOCITY_SAVE_SCALAR

        subroutine VELOCITY_SAVE_VECTOR(this,name,path,fname)
            Class(Velocity):: this
            integer:: i,N
            character(len=*),intent(in):: path,name,fname
            real(kind=c_double),dimension(:,:),allocatable:: out

            out=VELOCITY_GET_VECTOR(this,VELOCITY_FIND_VECTOR(this,name),2,1)
            N=SIZE(out,dim=1)
            OPEN(UNIT=13, FILE=TRIM(ADJUSTL(path))//"x"//TRIM(ADJUSTL(fname))//".dat", ACTION="write", STATUS="replace")
            DO i=1,N
                WRITE(13,"(*(G0,:,','))") out(i,:)
            END DO
            CLOSE(UNIT=13)

            out=VELOCITY_GET_VECTOR(this,VELOCITY_FIND_VECTOR(this,name),2,2)
            N=SIZE(out,dim=1)
            OPEN(UNIT=14, FILE=TRIM(ADJUSTL(path))//"y"//TRIM(ADJUSTL(fname))//".dat", ACTION="write", STATUS="replace")
            DO i=1,N
                WRITE(14,"(*(G0,:,','))") out(i,:)
            END DO
            CLOSE(UNIT=14)

        end subroutine VELOCITY_SAVE_VECTOR

end module velocitymodule