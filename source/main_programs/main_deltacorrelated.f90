program main
    use, intrinsic :: iso_c_binding
    use passivetracer
    use deltacorrelated
    implicit none

    integer:: i,N,j
    character(len=:),allocatable:: path,path2,path3,path4
    character(len=100):: fname,num
    real(kind=c_double):: dt
    TYPE(DeltaCorrelatedField):: u1
    TYPE(PassiveTracers):: passive_tracers
    integer :: values(1:8),k
    integer,dimension(:),allocatable :: seed
    real(kind=c_double):: r

    call date_and_time(values=values)

    call random_seed(size=k)
    allocate(seed(1:k))
    seed(:)= values(8)
    call random_seed(put=seed)
    call random_number(r)
    print*, r  

    dt=REAL(0.01,kind=c_double)
    N=500000
    path='./vel_params.txt'
    path2="./data/delta"

    call u1%init(dt,path)
    print*, "initialising passive tracers"
    call passive_tracers%init(dt,500)



    do i=1,1000
        call u1%update()
    end do

    do i=1,N

        call u1%update()
        call passive_tracers%update(PASSIVE_UPDATE_FUNCTION,u1%vel)

        passive_tracers%particle_manager%attribute_values(:,passive_tracers%particle_manager%get_pos('x'))=&
        &MODULO(passive_tracers%particle_manager%attribute_values(:,passive_tracers%particle_manager%get_pos('x')),1.0)
        passive_tracers%particle_manager%attribute_values(:,passive_tracers%particle_manager%get_pos('y'))=&
        &MODULO(passive_tracers%particle_manager%attribute_values(:,passive_tracers%particle_manager%get_pos('y')),1.0)
        if (MODULO(i,1000)==0) then
            write(fname,'(A3,I0)') 'div',i
            call u1%vel%save_scalar('divergence',path3,fname)
            write(fname,'(A3,I0)') 'vel',i
            call u1%vel%save_vector('velocity',path3,fname)
            write(fname,'(A3,I0)') 'pos',i
            call passive_tracers%particle_manager%save(path3,fname)
            print*, i
        end if
    end do
end program main
