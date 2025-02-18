program QGtest

use, intrinsic :: iso_c_binding
use loadQG
use inertialslowmanifold
implicit none


integer:: N,i,j
real(kind=c_double):: dt,timescale,basinscale,H1,H2,R_d,beta,f,eps
type(QGfield):: QG1
type(InertialTracers):: inertial_tracers
character(len=100):: fname

N=4096
dt=0.01
basinscale=3600e5
timescale=basinscale/N
H1=1e5
H2=3e5
R_d=25e5
f=6.28
beta=6.22
eps=0.01
print*, "working"

call QG1%init(dt,timescale,N,basinscale,R_d,H1,H2)

print*, "working 2"

call inertial_tracers%init(dt,2000,f,eps,beta)

print*, "working 3"

do i=0,26200

    call QG1%update() 
    call inertial_tracers%update(PASSIVE_UPDATE_FUNCTION,QG1%vel)
    write(fname,'(A3,I0)') 'pos',i

    inertial_tracers%particle_manager%attribute_values(:,inertial_tracers%particle_manager%get_pos('x'))=&
    &MODULO(inertial_tracers%particle_manager%attribute_values(:,inertial_tracers%particle_manager%get_pos('x')),1.0)
    inertial_tracers%particle_manager%attribute_values(:,inertial_tracers%particle_manager%get_pos('y'))=&
    &MODULO(inertial_tracers%particle_manager%attribute_values(:,inertial_tracers%particle_manager%get_pos('y')),1.0)
    
    if (MODULO(i,1000)==0) then
        call inertial_tracers%particle_manager%save('./data/inertialQG/',fname)
        print*, "saved"
        print*, i
    end if
end do



end program QGtest
