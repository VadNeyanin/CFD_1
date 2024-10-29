program yacheika
use, intrinsic :: iso_fortran_env, only : dp => real64
include 'omp_lib.h'
integer :: NX, NT, i, j, m, time_count, SITUATION
real(dp), allocatable :: u(:), un(:), x(:), gauss_matrix(:, :), gauss_answer(:)
real(dp), allocatable :: u_potok(:), un_potok(:), u_th(:)
real(dp) :: L, h, CFL, dt, Time, c, C0, C1, pi, k, coeff
!Задание констант
L = 1.0_dp
m = 12
C0 = 2.0_dp
C1 = 1.0_dp
NX = 101
NT = 25
CFL = 0.0001_dp
SITUATION = 1
allocate(u(1:NX), un(1:NX), x(1:NX), gauss_matrix(1:NX, 1:NX+1), gauss_answer(1:NX))
allocate(u_potok(1:NX), un_potok(1:NX), u_th(1:NX))
!Подготовка
pi = 4.0_dp*atan(1.0_dp)
k = m*pi/L
h = L/real(NX-1)
if (SITUATION == 0) then
	dt = CFL*h/0.6_dp
else
	dt = CFL*h/0.6_dp
end if
Time = dt*NT
write(*,*) 'Time =', Time, 'dt =', dt
x(1) = 0.0_dp
do i=2, NX
	x(i) = x(i-1) + h
end do
u(:) = 0.0_dp
un(:) = 0.0_dp
t = 0.0_dp
!Рассчитаем теоретические значения и зададим начальные условия
if (SITUATION == 0) then
	u_th(1) = C0+C1*sin(k*(x(1)))
	do i=2, NX
		u_th(i) = C0+C1*sin(k*(x(i)))
		call InitValue_0(x, u, NX, k)
		call InitValue_0(x, u_potok, NX, k)
	end do
else
	call InitValue_1(x, u, NX)
	call InitValue_1(x, u_potok, NX)
	do i = 1, NX
		if (x(i)-0.5_dp*Time*(0.2_dp+0.6_dp) .lt. 0.5_dp) then
			u_th(i) = 0.6_dp
		else
			u_th(i) = 0.2_dp
		end if
	end do
end if
!Теперь расчет по индивидуальной сетке
if (SITUATION == 0) then
	call BoundValue_0(gauss_answer, NX)!Не нужны
	call BoundValue_0(un_potok, NX)
else
	call BoundValue_1(gauss_answer, NX)
	call BoundValue_1(un_potok, NX)
end if
coeff = dt/h
do time_count=1, NT
	gauss_matrix(1, 1) = 1.0_dp+coeff*u(1)
	gauss_matrix(1, NX) = 1.0_dp-coeff*u(NX)
	if (u_potok(1) >= 0.0_dp) then
			un_potok(1) = u_potok(1)-dt/h*0.5_dp*(u_potok(1)**2.0_dp-u_potok(NX)**2.0_dp)
		else
			un_potok(1) = u_potok(1)-dt/h*0.5_dp*(u_potok(2)**2.0_dp-u_potok(1)**2.0_dp)
		end if
	gauss_matrix(1, NX+1) = u(NX)+u(1)
	do i=2, NX
		gauss_matrix(i, i) = 1.0_dp+coeff*u(i)
		gauss_matrix(i, i-1) = 1.0_dp-coeff*u(i-1)
		gauss_matrix(i, NX+1) = u(i-1)+u(i)
		if (u_potok(i) >= 0.0_dp) then
			un_potok(i) = u_potok(i)-dt/h*0.5_dp*(u_potok(i)**2.0_dp-u_potok(i-1)**2.0_dp)
		else
			un_potok(i) = u_potok(i)-dt/h*0.5_dp*(u_potok(i+1)**2.0_dp-u_potok(i)**2.0_dp)
		end if
	end do
	call gauss(gauss_matrix, gauss_answer, NX)
	u = gauss_answer
	do i = 1, NX
	if (x(i) .lt. 0.45_dp+0.6*time_count*dt) then
		u(i) = 0.6_dp
	end if
	end do
	u_potok(2:NX-1) = un_potok(2:NX-1)
end do

open(1, file = 'yacheika.plt', status = 'old')
write(1,*) 'VARIABLES = "x, м", "U, м/с", "U_potok, м/с", "U_th, м/с"'
write(1,*) 'ZONE I=',NX
do i=1, NX
	write(1,'(100E25.16)') x(i), u(i), u_potok(i), u_th(i)
end do
close(1)
end program

subroutine InitValue_0(x, u, NX, k)
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
integer, intent(in) :: NX
real(dp), intent(in) :: k
integer :: i
real(dp), intent(inout) :: u(1:NX), x(1:NX)
u(1) = 2.0_dp+sin(k*x(1))
do i = 2, NX
	u(i) = 2.0_dp+sin(k*x(i))
end do
end subroutine

subroutine BoundValue_0(u, NX)
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
integer, intent(in) :: NX
real(dp) :: u(1:NX)
u(1) = u(NX)
end subroutine

subroutine InitValue_1(x, u, NX)
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
integer, intent(in) :: NX
integer :: i
real(dp), intent(inout) :: u(1:NX), x(1:NX)
do i = 1, NX
	if (x(i) .lt. 0.5_dp) then
		u(i) = 0.6_dp
	else
		u(i) = 0.2_dp
	end if
end do
u(NX+1) = 0.2_dp
end subroutine

subroutine BoundValue_1(u, NX)
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
integer, intent(in) :: NX
real(dp) :: u(1:NX)
u(1)=0.6_dp
u(NX)=0.2_dp
end subroutine

subroutine gauss(gauss_matrix, gauss_answer, N)
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
real(dp) :: temp, e
integer, intent(in) :: N
integer :: i, j, k
real(dp), intent(inout) :: gauss_answer(N)
real(dp) :: matrix(N, N+1)
real(dp), intent(in) :: gauss_matrix(N, N+1)
matrix = gauss_matrix
do i=1, N
    if (matrix(i,i) == 0) then
        print*,'input error (gauss_matrix(i,i)=0)'
        stop
    end if
end do

do i=1, N
	do j=i+1, N
	temp = (matrix(j,i))/(matrix(i,i))
		do k=1, N+1
			matrix(j,k) = matrix(j,k)-matrix(i,k)*temp
		end do
	end do
end do

!обратный ход

gauss_answer(N) = matrix(N, N+1)/matrix(N, N)
do i=N-1, 1, -1
	temp = 0.0_dp
	do j=i+1, N
    	temp = temp + matrix(i,j)*gauss_answer(j)
	end do
	temp = matrix(i,N+1) - temp
	gauss_answer(i) = temp/matrix(i,i)
end do
!ответ
!do i=1, n
!    print*,"gauss_answer(",i,")=",gauss_answer(i)
!end do
end subroutine