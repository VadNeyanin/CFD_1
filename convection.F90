program yacheika
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
integer :: NX, NT, i, j, m, time_count, SITUATION
real(dp), allocatable :: u(:), un(:), x(:), gauss_matrix(:, :), gauss_answer(:), yacheika_matrix(:)
real(dp), allocatable :: u_potok(:), un_potok(:), u_chisl(:), u_pr_chisl(:)
real(dp) :: x_th(0:1000), u_th(0:1000)
real(dp) :: L, h, CFL, dt, t, Time, p, c, C0, C1, pi, k, G, FE, alpha, beta

!Задание констант
L = 1.0_dp
m = 12
c = 0.5_dp
C0 = 0.0_dp
C1 = 1.0_dp
NX = 301
NT = 180
CFL = 0.75_dp
SITUATION = 1
allocate(u(1:NX), un(1:NX), x(1:NX), gauss_matrix(1:NX, 1:NX+1), gauss_answer(1:NX), yacheika_matrix(1:NX))
allocate(u_potok(1:NX), un_potok(1:NX), u_chisl(1:NX), u_pr_chisl(1:NX))
!Подготовка
pi = 4.0_dp*atan(1.0_dp)
k = m*pi/L
h = L/real(NX-1)
dt = CFL*h/c
Time = dt*NT
write(*,*) Time, dt, h
x(1) = 0.0_dp
do i=2, NX
	x(i) = x(i-1) + h
end do
u(:) = 0.0_dp
un(:) = 0.0_dp
t = 0.0_dp
!Рассчитаем теоретические значения и зададим начальные условия
if (SITUATION == 0) then
	u_th(1) = C0+C1*sin(k*(x(1)-c*Time))
	do i=2, NX
		u_th(i) = C0+C1*sin(k*(x(i)-c*Time))
		call InitValue_0(x, u, NX, k)
		call InitValue_0(x, u_potok, NX, k)
	end do
else
	do i=1, NX
		if ((x(i)-c*Time) .lt. 0.5_dp) then
			u_th(i) = 0.6_dp
		else
			u_th(i) = 0.2_dp
		end if
	end do
	call InitValue_1(x, u, NX)
	call InitValue_1(x, u_potok, NX)
end if
!Теперь расчет по индивидуальной сетке
alpha = (0.5_dp/dt)+(0.5_dp*c/h)
beta = (0.5_dp/dt)-(0.5_dp*c/h)
do i=2, NX
	gauss_matrix(i, i) = alpha
	gauss_matrix(i, i-1) = beta
end do
gauss_matrix(1, 1) = alpha
gauss_matrix(1, NX) = beta
if (SITUATION == 0) then
		call BoundValue_0(gauss_answer, NX)!Не нужны
		call BoundValue_0(un_potok, NX)
	else
		call BoundValue_1(gauss_answer, NX)
		call BoundValue_1(un_potok, NX)
	end if
do time_count=1, NT
	
	gauss_matrix(1, NX+1) = 0.5_dp*((u(1)+u(NX))/dt - (c*(u(1)-u(NX))/h))
	do i=2, NX
		gauss_matrix(i, NX+1) = 0.5_dp*((u(i)+u(i-1))/dt - (c*(u(i)-u(i-1))/h))
		un_potok(i) = u_potok(i)-CFL*(u_potok(i)-u_potok(i-1))
	end do
	call gauss(gauss_matrix, gauss_answer, NX)
	u_potok = un_potok
	u = gauss_answer
	do i = 1, NX
		if (x(i) .lt. 0.3_dp+c*time_count*dt) then
			u(i) = 0.6_dp
		end if
	end do
end do
do i = 1, NX
	u_chisl(i) = sin(12.0_dp*pi/L*(x(i)-0.025_dp)-0.00054_dp)
	u_pr_chisl(i) = 0.98539868*sin(12.0_dp*pi/L*(x(i)-0.025_dp)-0.00031_dp)
end do
open(1, file = 'yacheika.plt', status = 'old')
write(1,*) 'VARIABLES = "x, м", "U, м/c", "U_prot, м/с", "U_th, м/с", "U_pr_chisl", "U_chisl"'
write(1,*) 'ZONE I=',NX
do i=1, NX
	write(1,'(100E25.16)') x(i), u(i), u_potok(i), u_th(i), u_pr_chisl(i), u_chisl(i)
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
u(1) = sin(k*x(1))
do i = 2, NX
	u(i) = sin(k*x(i))
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