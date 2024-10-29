program power
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
real(dp), ALLOCATABLE :: u(:), un(:), u_n_1_2(:), y(:), A(:), B(:), C(:), D(:)
real(dp) :: u0, y_max, time_max, dltY, dltTime, v, VNM, approx, new_time, buffer, etha, time
integer :: n_time, n_y, i, n, scheme, j

!Задание констант
u0 = 5.0_dp
y_max = 4.0_dp
time_max = 1.0_dp
n_y = 301
v = 1.0_dp
VNM = 0.01_dp!0.5, 0.01, 0.17
scheme = 3
allocate(A(2:n_y-1), B(2:n_y-1), C(2:n_y-1), D(2:n_y-1))
allocate(u(1:n_y), un(1:n_y), u_n_1_2(0:n_y+1), y(n_y))
dltY = y_max/(real(n_y)-1.0_dp)
if (scheme == 1) then
    dltTime = VNM*dltY**2.0_dp/v
else if (scheme == 2) then
    dltTime = dltY**2.0_dp/(2.0_dp*v)
    VNM = v*dltTime/dltY**2.0_dp
else
    dltTime = dltY**2.0_dp/(6.0_dp*v)
    VNM = v*dltTime/dltY**2.0_dp
end if
n_time = time_max/dltTime
write(*,*) n_time, dltTime, VNM
y(1) = 0.0_dp
do i = 2, n_y-1
    y(i) = y(i-1)+dltY
end do
y(n_y) = y_max
!НУ
u(:) = 0.0_dp
un(:) = 0.0_dp
u(1) = u0
u_n_1_2 = 0.0_dp
!Подготовка к прогонке
call InitValue(u, n_y)
call BoundValue(u, n_y, u0, 1)
if (scheme == 1) then
do j = 1, n_time
    do i = 2, n_y - 1
        A(i) = -VNM/2.0_dp
        C(i) = A(i)
        B(i) = 1.0_dp + VNM
        D(i) = -u(i)
    end do
    !Прогонкой найдем слой 1/2
    call Progonka(n_y, A, B, C, D, u_n_1_2)
    !Граничные условия для вспом временного слоя
    u_n_1_2(1) = u0
    u_n_1_2(n_y) = approx(u_n_1_2(n_y-1), u_n_1_2(n_y-2))
    u_n_1_2(0) = approx(u_n_1_2(1), u_n_1_2(2))
    u_n_1_2(n_y+1) = approx(u_n_1_2(n_y), u_n_1_2(n_y-1))
    !Следующий временной слой
    do i = 2, n_y-1
        un(i) = new_time(VNM, u(i), u_n_1_2(i+2), u_n_1_2(i+1), u_n_1_2(i), u_n_1_2(i-1), u_n_1_2(i-2))
    end do
    !Граничные условия для нового временного слоя
    call BoundValue(un, n_y, u0, 1)
    !Обновляем массив
    u(:) = un(:)
end do
open(1, file = 'U_y_auto_predict.plt', status = 'old')
write(1,*) 'VARIABLES = "η", "φ_ind", "φ_th"'
write(1,*) 'ZONE I=',n_y
    do i=1, n_y
        etha = y(i)/(2.*sqrt(v*time_max))
        write(1,'(100E25.16)') etha, u(i)/u0, erfc(y(i)/(2.*sqrt(v*time_max)))
    end do
close(1)
open(1, file = 'U_y_physical_predict.plt', status = 'old')
write(1,*) 'VARIABLES = "y, м", "U_ind, м/с", "U_th, м/с"'
write(1,*) 'ZONE I=',n_y
    do i=1, n_y
        write(1,'(100E25.16)') y(i), u(i), u0*erfc(y(i)/(2.*sqrt(v*time_max)))
    end do
close(1)
else if (scheme > 1) then
call BoundValue(u, n_y, u0, 1)
do n = 1, n_time
    do i = 2, n_y - 1
        un(i) = (1.0_dp-2.0_dp*VNM)*u(i)+VNM*(u(i-1) + u(i+1))
    end do
    call BoundValue(un, n_y, u0, 1)
    u = un
end do
open(1, file = 'U_y_auto_central.plt', status = 'old')
write(1,*) 'VARIABLES = "η", "φ_central", "φ_th"'
write(1,*) 'ZONE I=',n_y
    do i=1, n_y
        etha = y(i)/(2.*sqrt(v*time_max))
        write(1,'(100E25.16)') etha, un(i)/u0, erfc(y(i)/(2.*sqrt(v*time_max)))
    end do
close(1)
open(1, file = 'U_y_physical_central.plt', status = 'old')
write(1,*) 'VARIABLES = "y, м", "U_central, м/с", "U_th, м/с"'
write(1,*) 'ZONE I=',n_y
    do i=1, n_y
        write(1,'(100E25.16)') y(i), un(i), u0*erfc(y(i)/(2.*sqrt(v*time_max)))
    end do
close(1)
end if
end program power

SUBROUTINE InitValue(U, NX)
IMPLICIT NONE
integer i
integer, intent(in) :: NX
real(8), dimension(1:NX), intent(out) :: U
do i = 1, NX
U(i) = 0.0D0
end do
END SUBROUTINE

subroutine BoundValue(u, n_y, u0, k)
use, intrinsic :: iso_fortran_env, only : dp => real64
integer, intent(in) :: n_y, k
real(dp), intent(in) :: u0
real(dp), dimension(1:NX), intent(out) ::u
u(1) = u0
if (k == 2) then
    u(n_y) = (4.0_dp*u(n_y-1) - u(n_y-2))/3.0_dp
else if (k == 1) then
    u(n_y) = (48.0_dp*u(n_y-1)-36.0_dp*u(n_y-2)+16.0_dp*u(n_y-3)-3.0_dp*u(n_y-4))/25.0_dp
end if
END SUBROUTINE

!Метод прогонки; F - вектор решения
    !A, B, C, D – векторы прогоночных коэффициентов
    !Im – число узлов сетки (нумерация начинается с 1)
    subroutine Progonka(Im,A,B,C,D,F)
    use, intrinsic :: iso_fortran_env, only : dp => real64
    implicit none
    integer, intent(in):: Im
    real(dp), dimension(2:Im-1), intent(in) :: A, B, C, D
    real(dp), dimension(0:Im+1), intent(out) :: F
    real(dp), dimension(2:Im-1):: alpha, beta
    real(dp) :: k0
    integer :: i
!Прямой ход
    alpha(2) = -A(2) / B(2)
    beta(2) = -D(2) / B(2)
    do i = 3, Im-2
        k0 = B(i) + C(i)*alpha(i-1)
        alpha(i) = -A(i)/k0
        beta(i) = -(D(i) + C(i)*beta(i-1))/k0
    end do
    F(Im-1) = -(D(Im-1) + C(Im-1)*beta(Im-2))/(B(Im-1) + C(Im-1)*alpha(Im-2))
    !write(*,*) C(Im-1)*beta(Im-2)
    do i = (Im-2), 2, -1
        F(i) = alpha(i)*F(i+1) + beta(i)
    end do
    !write(*,*) F
    end subroutine

real(dp) function approx(u_2, u_3)
use, intrinsic :: iso_fortran_env, only : dp => real64
real(dp) :: u_2, u_3
approx = (2.0_dp*u_2-u_3)
end function

real(dp) function new_time(VNM, u, u_2, u_1, ui, u__1, u__2)
use, intrinsic :: iso_fortran_env, only : dp => real64
real(dp) :: u, ui, u_1, u_2, u__1, u__2, VNM, a
new_time = u + (VNM/12.0_dp)*(-u_2+16.0_dp*u_1-30.0_dp*ui + 16.0_dp*u__1-u__2)
end function