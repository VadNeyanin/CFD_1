program kurs
use, intrinsic :: iso_fortran_env, only : dp => real64
integer :: NI, NJ, SMAX, i, j, j_bort_chanel, iter, max_iter, i_h
real(dp) :: radius_coeff, r_bort, max_diff, diff, diff_first, sum_s_1, Force_1, approx_1, approx_2, degree
real(dp) :: dr, dphi, dh, h_coeff, EPS, pi, LAMBDA, h_1, old_p, sum_s, s, Force, h_b, square, Force_func, p_0
real(dp), allocatable :: X(:,:), Y(:,:), delta(:), rad(:,:), phi(:, :), p_theory(:)
real(dp), allocatable :: P(:,:), h(:, :), r(:), A(:, :), B(:, :), C(:, :), D(:, :), E(:, :), F(:, :)

open(1,FILE='input.txt')
   read(1,*) NI
   read(1,*) NJ
   read(1,*) h_coeff
   read(1,*) LAMBDA
   read(1,*) radius_coeff
   read(1,*) r_bort
   read(1,*) max_iter
   read(1,*) eps
   read(1,*) degree
   read(1,*) i_h
   read(1,*) p_0
close(IO)

h_2 = 1.0_dp
h_1 = h_coeff
dr = (1.0_dp - radius_coeff)/real(NJ-1)
pi = 4.0_dp*atan(1.0_dp)

allocate(X(NI, NJ), Y(NI, NJ), r(NJ), h(NI, NJ), A(NI, NJ), B(NI, NJ), C(NI, NJ), D(NI, NJ), E(NI, NJ), F(NI, NJ))
allocate(p(NI, NJ), rad(NI, NJ), phi(NI, NJ), p_theory(NJ))
!Поиск i, при котором попадаем на границу между бортиком и каналом
j_bort_chanel = nint(r_bort/dr)+1
dh = (h_coeff-1.0_dp)/(NI-1)
h = h_2
do i=1, NI
   do j = j_bort_chanel+1, NJ - j_bort_chanel
      h(i, j) = h_1 - (i-1)*dh
   end do
end do

!Построение сетки
do j=1, NJ
   r(j) = radius_coeff + real(j-1)*dr
   dphi = degree*pi/(180.0_dp*real(NI-1))
   do i=1, NI
      X(I, J) = r(j)*cos(real(i-1)*dphi)
      Y(I, J) = r(j)*sin(real(i-1)*dphi)
      phi(I, :) = real(i-1)*dphi
   end do
end do

write(*,*) j_bort_chanel, dr, dphi
!Задание начальных условий
P = p_0
call init_coeff(A, B, C, D, E, F, LAMBDA, h, dh, r, dr, dphi, j_bort_chanel, NI, NJ)
open(1, file="Res_history.plt", status="old")
write(1,*) 'VARIABLES = "Iterattion", "Res_p", "Force"'
write(1,*) 'ZONE I=', max_iter
write(*,*) "   P_res_i/P_res_1             ", "   Force                  ", "Iterattion"
do iter = 1, max_iter

   !-------------------------------------------------------------------------------------------------------------------!
   sum_s_1 = 0.0_dp
   Force_1 = 0.0_dp
   do j=1, NJ-1
      s = square(r(j+1), r(j), dphi)
      sum_s_1 = sum_s_1+s
      Force_1 = Force_1 + 0.25_dp*s*(p(1, j)+p(2, j+1)+p(1,j+1)+p(2,j))
   end do
   do i=2, NI-1
      s = square(r(2), r(1), dphi)
      sum_s_1 = sum_s_1+s
      Force_1 = Force_1 + 0.25_dp*s*(p(i, 1)+p(i+1, 2)+p(i,2)+p(i+1,1))
   end do
   sum_s = sum_s_1
   Force = Force_1

   max_diff = 0.0_dp
   do i = 2, NI-1
      do j = 2, j_bort_chanel-1
         old_p = p(i, j)
         p(i, j) = (A(i, j) * p(i-1, j) + B(i, j) * p(i+1, j) &
         + D(i, j) * p(i, j-1) + E(i, j) * p(i, j+1) - F(i, j)) / C(i, j)
         diff = abs(p(i, j) - old_p)
         Force = Force + Force_func(p, i, j, square(r(j+1), r(j), dphi), NI, NJ, p_0)
         if (diff > max_diff) then
            max_diff = diff
         end if
      end do
      old_p = p(i, j_bort_chanel)
      p(i, j_bort_chanel) = approx_1(p, h, i, j_bort_chanel, NI, NJ)
      Force = Force + Force_func(p, i, j_bort_chanel, square(r(j_bort_chanel+1), r(j_bort_chanel), dphi), NI, NJ, p_0)
      diff = abs(p(i, j_bort_chanel) - old_p)
      if (diff > max_diff) then
            max_diff = diff
         end if
   end do
   do i = 2, NI-1
      do j = j_bort_chanel+1, NJ-j_bort_chanel
         old_p = p(i, j)
         p(i, j) = (A(i, j) * p(i-1, j) + B(i, j) * p(i+1, j) &
         + D(i, j) * p(i, j-1) + E(i, j) * p(i, j+1) - F(i, j)) / C(i, j)
         diff = abs(p(i, j) - old_p)
         Force = Force + Force_func(p, i, j, square(r(j+1), r(j), dphi), NI, NJ, p_0)
         if (diff > max_diff) then
            max_diff = diff
         end if
      end do
      old_p = p(i, NJ-j_bort_chanel+1)
      p(i, NJ-j_bort_chanel+1) = approx_1(p, h, i, NJ-j_bort_chanel+1, NI, NJ)
      Force = Force+Force_func(p, i, NJ-j_bort_chanel, square(r(NJ-j_bort_chanel+1), r(NJ-j_bort_chanel), dphi), NI, NJ, p_0)
      diff = abs(p(i, NJ-j_bort_chanel+1) - old_p)
      if (diff > max_diff) then
         max_diff = diff
      end if
   end do
   do i = 2, NI-1
      do j = NJ-j_bort_chanel+2, NJ-1
         old_p = p(i, j)
         p(i, j) = (A(i, j) * p(i-1, j) + B(i, j) * p(i+1, j) &
         + D(i, j) * p(i, j-1) + E(i, j) * p(i, j+1) - F(i, j)) / C(i, j)
         diff = abs(p(i, j) - old_p)

         Force = Force + Force_func(p, i, j, square(r(j+1), r(j), dphi), NI, NJ, p_0)

         if (diff > max_diff) then
            max_diff = diff
         end if
      end do
   end do
   ! Проверка сходимости
   if (iter == 1) then
      diff_first = max_diff
   end if
   write(*,*) max_diff/diff_first, Force, iter
   write(1, '(I6,100E25.16,100E25.16)') iter, max_diff/diff_first, Force
   if (max_diff/diff_first < eps) then
      exit
   end if 
end do
close(1)

sum_s = 0.0_dp
do i=1, NI-1
   do j=1, NJ-1
      sum_s = sum_s+square(r(j+1), r(j), dphi)
   end do
end do
write(*,*) 'Square =', sum_s
do i=1, NI
   do j=1, NJ
      rad(i, j) = r(j)
   end do
end do
do j=j_bort_chanel+1, NJ-j_bort_chanel
   p_theory(j) = 1.0_dp + LAMBDA*sin(degree*pi/180_dp)*r(j)**2.0_dp*(h_1-h(i_h, j))*(h(i_h, j)-h_2)&
   /(h(i_h, j)**2.0_dp*(h_1**2.0_dp-h_2**2.0_dp))
end do
open(1, file='results.plt', status='old')
   write(1,*) 'VARIABLES = "X", "Y", "p", "h", "r", "phi"'
   write(1,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
   write(1,'(100E25.16)') X(1:NI,1:NJ)
   write(1,'(100E25.16)') Y(1:NI,1:NJ)
   write(1,'(100E25.16)') p(1:NI,1:NJ)
   write(1,'(100E25.16)') h(1:NI,1:NJ)
   write(1,'(100E25.16)') rad(1:NI,1:NJ)
   write(1,'(100E25.16)') phi(1:NI,1:NJ)*180.0_dp/pi
close(1)
open(1, file='p_theory.plt', status='old')
   write(1,*) 'VARIABLES = "r", "p_prog"'
   write(1,*) 'ZONE J=',NJ-2*j_bort_chanel
   do j=j_bort_chanel+1, NJ-j_bort_chanel
      write(1, '(100E25.16, 100E25.16)') r(j), p_theory(j)
   end do
close(1)
open(1, file='p_for_theory.plt', status='old')
   write(1,*) 'VARIABLES = "r", "p_prog"'
   write(1,*) 'ZONE J=',NJ
   do j=1, NJ
      write(1, '(100E25.16, 100E25.16)') r(j), p(i_h, j)
   end do
close(1)
end program

real(dp) function approx_1(p, h, i, j_bort_chanel, NI, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
integer, intent(in) :: NI, NJ, i, j_bort_chanel
real(dp), intent(in) :: p(NI, NJ), h(NI, NJ)
approx_1 = (h(i, j_bort_chanel-1)**3.0_dp*p(i, j_bort_chanel-1)+h(i, j_bort_chanel+1)**3.0_dp &
*p(i, j_bort_chanel+1))/(h(i, j_bort_chanel-1)**3.0_dp+h(i, j_bort_chanel+1)**3.0_dp)
end function

real(dp) function approx_2(p, h, i, j_bort_chanel, NI, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
integer, intent(in) :: NI, NJ, i, j_bort_chanel
real(dp), intent(in) :: p(NI, NJ), h(NI, NJ)
real(dp) :: h_b
h_b = (h(i, j_bort_chanel+1)/h(i, j_bort_chanel-1))**3.0_dp
approx_2 = (4.0_dp*(p(i, j_bort_chanel-1)+h_b*p(i, j_bort_chanel+1))&
-p(i, j_bort_chanel-2)-h_b*p(i, j_bort_chanel+2))/(3.0_dp*(1.0_dp+h_b))
end function

real(dp) function square(r_2, r_1, dphi)
use, intrinsic :: iso_fortran_env, only : dp => real64
real(dp), intent(in) :: r_2, r_1, dphi
square = 0.5_dp*(r_2**2.0-r_1**2.0)*dphi 
end function

real(dp) function Force_func(p, i, j, s, NI, NJ, p_0)
use, intrinsic :: iso_fortran_env, only : dp => real64
integer, intent(in) :: NI, NJ, i, j
real(dp), intent(in) :: p(NI, NJ)
real(dp), intent(in) :: s, p_0
Force_func = (0.25_dp*(p(i, j)+p(i, j+1)+p(i+1, j+1)+p(i+1, j))-p_0)*s
end function

subroutine init_coeff(A, B, C, D, E, F, LAMBDA, h, dh, r, dr, dphi, j_bort_chanel, NI, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
integer, intent(in) :: NI, NJ, j_bort_chanel
real(dp), intent(in) :: LAMBDA, dh, dr, dphi
real(dp), intent(inout) :: A(NI, NJ), B(NI, NJ), C(NI, NJ), D(NI, NJ), E(NI, NJ), F(NI, NJ)
real(dp), intent(inout) :: h(NI, NJ), r(NJ)
real(dp) :: r_j_plus_half, r_j_minus_half, h_i_minus_half, h_i_plus_half
integer :: i, j


A = 1.0_dp
B = 1.0_dp
C = 1.0_dp
D = 1.0_dp
E = 1.0_dp
F = 1.0_dp


do j=2, NJ-1
   do i=2, NI-1
      h_i_minus_half = (h(i, j)+h(i-1, j))/2.0_dp
      h_i_plus_half = (h(i, j)+h(i+1, j))/2.0_dp
      r_j_minus_half = (r(j)+r(j-1))/2.0_dp
      r_j_plus_half = (r(j)+r(j+1))/2.0_dp

      A(i, j) = h_i_minus_half**3.0_dp/(r(j)*dphi)**2.0_dp
      B(i, j) = h_i_plus_half**3.0_dp/(r(j)*dphi)**2.0_dp
      D(i, j) = h(i, j)**3.0_dp*r_j_minus_half/(r(j)*dr**2.0_dp)
      E(i, j) = h(i, j)**3.0_dp*r_j_plus_half/(r(j)*dr**2.0_dp)
      F(i, j) = LAMBDA*(h(i+1, j)-h(i-1, j))/(2.0_dp*dphi)
      C(i, j) = A(i, j) + B(i, j) + D(i, j) + E(i, j)
   end do
end do

end subroutine