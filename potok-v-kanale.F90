program PrandtlEquation
use, intrinsic :: iso_fortran_env, only : dp => real64
   integer, parameter:: IO = 12 ! input-output unit
   integer :: NI, NJ, SMAX, i, j
   real(dp) :: L, H, U0, MU, Nu, R0, P0, G
   real(dp) :: dx, dy, EPS, pi
   real(dp), allocatable :: X_Node(:,:),Y_Node(:,:), dU_dy(:), delta(:)
   real(dp), allocatable :: U_n(:,:), V_n(:,:), P_n(:,:)
   pi = 4.0_dp*atan(1.0_dp)
   write(*,*) 'Read input file'
   L = 2.0_dp
   H = 0.1_dp
   NI = 51
   NJ = 101
   EPS = 1.0e-8_dp
   SMAX = 50
   U0 = 1.0_dp
   MU = 1.0_dp
   R0 = 1000.0_dp
   P0 = 1.0_dp
   
   allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
   allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

   !----------------- Node variables -----------------------------
   allocate(U_n(NI,NJ))  ! Velocity U
   allocate(V_n(NI,NJ))  ! Velocity V
   allocate(P_n(NI,NJ))  ! Pressure

   allocate(dU_dy(NI))
   allocate(delta(NI))

   !----------------- Coordinate of nodes ------------------------
   dx=L/(NI-1)
   dy=H/(NJ-1)

   do i=1,NI
      do j=1,NJ
         X_Node(I,J)=(I-1)*dx
         Y_Node(I,J)=(J-1)*dy
      end do
   end do

   !----------------- Parameters ------------------------

   NU=MU/R0

   write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
   write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
   write(*,*)'ReL= ', U0*L/NU

   !----------------- Initial fields -----------------------------

   do I=1,NI
      do J=1,NJ
         U_n(I,J)=U0
         V_n(I,J)=1.0e-5
         P_n(I,J)=P0
      end do
   end do

   !---------------- Solve Prandtl equations ---------------------

   write(*,*) 'Solve Prandtl equations'
   call prandtl(NI, SMAX, NJ, U_n, V_n, P_n, dy, dx, NU, EPS, P0)

   !----------------- Output data ------------------------------

   write(*,*) 'Output data'
   open(IO,FILE='Results.plt')
   call output_fields(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
   close(IO)

   write(*, *) "Post processed data"
   call post_process(U_n, NI, NJ, dU_dy, delta, U0, dy, dx, NU, MU, R0, IO, H)
   close(IO)
end program

!************************************************************************************************

subroutine prandtl(NI, SMAX, NJ, U_n, V_n, P_n, dy, dx, NU, EPS, P0)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer, intent(in) :: NI, SMAX, NJ
   real(dp), intent(in) :: NU, dy, dx, EPS, P0
   real(dp), dimension(NI, NJ) :: U_n, V_n, P_n
   real(dp), dimension(NJ) :: U_s, U_s1, U_ni, V_s, V_s1, D, F
   real(dp) :: U_err, V_err, P_err, P_s, P_s1, P_ni
   integer :: i, s
   P_s = P0
   do i=2, NI, 1
      U_ni = U_n(I-1, :)
      U_s = U_ni
      V_s = V_n(I-1, :)
      P_ni = P_s
      do s=1, SMAX, 1
         !write(*,*) 'Hello'
         call solver_type_matrix(U_s, U_ni, V_s, NJ, dy, dx, NU, D, F, P_ni)
         !write(*,*) 'Hello_1'
         call find_p(P_s1, V_s, U_ni, D, F, dx, dy, NJ)
         !write(*,*) P_s1
         call find_u_s1(U_s1, D, F, P_s1, NJ)
         call find_v_s1(V_s1, U_s1, U_ni, dx, dy, NJ)
         U_err = MAXVAL(ABS(U_s1 - U_s)) / (MAXVAL(U_s1))
         V_err = MAXVAL(ABS(V_s1 - V_s)) / (MAXVAL(V_s1))
         !P_err = abs(P_s1 - P_s)/P_s1
         U_s = U_s1
         V_s = V_s1
         P_s = P_s1
         if ((U_err < EPS).and.(V_err < EPS)) then
            exit
         end if
      end do
      !write(*,*) '--------------------------------------------------------------'
      print *, U_err, V_err, P_s1, i, s
      !write(*,*) '--------------------------------------------------------------'
      U_n(i, :) = U_s
      V_n(i, :) = V_s
      P_n(i, :) = P_s1
   end do

end subroutine

!************************************************************************************************

subroutine solver_type_matrix(U_s, U_ni, V_s, NJ, dy, dx, NU, D, F, P_ni)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer, intent(in) :: NJ
   integer :: j
   real(dp), dimension(NJ), intent(inout) :: D, F
   real(dp), dimension(NJ), intent(in) :: U_s, U_ni, V_s
   real(dp) :: A_obr(1:NJ, 1:NJ), gauss_matrix(1:NJ, 1:NJ+1)
   real(dp), intent(in) :: dy, dx, NU, P_ni
   !write(*,*) "Hello_1"
   gauss_matrix = 0.0_dp
   call BC_U(gauss_matrix, D, F, NJ)
   
   do j = 2, NJ-1, 1
      gauss_matrix(j, j-1) = -(V_s(j-1)/(2.0_dp*dy)) - (NU/dy**2.0_dp) !A
      gauss_matrix(j, j) = (U_s(j)/dx) + (2.0_dp*NU/(dy**2.0_dp))!B
      gauss_matrix(j, j+1) = (V_s(j+1)/(2.0_dp*dy)) - (NU/dy**2.0_dp)!C
      F(j) = 1.0_dp/dx
      D(j) = (U_ni(j)**2.0_dp+P_ni)/dx
   end do
   call obr_matrix_gauss(gauss_matrix, A_obr, NJ)
   D = matmul(A_obr, D)
   F = matmul(A_obr, F)
end subroutine

!************************************************************************************************

subroutine BC_U(gauss_matrix, D, F, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer, intent(in) :: NJ
   real(dp), dimension(NJ), intent(inout) :: D, F
   real(dp), intent(inout) :: gauss_matrix(1:NJ, 1:NJ+1)

   gauss_matrix(1, 1) = 1.0_dp
   gauss_matrix(1, 2) = 0.0_dp
   D(1) = 0.0_dp
   F(1) = 0.0_dp
   gauss_matrix(NJ, NJ) = 1.0_dp
   gauss_matrix(NJ, NJ-1) = 0.0_dp !A
   D(NJ) = 0.0_dp
   F(NJ) = 0.0_dp
end subroutine

subroutine obr_matrix_gauss(gauss_matrix, A_obr, N)
use, intrinsic :: iso_fortran_env, only : dp => real64
implicit none
real(dp) :: temp
integer, intent(in) :: N
integer :: i, j, k, element
real(dp) :: gauss_answer(1:N), matrix(1:N, 1:N+1)
real(dp), intent(inout) :: A_obr(1:N, 1:N), gauss_matrix(1:N, 1:N+1)
A_obr = 0.0_dp
do element=1, N
   matrix = gauss_matrix
   matrix(element, N+1) = 1.0_dp
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

   gauss_answer(N) = matrix(N, N+1)/matrix(N, N)
   do i=N-1, 1, -1
      temp = 0.0_dp
      do j=i+1, N
         temp = temp + matrix(i,j)*gauss_answer(j)
      end do
      temp = matrix(i,N+1) - temp
      gauss_answer(i) = temp/matrix(i,i)
   end do
   A_obr(:, element) = gauss_answer(:)
end do
!A_obr = transpose(A_obr)
end subroutine

!************************************************************************************************

subroutine find_p(P_s1, V_s, U_ni, D, F, dx, dy, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer :: j
   integer, intent(in) :: NJ
   real(dp), intent(in) :: dx, dy
   real(dp), intent(inout) :: P_s1
   real(dp), dimension(NJ), intent(in) :: D, F, V_s, U_ni
   real(dp) :: coeff_D, coeff_F, coeff_V, G
   coeff_D = 0.0_dp
   coeff_F = 0.0_dp
   do j = 2, NJ-1
      coeff_D = coeff_D + D(j)*dy
      coeff_F = coeff_F + F(j)*dy
   end do
   coeff_D = (coeff_D + 0.5_dp*dy*(D(1)+D(NJ)))/dx
   coeff_F = (coeff_F + 0.5_dp*dy*(F(1)+F(NJ)))/dx
   G = 0.0_dp
   do j = 2, NJ
      G = G + 0.5_dp*(U_ni(j)+U_ni(j-1))*dy
   end do
   coeff_V = -G/dx + V_s(1) - V_s(NJ) 
   P_s1 = (coeff_D + coeff_V)/coeff_F
end subroutine

!************************************************************************************************

subroutine find_u_s1(U_s1, D, F, P_s1, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer, intent(in) :: NJ
   real(dp), intent(in) :: P_s1
   real(dp), dimension(NJ), intent(in) :: D, F
   integer :: j
   real(dp), dimension(NJ), intent(inout) :: U_s1
   write(*,*) D(NJ-2), F(NJ-2), D(NJ-5), F(NJ-5), D(NJ-8), F(NJ-8), P_s1
   do j = 1, NJ
      U_s1(j) = D(j) - F(j)*P_s1
   end do
end subroutine

!************************************************************************************************

subroutine find_v_s1(V_s1, U_s1, U_ni, dx, dy, NJ)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer, intent(in) :: NJ
   real(dp), intent(in) :: dx, dy
   real(dp), dimension(NJ), intent(in) :: U_s1, U_ni
   integer :: j = 0
   real(dp), dimension(NJ), intent(inout) :: V_s1
   V_s1(1) = 0.0_dp
   do j = 2, NJ, 1
      V_s1(j) = V_s1(j-1) - (dy/2.0_dp)*((U_s1(j) - U_ni(j))/dx + (U_s1(j-1) - U_ni(j-1))/dx)
   end do

end subroutine

!************************************************************************************************

subroutine output_fields(IO,NI,NJ,X,Y,U,V,P)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none

   integer NI,NJ,IO
   real(dp), dimension(NI,NJ):: X,Y
   real(dp), dimension(NI,NJ):: U,V,P

   write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"'
   write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
   write(IO,'(100E25.16)') X(1:NI,1:NJ)
   write(IO,'(100E25.16)') Y(1:NI,1:NJ)
   write(IO,'(100E25.16)') U(1:NI,1:NJ)
   write(IO,'(100E25.16)') V(1:NI,1:NJ)
   write(IO,'(100E25.16)') P(1:NI,1:NJ)

end subroutine

subroutine post_process(U_n, NI, NJ, dU_dy, delta, U0, dy, dx, NU, MU, R0, IO, H)
use, intrinsic :: iso_fortran_env, only : dp => real64
   implicit none
   integer, intent(in) :: NI, NJ, IO
   real(dp), intent(in) :: U0, dy, dx, NU, R0, H
   real(dp), dimension(NI, NJ), intent(in) :: U_n
   integer :: i, j, k = 1, max_k = 1
   real(dp) :: Re_x, Cf, Cfs, MU
   real(dp), dimension(NI), intent(inout) :: dU_dy, delta
   character(len=1024) :: filename

   do i = 1, NI, 1
      do j = 1, NJ, 1
         if (abs(U_n(i, j) - U0) / U0 < 0.01_dp) then
            delta(i) = j * dy
            exit
         end if
      end do
      if (mod(i, 10) == 1) then
         write (filename, "(a, f3.1)") "U_by_y_", (i * dx)
         open(IO, FILE=filename)
         write(IO, "(a)") "U/U0, y/delta"
         do j = 1, NJ
            write(IO, "(f15.13,a,f15.13)")U_n(i, j) / U0,",",dy * j / delta(i)
         end do
         k = k + 1
      end if
      dU_dy(i) = (3.0_dp*U_n(i, 1) - 4.0_dp*U_n(i, 2) + U_n(i, 3))/(2.0_dp*dy)
   end do

   max_k = k

   open(IO, FILE="Post.plt")

   write(IO,*) 'VARIABLES = "dU_dy", "delta", "X", "Re_x", "Cf", "cf_analiz"'
   write(IO,*) 'ZONE I=',NI-1
   do i =  2, NI, 1
      Re_x = U0 * dx * i / NU
      Cf = -MU*dU_dy(i)/(0.5_dp*R0*U0**2.0_dp)
      Cfs = 12.0_dp/100_dp
      write(IO,'(100E25.16)') dU_dy(i), delta(i), dx * i, Re_x, Cf, Cfs
   end do

   close(IO)
end subroutine
!************************************************************************************************