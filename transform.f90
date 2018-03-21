module transform_module
implicit none
private
public oned_fourier_transform1,oned_fourier_transform2, twod_fourier_transform,twod_backward_fourier_transform,&
     & oned_backward_fourier_transform1, oned_backward_fourier_transform2, oned_sine_transform2, oned_inverse_sine_transform2,&
     & mode_evolution_analysis3
contains

  subroutine oned_fourier_transform1(s,s_fft,m,n) !calculating 1d DFT of s(:,:) along the first dimension
    use, intrinsic :: iso_c_binding
    use precision,only:p_
    implicit none
    include 'fftw3.f03'  
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    type(C_PTR) :: plan
    integer:: j


    in=s
    plan = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
    do j=0,n-1
       call fftw_execute_dft(plan, in(:,j), out(:,j))     !Fourier transformation along the first dimension
    enddo
    call fftw_destroy_plan(plan)  
    s_fft=out
  end subroutine oned_fourier_transform1

  subroutine oned_fourier_transform2(s,s_fft,m,n) !calculating 1d DFT of s(:,:) along the second dimension
    use, intrinsic :: iso_c_binding
    use precision,only:p_
    implicit none
    include 'fftw3.f03'  
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    type(C_PTR) :: plan
    integer:: i

    !Fourier transformation along the second dimension
    in=s
    plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan)  
    s_fft=out
  end subroutine oned_fourier_transform2


  subroutine twod_fourier_transform(s,s_fft,m,n) !calculating 2d DFT of source, tested
    use, intrinsic :: iso_c_binding
    use precision,only:p_
    implicit none
    include 'fftw3.f03'  
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    type(C_PTR) :: plan
    integer:: i,j

    !Fourier transformation along the first dimension
    in=s
    plan = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
    do j=0,n-1
       call fftw_execute_dft(plan, in(:,j), out(:,j))
    enddo
    call fftw_destroy_plan(plan)  

    !  in=out

    !  call toroidal_filter(in,out,m,n)

!!$  write(*,*) 'given by fftw'
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j)),i=0,m-1)
!!$  enddo
!!$do j=0,n-1
!!$call my_fft(in(:,j),out(:,j),m)
!!$enddo
!!$ write(*,*) 'given by my_fft'
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j)),i=0,m-1)
!!$  enddo

    !Fourier transformation along the second dimension
    in=out
    plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan)  
    s_fft=out

  end subroutine twod_fourier_transform

  subroutine twod_backward_fourier_transform(field_fft,field,m,n) !tested
    use, intrinsic :: iso_c_binding
    use precision,only:p_
    implicit none
    include 'fftw3.f03'  
    integer,intent(in):: m,n
    !  complex(C_DOUBLE_COMPLEX),intent(in):: field_fft(0:m-1,0:n-1)
    complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
    real(p_),intent(out):: field(0:m-1,0:n-1)
    type(C_PTR) :: plan
    !  complex(C_DOUBLE_COMPLEX) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j

    in=field_fft
    plan = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_backward,FFTW_ESTIMATE)
    do j=0,n-1 !fourier transform along the first direction
       call fftw_execute_dft(plan, in(:,j), out(:,j))
    enddo
    call fftw_destroy_plan(plan)  

    in=out
    plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_backward,FFTW_ESTIMATE)
    do i=0,m-1 !fourier transform along the second direction
       call fftw_execute_dft(plan, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan)  

    field=real(out)/(m*n)
!!$  !check whether the results are the same
!!$  write(*,*) 'original data='
!!$  do j=0,n-1
!!$     write(*,*) (s(i,j),i=0,m-1)
!!$  enddo
!!$  write(*,*) 'DFT and reverse-DFT data='
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j))/(m*n),i=0,m-1)
!!$  enddo
  end subroutine twod_backward_fourier_transform


subroutine my_fft(in,out,m) !for test
  use precision,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in):: m
  complex(p_),intent(in)::in(0:m-1)
  complex(p_),intent(out)::out(0:m-1)
  complex(p_),parameter::ii=(0.0_p_,1._p_)
  real(p_):: sum
  integer:: i,ip

  do i=0,m-1
     sum=0._p_
     do ip=0,m-1
        sum=sum+in(ip)*exp(-twopi*ii/m*ip*i)
     enddo
     out(i)=sum
  enddo

end subroutine my_fft


subroutine oned_backward_fourier_transform1(field_fft,field,m,n) 
  use, intrinsic :: iso_c_binding
  use precision,only:p_
  implicit none
  include 'fftw3.f03'  
  integer,intent(in):: m,n
  complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
  real(p_),intent(out):: field(0:m-1,0:n-1)
  type(C_PTR) :: plan
  complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
  integer:: j

  in=field_fft
  plan = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_backward,FFTW_ESTIMATE)
  do j=0,n-1 !fourier transform along the first direction
     call fftw_execute_dft(plan, in(:,j), out(:,j))
  enddo
  call fftw_destroy_plan(plan)  

  field=real(out)/m
end subroutine oned_backward_fourier_transform1


subroutine oned_backward_fourier_transform2(field_fft,field,m,n) 
  use, intrinsic :: iso_c_binding
  use precision,only:p_
  implicit none
  include 'fftw3.f03'  
  integer,intent(in):: m,n
  complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
  real(p_),intent(out):: field(0:m-1,0:n-1)
  type(C_PTR) :: plan
  complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
  integer:: i

  in=field_fft
  plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_backward,FFTW_ESTIMATE)
  do i=0,m-1 !fourier transform along the second direction
     call fftw_execute_dft(plan, in(i,:), out(i,:))
  enddo
  call fftw_destroy_plan(plan)  
  field=real(out)/n
end subroutine oned_backward_fourier_transform2

subroutine oned_sine_transform2(s,s_fft,m,n) !calculating 1d DST of s(:,:) along the second dimension
  use, intrinsic :: iso_c_binding
  use precision,only:p_
  implicit none
  include 'fftw3.f03'  
  integer,intent(in):: m,n
  real(p_),intent(in):: s(0:m-1,0:n-1)
  real(p_),intent(out):: s_fft(0:m-1,0:n-1)
  real(p_) :: in(0:m-1,0:n-1)
  type(C_PTR) :: plan !a pointer, which needs to be distoried, otherwise may cause memory leak
  integer:: i

  in=s
  plan =fftw_plan_r2r_1d(n, in(1,:), s_fft(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     call dfftw_execute_r2r(plan,s(i,:),s_fft(i,:)) !DST along the second dimension
  enddo
  call fftw_destroy_plan(plan)  

end subroutine oned_sine_transform2

subroutine oned_inverse_sine_transform2(s_dst,s,m,n) !computing 1d inverse DST of s(:,:) along the second dimension
  use, intrinsic :: iso_c_binding
  use precision,only:p_
  implicit none
  include 'fftw3.f03'  
  integer,intent(in):: m,n
  real(p_),intent(in):: s_dst(0:m-1,0:n-1)
  real(p_),intent(out):: s(0:m-1,0:n-1)
  real(p_) :: in(0:m-1,0:n-1)
  type(C_PTR) :: plan
  integer:: i

  in=s_dst
  plan =fftw_plan_r2r_1d(n, in(1,:), s(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     call dfftw_execute_r2r(plan,s_dst(i,:),s(i,:))   !DST along the second dimension
  enddo
  call fftw_destroy_plan(plan)  
  s=s/(2*(n+1))
end subroutine oned_inverse_sine_transform2


subroutine dst_dft(s,s_spectrum,m,n)
  use, intrinsic :: iso_c_binding
  use precision,only:p_
  implicit none
  include 'fftw3.f03'  
  integer,intent(in):: m,n
  real(p_),intent(in):: s(0:m-1,0:n-1)
  complex(p_),intent(out):: s_spectrum(0:m-1,0:n-1)
  real(p_):: in1(0:m-1,0:n-1),s_dst(0:m-1,0:n-1)
  complex(p_):: in2(0:m-1,0:n-1)
  type(C_PTR) :: plan
  integer:: i,j

  in1=s
  plan =fftw_plan_r2r_1d(n,   in1(1,:), s_dst(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     call dfftw_execute_r2r(plan,in1(i,:), s_dst(i,:)) !DST along the second dimension
  enddo
  call fftw_destroy_plan(plan)  

  in2=s_dst
  plan = fftw_plan_dft_1d(m,  in2(:,1), s_spectrum(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
  do j=0,n-1
     call fftw_execute_dft(plan, in2(:,j), s_spectrum(:,j))   !DFT along the first dimension
  enddo
  call fftw_destroy_plan(plan)

end subroutine dst_dft

subroutine mode_evolution_analysis3(t,a,m,n,file_unit) !using DFT for toroidal direction and DST for radial direction
  use constants,only: one
  use precision,only:p_
  use perturbation_field_matrix,only: toroidal_mode_number_included
  use magnetic_coordinates,only:nsegment
  implicit none
  real(p_),intent(in):: t
  integer,intent(in):: m,n,file_unit
  real(p_),intent(in):: a(0:m-1,0:n-1)
  complex(p_):: a_spectrum(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

  call dst_dft(a,a_spectrum,m,n)

  ipositive=toroidal_mode_number_included/nsegment
  !inegative=m-ipositive

  write(file_unit,'(20(1pe14.4))') t,(real(a_spectrum(ipositive,j)),imag(a_spectrum(ipositive,j)),j=0,4)

end subroutine mode_evolution_analysis3
end module transform_module
