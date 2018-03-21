module filter_module
implicit none
private
public toroidal_filter,radial_fourier_filter,radial_sine_filter
contains
subroutine toroidal_filter(in,out,m,n)
  use precision,only:p_
  use perturbation_field_matrix,only: toroidal_mode_number_included !is postive and less than (m/2)*nsegment
  use magnetic_coordinates,only:nsegment !computational toroidal region is 1/nsegment of the full torus
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in)::   in(0:m-1,0:n-1)
  complex(p_),intent(out):: out(0:m-1,0:n-1)
  integer:: nharmonic
  integer:: i,j,i_negative_wn

  out=(0._p_,0._p_) !initialized to zero
  nharmonic=toroidal_mode_number_included/nsegment !harmonic that needs to be kept in the computational toroidal segment
  do j=0,n-1
     do i=0,m/2 !scanning over the positive wavenumber
        if(i.eq.nharmonic) then 
           out(i,j)=in(i,j) !positive wavenumber
           i_negative_wn=m-i
           if(i.eq.0) i_negative_wn=0
           out(i_negative_wn,j)=in(i_negative_wn,j) !negative or zero wavenumber
        endif
     enddo
  enddo
end subroutine toroidal_filter


subroutine radial_fourier_filter(in,out,m,n)
  use precision,only:p_
  use perturbation_field_matrix,only:  radial_harmonics_included
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in)::   in(0:m-1,0:n-1)
  complex(p_),intent(out):: out(0:m-1,0:n-1)
  integer:: i,j,j_negative_wn

  out=(0._p_,0._p_) !initialized to zero

  do i=0,m-1
     do j=0,n/2  !scanning over the positive wavenumber
        if(j.le.radial_harmonics_included) then 
           out(i,j)=in(i,j) !positive wavenumber
           j_negative_wn=n-j
           if(j.eq.0) j_negative_wn=0
           out(i,j_negative_wn)=in(i,j_negative_wn) !negative or zero wavenumber
        endif
     enddo
  enddo
end subroutine radial_fourier_filter

subroutine radial_sine_filter(s,m,n)
  use precision,only:p_
  use perturbation_field_matrix,only:  radial_harmonics_included
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(inout):: s(0:m-1,0:n-1)
  integer:: i,j

  do j=0,n-1
     do i=0,m-1
        if(j.gt.radial_harmonics_included) s(i,j)=0._p_
     enddo
  enddo

end subroutine radial_sine_filter
end module filter_module
