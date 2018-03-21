subroutine report(t)
  use magnetic_coordinates,only: mtoroidal,nflux2
  use perturbation_field_matrix, only:ex=>ex_left,ey=>ey_left,epar=>epar_left
  use perturbation_field_matrix, only:mf_par=>mf_par_left,mf_x=>mf_x_left,mf_y=>mf_y_left
  use perturbation_field_matrix,only:source_e1,source_e2 !,jper_x_i_left,jper_y_i_left
  use perturbation_field_matrix,only:source1,source2,source3,source_faraday1, source_faraday2,source_faraday3
  use precision,only:p_
  implicit none
  real(p_),intent(in):: t
  integer:: i,j

i=mtoroidal/2
j=nflux2/2

!write(*,'(7(1pe14.4))') t,ex(i,j),ey(i,j),epar(i,j),mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) jper_x_i_left(mtoroidal/2,nflux2/2),jper_y_i_left(mtoroidal/2,nflux2/2)!,source_e1(mtoroidal/2,nflux2/2),source_e2(mtoroidal/2,nflux2/2)
!write(*,*) source1(mtoroidal/3,nflux2/3),source2(mtoroidal/3,nflux2/3),source3(mtoroidal/3,nflux2/3)
!write(*,*) source_faraday1(mtoroidal/3,nflux2/3),source_faraday2(mtoroidal/3,nflux2/3),source_faraday3(mtoroidal/3,nflux2/3)
end subroutine 


subroutine mode_evolution_analysis(t)
  use constants,only: one
  use precision,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtoroidal,n=>nflux2
  use perturbation_field_matrix,only: ex_left
  use perturbation_field_matrix,only: toroidal_mode_number_included
  use domain_decomposition,only:myid,numprocs
  use magnetic_coordinates,only:nsegment
  use transform_module,only: twod_fourier_transform
  implicit none
  real(p_),intent(in):: t
  real(p_):: a(0:m-1,0:n-1)
  complex(p_):: a_fft(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative,jp,jn
  do i=0,m-1
     do j=0,n-1
        a(i,j)=ex_left(i+1,j+1)
     enddo
  enddo
  call twod_fourier_transform(a,a_fft,m,n)


ipositive=toroidal_mode_number_included/nsegment
inegative=m-ipositive
!!$jp=9
!!$jn=n-jp
!!$if(myid.eq.1) write(*,*) t,real(a_fft(ipositive,jp)),imag(a_fft(ipositive,jp)),real(a_fft(inegative,jn)),imag(a_fft(inegative,jn))


 write(*,'(13(1pe14.4))') t,(real(a_fft(ipositive,j)),imag(a_fft(ipositive,j)),j=0,5)


end subroutine mode_evolution_analysis


subroutine mode_evolution_analysis2(t) !for adiabatic electron model, using DFT for both toroidal and raial direction
  use constants,only: one
  use precision,only:p_
  use magnetic_coordinates,only: m=>mtoroidal,n=>nflux2
  use perturbation_field_matrix,only: ef_cyl_phi_left,ef_cyl_r_left,ef_cyl_z_left,den_left
  use perturbation_field_matrix,only: ef_cyl_phi_right,ef_cyl_r_right,ef_cyl_z_right
  use perturbation_field_matrix,only: toroidal_mode_number_included
  use magnetic_coordinates,only:nsegment
  use transform_module,only: twod_fourier_transform
  implicit none
  real(p_),intent(in):: t
  real(p_):: a(0:m-1,0:n-1)
  complex(p_):: a_fft(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

  do i=0,m-1
     do j=0,n-1
        a(i,j)=ef_cyl_phi_left(i+1,j+1)
     enddo
  enddo
  call twod_fourier_transform(a,a_fft,m,n)

  ipositive=toroidal_mode_number_included/nsegment
  inegative=m-ipositive

  write(*,'(20(1pe14.4))') t,(real(a_fft(ipositive,j)),imag(a_fft(ipositive,j)),j=0,4),&
       & den_left(m/2,n/2),ef_cyl_phi_left(m/2,n/2)

end subroutine mode_evolution_analysis2





