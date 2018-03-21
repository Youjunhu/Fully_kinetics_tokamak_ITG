subroutine field_solver(myid,numprocs,solve_id)
  use constants,only: one,one_half
  use precision,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtoroidal,n=>nflux2,dtheta
  use perturbation_field_matrix,only:coeff1_ex, coeff1_ey, coeff1_ex_yy, coeff1_ex_xy, coeff1_ey_xx, coeff1_ey_xy !coefficients of the 1st perpendicular equation
  use perturbation_field_matrix,only:coeff2_ex, coeff2_ey, coeff2_ex_yy, coeff2_ex_xy, coeff2_ey_xx, coeff2_ey_xy !coefficients of the 2nd perpendicular equation
  use perturbation_field_matrix,only:coeff3_epar_xx,coeff3_epar_yy,coeff3_epar_xy,coeff3_epar,coeff3_epar_x,coeff3_epar_y
  use perturbation_field_matrix,only:coeff1_ey_implicit, coeff2_ex_implicit
  use perturbation_field_matrix,only: source1,source2,source3 
  use perturbation_field_matrix,only: ex_left,ey_left,epar_left !as output
  use perturbation_field_matrix,only: filter_toroidal
  use transform_module,only: twod_fourier_transform,twod_backward_fourier_transform
  use filter_module,only: toroidal_filter
  implicit none
  integer,intent(in):: myid,numprocs,solve_id
  complex(p_):: source1_fft(m,n),source2_fft(m,n),source3_fft(m,n)
  complex(p_):: ex_fft(m,n),ey_fft(m,n),epar_fft(m,n), in(m,n),out(m,n)
  complex(p_),parameter:: ii=(0._p_,1._p_)
  real(p_)::    ex(m,n),    ey(m,n),    epar(m,n)
  integer,parameter:: nvar=3 !E_para,E_x,E_y
  complex(p_)::  a(nvar,nvar),b(nvar)
  INTEGER:: ipiv(nvar),info
  real(p_):: ktor,krad  !toroidal wavenumber and radial wavenumber
  real(p_):: cdt !,cimp
  integer:: i,j
  !  real(p_),parameter:: eps=1.0d-8

  call prepare_total_source_terms()

  if(solve_id.eq. 1) cdt=one_half
  if(solve_id.eq. 2) cdt=one
  !  if(solve_id.eq. 1) cimp=one_half
  ! if(solve_id.eq. 2) cimp=one
!!$  if(myid.eq.0) write(*,*) 'field=, before solving',ex_left(m/3,n/3), ey_left(m/3,n/3), 'source=',source1(m/3,n/3),&
!!$       & source2(m/3,n/3), source3(m/3,n/3)

!!$do i=1,m
!!$do j=1,n
!!$if(abs(source1(i,j)).ge.eps) write(*,*) 'source term1 is nonzero', source1(i,j)
!!$if(abs(source2(i,j)).ge.eps) write(*,*) 'source term2 is nonzero',source2(i,j)
!!$if(abs(source3(i,j)).ge.eps) write(*,*) 'source term3 is nonzero',source3(i,j)
!!$enddo
!!$enddo 

  call twod_fourier_transform(source1,source1_fft,m,n) !Fourier transform the source term
  call twod_fourier_transform(source2,source2_fft,m,n) !Fourier transform the source term
  call twod_fourier_transform(source3,source3_fft,m,n) !Fourier transform the source term

  ! if(myid.eq.0) write(*,*) coeff1_ex, coeff1_ey, coeff1_ex_yy, coeff1_ex_xy, coeff1_ey_xx, coeff1_ey_xy, 'first eq'
  !if(myid.eq.0) write(*,*) coeff2_ex, coeff2_ey, coeff2_ex_yy, coeff2_ex_xy, coeff2_ey_xx, coeff2_ey_xy, 'second eq'
  !solve field equations in Fourier-space

  do j=1,n !in Fourier space
     do i=1,m !in Fourier space
        krad=(j-1)/(n*dradcor)    !the corresponding wavenumber along radial direction (1/lambda, not twopi/lambda)
        if((j-1).gt.n/2) krad=(j-1-n)/(n*dradcor)      !corresponding to the negative wavenumber
        ktor=(i-1)/(m*dtor) !the corresponding wavenumber along toroidal direction (1/lambda, not twopi/lambda)
        if((i-1).gt.m/2) ktor=(i-1-m)/(m*dtor)   !corresponding to the negative wavenumber

        b(1)=source1_fft(i,j) !right-hand-side
        b(2)=source2_fft(i,j)
        b(3)=source3_fft(i,j)

        a(1,1)=coeff1_ex+(coeff1_ex_xy*(-krad*ktor)+coeff1_ex_yy*(-ktor**2))*cdt
        a(1,2)=coeff1_ey+coeff1_ey_implicit*cdt+(coeff1_ey_xy*(-krad*ktor)+coeff1_ey_xx*(-krad**2))*cdt
        a(1,3)=0._p_

        a(2,1)=coeff2_ex+coeff2_ex_implicit*cdt+(coeff2_ex_xy*(-krad*ktor)+coeff2_ex_yy*(-ktor**2))*cdt
        a(2,2)=coeff2_ey+(coeff2_ey_xy*(-krad*ktor)+coeff2_ey_xx*(-krad**2))*cdt
        a(2,3)=0._p_

        a(3,1)=0.0_p_
        a(3,2)=0.0_p_
        a(3,3)=coeff3_epar_xx*(-krad**2)+coeff3_epar_xy*(-krad*ktor)+coeff3_epar_yy*(-ktor*ktor)&
             & +coeff3_epar+coeff3_epar_x*(ii*krad)+coeff3_epar_y*(ii*ktor)
        call zgesv(nvar,1,a,nvar,ipiv,b,nvar,info) !lapack routine to slove linear system of equations, ref: http://www.netlib.org/lapack/explore-3.1.1-html/zgesv.f.html

        ex_fft(i,j)=b(1)
        ey_fft(i,j)=b(2)
        epar_fft(i,j)=b(3)
     enddo
  enddo

  if(filter_toroidal.eqv..true.) then
     in=ex_fft
     call toroidal_filter(in,out,m,n)
     ex_fft=out

     in=ey_fft
     call toroidal_filter(in,out,m,n)
     ey_fft=out

     in=epar_fft
     call toroidal_filter(in,out,m,n)
     epar_fft=out
  endif

  call twod_backward_fourier_transform(ex_fft,  ex,m,n)
  call twod_backward_fourier_transform(ey_fft,  ey,m,n)
  call twod_backward_fourier_transform(epar_fft,epar,m,n)

  do i=1,m
     do j=1,n
        ex_left(i,j)=ex(i,j)
        ey_left(i,j)=ey(i,j)
        epar_left(i,j)=epar(i,j)
     enddo
  enddo

  do j=1,n
     ex_left(m+1,j)=ex(1,j) !peroidic boundary condition
     ey_left(m+1,j)=ey(1,j) !peroidic boundary condition
     epar_left(m+1,j)=epar(1,j) !peroidic boundary condition
  enddo


end subroutine field_solver


subroutine field_solver_electromagnetic_case(dtao) 
  use constants,only: one,one_half
  use precision,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtoroidal,n=>nflux2,dtheta

  use perturbation_field_matrix,only:coeff1_ex, coeff1_ey, coeff1_ex_yy, coeff1_ex_xy, coeff1_ey_xx, coeff1_ey_xy !coefficients of the 1st perpendicular equation
  use perturbation_field_matrix,only:coeff2_ex, coeff2_ey, coeff2_ex_yy, coeff2_ex_xy, coeff2_ey_xx, coeff2_ey_xy !coefficients of the 2nd perpendicular equation
  use perturbation_field_matrix,only:coeff3_epar_xx,coeff3_epar_yy,coeff3_epar_xy,coeff3_epar,coeff3_epar_x,coeff3_epar_y
  use perturbation_field_matrix,only:coeff1_ey_implicit, coeff2_ex_implicit
  use perturbation_field_matrix,only: source1,source2,source3 
  use perturbation_field_matrix,only: ex_left,ey_left,epar_left !as output
  use perturbation_field_matrix,only: filter_toroidal,filter_radial
  use domain_decomposition,only:myid,numprocs
  use ions_module,only: nmarker_i,w_i
  use transform_module,only: twod_fourier_transform,twod_backward_fourier_transform
  use filter_module,only: toroidal_filter,radial_fourier_filter
  implicit none
  real(p_),intent(in):: dtao
  !  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  ! real(p_),intent(in):: vr_i(nmarker_i),vphi_i(nmarker_i),vz_i(nmarker_i)
  real(p_)::   epar(m,n), ex_new(m,n),  ey_new(m,n)
  real(p_)::   ex_old(m,n),  ey_old(m,n)
  complex(p_):: ex_fft(m,n),ey_fft(m,n),epar_fft(m,n), in(m,n),out(m,n)
  complex(p_):: source1_fft(m,n),source2_fft(m,n),source3_fft(m,n)
  complex(p_),parameter:: ii=(0._p_,1._p_)
  real(p_):: ktor,krad  !toroidal wavenumber and radial wavenumber
  integer,parameter:: nvar=3 !E_para,E_x,E_y
  complex(p_)::  a(nvar,nvar),b(nvar)
  INTEGER:: ipiv(nvar),info
  integer:: i,j,it
  !integer,parameter::num_iteration=6
  integer,parameter::num_iteration=1
  real(p_),parameter:: eps=1.0d-3
  real(p_):: w_i_tmp(nmarker_i)

  !  do it=1,num_iteration !iteration over the perpendicular electric field
  !write(*,*) 'iteration=',it, 'myid=',myid
!!$     do i=1,m !record the old value of the perpendicular electric field
!!$        do j=1,n
!!$           ex_old(i,j)=ex_left(i,j)
!!$           ey_old(i,j)=ey_left(i,j)
!!$        enddo
!!$     enddo

  !call push_ion_weight0(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_tmp,nmarker_i) !implicit scheme for ion weight evolution, requiring iteration

  !     call push_ion_weight_without_eper(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_tmp,nmarker_i) !implicit scheme for ion weight evolution, not requiring iteration

  !     call deposit_ions(radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i_tmp)

  call prepare_ion_source_terms()
  call prepare_total_source_terms()

  call twod_fourier_transform(source1,source1_fft,m,n) !Fourier transform the source term
  call twod_fourier_transform(source2,source2_fft,m,n) !Fourier transform the source term
  call twod_fourier_transform(source3,source3_fft,m,n) !Fourier transform the source term

  do j=1,n 
     do i=1,m
        krad=(j-1)/(n*dradcor)    !the corresponding wavenumber along radial direction (1/lambda, not twopi/lambda)
        if((j-1).gt.n/2) krad=(j-1-n)/(n*dradcor)      !corresponding to the negative wavenumber
        ktor=(i-1)/(m*dtor) !the corresponding wavenumber along toroidal direction (1/lambda, not twopi/lambda)
        if((i-1).gt.m/2) ktor=(i-1-m)/(m*dtor)   !corresponding to the negative wavenumber

        b(1)=source1_fft(i,j) !right-hand-side
        b(2)=source2_fft(i,j)
        b(3)=source3_fft(i,j)

        a(1,1)=coeff1_ex +(coeff1_ex_xy*(-krad*ktor)+coeff1_ex_yy*(-ktor**2))*dtao
        a(1,2)=coeff1_ey +coeff1_ey_implicit*dtao  +(coeff1_ey_xy*(-krad*ktor)+coeff1_ey_xx*(-krad**2))*dtao
        a(1,3)=0._p_

        a(2,1)=coeff2_ex +coeff2_ex_implicit*dtao +(coeff2_ex_xy*(-krad*ktor)+coeff2_ex_yy*(-ktor**2))*dtao
        a(2,2)=coeff2_ey +(coeff2_ey_xy*(-krad*ktor)+coeff2_ey_xx*(-krad**2))*dtao
        a(2,3)=0._p_

        a(3,1)=0.0_p_
        a(3,2)=0.0_p_
        a(3,3)=coeff3_epar +coeff3_epar_xx*(-krad**2)+coeff3_epar_xy*(-krad*ktor)+coeff3_epar_yy*(-ktor*ktor)&
             & +coeff3_epar_x*dtao*(ii*krad)+coeff3_epar_y*dtao*(ii*ktor)
        call zgesv(nvar,1,a,nvar,ipiv,b,nvar,info) !lapack routine to slove linear system of equations, ref: http://www.netlib.org/lapack/explore-3.1.1-html/zgesv.f.html

        ex_fft(i,j)=b(1)
        ey_fft(i,j)=b(2)
        epar_fft(i,j)=b(3)

     enddo
  enddo

  if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
     in=ex_fft
     call toroidal_filter(in,out,m,n)
     ex_fft=out

     in=ey_fft
     call toroidal_filter(in,out,m,n)
     ey_fft=out

     in=epar_fft
     call toroidal_filter(in,out,m,n)
     epar_fft=out
  endif

  if(filter_radial.eqv..true.) then !filter over the radial mode number, keeping only low-radial-harmonics of the perturbation
     in=ex_fft
     call radial_fourier_filter(in,out,m,n)
     ex_fft=out

     in=ey_fft
     call radial_fourier_filter(in,out,m,n)
     ey_fft=out

     in=epar_fft
     call radial_fourier_filter(in,out,m,n)
     epar_fft=out
  endif

  call twod_backward_fourier_transform(ex_fft,  ex_new,m,n)
  call twod_backward_fourier_transform(ey_fft,  ey_new,m,n)
  call twod_backward_fourier_transform(epar_fft,epar,m,n)

  do i=1,m
     do j=1,n
        ex_left(i,j)=ex_new(i,j)
        ey_left(i,j)=ey_new(i,j)
        epar_left(i,j)=epar(i,j)
     enddo
  enddo

  do j=1,n
     ex_left(m+1,j)=ex_new(1,j) !peroidic boundary condition
     ey_left(m+1,j)=ey_new(1,j) !peroidic boundary condition
     epar_left(m+1,j)=epar(1,j) !peroidic boundary condition
  enddo


  call communicate_field_value_between_neighbour_cells()

!!$     i=m/2 !choose an arbitrary location
!!$     j=n/2
!!$     if(abs(ex_old(i,j)-ex_new(i,j)).le.eps .and. abs(ey_old(i,j)-ey_new(i,j)).le.eps) exit !do not use this because it cause dead lock between mpi procs
  ! enddo
  !  if(it.eq.num_iteration+1) write(*,*) 'warning****, exceed maximal number of iteration.'

end subroutine field_solver_electromagnetic_case




subroutine field_solver_electromagnetic_case_iteration(dtao,radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i) !electrostatic case, field equations become algebra equations, iteration for implicit ion weight pusher
  use constants,only: one,one_half
  use precision,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtoroidal,n=>nflux2,dtheta

  use perturbation_field_matrix,only:coeff1_ex, coeff1_ey, coeff1_ex_yy, coeff1_ex_xy, coeff1_ey_xx, coeff1_ey_xy !coefficients of the 1st perpendicular equation
  use perturbation_field_matrix,only:coeff2_ex, coeff2_ey, coeff2_ex_yy, coeff2_ex_xy, coeff2_ey_xx, coeff2_ey_xy !coefficients of the 2nd perpendicular equation
  use perturbation_field_matrix,only:coeff3_epar_xx,coeff3_epar_yy,coeff3_epar_xy,coeff3_epar,coeff3_epar_x,coeff3_epar_y
  use perturbation_field_matrix,only:coeff1_ey_implicit, coeff2_ex_implicit
  use perturbation_field_matrix,only: source1,source2,source3 
  use perturbation_field_matrix,only: ex_left,ey_left,epar_left !as output
  use perturbation_field_matrix,only: filter_toroidal
  use domain_decomposition,only:myid,numprocs
  use ions_module,only: nmarker_i,w_i,active_i
  use transform_module,only: twod_fourier_transform,twod_backward_fourier_transform
  use filter_module,only:toroidal_filter
  use deposit_ions_module,only: deposit_ions
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  real(p_),intent(in):: vr_i(nmarker_i),vphi_i(nmarker_i),vz_i(nmarker_i)
  real(p_)::   epar(m,n), ex_new(m,n),  ey_new(m,n)
  real(p_)::   ex_old(m,n),  ey_old(m,n)
  complex(p_):: ex_fft(m,n),ey_fft(m,n),epar_fft(m,n), in(m,n),out(m,n)
  complex(p_):: source1_fft(m,n),source2_fft(m,n),source3_fft(m,n)
  complex(p_),parameter:: ii=(0._p_,1._p_)
  real(p_):: ktor,krad  !toroidal wavenumber and radial wavenumber
  integer,parameter:: nvar=3 !E_para,E_x,E_y
  complex(p_)::  a(nvar,nvar),b(nvar)
  INTEGER:: ipiv(nvar),info
  integer:: i,j,it
  integer,parameter::num_iteration=6
  !integer,parameter::num_iteration=1
  real(p_),parameter:: eps=1.0d-3
  real(p_):: w_i_tmp(nmarker_i)

  do it=1,num_iteration !iteration over the perpendicular electric field
     !write(*,*) 'iteration=',it, 'myid=',myid
!!$     do i=1,m !record the old value of the perpendicular electric field
!!$        do j=1,n
!!$           ex_old(i,j)=ex_left(i,j)
!!$           ey_old(i,j)=ey_left(i,j)
!!$        enddo
!!$     enddo

     call push_ion_weight0(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_tmp,nmarker_i) !implicit scheme for ion weight evolution, requiring iteration
!     call push_ion_weight_without_eper(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_tmp,nmarker_i) !implicit scheme for ion weight evolution, not requiring iteration

!     call deposit_ions(                   radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i_tmp)
     call deposit_ions (nmarker_i,active_i,radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i_tmp)
     call prepare_ion_source_terms()
     call prepare_total_source_terms()

     call twod_fourier_transform(source1,source1_fft,m,n) !Fourier transform the source term
     call twod_fourier_transform(source2,source2_fft,m,n) !Fourier transform the source term
     call twod_fourier_transform(source3,source3_fft,m,n) !Fourier transform the source term

     do j=1,n 
        do i=1,m
           krad=(j-1)/(n*dradcor)    !the corresponding wavenumber along radial direction (1/lambda, not twopi/lambda)
           if((j-1).gt.n/2) krad=(j-1-n)/(n*dradcor)      !corresponding to the negative wavenumber
           ktor=(i-1)/(m*dtor) !the corresponding wavenumber along toroidal direction (1/lambda, not twopi/lambda)
           if((i-1).gt.m/2) ktor=(i-1-m)/(m*dtor)   !corresponding to the negative wavenumber

           b(1)=source1_fft(i,j) !right-hand-side
           b(2)=source2_fft(i,j)
           b(3)=source3_fft(i,j)

           a(1,1)=coeff1_ex +(coeff1_ex_xy*(-krad*ktor)+coeff1_ex_yy*(-ktor**2))*dtao
           a(1,2)=coeff1_ey +coeff1_ey_implicit*dtao*0._p_  +(coeff1_ey_xy*(-krad*ktor)+coeff1_ey_xx*(-krad**2))*dtao
           a(1,3)=0._p_

           a(2,1)=coeff2_ex +coeff2_ex_implicit*dtao*0._p_ +(coeff2_ex_xy*(-krad*ktor)+coeff2_ex_yy*(-ktor**2))*dtao
           a(2,2)=coeff2_ey +(coeff2_ey_xy*(-krad*ktor)+coeff2_ey_xx*(-krad**2))*dtao
           a(2,3)=0._p_

           a(3,1)=0.0_p_
           a(3,2)=0.0_p_
           a(3,3)=coeff3_epar +coeff3_epar_xx*(-krad**2)+coeff3_epar_xy*(-krad*ktor)+coeff3_epar_yy*(-ktor*ktor)&
                & +coeff3_epar_x*dtao*(ii*krad)+coeff3_epar_y*dtao*(ii*ktor)
           call zgesv(nvar,1,a,nvar,ipiv,b,nvar,info) !lapack routine to slove linear system of equations, ref: http://www.netlib.org/lapack/explore-3.1.1-html/zgesv.f.html

           ex_fft(i,j)=b(1)
           ey_fft(i,j)=b(2)
           epar_fft(i,j)=b(3)

        enddo
     enddo

     if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
        in=ex_fft
        call toroidal_filter(in,out,m,n)
        ex_fft=out

        in=ey_fft
        call toroidal_filter(in,out,m,n)
        ey_fft=out

        in=epar_fft
        call toroidal_filter(in,out,m,n)
        epar_fft=out
     endif

     call twod_backward_fourier_transform(ex_fft,  ex_new,m,n)
     call twod_backward_fourier_transform(ey_fft,  ey_new,m,n)
     call twod_backward_fourier_transform(epar_fft,epar,m,n)

     do i=1,m
        do j=1,n
           ex_left(i,j)=ex_new(i,j)
           ey_left(i,j)=ey_new(i,j)
           epar_left(i,j)=epar(i,j)
        enddo
     enddo

     do j=1,n
        ex_left(m+1,j)=ex_new(1,j) !peroidic boundary condition
        ey_left(m+1,j)=ey_new(1,j) !peroidic boundary condition
        epar_left(m+1,j)=epar(1,j) !peroidic boundary condition
     enddo



     call communicate_field_value_between_neighbour_cells()

!!$     i=m/2 !choose an arbitrary location
!!$     j=n/2
!!$     if(abs(ex_old(i,j)-ex_new(i,j)).le.eps .and. abs(ey_old(i,j)-ey_new(i,j)).le.eps) exit !do not use this because it cause dead lock between mpi procs
  enddo
  !  if(it.eq.num_iteration+1) write(*,*) 'warning****, exceed maximal number of iteration.'

end subroutine field_solver_electromagnetic_case_iteration


subroutine field_solver_electrostatic_case(dtao,radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i) !electrostatic case, field equations become algebra equations, iteration for implicit ion weight pusher
  use constants,only: one,one_half
  use precision,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtoroidal,n=>nflux2,dtheta
  use perturbation_field_matrix,only:coeff1_ex,coeff1_ey,coeff1_ey_implicit !coefficients of the 1st perpendicular equation
  use perturbation_field_matrix,only:coeff2_ex, coeff2_ey,coeff2_ex_implicit !coefficients of the 2nd perpendicular equation
  use perturbation_field_matrix,only:coeff3_epar
  use perturbation_field_matrix,only: source1,source2,source3 
  use perturbation_field_matrix,only: ex_left,ey_left,epar_left !as output
  use perturbation_field_matrix,only: filter_toroidal
  use domain_decomposition,only:myid,numprocs
  use ions_module,only: nmarker_i,w_i,active_i
  use transform_module,only: twod_fourier_transform,twod_backward_fourier_transform
  use filter_module,only:toroidal_filter
  use deposit_ions_module,only: deposit_ions
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  real(p_),intent(in):: vr_i(nmarker_i),vphi_i(nmarker_i),vz_i(nmarker_i)
  real(p_)::   epar(m,n), ex_new(m,n),  ey_new(m,n)
  real(p_)::   ex_old(m,n),  ey_old(m,n)
  complex(p_):: ex_fft(m,n),ey_fft(m,n),epar_fft(m,n), in(m,n),out(m,n)
  integer,parameter:: nvar=3 !E_para,E_x,E_y
  real(p_)::  a(nvar,nvar),b(nvar)
  INTEGER:: ipiv(nvar),info
  integer:: i,j,it
!  integer,parameter::num_iteration=6
  integer,parameter::num_iteration=1
  real(p_),parameter:: eps=1.0d-3
  real(p_):: w_i_tmp(nmarker_i)

  do it=1,num_iteration !iteration over the perpendicular electric field
     !write(*,*) 'iteration=',it, 'myid=',myid
!!$     do i=1,m !record the old value of the perpendicular electric field
!!$        do j=1,n
!!$           ex_old(i,j)=ex_left(i,j)
!!$           ey_old(i,j)=ey_left(i,j)
!!$        enddo
!!$     enddo

!     call push_ion_weight0(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_tmp,nmarker_i) !implicit scheme for ion weight evolution, requiring iteration
     call push_ion_weight_without_eper(dtao,radcor_i,theta_i,alpha_i,w_i,w_i_tmp,nmarker_i) !implicit scheme for ion weight evolution, not requiring iteration
!     call deposit_ions(                  radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i_tmp)
    call deposit_ions(nmarker_i,active_i,radcor_i,theta_i,alpha_i,vr_i,vphi_i,vz_i,w_i_tmp)
     call prepare_ion_source_terms()
     call prepare_total_source_terms()

     do j=1,n !solve the algebra equation systems at each spatial grid
        do i=1,m

           a(1,1)=coeff1_ex !coefficients of the linear system
           a(1,2)=coeff1_ey +coeff1_ey_implicit*dtao
           a(1,3)=0._p_

           a(2,1)=coeff2_ex +coeff2_ex_implicit*dtao
           a(2,2)=coeff2_ey
           a(2,3)=0._p_

           a(3,1)=0.0_p_
           a(3,2)=0.0_p_
           a(3,3)=coeff3_epar

           b(1)=source1(i,j) !right-hand-side of the linear system
           b(2)=source2(i,j)
           b(3)=source3(i,j)
           call dgesv(nvar,1,a,nvar,ipiv,b,nvar,info) !lapack routine to solve linear system of equations, ref: http://www.netlib.org/lapack/explore-3.1.1-html/zgesv.f.html
           ex_new(i,j)=b(1) !the solution of the linear system
           ey_new(i,j)=b(2)
           epar(i,j)=b(3)
        enddo
     enddo

     if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
        call twod_fourier_transform(ex_new,ex_fft,m,n) !Fourier transform of ex
        call twod_fourier_transform(ey_new,ey_fft,m,n) !Fourier transform of ey
        call twod_fourier_transform(epar,epar_fft,m,n) !Fourier transform of epar
        in=ex_fft
        call toroidal_filter(in,out,m,n)
        ex_fft=out

        in=ey_fft
        call toroidal_filter(in,out,m,n)
        ey_fft=out

        in=epar_fft
        call toroidal_filter(in,out,m,n)
        epar_fft=out
        call twod_backward_fourier_transform(ex_fft,  ex_new,m,n)
        call twod_backward_fourier_transform(ey_fft,  ey_new,m,n)
        call twod_backward_fourier_transform(epar_fft,epar,m,n)
     endif

     do i=1,m
        do j=1,n
           ex_left(i,j)=ex_new(i,j)
           ey_left(i,j)=ey_new(i,j)
           epar_left(i,j)=epar(i,j)
        enddo
     enddo

     do j=1,n
        ex_left(m+1,j)=ex_new(1,j) !peroidic boundary condition
        ey_left(m+1,j)=ey_new(1,j) !peroidic boundary condition
        epar_left(m+1,j)=epar(1,j) !peroidic boundary condition
     enddo

     call communicate_field_value_between_neighbour_cells()

!!$     i=m/2 !choose an arbitrary location
!!$     j=n/2
!!$     if(abs(ex_old(i,j)-ex_new(i,j)).le.eps .and. abs(ey_old(i,j)-ey_new(i,j)).le.eps) exit !do not use this because it cause dead lock between mpi procs
  enddo
  !  if(it.eq.num_iteration+1) write(*,*) 'warning****, exceed maximal number of iteration.'

end subroutine field_solver_electrostatic_case








