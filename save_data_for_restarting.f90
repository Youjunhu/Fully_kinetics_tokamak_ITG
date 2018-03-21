subroutine save_data_for_restarting()
  use ions_module,only:nmarker_i,active_i,radcor_i,theta_i,alpha_i,w_i,&
       & r_i,phi_i,z_i,r_i_mid,phi_i_mid,z_i_mid,&
       & vr_i_integer,vphi_i_integer,vz_i_integer,vr_i,vphi_i,vz_i
  use domain_decomposition,only:myid
  implicit none
  character(len=64)::cfile
  integer:: ufile !file nunit

  cfile = 'myidxxxxx.pd'
  write(cfile(5:9),'(i5.5)') myid
  open(newunit=ufile,file=cfile,form='unformatted') !newunit is a keyword argument (as output) of open()
  write(ufile) nmarker_i
  write(ufile) 0

end subroutine save_data_for_restarting
