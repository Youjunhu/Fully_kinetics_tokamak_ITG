subroutine read_gfile(gfile_name,nwmax,psi,nx,nz,rleft,zmid,xdim,zdim,psi_axis,psi_lcfs,xmaxis,zmaxis, baxis,&
     & fpsi,qpsi,presspsi,pprime,ffprime,myid)
  use precision,only:p_
  use boundary,only: nlim,rlim,zlim,np_lcfs, x_lcfs,z_lcfs !as output
  use radial_module,only: sign_bphi

  implicit none
  character(*),intent(in):: gfile_name
  integer,intent(in)::nwmax,myid
  real(p_),intent(out):: psi(nwmax,nwmax)
  integer,intent(out):: nx,nz
  character(len=8):: ntitle(5),vid
  integer:: neq,ipestg
  real(p_)::xdim, zdim, rmajor_mk, rleft, zmid
  real(p_)::xmaxis, zmaxis, psi_axis, psi_lcfs, btorus,baxis
  real(p_):: dumaraya5(5),dumarayb5(5)
  real(p_),intent(out):: fpsi(nwmax),qpsi(nwmax),presspsi(nwmax)
  real(p_),intent(out):: ffprime(nwmax),pprime(nwmax)
  !  integer::np_lcfs !, nlim
  ! real(p_):: x_lcfs(nwmax),z_lcfs(nwmax) !,rlim(nwmax), zlim(nwmax)

  integer:: i,j,ierr
  logical:: reverse_tf,reverse_ip
  namelist /modify_magnetic_configuration/reverse_tf,reverse_ip

!!$  namelist/gfile_namelist/gfile_name
!!$  open(11,file='input.nmlt')
!!$  read(11,gfile_namelist)
!!$  close(11)
!!$  write(*,gfile_namelist)


  neq=123
  !open and read in eqdsk file
  open(neq,file=gfile_name,status='old')

  ipestg = 4
  read (neq, '(6a8, 3i4)') (ntitle(i), i=1,5), vid, ipestg, nx, nz
  read (neq,300) xdim, zdim, rmajor_mk, rleft, zmid
  read (neq,300) xmaxis, zmaxis, psi_axis, psi_lcfs, btorus

  read (neq,300) dumaraya5
  read (neq,300) dumarayb5

  read (neq ,300) (fpsi(i), i=1,nx)
  read (neq ,300) (presspsi(i), i=1,nx)
  read (neq ,300) (ffprime(i), i=1,nx)
  read (neq ,300) (pprime(i), i=1,nx)
  read (neq ,300) ((psi(i,j), i=1,nx), j=1,nz)
  read (neq ,300) (qpsi(i), i=1,nx)
  read (neq ,'(2i5)') np_lcfs, nlim
  allocate(x_lcfs(np_lcfs))
  allocate(z_lcfs(np_lcfs))
  allocate(rlim(nlim))
  allocate(zlim(nlim))
  read (neq ,300) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  read (neq ,300) (rlim(i), zlim(i), i=1,nlim)
  close(neq)



  !to verify that I have read the eqdsk file correctly, I write the data read to a new file called 'tmp.gfile'.
  !After the program finished, I compare the file 'tmp.gfile' with the original file using diff command
  !the output of diff command indicates that the two files are idential, which shows I have read the eqdsk file correctly
  neq=111
  open(neq,file='tmp.gfile')
  write (neq, '(6a8, 3i4)') (ntitle(i), i=1,5), vid, ipestg, nx, nz
  write (neq,300) xdim, zdim, rmajor_mk, rleft, zmid
  write (neq,300) xmaxis, zmaxis, psi_axis, psi_lcfs, btorus
  write (neq,300) dumaraya5
  write (neq,300) dumarayb5
  write (neq ,300) (fpsi(i), i=1,nx)
  write (neq ,300) (presspsi(i), i=1,nx)
  write (neq ,300) (ffprime(i), i=1,nx)
  write (neq ,300) (pprime(i), i=1,nx)
  write (neq ,300) ((psi(i,j), i=1,nx), j=1,nz)
  write (neq ,300) (qpsi(i), i=1,nx)
  write (neq ,'(2i5)') np_lcfs, nlim
  write (neq ,300) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  write (neq ,300) (rlim(i), zlim(i), i=1,nlim)
  close(neq)

300 format (5e16.9)


  open(31,file='input.nmlt')
  read(31,modify_magnetic_configuration)
  close(31)
  if(myid.eq.0) write(*,modify_magnetic_configuration)


  !Somtimes, I alter some quantities (e.g. reverse the toroidal magnetic field, or increase the pressure by a constant),in this case, the out gfile is different from the original one
!!$  do i=1,nx
!!$     presspsi(i)=presspsi(i)+0.5*presspsi(1) !increase the presure
!!$  enddo


  if(reverse_tf.eqv..true.)  fpsi=-fpsi !revert the toroidal magnetic field

  if(reverse_ip.eqv..true.) then
     psi=-psi !revert direction of the torodial current
     psi_axis=-psi_axis
     psi_lcfs=-psi_lcfs
     ffprime=-ffprime
     pprime=-pprime
  endif

  baxis=fpsi(1)/xmaxis


  !  write(*,*) 'xdim=',xdim, 'zdim=',zdim, 'rleft=',rleft, 'zmid=',zmid
  if(myid.eq.0) write(*,*) 'rleft=',rleft, 'rright=', rleft+xdim,'zlow=',zmid-zdim/2._p_,'zupp=',zmid+zdim/2._p_
  if(myid.eq.0)  write(*,*) 'magnetic location ', 'r_axis=',xmaxis, 'z_axis=',zmaxis, 'baxis=',baxis
  if(myid.eq.0) write(*,*) 'vacuum magnetic field at rcenter=',btorus, 'rcenter=',rmajor_mk
  if(myid.eq.0) write(*,*) 'psi_axis=',psi_axis,'psi_lcfs=',psi_lcfs
  if(myid.eq.0)   write(*,*) 'np_lcfs=',np_lcfs
  !write(*,*) 'x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)=', x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  if(myid.eq.0)   write(*,*) 'Cyclotron angular frequency of Deuterium ion at magnetic axis (10^6 rad/s)=', &
       & fpsi(1)/xmaxis*1.6022d-19/3.3452d-27/1.d6

  if(fpsi(1)>0) then 
     sign_bphi=1._p_
     if(myid.eq.0)    write(*,*) 'bphi>0'
  else
     sign_bphi=-1._p_
     if(myid.eq.0)      write(*,*) 'bphi<0'
  endif

  if(psi_lcfs>psi_axis) then
     if(myid.eq.0)      write(*,*) 'Iphi<0'
  else
     if(myid.eq.0)      write(*,*) 'Iphi>0'
  endif

  if(myid.eq.0) then
     open(321,file='lcfs.txt')
     do i=1,np_lcfs
        write(321,*) x_lcfs(i), z_lcfs(i)
     enddo
     close(321)

     open(321,file='limiter.txt')
     do i=1,nlim
        write(321,*) rlim(i), zlim(i)
     enddo
     close(321)
     write(*,*) 'wall r_min=',minval(rlim), 'wall r_max=',maxval(rlim), 'wall z_min=',minval(zlim), 'wall z_max=',maxval(zlim)
  endif
end subroutine read_gfile
