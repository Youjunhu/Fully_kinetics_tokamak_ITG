subroutine arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis,myid)
  !replacing one point near the high-field side of the midplane by a point that is exactly on the high-field side of the midplane
  !then arrange the arrays x_lcfs_new and z_lcfs_new so that the starting point (x_lcfs_new(1),z_lcfs_new(1)) is on the high-field-side of the midplane
  !midplane is defined as the z=z_axis plane, where z_axis the z coordinate of the magnetic axis
  use precision,only:p_
  implicit none
  integer,intent(in):: np_lcfs,myid
  real(p_),intent(inout):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  real(p_),intent(in):: x_axis,z_axis
  real(p_):: x_lcfs_new(np_lcfs),z_lcfs_new(np_lcfs) 
  real(p_):: Rout
  real:: r_major,r_minor,eps,direction
  integer:: kk(1),k,i


  !write(*,*) ' Z values of the uppermost point of LCFS: ', maxval(z_lcfs)
  !write(*,*) ' Z values of the lowest point of LCFS: ' ,minval(z_lcfs)

  !set the starting point of LCFS to be at the low-field-side of the midplane
  !maxval(x_lcfs)
  !kk=maxloc(x_lcfs) !return the index of the array for which R is the largest,in order to determine the low-field side of the midplane
kk=minloc(x_lcfs) !return the index of the array for which R is the smallest,in order to determine the high-field side of the midplane
  k=kk(1)
  !k=10
  !write(*,*) 'index of the point on the lcfs that have the largest R, k=',k
  !if((z_lcfs(k+1)-z_axis)*(z_lcfs(k-1)-z_axis)>0) stop 'error in selecting the point on LCFS'

  call major_radius_on_midplane(np_lcfs,x_lcfs,z_lcfs,x_axis,z_axis,Rout)

  !  r_major=(maxval(x_lcfs)+minval(x_lcfs))/2._p_ !by definition
  !  r_minor=(maxval(x_lcfs)-minval(x_lcfs))/2._p_ !by definition

  ! write(*,*) 'inverse aspect ratio of LCFS is (definied at low-field side) ', (Rout-x_axis)/x_axis
  !  write(*,*) 'r_axis=',x_axis, 'r_major=', r_major, 'r_minor=',r_minor
  !  eps= r_minor/r_major !standard definition
  ! write(*,*) 'inverse aspect ratio of LCFS (i.e., r_minor/r_major) is ', eps
  ! write(*,*) 'ellipticity (elongation) of LCFS is ', (maxval(z_lcfs)-minval(z_lcfs))/2._p_/r_minor
  !write(*,*) 'upper triangularity of LCFS is ', (r_major-x_lcfs(maxloc(z_lcfs)))/r_minor, &
  !          & 'lower triangularity of LCFS is ', (r_major-x_lcfs(minloc(z_lcfs)))/r_minor
  !replace one point of LCFS with the new point
  x_lcfs(k)=Rout
  z_lcfs(k)=z_axis

  !arrange the arrays so that x_lcfs_new and z_lcfs_new start from the low/high-field-side of the midplane
  do i=1,np_lcfs
     if(k+i-1.le.np_lcfs) then
        x_lcfs_new(i)=x_lcfs(k+i-1)
        z_lcfs_new(i)=z_lcfs(k+i-1)
     else
        x_lcfs_new(i)=x_lcfs(k+i-np_lcfs)
        z_lcfs_new(i)=z_lcfs(k+i-np_lcfs)
     endif
  enddo

  !use x_lcfs and z_lcfs to store the new data
  x_lcfs=x_lcfs_new
  z_lcfs=z_lcfs_new

  !check wheter the direction of the sequecne (r(i),z(i)) with i increasing is clockwise or anticlockwise when viewed along grad_phi direction, if clockwise, switch it to anticlockwise
  !This is achieved by using the determination of the direction matrix (a well known method in graphic theory).
  !Because the contours of Psi considered here are always convex polygons (instead of concave polygons), we can select any vertex on the curve to calculate the direction matrix. (refer to wikipedia about the direction matrix)
  direction=(x_lcfs(2)-x_lcfs(1))*(z_lcfs(3)-z_lcfs(1))-(x_lcfs(3)-x_lcfs(1))*(z_lcfs(2)-z_lcfs(1))
  if(direction .lt. 0.) then
     if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) &
          & with i increasing is clockwise, switch it to anticlockwise'
           
     do i=1,np_lcfs !switch it to anticlockwise
        x_lcfs(i)=x_lcfs_new(np_lcfs+1-i)
        z_lcfs(i)=z_lcfs_new(np_lcfs+1-i)
     enddo
  else if (direction .gt. 0.) then
     if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) with i increasing is anticlockwise'
  else
     stop 'the three vertex (points) used in calculating the direction matrix is collinear'
  endif


  if (myid.eq.0) then
     open(123,file='lcfs2.txt')
     do i=1,np_lcfs
        write(123,*) x_lcfs(i),z_lcfs(i)
     enddo
     close(123)
  endif

end subroutine arrange_lcfs



subroutine major_radius_on_midplane(mpoloidal,rs,zs,r_axis,z_axis,Rout)
  !given a flux surface, this subroutine determines the major radius of the point on the low/high-field side of the middle plane
  use precision,only:p_
  use interpolate_module,only: linear_1d_interpolation_tmp

  implicit none
  integer,intent(in):: mpoloidal
  real(p_),intent(in):: rs(mpoloidal), zs(mpoloidal),r_axis,z_axis
  real(p_),intent(out):: Rout !the major radius of the point on the low/high-field-side of the mid-plane
  integer:: i,j,k1,k2,n
  real(p_):: r_select(mpoloidal),z_select(mpoloidal)
  !real(p_),dimension(:),allocatable:: x,z,tmp_y2
  real(p_):: tmp
 
  n=1
  do i=1,mpoloidal-1 
!     if(rs(i).gt.r_axis) then !select the low-field side (i.e., the part with r larger than r_axis) of a flux surface
     if(rs(i).lt.r_axis) then !select the high-field side (i.e., the part with r less than r_axis) of a flux surface
        r_select(n)=rs(i)
        z_select(n)=zs(i)
        n=n+1
     endif
  enddo
!  write(*,*) 'n-1= ', n-1
  !order the array according to the value of z_select
  do i=1,n-1
     do j=i+1,n-1
        if(z_select(j).le.z_select(i)) then
           !exchange the z value 
           tmp=z_select(i)
           z_select(i)=z_select(j)
           z_select(j)=tmp
           !also exchange the r value (I forgot this step in the older version, which cause a serious error)
           tmp=r_select(i)
           r_select(i)=r_select(j)
           r_select(j)=tmp
        endif
     enddo
  end do

  call linear_1d_interpolation_tmp(n-1,z_select,r_select,z_axis,Rout)  

end subroutine 







