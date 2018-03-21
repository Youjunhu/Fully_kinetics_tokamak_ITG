!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE pputil

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ppinit, ppexit, init_pmove, end_pmove, pmove, pmove2
!
  INTEGER, SAVE :: me, nvp,npp,GCLR,TCLR
  INTEGER, SAVE :: pmove_tag=0
  INTEGER, SAVE :: TUBE_COMM,GRID_COMM
  REAL(8), DIMENSION(:), ALLOCATABLE, SAVE :: s_buf, r_buf
  logical, DIMENSION(:), ALLOCATABLE, SAVE :: s_buf2, r_buf2 !yj
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: s_counts, s_displ
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: r_counts, r_displ
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ipsend, iphole
!
CONTAINS
!
!===========================================================================
  SUBROUTINE init_pmove(xp, np, lz, ierr)
!
    INCLUDE 'mpif.h'
!
    REAL(8), DIMENSION(:), INTENT(in) :: xp
    INTEGER, INTENT(in) :: np
    REAL(8), INTENT(in) :: lz
    INTEGER, INTENT(out) :: ierr
!
!  Local vars
    INTEGER :: nsize, ksize
    INTEGER :: i, ip, iz, ih, iwork
    REAL(8) :: dzz, xt
    INTEGER, DIMENSION(0:nvp-1) :: isb
!
!----------------------------------------------------------------------
!              0.   Allocate fixed size arrays
!
    IF( .not. ALLOCATED(s_counts) ) ALLOCATE(s_counts(0:nvp-1))
    IF( .not. ALLOCATED(s_displ) ) ALLOCATE(s_displ(0:nvp-1))
    IF( .not. ALLOCATED(r_counts) ) ALLOCATE(r_counts(0:nvp-1))
    IF( .not. ALLOCATED(r_displ) ) ALLOCATE(r_displ(0:nvp-1))
!
!----------------------------------------------------------------------
!              1.  Construct send buffer
!
    dzz = lz / nvp
    s_counts = 0
    DO ip = 1,np
       xt = MODULO(xp(ip), lz)            !!! Assume periodicity
       iz = INT(xt/dzz)
       IF( iz .ne. GCLR )THEN
          s_counts(iz) = s_counts(iz)+1
       END IF
    END DO
    s_displ(0) = 0
    DO i=1,nvp-1
       s_displ(i) = s_displ(i-1) + s_counts(i-1)
    END DO
!
    nsize = sum(s_counts)
    IF( .not. ALLOCATED(s_buf) ) THEN
       ksize=2*nsize         ! To prevent too much futur reallocations
       ALLOCATE(s_buf(1:ksize))
       ALLOCATE(s_buf2(1:ksize))
       ALLOCATE(ipsend(1:ksize))
       ALLOCATE(iphole(1:ksize))
    ELSE IF ( SIZE(s_buf) .LT. nsize ) THEN
       DEALLOCATE(s_buf)
       DEALLOCATE(s_buf2)
       DEALLOCATE(ipsend)
       DEALLOCATE(iphole)
       ALLOCATE(s_buf(1:nsize))
       ALLOCATE(s_buf2(1:nsize))
       ALLOCATE(ipsend(1:nsize))
       ALLOCATE(iphole(1:nsize))
    END IF
!
!----------------------------------------------------------------------
!              2.  Construct (sorted) pointers to holes
!
    isb(0:nvp-1) = s_displ(0:nvp-1)
    ih = 0
    DO ip=1,np
       xt = MODULO(xp(ip), lz)            !!! Assume periodicity
       iz = INT(xt/dzz)
       IF( iz .ne. GCLR ) THEN
          isb(iz) = isb(iz)+1
          ipsend(isb(iz)) = ip
          ih = ih+1
          iphole(ih) = ip
       END IF
    END DO
!
!----------------------------------------------------------------------
!              3.  Construct receive buffer
!
    CALL MPI_ALLTOALL(s_counts, 1, MPI_INTEGER, &
         & r_counts, 1, MPI_INTEGER, TUBE_COMM, ierr)
    r_displ(0) = 0
    DO i=1,nvp-1
       r_displ(i) = r_displ(i-1) + r_counts(i-1)
    END DO
!
    nsize = sum(r_counts)
    IF( .not. ALLOCATED(r_buf) ) THEN
       ksize=2*nsize         ! To prevent too much futur reallocations
       ALLOCATE(r_buf(1:ksize))
       ALLOCATE(r_buf2(1:ksize))
    ELSE IF ( SIZE(r_buf) .LT. nsize ) THEN
       DEALLOCATE(r_buf)
       DEALLOCATE(r_buf2)
       ALLOCATE(r_buf(1:nsize))
       ALLOCATE(r_buf2(1:nsize))
    END IF
!
!  Check for part. array overflow
    ierr = 0
    nsize = np - sum(s_counts) + sum(r_counts)
    if( nsize .gt. size(xp) ) then
       write(*,*) 'PE', me, 'Particle array overflow'
       ierr = 1
    end if
!    call ppsum(ierr)
!
!----------------------------------------------------------------------
!              4.  Initialize tag
!
    pmove_tag = 101
!
  END SUBROUTINE init_pmove
!===========================================================================
  SUBROUTINE pmove(xp, np_old, np_new, ierr)
!
    INCLUDE 'mpif.h'
!
    REAL(8), DIMENSION(:), INTENT(inout) :: xp
    INTEGER, INTENT(in) :: np_old
    INTEGER, INTENT(out) :: np_new
    INTEGER, INTENT(out) :: ierr
!
!  Local vars
    INTEGER :: nsize
    INTEGER :: i, ip, iz, ih, isrc
    INTEGER :: nhole, mhole, nrrecv, nrsend, nptot_old, nptot
    INTEGER :: ind, count, tot_count, iwork
!
!  Local arrays
    INTEGER, DIMENSION(1:nvp) :: s_requ, r_requ, id_source
    INTEGER :: isrt, iend
    INTEGER :: status(MPI_STATUS_SIZE), arr_status(MPI_STATUS_SIZE,nvp)
!
!----------------------------------------------------------------------
!              1.  Fill send buffer
!
    DO i=0,nvp-1
       IF( s_counts(i) .GT. 0 ) THEN
          isrt = s_displ(i)+1
          iend = s_displ(i)+s_counts(i)
          s_buf(isrt:iend) = xp(ipsend(isrt:iend))
       END IF
    END DO
!----------------------------------------------------------------------
!              2.   Initiate non-blocking send/receive
!
    pmove_tag = pmove_tag+1
    nrrecv=0             !......... Start non-blocking receive
    DO i=0,nvp-1
       IF(r_counts(i) .GT. 0 ) THEN
          nrrecv=nrrecv+1
          id_source(nrrecv) = i
          isrt = r_displ(i)+1
          CALL MPI_IRECV(r_buf(isrt), r_counts(i), MPI_REAL8,&
               & i, pmove_tag, TUBE_COMM, r_requ(nrrecv), ierr)
       END IF
    END DO
    nrsend=0             !......... Start non-blocking SYNCHRONOUS send
    DO i=0,nvp-1
       IF(s_counts(i) .GT. 0 ) THEN
          nrsend=nrsend+1
          isrt = s_displ(i)+1
          CALL MPI_ISSEND(s_buf(isrt), s_counts(i), MPI_REAL8,&
               & i, pmove_tag, TUBE_COMM, s_requ(nrsend), ierr)
       END IF
    END DO
!
!----------------------------------------------------------------------
!              3.   Remove holes and compress part. arrays
!
    nhole = sum(s_counts)
    ip = np_old
    DO ih = nhole, 1, -1
       xp(iphole(ih)) = xp(ip)
       ip = ip-1
    END DO
    np_new = ip
!
!----------------------------------------------------------------------
!              4.   Store incoming part. to the part. arrays
!
    tot_count = 0
    DO i=1,nrrecv
       CALL MPI_WAITANY(nrrecv, r_requ, ind, status, ierr)
       isrc = id_source(ind)
       CALL MPI_GET_COUNT(status, MPI_REAL8, count, ierr)
       IF( count .ne. r_counts(isrc) ) THEN
          WRITE(*,*) 'PE',me, '  Counts mismatched from PE',isrc,&
               & count,  r_counts(isrc)
       END IF
       tot_count =  tot_count+count
    END DO
!
    IF( tot_count .GT. 0 ) THEN
       isrt = np_new + 1
       iend = np_new + tot_count
       xp(isrt:iend) = r_buf(1:tot_count)
       np_new = iend
    END IF
!
!----------------------------------------------------------------------
!              5.   Epilogue
!
!... Wait for any non-blocking comm. requests
    IF( nrsend.gt.0) CALL MPI_WAITALL(nrsend, s_requ, arr_status, ierr)
!
!... Check consistency
    ierr = 0
    CALL MPI_ALLREDUCE(np_old, nptot_old, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(np_new, nptot, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    IF( nptot.ne.nptot_old ) THEN
       IF(me.eq.0) WRITE(*,*) 'PMOVE: mismatch in total numbers:',&
            & nptot_old, nptot
       ierr = 1
    END IF
!    call ppsum(ierr)
!
!----------------------------------------------------------------------!
  END SUBROUTINE pmove


  SUBROUTINE pmove2(xp, np_old, np_new, ierr) !differ from pmove() in that it deal with logical array, instead of real array, by yjhu
!
    INCLUDE 'mpif.h'
!
    logical, DIMENSION(:), INTENT(inout) :: xp
    INTEGER, INTENT(in) :: np_old
    INTEGER, INTENT(out) :: np_new
    INTEGER, INTENT(out) :: ierr
!
!  Local vars
    INTEGER :: nsize
    INTEGER :: i, ip, iz, ih, isrc
    INTEGER :: nhole, mhole, nrrecv, nrsend, nptot_old, nptot
    INTEGER :: ind, count, tot_count, iwork
!
!  Local arrays
    INTEGER, DIMENSION(1:nvp) :: s_requ, r_requ, id_source
    INTEGER :: isrt, iend
    INTEGER :: status(MPI_STATUS_SIZE), arr_status(MPI_STATUS_SIZE,nvp)
!
!----------------------------------------------------------------------
!              1.  Fill send buffer
!
    DO i=0,nvp-1
       IF( s_counts(i) .GT. 0 ) THEN
          isrt = s_displ(i)+1
          iend = s_displ(i)+s_counts(i)
          s_buf2(isrt:iend) = xp(ipsend(isrt:iend))
       END IF
    END DO
!----------------------------------------------------------------------
!              2.   Initiate non-blocking send/receive
!
    pmove_tag = pmove_tag+1
    nrrecv=0             !......... Start non-blocking receive
    DO i=0,nvp-1
       IF(r_counts(i) .GT. 0 ) THEN
          nrrecv=nrrecv+1
          id_source(nrrecv) = i
          isrt = r_displ(i)+1
          CALL MPI_IRECV(r_buf2(isrt), r_counts(i), MPI_logical,&
               & i, pmove_tag, TUBE_COMM, r_requ(nrrecv), ierr)
       END IF
    END DO
    nrsend=0             !......... Start non-blocking SYNCHRONOUS send
    DO i=0,nvp-1
       IF(s_counts(i) .GT. 0 ) THEN
          nrsend=nrsend+1
          isrt = s_displ(i)+1
          CALL MPI_ISSEND(s_buf2(isrt), s_counts(i), MPI_logical,&
               & i, pmove_tag, TUBE_COMM, s_requ(nrsend), ierr)
       END IF
    END DO
!
!----------------------------------------------------------------------
!              3.   Remove holes and compress part. arrays
!
    nhole = sum(s_counts)
    ip = np_old
    DO ih = nhole, 1, -1
       xp(iphole(ih)) = xp(ip)
       ip = ip-1
    END DO
    np_new = ip
!
!----------------------------------------------------------------------
!              4.   Store incoming part. to the part. arrays
!
    tot_count = 0
    DO i=1,nrrecv
       CALL MPI_WAITANY(nrrecv, r_requ, ind, status, ierr)
       isrc = id_source(ind)
       CALL MPI_GET_COUNT(status, MPI_logical, count, ierr)
       IF( count .ne. r_counts(isrc) ) THEN
          WRITE(*,*) 'PE',me, '  Counts mismatched from PE',isrc,&
               & count,  r_counts(isrc)
       END IF
       tot_count =  tot_count+count
    END DO
!
    IF( tot_count .GT. 0 ) THEN
       isrt = np_new + 1
       iend = np_new + tot_count
       xp(isrt:iend) = r_buf2(1:tot_count)
       np_new = iend
    END IF
!
!----------------------------------------------------------------------
!              5.   Epilogue
!
!... Wait for any non-blocking comm. requests
    IF( nrsend.gt.0) CALL MPI_WAITALL(nrsend, s_requ, arr_status, ierr)
!
!... Check consistency
    ierr = 0
    CALL MPI_ALLREDUCE(np_old, nptot_old, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(np_new, nptot, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    IF( nptot.ne.nptot_old ) THEN
       IF(me.eq.0) WRITE(*,*) 'PMOVE: mismatch in total numbers:',&
            & nptot_old, nptot
       ierr = 1
    END IF
!    call ppsum(ierr)
!
!----------------------------------------------------------------------!
  END SUBROUTINE pmove2

!===========================================================================
  SUBROUTINE end_pmove(ierr)
!
    INCLUDE 'mpif.h'
    INTEGER, INTENT(OUT) :: ierr
!
!   Local vars
!----------------------------------------------------------------------!
  END SUBROUTINE end_pmove
!
!
 	SUBROUTINE ppinit(idproc,nproc,ntube,com1,com2)
!
     INCLUDE 'mpif.h'
     INTEGER, INTENT(IN) :: ntube
     INTEGER, INTENT(OUT) :: nproc
	INTEGER, INTENT(OUT) :: idproc,com1,com2
	INTEGER :: ierr,npp
!
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npp, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)
	nproc = npp
	IDPROC = me
! 
	GCLR=INT(me/ntube)
	TCLR=MOD(me,ntube)
!
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,GCLR,&
     	 &	TCLR,GRID_COMM,ierr)
	CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,TCLR,&
         &	GCLR,TUBE_COMM,ierr)
!
	com1=TUBE_COMM
	com2=GRID_COMM
	nvp=npp/ntube
!
	END SUBROUTINE ppinit
!
!===========================================================================
  SUBROUTINE ppexit
    INTEGER :: ierr
    CALL MPI_FINALIZE(ierr)
    STOP
  END SUBROUTINE ppexit
!
END MODULE pputil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
