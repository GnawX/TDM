PROGRAM TDM
  IMPLICIT NONE
  INTEGER,PARAMETER       :: q = SELECTED_REAL_KIND(16)
  INTEGER :: nat,n1,n2,n3,i,j,k,ii
  REAL(q) :: vec(3,3),mu(3),dv
  REAL(q), ALLOCATABLE :: vdata1(:),vdata2(:), r(:,:)
  
  CALL READ_CUBE_HEADER('homo.cube',nat,n1,n2,n3,vec)
  
  dv = ABS((vec(1,2)*vec(2,3)-vec(1,3)*vec(2,2))*vec(3,1) + &
       (vec(1,3)*vec(2,1)-vec(1,1)*vec(2,3))*vec(3,2) + &
       (vec(1,1)*vec(2,2)-vec(1,2)*vec(2,1))*vec(3,3))
       
  ALLOCATE(vdata1(n1*n2*n3),vdata2(n1*n2*n3))
  ALLOCATE(r(n1*n2*n3,3))
  
  CALL READ_CUBE_DATA('homo.cube',nat,n1,n2,n3,vdata1)
  CALL READ_CUBE_DATA('lumo.cube',nat,n1,n2,n3,vdata2)
  
  ii=0
  DO i=1,n1
     DO j=1,n2
        DO k=1,n3
           ii=ii+1
           r(ii,1)=i*vec(1,1)+j*vec(2,1)+k*vec(3,1)
           r(ii,2)=i*vec(1,2)+j*vec(2,2)+k*vec(3,2)
           r(ii,3)=i*vec(1,3)+j*vec(2,3)+k*vec(3,3)
        ENDDO
     ENDDO
  ENDDO
  
  mu(1)=SUM(vdata1*vdata2*r(:,1))*dv
  mu(2)=SUM(vdata1*vdata2*r(:,2))*dv
  mu(3)=SUM(vdata1*vdata2*r(:,3))*dv
  
  
  WRITE(*,FMT='(3F15.5)') (mu(i),i=1,3)

END PROGRAM

SUBROUTINE READ_CUBE_HEADER(FIN,NAT,N1,N2,N3,VEC)
  IMPLICIT NONE
  INTEGER,PARAMETER       :: q = SELECTED_REAL_KIND(16)
  CHARACTER(20), INTENT(IN)   :: FIN
  INTEGER, INTENT(OUT)    :: NAT,N1,N2,N3
  REAL(q), INTENT(OUT)    :: VEC(3,3)
  
  INTEGER :: i,iu 
  REAL(q) :: tmp
  
  iu = 10
  OPEN(UNIT=iu,FILE=FIN)
  READ(iu,*)
  READ(iu,*)
  READ(iu,*) NAT, tmp, tmp, tmp
  READ(iu,*) N1, (VEC(1,i), i=1,3)
  READ(iu,*) N2, (VEC(2,i), i=1,3)
  READ(iu,*) N3, (VEC(3,i), i=1,3)
  CLOSE(iu)
END SUBROUTINE READ_CUBE_HEADER

SUBROUTINE READ_CUBE_DATA(FIN,NAT,N1,N2,N3,VDATA)
  IMPLICIT NONE
  INTEGER,PARAMETER       :: q = SELECTED_REAL_KIND(16)
  CHARACTER(20), INTENT(IN)   :: FIN
  INTEGER, INTENT(IN)     :: NAT,N1,N2,N3
  REAL(q), INTENT(OUT)    :: VDATA(N1*N2*N3)
  
  INTEGER :: iu,i,j,k,m,mm,ii
  REAL(q) :: vtmp(N1,N2,N3)
  
  iu = 100
  
  OPEN(UNIT=iu,FILE=FIN)
  DO i=1,NAT+6
     READ(iu,*)
  ENDDO
  ii = 0
  DO i=1,N1
     DO j=1,N2
        m=N3/6
        mm=N3-m*6
        IF (mm==0) THEN
           DO k=1,m
              READ(iu,*) vtmp(i,j,6*k-5:6*k)
           ENDDO 
        ELSE
           DO k=1,m
              READ(iu,*) vtmp(i,j,6*k-5:6*k)
           ENDDO
           READ(iu,*) vtmp(i,j,m*6+1:N3)
        ENDIF
        ii = ii + 1
        VDATA((ii-1)*N3+1:ii*N3)=vtmp(i,j,:)
     ENDDO
  ENDDO
  CLOSE(iu)

END SUBROUTINE READ_CUBE_DATA
