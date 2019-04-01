PROGRAM TDM
IMPLICIT NONE
INTEGER,PARAMETER       :: q = SELECTED_REAL_KIND(10)
INTEGER :: nat,n1,n2,n3
REAL(q) :: vec(3,3)
REAL(q), ALLOCATABLE :: vdata(:)

CALL READ_CUBE_HEADER('homo.cube',nat,n1,n2,n3,vec)
ALLOCATE(vdata(n1*n2*n3))
CALL READ_CUBE_DATA('homo.cube',nat,n1,n2,n3,vdata)

END PROGRAM

SUBROUTINE READ_CUBE_HEADER(FIN,NAT,N1,N2,N3,VEC)
IMPLICIT NONE
INTEGER,PARAMETER       :: q = SELECTED_REAL_KIND(10)
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
INTEGER,PARAMETER       :: q = SELECTED_REAL_KIND(10)
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
         READ(iu,*) vtmp(i,j,m*6:3
      ENDIF
      ii = ii + 1
      VDATA((ii-1)*N3+1:ii*N3)=vtmp(i,j,:)
   ENDDO
ENDDO
CLOSE(iu)

END SUBROUTINE READ_CUBE_DATA
