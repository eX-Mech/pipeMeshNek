program pipeMeshNek

!==============================================================================
! geometry definitions
!
! (X) : border number
!
!      ____________________________L_(2)_______________________________
!      |                                                               |
!      |                                                               |
!      R (1)                                                           | (3)
!      |                                                               |
!  __ .| __ . __ . __ . __ . __ . __ . __ . __ . __ . __ . __ . __ . __| . __ 
!    (0,0,0)            axis of symmetry (4)
!
!==============================================================================


!==============================================================================
! variables definition

  IMPLICIT NONE
   ! geometry variables
   REAL(KIND=8)   :: R   ! pipe radius
   REAL(KIND=8)   :: L   ! pipe length
   REAL(KIND=8)   :: lSq ! side of the inner square
   REAL(KIND=8)   :: iR  ! inner element radius
   REAL(KIND=8)   :: oR  ! outer element radius
   ! mesh variables
   INTEGER        :: nR  ! number of elements on the radius (has to be > nTh/8)
   INTEGER        :: nL  ! number of elements along the pipe
   INTEGER        :: nTh ! number of circular sweeps (has to be >= 8 and a power of 2)
   REAL(KIND=8)   :: rR  ! element ratio along the radius
   REAL(KIND=8)   :: rL  ! element ratio along the pipe
   REAL(KIND=8)   :: sR  ! "summation"
   REAL(KIND=8)   :: dR  ! delta element
   REAL(KIND=8)   :: dL  ! delta element
   REAL(KIND=8)   :: de  ! delta element
   REAL(KIND=8)   :: rP0 ! radius "penalization", inner part
   REAL(KIND=8)   :: rP1 ! radius "penalization", outer part
   TYPE element
      INTEGER          :: num
      INTEGER          :: group
      CHARACTER(LEN=1) :: groupL = 'a'
      ! vertices
      REAL(KIND=8), DIMENSION(8) :: x
      REAL(KIND=8), DIMENSION(8) :: y
      REAL(KIND=8), DIMENSION(8) :: z
      ! faces
      LOGICAL,      DIMENSION(8) :: curvedEdge
      REAL(KIND=8), DIMENSION(8) :: curvedEdgeR = 0
      ! boundary conditions
      CHARACTER(LEN=2), DIMENSION(6)   :: bcType
      REAL(KIND=8),     DIMENSION(6,5) :: bcParameters = 0
   END TYPE element
   TYPE(element), DIMENSION(:), ALLOCATABLE  :: elem ! elements constituting the mesh
   INTEGER        :: nEl ! number of elements (total)
   INTEGER        :: nSq ! number of elements on the side of the inner square
   INTEGER        :: nPp ! number of elements in the central tube
   REAL(KIND=8)   :: alpha ! circular sweep angle
   INTEGER        :: nCurvedEdges ! number of curved edges
   INTEGER        :: nFpp  ! number of elements on one face of the            pipe
   INTEGER        :: nFpp4 ! number of elements on a quarter of a face of the pipe
   INTEGER        :: nPolynom ! Polynomial degree (-> nPolynom+1 Grid points)
   REAL(KIND=8), DIMENSION(20) :: xGridNodes ! (Local) Lobatto node-positions
   INTEGER        :: nsteps, iostep
   REAL(KIND=8)   :: Re, dt
   REAL(KIND=8)   :: reTau ! Friction-Reynolds Number (Re_tau)
   REAL(KIND=8)   :: deltaR, deltaRmin, deltaRmax, thMin, thMax ! wall-unit mesh-quality measures
   REAL(KIND=8)   :: deltaZmin, deltaZmax, deltaRThmin, deltaRThmax ! 
   INTEGER        :: np1, np10 ! number of gridPoints below 1 (10) wall units
   ! miscellaneous
   LOGICAL        :: debugFlag
   INTEGER        :: i, j, k, m, row, col
   REAL(KIND=8)   :: xTmp1, xTmp2
   REAL(KIND=8)   :: yTmp1, yTmp2
   INTEGER        :: fid1 = 100, fid3d = 103, fid2d = 102
   LOGICAL        :: existFlag
   CHARACTER(LEN=24) :: nameRea
   REAL(KIND=8)   :: PI = 4d0*DATAN(1d0)
   REAL(KIND=8)   :: SQ22 = SQRT(2d0)/2d0


!==============================================================================
! read input file 'INPUTgeometry'

   INQUIRE (FILE='INPUTgeometry', EXIST=existFlag)
   IF (.NOT.existFlag) THEN
      WRITE(*,*) '*************************************'
      WRITE(*,*) '*** ERROR:                        ***'
      WRITE(*,*) '*** File not found                ***'
      WRITE(*,*) '*** INPUTgeometry                 ***'
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ELSE
      OPEN(UNIT=fid1, FILE='INPUTgeometry', STATUS='old', FORM='formatted', ACTION='read')
   ENDIF

   READ(fid1,*) R
   READ(fid1,*) L
   READ(fid1,*) ! jump one line
   READ(fid1,*) nR
   READ(fid1,*) rP0
   READ(fid1,*) rP1
   READ(fid1,*) nL
   READ(fid1,*) nTh
   READ(fid1,*) ! jump one line
   READ(fid1,*) rR
   READ(fid1,*) rL
   READ(fid1,*) ! jump one line
   READ(fid1,*) Re
   READ(fid1,*) dt
   READ(fid1,*) nsteps
   READ(fid1,*) iostep
   READ(fid1,*) ! jump one line
   READ(fid1,*) nPolynom 
   READ(fid1,*) reTau 
   READ(fid1,*) ! jump one line
   READ(fid1,*) debugFlag

   CLOSE(fid1)


!==============================================================================
! check the input

   ! check if the number of circular sweeps is >= 8 and a power of 2
   !
   IF ( nTh.LT.8 .OR. MODULO(LOG(nTh*1d0)/LOG(2d0),1d0).GT.0d0 ) THEN
      WRITE(*,*) '******************************************'
      WRITE(*,*) '*** ERROR:                             ***'
      WRITE(*,*) '*** The number of circular sweeps      ***'
      WRITE(*,*) '*** (nTh) has to be >= 8               ***'
      WRITE(*,*) '*** and a power of 2.                  ***'
      WRITE(*,*) '******************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF

   ! check if nR > nTh/8
   !
   IF ( nR .LE. nTh/8 ) THEN
      WRITE(*,*) '******************************************'
      WRITE(*,*) '*** ERROR:                             ***'
      WRITE(*,*) '*** The number of elements on the      ***'
      WRITE(*,*) '*** radius (nR) has to be >= nTh/8.    ***'
      WRITE(*,*) '******************************************'
      WRITE(*,*) 'STOP.'
      STOP
   ENDIF

!==============================================================================
! do some preliminary computations
   
   ! mesh
   !
   alpha = 2d0*PI/nTH

   nSq = nTh / 8

   nFpp  = 4*nSq**2 + (nR-nSq)*nTh
   nFpp4 = nFpp / 4

   nPp = nFpp * nL

   nEl = nPp

   ! ratios
   !
   sR = 0d0
   DO i = 1, (nR-nSq)
      sR = sR + rR**i
   ENDDO
   dR  = R/(nSq+sR)
   lSq = nSq*dR

   sR = 0d0
   DO i = 1, nL
      sR = sR + rL**(i-1)
   ENDDO
   dL = L/sR

   CALL lobatto_set(nPolynom, xGridNodes)
   np1=-1 ! outmost grid point in counting algorithm is wall!
   np10=-1 ! outmost grid point in counting algorithm is wall!
   deltaR = 0
   deltaRmin = R
   deltaRThmin = R*alpha
   deltaZmin = L
   thMin = PI/2.
   deltaRmax = 0
   deltaRThmax = 0
   deltaZmax = 0
   thMax = PI/2.

   ALLOCATE( elem(nEl) )

!==============================================================================
! write a little output to screen

   WRITE(*,*)
   WRITE(*,'(a)') '--> Geometry data:'
   WRITE(*,*)
   WRITE(*,'(4x,a4,1x,f7.3)') 'R =', R
   WRITE(*,'(4x,a4,1x,f7.3)') 'L =', L
   WRITE(*,*)
   WRITE(*,'(4x,a4,1x,f7.3,3x,a6,1x,f7.3,1x,a4,1x,f7.3)') 'rR =', rR, '| e0 =', dR, 'e1 =', rR**(nR  -1)*dR
   WRITE(*,'(4x,a4,1x,f7.3,3x,a6,1x,f7.3,1x,a4,1x,f7.3)') 'rL =', rL, '| e0 =', dL, 'e1 =', rL**(nL  -1)*dL
   WRITE(*,*)
   WRITE(*,*)
   WRITE(*,'(a)') '--> Element data:'
   WRITE(*,*)
   WRITE(*,'(4x,a6,1x,i10)') 'nFpp =', nFpp
   WRITE(*,*)
   !WRITE(*,'(4x,a6,1x,i10)') 'nPp  =', nPp
   !WRITE(*,*)
   WRITE(*,'(4x,a6,1x,i10)') 'nEl  =', nEl
   WRITE(*,*)

!==============================================================================
! create the mesh

   i = 0

   ! 1/4 of first face
   !
   ! SQUARE PART
   !
   DO row = 1, nSq
      !
      DO col = 1, nSq
         !
         i = i + 1
         !
         elem(i)%num = i
         elem(i)%group = 1

         elem(i)%curvedEdge(:) = .FALSE.

         elem(i)%bcType(:) = 'E'
         elem(i)%bcType(5) = 'v'

         ! front face
         elem(i)%bcParameters(5,1) = i-nFpp
         elem(i)%bcParameters(5,2) = 6
         ! back face
         elem(i)%bcParameters(6,1) = i+nFpp
         elem(i)%bcParameters(6,2) = 5
         ! right face
         IF ( col .NE. nSq ) THEN
            elem(i)%bcParameters(3,1) = i+1
            elem(i)%bcParameters(3,2) = 1
         ELSE
            elem(i)%bcParameters(3,1) = nSq**2+nSq
            elem(i)%bcParameters(3,2) = 1
         ENDIF
         ! left face
         IF ( col .NE. 1 ) THEN
            elem(i)%bcParameters(1,1) = i-1
            elem(i)%bcParameters(1,2) = 3
         ELSE
            elem(i)%bcParameters(1,1) = nFpp4+1
            elem(i)%bcParameters(1,2) = 4
         ENDIF
         ! bottom face
         IF ( row .NE. 1 ) THEN
            elem(i)%bcParameters(4,1) = i-nSq
            elem(i)%bcParameters(4,2) = 2
         ELSE
            elem(i)%bcParameters(4,1) = 1-nFpp4 + (col-1)*nSq
            elem(i)%bcParameters(4,2) = 1
         ENDIF
         ! top face
         IF ( row .NE. nSq ) THEN
            elem(i)%bcParameters(2,1) = i+nSq
            elem(i)%bcParameters(2,2) = 4
         ELSE
            elem(i)%bcParameters(2,1) = nSq**2+2*nSq-col+1
            elem(i)%bcParameters(2,2) = 1
         ENDIF

         IF ( col > row ) THEN
            ! lower right part

            iR = lSq/nSq * (col-1) * (rP0 + (rP1-rP0)/nSq*(col-1))
            oR = lSq/nSq *  col    * (rP0 + (rP1-rP0)/nSq* col   )

            elem(i)%curvedEdge(1) = .TRUE.
            elem(i)%curvedEdge(3) = .TRUE.
            elem(i)%curvedEdge(5) = .TRUE.
            elem(i)%curvedEdge(7) = .TRUE.
            elem(i)%curvedEdgeR(1) = -iR
            elem(i)%curvedEdgeR(3) = oR
            elem(i)%curvedEdgeR(5) = -iR
            elem(i)%curvedEdgeR(7) = oR

            elem(i)%y(1) = elem(i-1)%y(4)
            elem(i)%y(2) = elem(i-1)%y(3)
!            elem(i)%y(3) = oR * DCOS( DASIN((elem(i-1)%x(3)+(oR-lSq/nSq*col)*SQ22)/oR) ) &
!                           - (oR - lSq/nSq*col)*SQ22 * 1.5*(nSq-col+2)
!            elem(i)%y(4) = oR * DCOS( DASIN((elem(i-1)%x(4)+(oR-lSq/nSq*col)*SQ22)/oR) ) &
!                           - (oR - lSq/nSq*col)*SQ22 * 1.5*(nSq-col+2)
            elem(i)%y(3) = oR * DCOS( DASIN(elem(i-1)%x(3)/oR) ) &
                           - (oR*DCOS(DASIN(lSq/nSq*col*SQ22/oR)) - lSq/nSq*col*SQ22)
            elem(i)%y(4) = oR * DCOS( DASIN(elem(i-1)%x(4)/oR) ) &
                           - (oR*DCOS(DASIN(lSq/nSq*col*SQ22/oR)) - lSq/nSq*col*SQ22)

            elem(i)%x(1) = elem(i-1)%x(4)
            elem(i)%x(2) = elem(i-1)%x(3)
            elem(i)%x(3) = elem(i)%x(2)
            elem(i)%x(4) = elem(i)%x(1)

            if ( debugFlag ) then
               write(*,*) 'lr ', row, col, rP0+(rP1-rP0)/nSq*col
            endif

         ELSEIF ( row > col ) THEN
            ! upper left part

            iR = lSq/nSq * (row-1) * (rP0 + (rP1-rP0)/nSq*(row-1))
            oR = lSq/nSq *  row    * (rP0 + (rP1-rP0)/nSq* row   )

            elem(i)%curvedEdge(2) = .TRUE.
            elem(i)%curvedEdge(4) = .TRUE.
            elem(i)%curvedEdge(6) = .TRUE.
            elem(i)%curvedEdge(8) = .TRUE.
            elem(i)%curvedEdgeR(2) = oR
            elem(i)%curvedEdgeR(4) = -iR
            elem(i)%curvedEdgeR(6) = oR
            elem(i)%curvedEdgeR(8) = -iR

            elem(i)%y(1) = elem(i-nSq)%y(2)
            elem(i)%y(2) = elem(i)%y(1)
            elem(i)%y(3) = elem(i-nSq)%y(3)
            elem(i)%y(4) = elem(i)%y(3)

            elem(i)%x(1) = elem(i-nSq)%x(2)
!            elem(i)%x(2) = oR * DSIN( DACOS((elem(i)%y(2)+(oR-lSq/nSq*row)*SQ22)/oR) ) &
!                           - (oR - lSq/nSq*row)*SQ22 * 1.5*(nSq-row+2)
!            elem(i)%x(3) = oR * DSIN( DACOS((elem(i)%y(4)+(oR-lSq/nSq*row)*SQ22)/oR) ) &
!                           - (oR - lSq/nSq*row)*SQ22 * 1.5*(nSq-row+2)
            elem(i)%x(2) = oR * DSIN( DACOS(elem(i)%y(2)/oR) ) &
                           - (oR*DCOS(DASIN(lSq/nSq*row*SQ22/oR)) - lSq/nSq*row*SQ22)
            elem(i)%x(3) = oR * DSIN( DACOS(elem(i)%y(4)/oR) ) &
                           - (oR*DCOS(DASIN(lSq/nSq*row*SQ22/oR)) - lSq/nSq*row*SQ22)
            elem(i)%x(4) = elem(i-nSq)%x(3)

            if ( debugFlag ) then
               write(*,*) 'ul ', row, col, rP0+(rP1-rP0)/nSq*row
            endif

         ELSE
            ! diagonal

            oR = lSq/nSq * row * (rP0 + (rP1-rP0)/nSq*row)

            elem(i)%curvedEdge(2) = .TRUE.
            elem(i)%curvedEdge(3) = .TRUE.
            elem(i)%curvedEdge(6) = .TRUE.
            elem(i)%curvedEdge(7) = .TRUE.
            elem(i)%curvedEdgeR(2) = oR
            elem(i)%curvedEdgeR(3) = oR
            elem(i)%curvedEdgeR(6) = oR
            elem(i)%curvedEdgeR(7) = oR

            elem(i)%y(1) = lSq/nSq * (row-1) * SQ22
            elem(i)%y(2) = elem(i)%y(1)
            elem(i)%y(3) = lSq/nSq * row     * SQ22
!            elem(i)%y(4) = oR * DCOS( DASIN((lSq/nSq*(row-1)*SQ22+(oR-lSq/nSq*col)*SQ22)/oR) ) &
!                           - (oR - lSq/nSq*row)*SQ22 * 1.5*(nSq-row+2)
            elem(i)%y(4) = oR * DCOS( DASIN(lSq/nSq*(row-1)*SQ22/oR) ) &
                           - (oR*DCOS(DASIN(lSq/nSq*col*SQ22/oR)) - lSq/nSq*col*SQ22)

            elem(i)%x(1) = lSq/nSq * (row-1) * SQ22
!            elem(i)%x(2) = oR * DSIN( DACOS((elem(i)%y(2)+(oR-lSq/nSq*col)*SQ22)/oR) ) &
!                           - (oR - lSq/nSq*row)*SQ22 * 1.5*(nSq-row+2)
            elem(i)%x(2) = oR * DSIN( DACOS(elem(i)%y(2)/oR) ) &
                           - (oR*DCOS(DASIN(lSq/nSq*row*SQ22/oR)) - lSq/nSq*row*SQ22)
            elem(i)%x(3) = lSq/nSq * row     * SQ22
            elem(i)%x(4) = elem(i)%x(1)

            if ( debugFlag ) then
               write(*,*) 'd  ', row, oR
            endif

         ENDIF

         elem(i)%y(5) = elem(i)%y(1)
         elem(i)%y(6) = elem(i)%y(2)
         elem(i)%y(7) = elem(i)%y(3)
         elem(i)%y(8) = elem(i)%y(4)

         elem(i)%x(5) = elem(i)%x(1)
         elem(i)%x(6) = elem(i)%x(2)
         elem(i)%x(7) = elem(i)%x(3)
         elem(i)%x(8) = elem(i)%x(4)

         elem(i)%z(1) = 0d0
         elem(i)%z(2) = 0d0
         elem(i)%z(3) = 0d0
         elem(i)%z(4) = 0d0
         elem(i)%z(5) = 0d0 + dL
         elem(i)%z(6) = 0d0 + dL
         elem(i)%z(7) = 0d0 + dL
         elem(i)%z(8) = 0d0 + dL

      ENDDO

   ENDDO
   !
   ! CIRCULAR PART
   !
   iR  = lSq*rP1
   oR  = lSq + dR*rR
   col = 1
   DO row = 1, nTh/8
      !
      i = i + 1
      !
      elem(i)%num = i
      elem(i)%group = 1

      elem(i)%curvedEdge(:) = .FALSE.
      elem(i)%curvedEdge(1) = .TRUE.
      elem(i)%curvedEdge(3) = .TRUE.
      elem(i)%curvedEdge(5) = .TRUE.
      elem(i)%curvedEdge(7) = .TRUE.
      elem(i)%curvedEdgeR(1) = -iR
      elem(i)%curvedEdgeR(3) = oR
      elem(i)%curvedEdgeR(5) = -iR
      elem(i)%curvedEdgeR(7) = oR

      elem(i)%bcType(:) = 'E'
      elem(i)%bcType(5) = 'v'

      ! front face
      elem(i)%bcParameters(5,1) = i-nFpp
      elem(i)%bcParameters(5,2) = 6
      ! back face
      elem(i)%bcParameters(6,1) = i+nFpp
      elem(i)%bcParameters(6,2) = 5
      ! right face
      elem(i)%bcParameters(3,1) = i+nTh/4
      elem(i)%bcParameters(3,2) = 1
      ! left face
      elem(i)%bcParameters(1,1) = row*nSq
      elem(i)%bcParameters(1,2) = 3
      ! bottom face
      IF ( row .NE. 1 ) THEN
         elem(i)%bcParameters(4,1) = i-1
         elem(i)%bcParameters(4,2) = 2
      ELSE
         elem(i)%bcParameters(4,1) = 1-nFpp4+nSq**2+nTh/8
         elem(i)%bcParameters(4,2) = 2
      ENDIF
      ! top face
      elem(i)%bcParameters(2,1) = i+1
      elem(i)%bcParameters(2,2) = 4

      elem(i)%y(1) = elem(nSq*row)%y(4)
      elem(i)%y(2) = elem(nSq*row)%y(3)
      elem(i)%y(3) = oR * COS(alpha*row)
      elem(i)%y(4) = oR * COS(alpha*(row-1))
      elem(i)%y(5) = elem(i)%y(1)
      elem(i)%y(6) = elem(i)%y(2)
      elem(i)%y(7) = elem(i)%y(3)
      elem(i)%y(8) = elem(i)%y(4)

      elem(i)%x(1) = elem(nSq*row)%x(4)
      elem(i)%x(2) = elem(nSq*row)%x(3)
      elem(i)%x(3) = oR * SIN(alpha*row)
      elem(i)%x(4) = oR * SIN(alpha*(row-1))
      elem(i)%x(5) = elem(i)%x(1)
      elem(i)%x(6) = elem(i)%x(2)
      elem(i)%x(7) = elem(i)%x(3)
      elem(i)%x(8) = elem(i)%x(4)

      elem(i)%z(1) = 0d0
      elem(i)%z(2) = 0d0
      elem(i)%z(3) = 0d0
      elem(i)%z(4) = 0d0
      elem(i)%z(5) = 0d0 + dL
      elem(i)%z(6) = 0d0 + dL
      elem(i)%z(7) = 0d0 + dL
      elem(i)%z(8) = 0d0 + dL

   ENDDO
   row = 1
   DO col = 1, nTh/8
      !
      i = i + 1
      !
      elem(i)%num = i
      elem(i)%group = 1

      elem(i)%curvedEdge(:) = .FALSE.
      elem(i)%curvedEdge(1) = .TRUE.
      elem(i)%curvedEdge(3) = .TRUE.
      elem(i)%curvedEdge(5) = .TRUE.
      elem(i)%curvedEdge(7) = .TRUE.
      elem(i)%curvedEdgeR(1) = -iR
      elem(i)%curvedEdgeR(3) = oR
      elem(i)%curvedEdgeR(5) = -iR
      elem(i)%curvedEdgeR(7) = oR

      elem(i)%bcType(:) = 'E'
      elem(i)%bcType(5) = 'v'

      ! front face
      elem(i)%bcParameters(5,1) = i-nFpp
      elem(i)%bcParameters(5,2) = 6
      ! back face
      elem(i)%bcParameters(6,1) = i+nFpp
      elem(i)%bcParameters(6,2) = 5
      ! right face
      elem(i)%bcParameters(3,1) = i+nTh/4
      elem(i)%bcParameters(3,2) = 1
      ! left face
      elem(i)%bcParameters(1,1) = nSq**2-col+1
      elem(i)%bcParameters(1,2) = 2
      ! bottom face
      elem(i)%bcParameters(4,1) = i-1
      elem(i)%bcParameters(4,2) = 2
      ! top face
      IF ( col .NE. nSq ) THEN
         elem(i)%bcParameters(2,1) = i+1
         elem(i)%bcParameters(2,2) = 4
      ELSE
         elem(i)%bcParameters(2,1) = nFpp4+nSq**2+1
         elem(i)%bcParameters(2,2) = 4
      ENDIF

      elem(i)%y(1) = elem(nSq**2+1-col)%y(3)
      elem(i)%y(2) = elem(nSq**2+1-col)%y(2)
      elem(i)%y(3) = oR * COS(PI/4+alpha*col)
      elem(i)%y(4) = oR * COS(PI/4+alpha*(col-1))
      elem(i)%y(5) = elem(i)%y(1)
      elem(i)%y(6) = elem(i)%y(2)
      elem(i)%y(7) = elem(i)%y(3)
      elem(i)%y(8) = elem(i)%y(4)

      elem(i)%x(1) = elem(nSq**2+1-col)%x(3)
      elem(i)%x(2) = elem(nSq**2+1-col)%x(2)
      elem(i)%x(3) = oR * SIN(PI/4+alpha*col)
      elem(i)%x(4) = oR * SIN(PI/4+alpha*(col-1))
      elem(i)%x(5) = elem(i)%x(1)
      elem(i)%x(6) = elem(i)%x(2)
      elem(i)%x(7) = elem(i)%x(3)
      elem(i)%x(8) = elem(i)%x(4)

      elem(i)%z(1) = 0d0
      elem(i)%z(2) = 0d0
      elem(i)%z(3) = 0d0
      elem(i)%z(4) = 0d0
      elem(i)%z(5) = 0d0 + dL
      elem(i)%z(6) = 0d0 + dL
      elem(i)%z(7) = 0d0 + dL
      elem(i)%z(8) = 0d0 + dL

   ENDDO
   !
   ! other columns of the circular part
   !
   de = dR*rR**2
   DO row = 1, (nR-nSq-1)

      !iR  = lSq +      dR
      !oR  = lSq + de + dR
      iR  = oR
      oR  = iR + de

      !dR  = de + dR
      de  = rR*de

      DO col = 1, nTh/4
         !
         i = i + 1
         !
         elem(i)%num = i
         elem(i)%group = 1
      
         elem(i)%curvedEdge(:) = .FALSE.
         elem(i)%curvedEdge(3) = .TRUE.
         elem(i)%curvedEdge(7) = .TRUE.
         elem(i)%curvedEdge(1) = .TRUE.
         elem(i)%curvedEdge(5) = .TRUE.
         elem(i)%curvedEdgeR(3) = oR
         elem(i)%curvedEdgeR(7) = oR
         elem(i)%curvedEdgeR(1) = -iR
         elem(i)%curvedEdgeR(5) = -iR
      
         elem(i)%bcType(:) = 'E'
         elem(i)%bcType(5) = 'v'

         ! front face
         elem(i)%bcParameters(5,1) = i-nFpp
         elem(i)%bcParameters(5,2) = 6
         ! back face
         elem(i)%bcParameters(6,1) = i+nFpp
         elem(i)%bcParameters(6,2) = 5
         ! right face
         elem(i)%bcParameters(3,1) = i+nTh/4
         elem(i)%bcParameters(3,2) = 1
         ! left face
         elem(i)%bcParameters(1,1) = i-nTh/4
         elem(i)%bcParameters(1,2) = 3
         ! bottom face
         IF ( col .NE. 1 ) THEN
            elem(i)%bcParameters(4,1) = i-1
            elem(i)%bcParameters(4,2) = 2
         ELSE
            elem(i)%bcParameters(4,1) = -nFpp4+nSq**2+nTh/4*row
            elem(i)%bcParameters(4,2) = 2
         ENDIF
         ! top face
         IF ( col .NE. nTh/4 ) THEN
            elem(i)%bcParameters(2,1) = i+1
            elem(i)%bcParameters(2,2) = 4
         ELSE
            elem(i)%bcParameters(2,1) = nFpp4+nSq**2+nTh/4*row+1
            elem(i)%bcParameters(2,2) = 4
         ENDIF
      
         elem(i)%y(1) = elem(i-nTh/4)%y(4)
         elem(i)%y(2) = elem(i-nTh/4)%y(3)
         elem(i)%y(3) = oR * COS(alpha*col)
         elem(i)%y(4) = oR * COS(alpha*(col-1))
         elem(i)%y(5) = elem(i)%y(1)
         elem(i)%y(6) = elem(i)%y(2)
         elem(i)%y(7) = elem(i)%y(3)
         elem(i)%y(8) = elem(i)%y(4)
      
         elem(i)%x(1) = elem(i-nTh/4)%x(4)
         elem(i)%x(2) = elem(i-nTh/4)%x(3)
         elem(i)%x(3) = oR * SIN(alpha*col)
         elem(i)%x(4) = oR * SIN(alpha*(col-1))
         elem(i)%x(5) = elem(i)%x(1)
         elem(i)%x(6) = elem(i)%x(2)
         elem(i)%x(7) = elem(i)%x(3)
         elem(i)%x(8) = elem(i)%x(4)
      
         elem(i)%z(1) = 0d0
         elem(i)%z(2) = 0d0
         elem(i)%z(3) = 0d0
         elem(i)%z(4) = 0d0
         elem(i)%z(5) = 0d0 + dL
         elem(i)%z(6) = 0d0 + dL
         elem(i)%z(7) = 0d0 + dL
         elem(i)%z(8) = 0d0 + dL
      
      ENDDO

   ENDDO
   !
   ! add "Wall" boundary condition on the external elements
   !
   DO j = nFpp4-nTh/4+1, nFpp4

      elem(j)%bcType(3) = 'W'
   ENDDO

   ! calculate a few grid measures in wall-units
   ! Cf. mesh description in [El Khoury 13]
   ! delta-r^+ : distance between (pipe) wall and outmost gridpoint
   ! delta-R^+ : min and max
   ! delta-R^Theta+ : min and max
   ! delta-Z^+ : min and max
   ! And: element corner-angle (in degrees)
   ! theta : min and max
   
   ! OBS: accounts (only) for linear (!) element deformation!
   ! BUT the curved mesh faces are concentric for the outer element-rings!
   ! --> should be correct

   ! Loop over all elements
   DO j = 1, nFpp4

     ! Minimal and maximal element-corner angle
     xTmp1 = (elem(j)%x(4)-elem(j)%x(1))
     yTmp1 = (elem(j)%y(4)-elem(j)%y(1))

     xTmp2 = (elem(j)%x(2)-elem(j)%x(1))
     yTmp2 = (elem(j)%y(2)-elem(j)%y(1))

     thMin = min(thMin, DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     thMax = max(thMax, DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     xTmp1 = (elem(j)%x(3)-elem(j)%x(2))
     yTmp1 = (elem(j)%y(3)-elem(j)%y(2))

     thMin = min(thMin, PI-DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     thMax = max(thMax, PI-DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     xTmp2 = (elem(j)%x(4)-elem(j)%x(3))
     yTmp2 = (elem(j)%y(4)-elem(j)%y(3))

     thMin = min(thMin, PI-DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     thMax = max(thMax, PI-DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     xTmp1 = (elem(j)%x(1)-elem(j)%x(4))
     yTmp1 = (elem(j)%y(1)-elem(j)%y(4))

     thMin = min(thMin, DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ & 
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))

     thMax = max(thMax, DACOS((xTmp1*xTmp2 + yTmp1*yTmp2)/ &
       SQRT((xTmp1**2+yTmp1**2)*(xTmp2**2+yTmp2**2))))


    ! Iterate over all GLL-Nodes (only of interior elements)
    ! (Otherwise, delta-r-max would immediately be on the outer ring)
      IF (J .GT. nFpp4-nTh/4) THEN
        CYCLE
      ENDIF
      DO k = 2, nPolynom+1
        DO m = 2, nPolynom+1
        ! xTmp1 and xTmp2, respectively, 
        ! correspond to temporary points G and F in (GeoGebra-file) sketch.ggb

          xTmp1 = elem(j)%x(1) + (1.+xGridNodes(k))*.5*(elem(j)%x(4)-elem(j)%x(1))
          yTmp1 = elem(j)%y(1) + (1.+xGridNodes(k))*.5*(elem(j)%y(4)-elem(j)%y(1))

          xTmp2 = elem(j)%x(2) + (1.+xGridNodes(k))*.5*(elem(j)%x(3)-elem(j)%x(2))
          yTmp2 = elem(j)%y(2) + (1.+xGridNodes(k))*.5*(elem(j)%y(3)-elem(j)%y(2))

          deltaRmin = min(deltaRmin,SQRT( &
          ((xGridNodes(m)-xGridNodes(m-1))*.5*(xTmp2-xTmp1))**2 + &
           ((xGridNodes(m)-xGridNodes(m-1))*.5*(yTmp2-yTmp1))**2))

          deltaRmax = max(deltaRmax,SQRT( &
          ((xGridNodes(m)-xGridNodes(m-1))*.5*(xTmp2-xTmp1))**2 + &
           ((xGridNodes(m)-xGridNodes(m-1))*.5*(yTmp2-yTmp1))**2))
          

          ! Minimize/Maximize along both dimensions ..
          xTmp1 = elem(j)%x(1) + (1.+xGridNodes(k))*.5*(elem(j)%x(2)-elem(j)%x(1))
          yTmp1 = elem(j)%y(1) + (1.+xGridNodes(k))*.5*(elem(j)%y(2)-elem(j)%y(1))

          xTmp2 = elem(j)%x(4) + (1.+xGridNodes(k))*.5*(elem(j)%x(3)-elem(j)%x(4))
          yTmp2 = elem(j)%y(4) + (1.+xGridNodes(k))*.5*(elem(j)%y(3)-elem(j)%y(4))

          deltaRmin = min(deltaRmin,SQRT( &
          ((xGridNodes(m)-xGridNodes(m-1))*.5*(xTmp2-xTmp1))**2 + &
           ((xGridNodes(m)-xGridNodes(m-1))*.5*(yTmp2-yTmp1))**2))

          deltaRmax = max(deltaRmax,SQRT( &
          ((xGridNodes(m)-xGridNodes(m-1))*.5*(xTmp2-xTmp1))**2 + &
           ((xGridNodes(m)-xGridNodes(m-1))*.5*(yTmp2-yTmp1))**2))

        ENDDO
      ENDDO
   ENDDO

   ! Iterate one elements-row inwards from wall
   j = nFpp4-nTh/4+1
   DO WHILE ( (deltaR*reTau/R) .LT. 10. )
      ! NOTE: that the element vertices are double-counted (intentionally)!
      DO k = 1, nPolynom+1
        xTmp1 = elem(j)%x(4) + (1.+xGridNodes(k))*.5*(elem(j)%x(1)-elem(j)%x(4))
        yTmp1 = elem(j)%y(4) + (1.+xGridNodes(k))*.5*(elem(j)%y(1)-elem(j)%y(4))
        deltaR = R - SQRT(xTmp1**2 + yTmp1**2)

        IF ( (deltaR*reTau/R) .LT. 10. ) THEN
          nP10 = nP10 +1
          IF ( (deltaR*reTau/R) .LT. 1. ) THEN
            nP1 = nP1 +1
          ENDIF
        ENDIF
      ENDDO
      j = INT(elem(j)%bcParameters(1,1))

      ! For "interior face" BC-Type ('E') the connectivity seems to be not guaranteed
      IF (nP10 .GT. 1000) THEN
        WRITE(*,*) 'pipeMeshNek seems stuck while counting grid points'
        WRITE(*,*) 'STOP.'
        STOP
      ENDIF
   ENDDO


   j = nFpp4-nTh/4+2
   xTmp1 = elem(j)%x(1) + (1.+xGridNodes(nPolynom))*.5*(elem(j)%x(4)-elem(j)%x(1))
   yTmp1 = elem(j)%y(1) + (1.+xGridNodes(nPolynom))*.5*(elem(j)%y(4)-elem(j)%y(1))
   deltaR = R-SQRT(xTmp1**2 + yTmp1**2)

   DO k = 2, nPolynom+1
   xTmp1 = elem(j)%x(3) + (1.+xGridNodes(k))*.5*(elem(j)%x(4)-elem(j)%x(3))
   yTmp1 = elem(j)%y(3) + (1.+xGridNodes(k))*.5*(elem(j)%y(4)-elem(j)%y(3))

   xTmp2 = elem(j)%x(3) + (1.+xGridNodes(k-1))*.5*(elem(j)%x(4)-elem(j)%x(3))
   yTmp2 = elem(j)%y(3) + (1.+xGridNodes(k-1))*.5*(elem(j)%y(4)-elem(j)%y(3))

   deltaRThmin = min(deltaRThmin, &
     R*ABS(DATAN(yTmp1/xTmp1) - DATAN(yTmp2/xTmp2)))
   deltaRThmax = max(deltaRThmax, &
     R*ABS(DATAN(yTmp1/xTmp1) - DATAN(yTmp2/xTmp2)))
   ENDDO
   WRITE(*,*) 'delta r Plus (outmost grid) ', deltaR*reTau/R
   WRITE(*,*) '# grid points < delta r+=1  ', np1
   WRITE(*,*) '# grid points < delta r+=10 ', np10
   WRITE(*,*) 'Theta min [°]               ', thMin/PI*180.
   WRITE(*,*) 'Theta max [°]               ', thMax/PI*180.
   WRITE(*,*) 'delta r Plus min            ', deltaRmin*reTau/R
   ! squareroot of two for diagonals
   WRITE(*,*) 'delta r Plus max            ', SQRT(2d0)*deltaRmax*reTau/R 
   WRITE(*,*) 'delta R Theta Plus min      ', deltaRThmin*reTau/R
   WRITE(*,*) 'delta R Theta Plus max      ', deltaRThmax*reTau/R 

   !
   ! "mirror" the first 1/4 face on the other quarters
   !
   DO i = i+1, 2*nFpp4


      elem(i)%num = i
      elem(i)%group = 1

      elem(i)%curvedEdge(:)  = elem(i-nFpp4)%curvedEdge(:)
      elem(i)%curvedEdgeR(:) = elem(i-nFpp4)%curvedEdgeR(:)

      elem(i)%bcType(:) = elem(i-nFpp4)%bcType(:)
      elem(i)%bcParameters(:,1) = elem(i-nFpp4)%bcParameters(:,1)+nFpp4
      elem(i)%bcParameters(:,2) = elem(i-nFpp4)%bcParameters(:,2)

      elem(i)%y(:) = - elem(i-nFpp4)%x(:)
      elem(i)%x(:) =   elem(i-nFpp4)%y(:)
      elem(i)%z(:) =   elem(i-nFpp4)%z(:)

   ENDDO
   DO i = i, 3*nFpp4

      elem(i)%num = i
      elem(i)%group = 1

      elem(i)%curvedEdge(:)  = elem(i-2*nFpp4)%curvedEdge(:)
      elem(i)%curvedEdgeR(:) = elem(i-2*nFpp4)%curvedEdgeR(:)

      elem(i)%bcType(:) = elem(i-2*nFpp4)%bcType(:)
      elem(i)%bcParameters(:,1) = elem(i-2*nFpp4)%bcParameters(:,1)+2*nFpp4
      elem(i)%bcParameters(:,2) = elem(i-2*nFpp4)%bcParameters(:,2)

      elem(i)%y(:) = - elem(i-2*nFpp4)%y(:)
      elem(i)%x(:) = - elem(i-2*nFpp4)%x(:)
      elem(i)%z(:) =   elem(i-2*nFpp4)%z(:)

   ENDDO
   DO i = i, 4*nFpp4

      elem(i)%num = i
      elem(i)%group = 1

      elem(i)%curvedEdge(:)  = elem(i-3*nFpp4)%curvedEdge(:)
      elem(i)%curvedEdgeR(:) = elem(i-3*nFpp4)%curvedEdgeR(:)

      elem(i)%bcType(:) = elem(i-3*nFpp4)%bcType(:)
      elem(i)%bcParameters(:,1) = elem(i-3*nFpp4)%bcParameters(:,1)+3*nFpp4
      elem(i)%bcParameters(:,2) = elem(i-3*nFpp4)%bcParameters(:,2)

      elem(i)%y(:) =   elem(i-3*nFpp4)%x(:)
      elem(i)%x(:) = - elem(i-3*nFpp4)%y(:)
      elem(i)%z(:) =   elem(i-3*nFpp4)%z(:)

   ENDDO
   i = i - 1
   !
   ! correct boundary conditions at "glue side"
   !
   DO j = 1, nSq
      elem(j)%bcParameters(4,1) = 3*nFpp4+(j-1)*nSq+1
      elem(3*nFpp4+(j-1)*nSq+1)%bcParameters(1,1) = j
   ENDDO
   DO j = 1, (nR-nSq)
      elem(nSq**2+1+(j-1)*nTh/4)%bcParameters(4,1) = 3*nFpp4+nSq**2+nTh/4*j
      elem(3*nFpp4+nSq**2+nTh/4*j)%bcParameters(2,1) = nSq**2+1+(j-1)*nTh/4
   ENDDO
   !
	 ! advance the face to the end of the pipe
   !
   de = rL * dL
   deltaZmin = min(deltaZmin,de)
   deltaZmax = max(deltaZmax,de)
   DO j = 1, nL-1

      deltaZmin = min(deltaZmin,de)
      deltaZmax = max(deltaZmax,de)

      CALL advanceFace (1 + nFpp*(j-1), nFpp*j, de, elem)

      de = rL * de

      i = i + nFpp

      IF ( j .EQ. 1 ) THEN

         elem(1+nFpp:nFpp*2)%bcType(5) = 'E'

      ENDIF

   ENDDO
   !
   ! add PERIODIC boundary conditions on first and last face
   !
   DO j = i-nFpp+1, i

      ! write last face
      elem(j)%bcType(6) = 'P' ! O
      elem(j)%bcParameters(6,1) = j-nPp+nFpp
      ! and correct periodicity on the first face
      elem(j-nPp+nFpp)%bcType(5) = 'P' ! v
      elem(j-nPp+nFpp)%bcParameters(5,1) = j

   ENDDO

   WRITE(*,*) 'delta z Plus min            ', 0.5*(xGridNodes(nPolynom+1)-xGridNodes(nPolynom)) &
     *deltaZmin*reTau/R
   WRITE(*,*) 'delta z Plus max            ', 0.5*(xGridNodes(nPolynom/2+2)-xGridNodes(nPolynom/2+1)) &
     *deltaZmax*reTau/R 


!==============================================================================
! clean up

   ! eliminate useless parameters in boundary conditions
   !
   ! now we also have periodic elements
   !!!DO i = 1, nEl
   !!!   DO j = 1, 6
   !!!      IF ( elem(i)%bcType(j) .NE. 'E' ) THEN
   !!!         elem(i)%bcParameters(j,:) = 0
   !!!      ENDIF
   !!!   ENDDO
   !!!ENDDO

   ! count curved edges
   !
   nCurvedEdges = 0
   DO i = 1, nEl
      DO j = 1, 8
         IF ( elem(i)%curvedEdge(j) )  nCurvedEdges = nCurvedEdges + 1
      ENDDO
   ENDDO

   ! set to zero any "approximate" zero
   !
   DO i = 1, nEl
      DO j = 1, 8
         IF ( abs(elem(i)%x(j)) .LE. 1d-15 ) THEN
            elem(i)%x(j) = 0d0
         ENDIF
         IF ( abs(elem(i)%y(j)) .LE. 1d-15 ) THEN
            elem(i)%y(j) = 0d0
         ENDIF
         IF ( abs(elem(i)%z(j)) .LE. 1d-15 ) THEN
            elem(i)%z(j) = 0d0
         ENDIF
      ENDDO
   ENDDO


!==============================================================================
! initialize nameRea

   WRITE(nameRea,'(a)') 'base.rea'

   INQUIRE (FILE=trim(nameRea), EXIST=existFlag)
   IF (existFlag) THEN
      WRITE(*,*) '*************************************'
      !WRITE(*,*) '*** ERROR:                        ***'
      WRITE(*,*) '*** WARNING:                      ***'
      WRITE(*,*) '*** File already present          ***'
      WRITE(*,*) '*** ', trim(nameRea), '                      ***'
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'OVERWRITING.'
      !WRITE(*,*) 'STOP.'
      !STOP
   ENDIF
   !ELSE
      !OPEN(UNIT=fid3d, FILE=trim(nameRea), STATUS='new', ACTION='write')
      OPEN(UNIT=fid3d, FILE=trim(nameRea), ACTION='write')
      CALL initializeMeshFile(fid3d, Re, dt, nsteps, iostep, debugFlag)
   !ENDIF

   WRITE(nameRea,'(a)') 'base2d.rea'

   INQUIRE (FILE=trim(nameRea), EXIST=existFlag)
   IF (existFlag) THEN
      WRITE(*,*) '*************************************'
      !WRITE(*,*) '*** ERROR:                        ***'
      WRITE(*,*) '*** WARNING:                      ***'
      WRITE(*,*) '*** File already present          ***'
      WRITE(*,*) '*** ', trim(nameRea), '                    ***'
      WRITE(*,*) '*************************************'
      WRITE(*,*) 'OVERWRITING.'
      !WRITE(*,*) 'STOP.'
      !STOP
   ENDIF
   !ELSE
      !OPEN(UNIT=fid2d, FILE=trim(nameRea), STATUS='new', ACTION='write')
      OPEN(UNIT=fid2d, FILE=trim(nameRea), ACTION='write')
      CALL initializeMeshFile(fid2d, Re, dt, nsteps, iostep, debugFlag)
   !ENDIF

!==============================================================================
! write element data

   WRITE(fid3d, '(a)') '  ***** MESH DATA *****  6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8'
   WRITE(fid3d, *)     nEl,  ' 3 ', nEl,  ' NEL,NDIM,NELV'

   WRITE(fid2d, '(a)') '  ***** MESH DATA ***** 1st line is X of corner 1,2,3,4. 2nd line is Y'
   WRITE(fid2d, *)     nFpp, ' 2 ', nFpp, ' NEL,NDIM,NELV'

   DO i = 1, nEl

      CALL writeElement ( fid3d, fid2d, elem(i), i .LE. nFpp )

   ENDDO

!==============================================================================
! write curved side data

   WRITE(fid3d, '(1x,a28)') '***** CURVED SIDE DATA *****'
   WRITE(fid3d, '(2x,i10,1x,a52)') nCurvedEdges, 'Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE'

   WRITE(fid2d, '(1x,a28)') '***** CURVED SIDE DATA *****'
   WRITE(fid2d, '(2x,i10,1x,a52)') nCurvedEdges/(2*nL), 'Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE'

   DO i = 1, nEl

      CALL writeCurvedEdges ( fid3d, fid2d, elem(i), i .LE. nFpp  )

   ENDDO

!==============================================================================
! write boundary conditions

   WRITE(fid3d, '(a)') '  ***** BOUNDARY CONDITIONS *****'
   WRITE(fid3d, '(a)') '  ***** FLUID BOUNDARY CONDITIONS *****'

   WRITE(fid2d, '(a)') '  ***** BOUNDARY CONDITIONS *****'
   WRITE(fid2d, '(a)') '  ***** FLUID BOUNDARY CONDITIONS *****'

   DO i = 1, nEl

      CALL writeBoundaryConditions ( fid3d, fid2d, elem(i),  i .LE. nFpp )

   ENDDO

!==============================================================================
! finalize nameRea

   CALL finalizeMeshFile(fid3d)
   CLOSE(fid3d)

   CALL finalizeMeshFile(fid2d)
   CLOSE(fid2d)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! debug: plot with Matlab

   IF ( debugFlag ) THEN
      OPEN(UNIT=fid3d, FILE='plotElements', ACTION='write')
      DO i = 1, nEl
         DO j = 1, 8
            WRITE(fid3d, *) elem(i)%x(j), elem(i)%y(j), elem(i)%z(j), i
         ENDDO
      ENDDO
      CLOSE(fid3d)
   ENDIF

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!==============================================================================
! end of program

   DEALLOCATE(elem)

   WRITE(*,*) 'ALL DONE.'


!==============================================================================
contains

!------------------------------------------------------------------------------

   subroutine advanceFace (firstElem, lastElem, dZ, elem)

      IMPLICIT NONE
      ! input variables
      INTEGER,      INTENT(IN) :: firstElem
      INTEGER,      INTENT(IN) :: lastElem
      REAL(KIND=8), INTENT(IN) :: dZ
      ! input/output variable
      TYPE(element), DIMENSION(:) :: elem
      ! local variables
      INTEGER :: j
      INTEGER :: faceElements


      faceElements = lastElem - firstElem + 1

      DO j = lastElem+1, lastElem+faceElements

         elem(j)%num = j
         elem(j)%group = elem(j-faceElements)%group
         
         elem(j)%curvedEdge(:)  = elem(j-faceElements)%curvedEdge(:)
         elem(j)%curvedEdgeR(:) = elem(j-faceElements)%curvedEdgeR(:)
         
         elem(j)%bcType(:) = elem(j-faceElements)%bcType(:)
         elem(j)%bcParameters(:,1) = elem(j-faceElements)%bcParameters(:,1) + faceElements
         elem(j)%bcParameters(:,2) = elem(j-faceElements)%bcParameters(:,2)
         
         elem(j)%x(:) = elem(j-faceElements)%x(:)
         elem(j)%y(:) = elem(j-faceElements)%y(:)

         elem(j)%z(1:4) = elem(j-faceElements)%z(5:8)
         elem(j)%z(5:8) = elem(j-faceElements)%z(5:8) + dZ

      ENDDO

   end subroutine advanceFace

!------------------------------------------------------------------------------

   subroutine initializeMeshFile (fid, Re, dt, nsteps, iostep, debugFlag)
   
      IMPLICIT NONE
      ! input variables
      INTEGER,      INTENT(IN) :: fid, nsteps, iostep
      REAL(KIND=8), INTENT(IN) :: Re, dt
      LOGICAL,      INTENT(IN) :: debugFlag

      WRITE(fid, '(a)') ' ****** PARAMETERS *****'
      WRITE(fid, '(a)') '   2.6000      NEKTON VERSION'

      ! super dirty hack to save code-lines
      WRITE(fid, '(a,i1, a)') '   ', mod(fid,100), ' DIMENSIONAL RUN'

      WRITE(fid, '(a)') '         118  PARAMETERS FOLLOW'
      WRITE(fid, '(a)') '   1.00000     P001: DENSITY'
      if ( debugFlag ) then
         WRITE(fid, '(a)') '  -10.         P002: VISCOS'
      else
         WRITE(fid, '(2x,f8.2,a17)') -Re, '     P002: VISCOS'
      endif
      WRITE(fid, '(a)') '   0.00000     P003: : : BETAG'
      WRITE(fid, '(a)') '   0.00000     P004: : : GTHETA'
      WRITE(fid, '(a)') '   0.00000     P005: : : PGRADX'
      WRITE(fid, '(a)') '   0.00000     P006:'
      WRITE(fid, '(a)') '   0.00000     P007: RHOCP'
      WRITE(fid, '(a)') '   0.00000     P008: CONDUCT'
      WRITE(fid, '(a)') '   0.00000     P009:'
      WRITE(fid, '(a)') '   0.00000     P010: FINTIME'
      if ( debugFlag ) then
         WRITE(fid, '(a)') '   3           P011: NSTEPS'
      else
         WRITE(fid, '(3x,i7,a17)') nsteps, '     P011: NSTEPS'
      endif
      WRITE(fid, '(3x,es8.2,a13)') dt, '    P012: DT'
      WRITE(fid, '(a)') '   0.00000     P013: IOCOMM'
      WRITE(fid, '(a)') '   0.00000     P014: IOTIME'
      if ( debugFlag ) then
         WRITE(fid, '(a)') '   1           P015: IOSTEP'
      else
         WRITE(fid, '(3x,i7,a17)') iostep, '     P015: IOSTEP'
      endif
      WRITE(fid, '(a)') '   0.00000     P016: PSSOLVER: 0=default'
      WRITE(fid, '(a)') '   1.00000     P017:'
      WRITE(fid, '(a)') '   0.00000     P018: GRID < 0 --> # cells on screen'
      WRITE(fid, '(a)') '   0.00000     P019: INTYPE'
      WRITE(fid, '(a)') '   0.00000     P020: NORDER'
      WRITE(fid, '(a)') '   1.0000E-09  P021: DIVERGENCE'
      WRITE(fid, '(a)') '   1.0000E-09  P022: HELMHOLTZ'
      WRITE(fid, '(a)') '   0.00000     P023: NPSCAL'
      WRITE(fid, '(a)') '   1.0000E-02  P024: TOLREL'
      WRITE(fid, '(a)') '   1.0000E-02  P025: TOLABS'
      WRITE(fid, '(a)') '   1.00000     P026: COURANT/NTAU'
      WRITE(fid, '(a)') '   3.00000     P027: TORDER'
      WRITE(fid, '(a)') '   0.00000     P028: TORDER: mesh velocity (0: p28=p27)'
      WRITE(fid, '(a)') '   0.00000     P029: = magnetic visc if > 0, = -1/Rm if < 0'
      WRITE(fid, '(a)') '   0.00000     P030: > 0 ==> properties set in uservp()'
      WRITE(fid, '(a)') '   0.00000     P031: NPERT: #perturbation modes'
      WRITE(fid, '(a)') '   0.00000     P032: #BCs in re2 file, if > 0'
      WRITE(fid, '(a)') '   0.00000     P033: : :'
      WRITE(fid, '(a)') '   0.00000     P034: : :'
      WRITE(fid, '(a)') '   0.00000     P035: : :'
      WRITE(fid, '(a)') '   0.00000     P036: : : XMAGNET'
      WRITE(fid, '(a)') '   0.00000     P037: : : NGRIDS'
      WRITE(fid, '(a)') '   0.00000     P038: : : NORDER2'
      WRITE(fid, '(a)') '   0.00000     P039: : : NORDER3'
      WRITE(fid, '(a)') '   0.00000     P040:'
      WRITE(fid, '(a)') '   0.00000     P041: 1-->multiplicative SEMG'
      WRITE(fid, '(a)') '   0.00000     P042: 0=gmres/1=pcg'
      WRITE(fid, '(a)') '   0.00000     P043: 0=semg/1=schwarz'
      WRITE(fid, '(a)') '   0.00000     P044: 0=E-based/1=A-based prec.'
      WRITE(fid, '(a)') '   0.00000     P045: Relaxation factor for DTFS'
      WRITE(fid, '(a)') '   0.00000     P046: reserved'
      WRITE(fid, '(a)') '   0.00000     P047: vnu: mesh matieral prop.'
      WRITE(fid, '(a)') '   0.00000     P048: : :'
      WRITE(fid, '(a)') '   0.00000     P049: : :'
      WRITE(fid, '(a)') '   0.00000     P050: : :'
      WRITE(fid, '(a)') '   0.00000     P051:'
      WRITE(fid, '(a)') '   0.00000     P052: IOHIS'
      WRITE(fid, '(a)') '   0.00000     P053:'
      WRITE(fid, '(a)') '  -3.00000     P054: fixed flow rate dir: |p54|=1,2,3=x,y,z'
      WRITE(fid, '(a)') '   1.00000     P055: vol.flow rate (p54>0) or Ubar (p54<0)'
      WRITE(fid, '(a)') '   0.00000     P056: : :'
      WRITE(fid, '(a)') '   0.00000     P057: : :'
      WRITE(fid, '(a)') '   0.00000     P058:'
      WRITE(fid, '(a)') '   0.00000     P059: !=0 --> full Jac. eval. for each el.'
      WRITE(fid, '(a)') '   0.00000     P060: !=0 --> init. velocity to small nonzero'
      WRITE(fid, '(a)') '   0.00000     P061:'
      WRITE(fid, '(a)') '   0.00000     P062: >0 --> force byte_swap for output'
      WRITE(fid, '(a)') '   8.00000     P063: =8 --> force 8-byte output'
      WRITE(fid, '(a)') '   0.00000     P064: =1 --> perturbation restart'
      WRITE(fid, '(a)') '   1.00000     P065: #iofiles (eg, 0 or 64); <0 --> sep. dirs'
      WRITE(fid, '(a)') '   6.00000     P066: output : <0=ascii, else binary'
      WRITE(fid, '(a)') '   6.00000     P067: restart: <0=ascii, else binary'
      WRITE(fid, '(a)') '   20000       P068: iastep: freq for avg_all (0=iostep)'
      WRITE(fid, '(a)') '   0.00000     P069: : :     frequency of srf dump'
      WRITE(fid, '(a)') '   0.00000     P070: : :'
      WRITE(fid, '(a)') '   0.00000     P071: : :'
      WRITE(fid, '(a)') '   0.00000     P072: : :'
      WRITE(fid, '(a)') '   0.00000     P073:'
      WRITE(fid, '(a)') '   0.00000     P074: verbose Helmholtz'
      WRITE(fid, '(a)') '   0.00000     P075: : :'
      WRITE(fid, '(a)') '   0.00000     P076: : :'
      WRITE(fid, '(a)') '   0.00000     P077: : :'
      WRITE(fid, '(a)') '   0.00000     P078: : :'
      WRITE(fid, '(a)') '   0.00000     P079: : :'
      WRITE(fid, '(a)') '   0.00000     P080: : :'
      WRITE(fid, '(a)') '   0.00000     P081: : :'
      WRITE(fid, '(a)') '   0.00000     P082: : :'
      WRITE(fid, '(a)') '   0.00000     P083:'
      WRITE(fid, '(a)') '   0.00000     P084: !=0 --> sets initial timestep if p12>0'
      WRITE(fid, '(a)') '   0.00000     P085: dt ratio if p84 !=0, for timesteps>0'
      WRITE(fid, '(a)') '   0.00000     P086: reserved'
      WRITE(fid, '(a)') '   0.00000     P087: : :'
      WRITE(fid, '(a)') '   0.00000     P088: : :'
      WRITE(fid, '(a)') '   0.00000     P089: : :   coarse grid weighting (default=10.'
      WRITE(fid, '(a)') '   0.00000     P090: : :'
      WRITE(fid, '(a)') '   0.00000     P091: : :'
      WRITE(fid, '(a)') '   0.00000     P092:'
      WRITE(fid, '(a)') '   20.0000     P093: Number of previous pressure solns saved'
      WRITE(fid, '(a)') '   9.00000     P094: start projecting velocity after p94 step'
      WRITE(fid, '(a)') '   9.00000     P095: start projecting pressure after p95 step'
      WRITE(fid, '(a)') '   0.00000     P096: : :   which saving algorithm 1 = discard'
      WRITE(fid, '(a)') '   0.00000     P097: : :   0 == > no iterative refinement'
      WRITE(fid, '(a)') '   0.00000     P098:'
      WRITE(fid, '(a)') '   3.00000     P099: dealiasing: <0--> off/3--> old/4--> new'
      WRITE(fid, '(a)') '   0.00000     P100:'
      WRITE(fid, '(a)') '   0.00000     P101: Number of additional modes to filter'
      WRITE(fid, '(a)') '   0.00000     P102: Dump out divergence at each time step'
      WRITE(fid, '(a)') '   0.01000     P103: weight of stabilizing filter (.01)'
      WRITE(fid, '(a)') '   0.00000     P104: : :'
      WRITE(fid, '(a)') '   0.00000     P105: : :'
      WRITE(fid, '(a)') '   0.00000     P106:'
      WRITE(fid, '(a)') '   0.00000     P107: !=0 --> add to h2 array in hlmhotz eqn'
      WRITE(fid, '(a)') '   0.00000     P108: : :'
      WRITE(fid, '(a)') '   0.00000     P109: : :'
      WRITE(fid, '(a)') '   0.00000     P110: : :'
      WRITE(fid, '(a)') '   0.00000     P111: : :'
      WRITE(fid, '(a)') '   0.00000     P112: : :'
      WRITE(fid, '(a)') '   0.00000     P113: : :'
      WRITE(fid, '(a)') '   0.00000     P114: : :'
      WRITE(fid, '(a)') '   0.00000     P115:'
      WRITE(fid, '(a)') '   0.00000     P116: !=0: x elements for fast tensor product'
      WRITE(fid, '(a)') '   0.00000     P117: !=0: y elements for fast tensor product'
      WRITE(fid, '(a)') '   0.00000     P118: !=0: z elements for fast tensor product'
      WRITE(fid, '(a)') '      4  Lines of passive scalar data follows2 CONDUCT; 2RHOCP'
      WRITE(fid, '(a)') '   1.00000       1.00000       1.00000       1.00000       1.00000'
      WRITE(fid, '(a)') '   1.00000       1.00000       1.00000       1.00000'
      WRITE(fid, '(a)') '   1.00000       1.00000       1.00000       1.00000       1.00000'
      WRITE(fid, '(a)') '   1.00000       1.00000       1.00000       1.00000'
      WRITE(fid, '(a)') '          13   LOGICAL SWITCHES FOLLOW'
      WRITE(fid, '(a)') ' T      IFFLOW'
      WRITE(fid, '(a)') ' F      IFHEAT'
      WRITE(fid, '(a)') ' T      IFTRAN'
      WRITE(fid, '(a)') ' T F F F F F F F F F F  IFNAV & IFADVC (convection in P.S. fields)'
      WRITE(fid, '(a)') ' F F T T T T T T T T T T  IFTMSH (IF mesh for this field is T mesh)'
      WRITE(fid, '(a)') ' F      IFAXIS'
      WRITE(fid, '(a)') ' F      IFSTRS'
      WRITE(fid, '(a)') ' F      IFSPLIT'
      WRITE(fid, '(a)') ' F      IFMGRID'
      WRITE(fid, '(a)') ' F      IFMODEL'
      WRITE(fid, '(a)') ' F      IFKEPS'
      WRITE(fid, '(a)') ' F      IFMVBD'
      WRITE(fid, '(a)') ' F      IFCHAR'
      WRITE(fid, '(a)') '   2.00000       2.00000      -1.00000      -1.00000     XFAC,YFAC,XZERO,YZERO'

   end subroutine initializeMeshFile

!------------------------------------------------------------------------------

   subroutine finalizeMeshFile (fid)
   
      IMPLICIT NONE
      ! input variables
      INTEGER, INTENT(IN) :: fid


      WRITE(fid, '(a)') '  ***** NO THERMAL BOUNDARY CONDITIONS *****'
      WRITE(fid, '(a)') '   0 PRESOLVE/RESTART OPTIONS  *****'
      WRITE(fid, '(a)') '   7         INITIAL CONDITIONS *****'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') 'C Default'
      WRITE(fid, '(a)') '  ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q'
      WRITE(fid, '(a)') '   4                 Lines of Drive force data follow'
      WRITE(fid, '(a)') 'C'
      WRITE(fid, '(a)') 'C'
      WRITE(fid, '(a)') 'C'
      WRITE(fid, '(a)') 'C'
      WRITE(fid, '(a)') '  ***** Variable Property Data ***** Overrrides Parameter data.'
      WRITE(fid, '(a)') '   1 Lines follow.'
      WRITE(fid, '(a)') '   0 PACKETS OF DATA FOLLOW'
      WRITE(fid, '(a)') '  ***** HISTORY AND INTEGRAL DATA *****'
      WRITE(fid, '(a)') '   0   POINTS.  Hcode, I,J,H,IEL'
      WRITE(fid, '(a)') '  ***** OUTPUT FIELD SPECIFICATION *****'
      WRITE(fid, '(a)') '   6 SPECIFICATIONS FOLLOW'
      WRITE(fid, '(a)') '   T      COORDINATES'
      WRITE(fid, '(a)') '   T      VELOCITY'
      WRITE(fid, '(a)') '   T      PRESSURE'
      WRITE(fid, '(a)') '   T      TEMPERATURE'
      WRITE(fid, '(a)') '   F      TEMPERATURE GRADIENT'
      WRITE(fid, '(a)') '   0      PASSIVE SCALARS'
      WRITE(fid, '(a)') '  ***** OBJECT SPECIFICATION *****'
      WRITE(fid, '(a)') '       0 Surface Objects'
      WRITE(fid, '(a)') '       0 Volume  Objects'
      WRITE(fid, '(a)') '       0 Edge    Objects'
      WRITE(fid, '(a)') '       0 Point   Objects'

   end subroutine finalizeMeshFile

!------------------------------------------------------------------------------

   subroutine writeElement (fid3d, fid2d, elem, isFirstFace)

      IMPLICIT NONE
      ! input variables
      INTEGER,       INTENT(IN) :: fid3d
      INTEGER,       INTENT(IN) :: fid2d
      TYPE(element), INTENT(IN) :: elem
      LOGICAL,       INTENT(IN) :: isFirstFace
      ! local variables

      WRITE(fid3d, '(a18,1x,i10,a4,i3,a1,a11,i5)') &
         '           ELEMENT', elem%num, ' [  ', 1, elem%groupL, ']    GROUP ', elem%group

      WRITE(fid3d, '(4(es14.6e2))') elem%x(1), elem%x(2), elem%x(3), elem%x(4)
      WRITE(fid3d, '(4(es14.6e2))') elem%y(1), elem%y(2), elem%y(3), elem%y(4)
      WRITE(fid3d, '(4(es14.6e2))') elem%z(1), elem%z(2), elem%z(3), elem%z(4)

      WRITE(fid3d, '(4(es14.6e2))') elem%x(5), elem%x(6), elem%x(7), elem%x(8)
      WRITE(fid3d, '(4(es14.6e2))') elem%y(5), elem%y(6), elem%y(7), elem%y(8)
      WRITE(fid3d, '(4(es14.6e2))') elem%z(5), elem%z(6), elem%z(7), elem%z(8)

      IF (isFirstFace) THEN

        WRITE(fid2d, '(a18,1x,i10,a4,i3,a1,a11,i5)') &
          '          ELEMENT', elem%num, ' [  ', 1, elem%groupL, ']    GROUP ', elem%group

        WRITE(fid2d, '(4(es14.6e2))') elem%x(1), elem%x(2), elem%x(3), elem%x(4)
        WRITE(fid2d, '(4(es14.6e2))') elem%y(1), elem%y(2), elem%y(3), elem%y(4)
      ENDIF

   end subroutine writeElement

!------------------------------------------------------------------------------

   subroutine writeCurvedEdges (fid3d, fid2d, elem, isFirstFace)

      IMPLICIT NONE
      ! input variables
      INTEGER,       INTENT(IN) :: fid3d
      INTEGER,       INTENT(IN) :: fid2d
      TYPE(element), INTENT(IN) :: elem
      LOGICAL,       INTENT(IN) :: isFirstFace
      ! local variables
      INTEGER :: j

      IF ( nEl < 1e3 ) THEN

         DO j = 1, 8 ! cycle on the eight edges
            IF ( elem%curvedEdge(j) ) THEN
               WRITE(fid3d, '(i3,i3,f10.5,4(f14.5),5x,a1)') &
                  j, elem%num, elem%curvedEdgeR(j), 0d0, 0d0, 0d0, 0d0, 'C'
            ENDIF
         ENDDO
      ELSEIF ( nEl < 1e6 ) THEN

         DO j = 1, 8 ! cycle on the eight edges
            IF ( elem%curvedEdge(j) ) THEN
               WRITE(fid3d, '(i2,i6,f10.5,4(f14.5),5x,a1)') &
                  j, elem%num, elem%curvedEdgeR(j), 0d0, 0d0, 0d0, 0d0, 'C'
            ENDIF
         ENDDO
      ELSE

         DO j = 1, 8 ! cycle on the eight edges
            IF ( elem%curvedEdge(j) ) THEN
               WRITE(fid3d, '(i2,i10,f10.5,4(f14.5),5x,a1)') &
                  j, elem%num, elem%curvedEdgeR(j), 0d0, 0d0, 0d0, 0d0, 'C'
            ENDIF
         ENDDO
      ENDIF

      IF (isFirstFace) THEN
        IF ( nFpp < 1e3 ) THEN

          DO j = 1, 4 ! cycle on the first four edges
          IF ( elem%curvedEdge(j) ) THEN
            WRITE(fid2d, '(i3,i3,f10.5,4(f14.5),5x,a1)') &
              j, elem%num, elem%curvedEdgeR(j), 0d0, 0d0, 0d0, 0d0, 'C'
          ENDIF
          ENDDO
        ELSEIF ( nFpp < 1e6 ) THEN

          DO j = 1, 4 ! cycle on the first four edges
          IF ( elem%curvedEdge(j) ) THEN
            WRITE(fid2d, '(i2,i6,f10.5,4(f14.5),5x,a1)') &
              j, elem%num, elem%curvedEdgeR(j), 0d0, 0d0, 0d0, 0d0, 'C'
          ENDIF
          ENDDO
        ELSE

          DO j = 1, 4 ! cycle on the first four edges
          IF ( elem%curvedEdge(j) ) THEN
            WRITE(fid2d, '(i2,i10,f10.5,4(f14.5),5x,a1)') &
              j, elem%num, elem%curvedEdgeR(j), 0d0, 0d0, 0d0, 0d0, 'C'
          ENDIF
          ENDDO
        ENDIF
      ENDIF

   end subroutine writeCurvedEdges

!------------------------------------------------------------------------------

   subroutine writeBoundaryConditions (fid3d, fid2d, elem, isFirstFace)

      IMPLICIT NONE
      ! input variables
      INTEGER,       INTENT(IN) :: fid3d
      INTEGER,       INTENT(IN) :: fid2d
      TYPE(element), INTENT(IN) :: elem
      LOGICAL,       INTENT(IN) :: isFirstFace
      ! local variables
      INTEGER :: j

      IF ( nEl < 1e3 ) THEN

         DO j = 1, 6 ! cycle on the six faces
            WRITE(fid3d, '(1x,a2,1x,i3,i3,5(g14.6))') &
               elem%bcType(j), elem%num, j, elem%bcParameters(j,:)
         ENDDO

      ELSEIF ( nEl < 1e5 ) THEN

         DO j = 1, 6 ! cycle on the six faces
            WRITE(fid3d, '(1x,a2,1x,i5,i1,5(g14.6))') &
               elem%bcType(j), elem%num, j, elem%bcParameters(j,:)
         ENDDO

      ELSEIF ( nEl < 1e6 ) THEN

         DO j = 1, 6 ! cycle on the six faces
            WRITE(fid3d, '(1x,a2,1x,i6,5(g14.6))') &
               elem%bcType(j), elem%num, elem%bcParameters(j,:)
         ENDDO
      ELSE

         DO j = 1, 6 ! cycle on the six faces
            WRITE(fid3d, '(1x,a2,1x,i11,i1,5(g18.11))') &
               elem%bcType(j), elem%num, j, elem%bcParameters(j,:)
         ENDDO
      ENDIF

      IF (isFirstFace) THEN
        IF ( nFpp < 1e3 ) THEN
          DO j = 1, 4 ! cycle on the four faces
            WRITE(fid2d, '(1x,a2,1x,i3,i3,5(g14.6))') &
              elem%bcType(j), elem%num, j, elem%bcParameters(j,:)
          ENDDO

        ELSEIF ( nFpp < 1e5 ) THEN
          DO j = 1, 4 ! cycle on the four faces
            WRITE(fid2d, '(1x,a2,1x,i5,i1,5(g14.6))') &
              elem%bcType(j), elem%num, j, elem%bcParameters(j,:)
          ENDDO

        ELSEIF ( nFpp < 1e6 ) THEN
          DO j = 1, 4 ! cycle on the four faces
            WRITE(fid2d, '(1x,a2,1x,i6,5(g14.6))') &
              elem%bcType(j), elem%num, elem%bcParameters(j,:)
          ENDDO

        ELSE

          DO j = 1, 4 ! cycle on the four faces
            WRITE(fid2d, '(1x,a2,1x,i11,i1,5(g18.11))') &
              elem%bcType(j), elem%num, j, elem%bcParameters(j,:)
          ENDDO

        ENDIF
      ENDIF

   end subroutine writeBoundaryConditions

!------------------------------------------------------------------------------


  !*****************************************************************************80
  !
  !! LOBATTO_SET sets abscissas Lobatto quadrature.
  !  taken from
  !   https://people.sc.fsu.edu/~jburkardt/f_src/quadrule/quadrule.html
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the polynomial order.
  !    N must be between 2 and 20.
  !
  !    Output, real ( kind = 8 ) X(N), the abscissas.
  !
  subroutine lobatto_set ( n, x)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n

    REAL ( KIND = 8 ) x(n+1)

    if ( n == 1 ) then
      x(1) =  - 1.0D+00
      x(2) =    1.0D+00
    else if ( n == 2 ) then

      x(1) =  - 1.0D+00
      x(2) =    0.0D+00
      x(3) =    1.0D+00
    else if ( n == 3 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.447213595499957939281834733746D+00
      x(3) =    0.447213595499957939281834733746D+00
      x(4) =    1.0D+00
    else if ( n == 4 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.654653670707977143798292456247D+00
      x(3) =    0.0D+00
      x(4) =    0.654653670707977143798292456247D+00
      x(5) =    1.0D+00
    else if ( n == 5 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.765055323929464692851002973959D+00
      x(3) =  - 0.285231516480645096314150994041D+00
      x(4) =    0.285231516480645096314150994041D+00
      x(5) =    0.765055323929464692851002973959D+00
      x(6) =    1.0D+00
    else if ( n == 6 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.830223896278566929872032213967D+00
      x(3) =  - 0.468848793470714213803771881909D+00
      x(4) =    0.0D+00
      x(5) =    0.468848793470714213803771881909D+00
      x(6) =    0.830223896278566929872032213967D+00
      x(7) =    1.0D+00
    else if ( n == 7 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.871740148509606615337445761221D+00
      x(3) =  - 0.591700181433142302144510731398D+00
      x(4) =  - 0.209299217902478868768657260345D+00
      x(5) =    0.209299217902478868768657260345D+00
      x(6) =    0.591700181433142302144510731398D+00
      x(7) =    0.871740148509606615337445761221D+00
      x(8) =    1.0D+00
    else if ( n == 8 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.899757995411460157312345244418D+00
      x(3) =  - 0.677186279510737753445885427091D+00
      x(4) =  - 0.363117463826178158710752068709D+00
      x(5) =    0.0D+00
      x(6) =    0.363117463826178158710752068709D+00
      x(7) =    0.677186279510737753445885427091D+00
      x(8) =    0.899757995411460157312345244418D+00
      x(9) =    1.0D+00
    else if ( n == 9 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.919533908166458813828932660822D+00
      x(3) =  - 0.738773865105505075003106174860D+00
      x(4) =  - 0.477924949810444495661175092731D+00
      x(5) =  - 0.165278957666387024626219765958D+00
      x(6) =    0.165278957666387024626219765958D+00
      x(7) =    0.477924949810444495661175092731D+00
      x(8) =    0.738773865105505075003106174860D+00
      x(9) =    0.919533908166458813828932660822D+00
      x(10) =   1.0D+00
    else if ( n == 10 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.934001430408059134332274136099D+00
      x(3) =  - 0.784483473663144418622417816108D+00
      x(4) =  - 0.565235326996205006470963969478D+00
      x(5) =  - 0.295758135586939391431911515559D+00
      x(6) =    0.0D+00
      x(7) =    0.295758135586939391431911515559D+00
      x(8) =    0.565235326996205006470963969478D+00
      x(9) =    0.784483473663144418622417816108D+00
      x(10) =   0.934001430408059134332274136099D+00
      x(11) =   1.0D+00
    else if ( n == 11 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.944899272222882223407580138303D+00
      x(3) =  - 0.819279321644006678348641581717D+00
      x(4) =  - 0.632876153031869677662404854444D+00
      x(5) =  - 0.399530940965348932264349791567D+00
      x(6) =  - 0.136552932854927554864061855740D+00
      x(7) =    0.136552932854927554864061855740D+00
      x(8) =    0.399530940965348932264349791567D+00
      x(9) =    0.632876153031869677662404854444D+00
      x(10) =   0.819279321644006678348641581717D+00
      x(11) =   0.944899272222882223407580138303D+00
      x(12) =   1.0D+00
    else if ( n == 12 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.953309846642163911896905464755D+00
      x(3) =  - 0.846347564651872316865925607099D+00
      x(4) =  - 0.686188469081757426072759039566D+00
      x(5) =  - 0.482909821091336201746937233637D+00
      x(6) =  - 0.249286930106239992568673700374D+00
      x(7) =    0.0D+00
      x(8) =    0.249286930106239992568673700374D+00
      x(9) =    0.482909821091336201746937233637D+00
      x(10) =   0.686188469081757426072759039566D+00
      x(11) =   0.846347564651872316865925607099D+00
      x(12) =   0.953309846642163911896905464755D+00
      x(13) =   1.0D+00
    else if ( n == 13 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.959935045267260901355100162015D+00
      x(3) =  - 0.867801053830347251000220202908D+00
      x(4) =  - 0.728868599091326140584672400521D+00
      x(5) =  - 0.550639402928647055316622705859D+00
      x(6) =  - 0.342724013342712845043903403642D+00
      x(7) =  - 0.116331868883703867658776709736D+00
      x(8) =    0.116331868883703867658776709736D+00
      x(9) =    0.342724013342712845043903403642D+00
      x(10) =   0.550639402928647055316622705859D+00
      x(11) =   0.728868599091326140584672400521D+00
      x(12) =   0.867801053830347251000220202908D+00
      x(13) =   0.959935045267260901355100162015D+00
      x(14) =   1.0D+00
    else if ( n == 14 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.965245926503838572795851392070D+00
      x(3) =  - 0.885082044222976298825401631482D+00
      x(4) =  - 0.763519689951815200704118475976D+00
      x(5) =  - 0.606253205469845711123529938637D+00
      x(6) =  - 0.420638054713672480921896938739D+00
      x(7) =  - 0.215353955363794238225679446273D+00
      x(8) =    0.0D+00
      x(9) =    0.215353955363794238225679446273D+00
      x(10) =   0.420638054713672480921896938739D+00
      x(11) =   0.606253205469845711123529938637D+00
      x(12) =   0.763519689951815200704118475976D+00
      x(13) =   0.885082044222976298825401631482D+00
      x(14) =   0.965245926503838572795851392070D+00
      x(15) =   1.0D+00
    else if ( n == 15 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.969568046270217932952242738367D+00
      x(3) =  - 0.899200533093472092994628261520D+00
      x(4) =  - 0.792008291861815063931088270963D+00
      x(5) =  - 0.652388702882493089467883219641D+00
      x(6) =  - 0.486059421887137611781890785847D+00
      x(7) =  - 0.299830468900763208098353454722D+00
      x(8) =  - 0.101326273521949447843033005046D+00
      x(9) =    0.101326273521949447843033005046D+00
      x(10) =   0.299830468900763208098353454722D+00
      x(11) =   0.486059421887137611781890785847D+00
      x(12) =   0.652388702882493089467883219641D+00
      x(13) =   0.792008291861815063931088270963D+00
      x(14) =   0.899200533093472092994628261520D+00
      x(15) =   0.969568046270217932952242738367D+00
      x(16) =   1.0D+00
    else if ( n == 16 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.973132176631418314156979501874D+00
      x(3) =  - 0.910879995915573595623802506398D+00
      x(4) =  - 0.815696251221770307106750553238D+00
      x(5) =  - 0.691028980627684705394919357372D+00
      x(6) =  - 0.541385399330101539123733407504D+00
      x(7) =  - 0.372174433565477041907234680735D+00
      x(8) =  - 0.189511973518317388304263014753D+00
      x(9) =    0.0D+00
      x(10) =   0.189511973518317388304263014753D+00
      x(11) =   0.372174433565477041907234680735D+00
      x(12) =   0.541385399330101539123733407504D+00
      x(13) =   0.691028980627684705394919357372D+00
      x(14) =   0.815696251221770307106750553238D+00
      x(15) =   0.910879995915573595623802506398D+00
      x(16) =   0.973132176631418314156979501874D+00
      x(17) =   1.0D+00
    else if ( n == 17 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.976105557412198542864518924342D+00
      x(3) =  - 0.920649185347533873837854625431D+00
      x(4) =  - 0.835593535218090213713646362328D+00
      x(5) =  - 0.723679329283242681306210365302D+00
      x(6) =  - 0.588504834318661761173535893194D+00
      x(7) =  - 0.434415036912123975342287136741D+00
      x(8) =  - 0.266362652878280984167665332026D+00
      x(9) =  - 0.897490934846521110226450100886D-01
      x(10) =   0.897490934846521110226450100886D-01
      x(11) =   0.266362652878280984167665332026D+00
      x(12) =   0.434415036912123975342287136741D+00
      x(13) =   0.588504834318661761173535893194D+00
      x(14) =   0.723679329283242681306210365302D+00
      x(15) =   0.835593535218090213713646362328D+00
      x(16) =   0.920649185347533873837854625431D+00
      x(17) =   0.976105557412198542864518924342D+00
      x(18) =   1.0D+00
    else if ( n == 18 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.978611766222080095152634063110D+00
      x(3) =  - 0.928901528152586243717940258797D+00
      x(4) =  - 0.852460577796646093085955970041D+00
      x(5) =  - 0.751494202552613014163637489634D+00
      x(6) =  - 0.628908137265220497766832306229D+00
      x(7) =  - 0.488229285680713502777909637625D+00
      x(8) =  - 0.333504847824498610298500103845D+00
      x(9) =  - 0.169186023409281571375154153445D+00
      x(10) =   0.0D+00
      x(11) =   0.169186023409281571375154153445D+00
      x(12) =   0.333504847824498610298500103845D+00
      x(13) =   0.488229285680713502777909637625D+00
      x(14) =   0.628908137265220497766832306229D+00
      x(15) =   0.751494202552613014163637489634D+00
      x(16) =   0.852460577796646093085955970041D+00
      x(17) =   0.928901528152586243717940258797D+00
      x(18) =   0.978611766222080095152634063110D+00
      x(19) =   1.0D+00
    else if ( n == 19 ) then

      x(1) =  - 1.0D+00
      x(2) =  - 0.980743704893914171925446438584D+00
      x(3) =  - 0.935934498812665435716181584931D+00
      x(4) =  - 0.866877978089950141309847214616D+00
      x(5) =  - 0.775368260952055870414317527595D+00
      x(6) =  - 0.663776402290311289846403322971D+00
      x(7) =  - 0.534992864031886261648135961829D+00
      x(8) =  - 0.392353183713909299386474703816D+00
      x(9) =  - 0.239551705922986495182401356927D+00
      x(10) = - 0.805459372388218379759445181596D-01
      x(11) =   0.805459372388218379759445181596D-01
      x(12) =   0.239551705922986495182401356927D+00
      x(13) =   0.392353183713909299386474703816D+00
      x(14) =   0.534992864031886261648135961829D+00
      x(15) =   0.663776402290311289846403322971D+00
      x(16) =   0.775368260952055870414317527595D+00
      x(17) =   0.866877978089950141309847214616D+00
      x(18) =   0.935934498812665435716181584931D+00
      x(19) =   0.980743704893914171925446438584D+00
      x(20) =   1.0D+00
    else

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LOBATTO_SET - Fatal error!'
      write ( *, '(a,i8)' ) '  Illegal value of polyonmial degree N = ', n
      write ( *, '(a)' ) '  Legal values are between 2 and 19.'
      stop 1

    end if

    return
  end subroutine lobatto_set 

!==============================================================================
end program pipeMeshNek
