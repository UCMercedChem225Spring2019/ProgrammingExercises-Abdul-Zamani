      program pgrm_03_03
!
!     This program computes components of the Hartree-Fock energy and the number
!     of electrons using contraction of various matrices loaded from
!     user-provided files (that presumably come from Gaussian calculations).
!
!     At run time, the program expects 5 command line arguments:
!       (1) the number of electrons;
!       (2) the number of atomic orbital basis functions;
!       (3) an input file containing core-Hamiltonian matrix elements (in
!           symmetric upper/column storage form);
!       (4) an input file containing the overlap matrix elements (in
!           symmetric upper/column storage form); and
!       (5) an input file containing Fock matrix elements (in symmetric
!           upper/column storage form).
!
!     At run time, the program outputs 6 items:
!       (1) the MO energies;
!       (2) the MO coefficients;
!       (3) the density matrix;
!       (4) the one-electron contribution to the total electronic energy;
!       (5) the two-electron contribution to the total electronic energy; and
!       (6) the tr(PS), which should equal the number of electrons.
!
!
!     A.Zamani 4/18/19 | CHEM 225 S19
!
!
      implicit none
      integer,parameter::IIn=10
      integer::i,iError,nElectrons,nOcc,nBasis,lenSym, j, k
      real::oneElectronEnergy,twoElectronEnergy,tracePS
      real,dimension(:),allocatable::symFock,symCoreHamiltonian, &
        symOverlap, moEnergies,tempSymMatrix
      real,dimension(:,:),allocatable::sqFock,fockTilde, &
        sqCoreHamiltonian,InvSqrtOverlap,moCoefficients,densityMatrix, &
        tempSqMatrix,  SSPEV_Scratch
      character(len=256)::cmdlineArg
      !
      !
      !Begin by reading the number of basis functions, allocating array memory,
      !and loading the symmetric Fock and overlap matrices from input files
      !provided on the command line.
      !
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nElectrons
      nOcc = nElectrons/2
      !
      call Get_Command_Argument(2,cmdlineArg)
      read(cmdlineArg,'(I)') nBasis
      lenSym = (nBasis*(nBasis+1))/2
      allocate(symFock(lenSym),symCoreHamiltonian(lenSym), &
        symOverlap(lenSym),tempSymMatrix(lenSym))
      allocate(moEnergies(nBasis))
      allocate(sqFock(nBasis,nBasis),fockTilde(nBasis,nBasis), &
        sqCoreHamiltonian(nBasis,nBasis), &
        invSqrtOverlap(nBasis,nBasis), &
        moCoefficients(nBasis,nBasis),densityMatrix(nBasis,nBasis), &
        tempSqMatrix(nBasis,nBasis),SSPEV_Scratch(nBasis,3))
      !
      !
      !.h
      Call Get_Command_Argument(3,cmdlineArg)
      Open(Unit=IIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.'
        STOP
      endIf
      DO i = 1, lenSym !Read triangle elements
         read(IIn,*) symCoreHamiltonian(i)
      END DO
      Close(Unit=IIn)
      !print*,1 !Uncomment to check what went wrong
      !
      !
      !.s
      Call Get_Command_Argument(4,cmdlineArg)
      Open(Unit=IIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.'
        STOP
      endIf
      DO i = 1, lenSym !Read triangle elements
         read(IIn,*) symOverlap(i)
      END DO
      Close(Unit=IIn)
      !print*,2 !Uncomment to check what went wrong
      !
      !
      !.f
      Call Get_Command_Argument(5,cmdlineArg)
      Open(Unit=IIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.'
        STOP
      endIf
      DO i = 1, lenSym !Read triangle elements
         read(IIn,*) symFock(i)
      END DO
      Close(Unit=IIn)
      !print*,3 !Uncomment to check what went wrong

      call SymmetricPacked2Matrix_UpperPacked(nBasis,symOverlap, &
        invSqrtOverlap)
      !OPERATIONS BEGIN     
      !
      !Form the square-root of the overlap matrix.
      !
      !WARNING!!! Change in the overlap matrix
      !Use temp array to switch back, then set 
      !symOverlap equal to tempSymMatrix.
      tempSymMatrix = symOverlap
      call InvSQRT_SymMatrix(nBasis,symOverlap,invSqrtOverlap)
      symOverlap = tempSymMatrix
      !
      !Form fTilde and solve for the MO energies and coefficients.
      !
      call SymmetricPacked2Matrix_UpperPacked(nBasis,symFock,sqFock)
      tempSqMatrix = MatMul(invSqrtOverlap,sqFock)
      fockTilde = MatMul(tempSqMatrix,invSqrtOverlap)
      call Sq2SymMatrix(nBasis,fockTilde,tempSymMatrix)
      call SSPEV('V','U',nBasis,tempSymMatrix,moEnergies,  &
        moCoefficients,nBasis,SSPEV_Scratch,iError)
      If(iError.ne.0) then
        write(*,*)' Failure in SSPEV.'
        STOP
      endIf
      write(*,*)' MO Energies:'
      call Print_Matrix_Full_Real(Reshape(moEnergies,[nBasis,1]), &
        nBasis,1)
      tempSqMatrix = moCoefficients
      moCoefficients = MatMul(invSqrtOverlap,tempSqMatrix)
      write(*,*)' MO Coefficients:'
      call Print_Matrix_Full_Real(moCoefficients,nBasis,nBasis)
      !
      !Form Density Matrix: P = c*cT...WHY DOESNT 2*(c*cT) WORK!!!?
      !
      !densityMatrix = MatMul(moCoefficients,transpose(moCoefficients))
!!!      
      densityMatrix = 0
      do i=1,nBasis
        do j=1,nBasis
          do k=1,nOcc
            densityMatrix(i,j) = densityMatrix(i,j) + &
             2*(moCoefficients(i,k)*moCoefficients(j,k))
          endDo
        endDo
      endDo
!!!
      write(*,*)' Density Matrix:' 
      call Print_Matrix_Full_Real(densityMatrix, nBasis, nBasis)
      ! 
      !Do tr(PH) and unpack symCoreHamiltonian
      !
      call SymmetricPacked2Matrix_UpperPacked(nBasis, &
           symCoreHamiltonian,sqCoreHamiltonian) 
      tempSqMatrix = MatMul(densityMatrix,sqCoreHamiltonian)
    
      oneElectronEnergy = 0  !Trace 
      do i=1, nBasis
        oneElectronEnergy = oneElectronEnergy + tempSqMatrix(i,i)
      enddo   
      write(*,*)' One electron energy contribution: ',oneElectronEnergy

      !Do 1/2tr(PG) and form G = F -H
      !Let G be the tempSqMatrix for the following:
      tempSqMatrix = sqFock - sqCoreHamiltonian  
      !let trace PG be the tempSqMatrix for the following:
      tempSqMatrix = (MatMul(densityMatrix,tempSqMatrix))/2
      !Divide trace PG by 2      
      twoElectronEnergy = 0 !Trace
      do i=1, nBasis
        twoElectronEnergy = twoElectronEnergy + tempSqMatrix(i,i)
      enddo
      write(*,*)' Two electron energy contribution: ', &
        twoElectronEnergy 
      !
      !Do tr(PS) and now unpack S overlap
      !
      call SymmetricPacked2Matrix_UpperPacked(nBasis, &
           symOverlap, tempSqMatrix)!read sqOverlap into temp
      tempSqMatrix = MatMul(densityMatrix,tempSqMatrix)
      !test below
      call Print_Matrix_Full_Real(tempSqMatrix, nBasis, nBasis)
      !PS is P(tempSqMatrix)
      oneElectronEnergy = 0  !Trace
      do i=1, nBasis !sum over nElec for tr(PS)
        !apparently summing over nBasis works...should it?
        tracePS = tracePS + tempSqMatrix(i,i)
      enddo
      write(*,*)' The number of electrons is: ',tracePS
      end program pgrm_03_03
      !
      ! 
      !
      Subroutine SymmetricPacked2Matrix_UpperPacked(N,ArrayIn,AMatOut)
      !
      !This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2
      !long.
      !It then converts that form to the N-by-N matrix AMatOut taking
      !ArrayIn to be in upper-packed storage form. Note: The storage mode
      !also assumes the upper-packed storage is packed by columns.
      !
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
      !
      Integer::i,j,k
      !
      !Loop through the elements of AMatOut and fill them appropriately
      !from Array_Input.
      !
      k = 1
      Do i = 1, N
       Do j = 1, i ! switch i and 1 to start from tail of upper
                   ! triangle
          AMatOut(i,j) = ArrayIn(k)
        If (j.ne.i) then
          AMatOUT(j,i) = ArrayIn(k)
        EndIf
          k = k + 1 ! iterates from ArrayIN to include each element
       EndDo
      EndDo
      Return
      End Subroutine SymmetricPacked2Matrix_UpperPacked
      !
      !This subroutine reads a symmetric matrix stored in a user-provided
      !input file and forms its inverse square-root. The product of the
      !inverse square-root with itself is then verified to be the
      !identity.
      !
      subroutine InvSQRT_SymMatrix(NDim, InputSymMatrix, &
        InvSqrtInputMatrix)

      implicit none
      integer:: i,NDim, Ierror
      real, dimension(NDim) :: EVal
      real,dimension(NDim) :: inputSymMatrix
      real,dimension(3*NDim) :: TempVec
      real,dimension(NDim,NDim) ::inputSqMatrix, &
        invSqrtInputMatrix, EVec, EValMatrix

      call sspev('V','U', NDim, inputSymMatrix, EVal, &
        EVec, NDim, TempVec, Ierror)
       If(Ierror.ne.0) then
         Write(*,*)' Failure in SSPEV.'
         STOP
       endIf

      Do i = 1, NDim
        EValMatrix(i,i) = 1/(sqrt(EVal(i)))
      EndDo

      invSqrtInputMatrix = MatMul(MatMul(EVec, EValMatrix) &
        , transpose(EVec))

      end subroutine InvSQRT_SymMatrix
        !
        !The following subroutine symmetrizes and linearizes to an upper/column
        !storage form of a symmetric square matrix.
        !
      subroutine Sq2SymMatrix(NDim,SqMatrix,SymMatrix)

      integer::i,j,k,NDim
      real,dimension(NDim) :: SymMatrix
      real,dimension(NDim,NDim) :: SqMatrix

       k = 0
       do j=1,NDim
         do i=1,j
             k = k+1
             SymMatrix(k)=SqMatrix(i,j)
         enddo
       enddo

      end subroutine Sq2SymMatrix
      !
      !

      Subroutine Print_Matrix_Full_Real(AMat,M,N)
      !
      !This subroutine prints a real matrix that is fully dimension -
      !i.e.,
      !not stored in packed form. AMat is the matrix, which is
      !dimensioned
      !(M,N).
      !
      !The output of this routine is sent to unit number 6 (set by the
      !local
      !parameter integer IOut).
      !
      !
      !Variable Declarations
      !
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
      !
      !Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
      !
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
      !
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
      !
      Return
      End Subroutine Print_Matrix_Full_Real

