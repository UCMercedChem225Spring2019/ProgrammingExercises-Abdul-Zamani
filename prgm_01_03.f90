      Program prgm_01_03
!
!     This program reads two 3x3 matrices from a user-defined input files. After the
!     files are opened and read, they are closed and then printed. The
!     matrix product of these two matrices is formed using the f90
!     intrinsic fxn MatMul. The final matrix is also printed.
!
!     CHEM 225 Spring 2019
!     A. Zamani, 1/25/2019.
!
      implicit none
      integer,parameter::inFileUnitA=10, inFileUnitB=11
      integer::errorFlag,i
      real,dimension(3,3)::matrixInA, matrixInB, matrixProduct
      character(len=128)::fileNameA, filenameB
!
!     Start by asking the user for the name of the data file.
!
      write(*,*)' What is the name of the first input data file?'
      read(*,*) fileNameA
      write(*,*)' What is the name of the second input data file?'
      read(*,*) filenameB
      write(*,*)
!
!     Open the data file and read matrixInA from that file.
!
      open(unit=inFileUnitA,file=TRIM(fileNameA),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitA,*) matrixInA(1,i),matrixInA(2,i),matrixInA(3,i)
      endDo
      close(inFileUnitA)
!
!     Open the data file and read matrixInB from that file.
!
      open(unit=inFileUnitB,file=TRIM(fileNameB),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitB,*) matrixInB(1,i),matrixInB(2,i),matrixInB(3,i)
      endDo
      close(inFileUnitB)
!
!     Call the subroutine PrintMatrix to print matrixInA & matrixInB.
!
      call PrintMatrix3x3(matrixInA)
      call PrintMatrix3x3(matrixInB)
!
!     Define matrix product using the intrinsic fxn
!     Call the subroutine PrintMatrix to print matrixProduct
!
      matrixProduct = MatMul(matrixInA,matrixInB)
      call PrintMatrix3x3(matrixProduct)

!
  999 continue
      End Program prgm_01_01


      Subroutine PrintMatrix3x3(matrix)
!
!     This subroutine prints a 3x3 real matrix. The output is written to StdOut.
!
      implicit none
      real,dimension(3,3),intent(in)::matrix ! Read into this
      real,dimension(3,3) :: StdOut !Print this
      integer::i
!
!     Format statements.
!
 1000 format(3(2x,f5.1))
!
!     Do the printing job for both matrices.
!
      write(*,*)' Printing Matrix... '
        do i=1,3                  !Iterate for all element positions in a 3x3
          StdOut(3,i)= matrix(3,i)
          StdOut(2,i)= matrix(2,i)!Transpose below or reverse elements
          StdOut(1,i)= matrix(1,i)
        enddo      
      write(*,*)
      write(*,1000)transpose(StdOut) !call format statement & transpose
      write(*,*)  
      return
      End Subroutine PrintMatrix3x3
