! This is a template code for converting some identified sections
! of data in a text file and keep the rest part unchanged for
! general purposes. In other words, an output file will be generated
! with some data converted and all the rest simply copied from
! the input file. Specifically, this code converts Cartesian
! coordinates into fractionals coordinates of particles in a cell
! with respect to the three period vectors of a periodic-structured
! system, as an example. Although not necessarily, it is assumed that
! each section to be converted contains all data needed for the
! converting operation, which means in this example all the period
! vectors are also available in each section. The identifier of each
! section should be in the first line of the section, but no data to
! be converted in that line. Sections of the same structure/format
! and the same way for conversion can be identified the same, then
! converted by the same routine(s) of the code, supplied by the user.
! Those sections are regarded as in the same group. Different sections
! of the same group may have different amount of data to be converted.
! If the data section does not has a unique identifier in the original
! input file, a line with such can be inserted into the input file,
! then the user has a choice to copy such a line into the output file
! or not.
!
!
! To use this template code, in most cases, users should edit
!     1) the KEY_WORD_GROUPS module;
!     2) SUBROUTINE DIRECTION()
!     3) SUBROUTINEs for specific conversion of each group of
!                    data sections.
!
! The KEY_WORD_GROUPS module is used to specify total number
! of groups of data sections to be converted (3 here);
!    the key words (identifier) for each group;
!    whether the key word line will be copied into the output file;
!    the number of lines after the key word line should be copied
!    into the output file directly, which do not contain any data
!    to be converted and/or needed.
!
! Following the example routines at the end
!    SUBROUTINE NICE_MOLECULE_PROCESSING()
!    SUBROUTINE XYZ_FORMAT
!    SUBROUTINE GEOMETRY()
! users should code their own routines of such to perform their
! specific and detailed conversion. For that purposes, many tools
! can be used from the FOR_READ_AND_WRITE module. Once they are
! coded, they should be registered in the routione
!    SUBROUTINE DIRECTION()
! based on the groups they are dealing with.
!
!
! Any comments please send to gang.liu@queensu.ca
!
! Copyright (c) CAC (HPCVL), Queen's University, Mar. 2017




DOUBLE PRECISION FUNCTION VECTOR_DOT(A, B)
    IMPLICIT NONE
    DOUBLE PRECISION :: A(3), B(3)
    VECTOR_DOT = A(1)*B(1) + A(2)*B(2) + A(3)*B(3)
    RETURN
END FUNCTION VECTOR_DOT




SUBROUTINE VECTOR_CROSS(A, B, C)
    IMPLICIT NONE
    DOUBLE PRECISION :: A(3), B(3), C(3)
    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)
    RETURN
END SUBROUTINE VECTOR_CROSS




MODULE VECTOR_ROUTINE_INTERFACE
    INTERFACE
         DOUBLE PRECISION FUNCTION VECTOR_DOT(A, B)
             IMPLICIT NONE
             DOUBLE PRECISION :: A(3), B(3)
         END FUNCTION VECTOR_DOT

         SUBROUTINE VECTOR_CROSS(A, B, C)
             IMPLICIT NONE
             DOUBLE PRECISION :: A(3), B(3), C(3)
         END SUBROUTINE VECTOR_CROSS
    END INTERFACE
END MODULE VECTOR_ROUTINE_INTERFACE




SUBROUTINE CARTESIAN_TO_FRACTIONAL(A, B, C, R, F, N)
    ! TO CALCULATE FRACTIONAL COORDINATES (F ARRAY)
    ! FROM CARTESIAN ONES (R ARRAY) .
    ! A, B, AND C ARE THE THREE KNOWN INDEPENDENT
    ! PERIOD VECTORS .
    USE VECTOR_ROUTINE_INTERFACE
    IMPLICIT NONE
    INTEGER          :: N, I
    DOUBLE PRECISION :: A(3), B(3), C(3), R(3,N), F(3,N)
    DOUBLE PRECISION :: ZERO = 0.1D-30,   SGM(3,3), OMG
    CALL  VECTOR_CROSS(B, C, SGM(1:3, 1))
    CALL  VECTOR_CROSS(C, A, SGM(1:3, 2))
    CALL  VECTOR_CROSS(A, B, SGM(1:3, 3))
    OMG = VECTOR_DOT  (A,    SGM(1:3, 1))
    IF(OMG.LT.ZERO) THEN
       WRITE(*,*) "THE CELL VOLUME IS NOT POSITIVE ENOUGH: ", OMG, ZERO
       WRITE(*,*) "A: ", A
       WRITE(*,*) "B: ", B
       WRITE(*,*) "C: ", C
       WRITE(*,*) "CODE STOPPED IN 'CARTESIAN_TO_FRACTIONAL()' SUBROUTINE."
       STOP
    END IF
    SGM = SGM / OMG
    DO I = 1, N, 1
       F(1, I) = VECTOR_DOT(R(1:3, I), SGM(1:3, 1))
       F(2, I) = VECTOR_DOT(R(1:3, I), SGM(1:3, 2))
       F(3, I) = VECTOR_DOT(R(1:3, I), SGM(1:3, 3))
    END DO
    RETURN
END SUBROUTINE CARTESIAN_TO_FRACTIONAL




SUBROUTINE CARTESIAN_TO_FRACTIONAL_CONVERT(ABC, R, F, N)
    IMPLICIT NONE
    INTEGER          :: N
    DOUBLE PRECISION :: ABC(3,3), R(3,N), F(3,N)
    CALL  CARTESIAN_TO_FRACTIONAL(ABC(1:3, 1), ABC(1:3, 2), ABC(1:3, 3), R, F, N)
    RETURN
END SUBROUTINE CARTESIAN_TO_FRACTIONAL_CONVERT




SUBROUTINE FRACTIONAL_TO_CARTESIAN(A, B, C, F, R_BACK, N)
    ! TO CALCULATE CARTESIAN COORDINATES BACK (R_BACK ARRAY)
    ! FROM FRACTIONAL ONES (F ARRAY).
    ! A, B, AND C ARE THE THREE KNOWN INDEPENDENT
    ! PERIOD VECTORS .
    IMPLICIT NONE
    INTEGER          :: N, I
    DOUBLE PRECISION :: A(3), B(3), C(3), F(3,N), R_BACK(3,N)
    DO I = 1, N, 1
       R_BACK(1,I) = F(1,I)*A(1) + F(2,I)*B(1) + F(3,I)*C(1)
       R_BACK(2,I) = F(1,I)*A(2) + F(2,I)*B(2) + F(3,I)*C(2)
       R_BACK(3,I) = F(1,I)*A(3) + F(2,I)*B(3) + F(3,I)*C(3)
    END DO
    RETURN
END SUBROUTINE FRACTIONAL_TO_CARTESIAN




SUBROUTINE FRACTIONAL_TO_CARTESIAN_CONVERT(ABC, F, R_BACK, N)
    IMPLICIT NONE
    INTEGER          :: N
    DOUBLE PRECISION :: ABC(3,3), F(3,N), R_BACK(3,N)
    CALL  FRACTIONAL_TO_CARTESIAN(ABC(1:3, 1), ABC(1:3, 2), ABC(1:3, 3), F, R_BACK, N)
    RETURN
END SUBROUTINE FRACTIONAL_TO_CARTESIAN_CONVERT




SUBROUTINE CARTESIAN_TO_FRACTIONAL_VERIFY(A, B, C, R, F, R_BACK, N, &
                                               RELATIVE_ERROR_ALLOWED)
    ! TO VERIFY THE RESULTS BY CALCULATING R_BACK
    ! ARRAY AND COMPARING IT WITH R ARRAY.
    ! A, B, AND C ARE THE THREE KNOWN INDEPENDENT
    ! PERIOD VECTORS .
    IMPLICIT NONE
    INTEGER          :: N, I, J, K
    DOUBLE PRECISION :: A(3), B(3), C(3), R(3,N), F(3,N), R_BACK(3,N)
    DOUBLE PRECISION :: ZERO = 0.1D-30, RELATIVE_ERROR_ALLOWED, RLTV, D1, D2

    CALL FRACTIONAL_TO_CARTESIAN(A, B, C, F, R_BACK, N)

    DO I = 1, N, 1
    DO J = 1, 3, 1
           D1 = ABS(R(J,I))
           D2 = ABS(R_BACK(J,I))
           IF(D1 .LT. ZERO) THEN
                 IF(D2 .LT. ZERO) THEN
                 ELSE
                       D1 = D2
                       CALL COMPARISON()
                 END IF
           ELSE
                 CALL COMPARISON()
           END IF
    END DO
    END DO
    RETURN

    CONTAINS
    SUBROUTINE COMPARISON()
               IMPLICIT NONE
               RLTV = ABS((R(J,I) - R_BACK(J,I)) / D1)
               IF(RLTV .GT. RELATIVE_ERROR_ALLOWED ) THEN
                   WRITE(*,*)
                   WRITE(*,*)
                   WRITE(*,*) "DIFFERENCE: ", R(J,I), R_BACK(J,I), " !!!"
                   WRITE(*,*) "A: ", A
                   WRITE(*,*) "B: ", B
                   WRITE(*,*) "C: ", C
                   WRITE(*,*) "R: ", R(1:3,I)
                   WRITE(*,*) "F: ", F(1:3,I)
                   WRITE(*,*) "R BACK: ", R_BACK(1:3,I)
                   WRITE(*,*) "J, I, N : ", J, I, N
                   WRITE(*,*) "STOPPED IN 'CARTESIAN_TO_FRACTIONAL_VERIFY()'."
                   STOP
               END IF
               RETURN
    END SUBROUTINE COMPARISON
END SUBROUTINE CARTESIAN_TO_FRACTIONAL_VERIFY




SUBROUTINE CARTESIAN_TO_FRACTIONAL_CHECK(ABC, R, F, R_BACK, N, &
                                          RELATIVE_ERROR_ALLOWED)
    IMPLICIT NONE
    INTEGER          :: N
    DOUBLE PRECISION :: ABC(3,3), R(3,N), F(3,N), R_BACK(3,N)
    DOUBLE PRECISION :: RELATIVE_ERROR_ALLOWED
    CALL CARTESIAN_TO_FRACTIONAL_VERIFY(ABC(1:3, 1), ABC(1:3, 2), ABC(1:3, 3), &
                                         R, F, R_BACK, N, RELATIVE_ERROR_ALLOWED)
    RETURN
END SUBROUTINE CARTESIAN_TO_FRACTIONAL_CHECK




SUBROUTINE FRACTIONAL_TO_CARTESIAN_VERIFY(A, B, C, F, R, F_BACK, N, &
                                               RELATIVE_ERROR_ALLOWED)
    ! TO VERIFY THE RESULTS BY CALCULATING F_BACK
    ! ARRAY AND COMPARING IT WITH F ARRAY.
    ! A, B, AND C ARE THE THREE KNOWN INDEPENDENT
    ! PERIOD VECTORS .
    IMPLICIT NONE
    INTEGER          :: N, I, J, K
    DOUBLE PRECISION :: A(3), B(3), C(3), F(3,N), R(3,N), F_BACK(3,N)
    DOUBLE PRECISION :: ZERO = 0.1D-30, RELATIVE_ERROR_ALLOWED, RLTV, D1, D2

    CALL CARTESIAN_TO_FRACTIONAL(A, B, C, R, F_BACK, N)

    DO I = 1, N, 1
    DO J = 1, 3, 1
           D1 = ABS(F(J,I))
           D2 = ABS(F_BACK(J,I))
           IF(D1 .LT. ZERO) THEN
                 IF(D2 .LT. ZERO) THEN
                 ELSE
                       D1 = D2
                       CALL COMPARISON()
                 END IF
           ELSE
                 CALL COMPARISON()
           END IF
    END DO
    END DO
    RETURN

    CONTAINS
    SUBROUTINE COMPARISON()
               IMPLICIT NONE
               RLTV = ABS((F(J,I) - F_BACK(J,I)) / D1)
               IF(RLTV .GT. RELATIVE_ERROR_ALLOWED ) THEN
                   WRITE(*,*)
                   WRITE(*,*)
                   WRITE(*,*) "DIFFERENCE: ", F(J,I), F_BACK(J,I), " !!!"
                   WRITE(*,*) "A: ", A
                   WRITE(*,*) "B: ", B
                   WRITE(*,*) "C: ", C
                   WRITE(*,*) "F: ", F(1:3,I)
                   WRITE(*,*) "R: ", R(1:3,I)
                   WRITE(*,*) "F BACK: ", F_BACK(1:3,I)
                   WRITE(*,*) "J, I, N : ", J, I, N
                   WRITE(*,*) "STOPPED IN 'FRACTIONAL_TO_CARTESIAN_VERIFY()'."
                   STOP
               END IF
               RETURN
    END SUBROUTINE COMPARISON
END SUBROUTINE FRACTIONAL_TO_CARTESIAN_VERIFY




SUBROUTINE FRACTIONAL_TO_CARTESIAN_CHECK(ABC, F, R, F_BACK, N, &
                                          RELATIVE_ERROR_ALLOWED)
    IMPLICIT NONE
    INTEGER          :: N
    DOUBLE PRECISION :: ABC(3,3), F(3,N), R(3,N), F_BACK(3,N)
    DOUBLE PRECISION :: RELATIVE_ERROR_ALLOWED
    CALL FRACTIONAL_TO_CARTESIAN_VERIFY(ABC(1:3, 1), ABC(1:3, 2), ABC(1:3, 3), &
                                         F, R, F_BACK, N, RELATIVE_ERROR_ALLOWED)
    RETURN
END SUBROUTINE FRACTIONAL_TO_CARTESIAN_CHECK




MODULE FOR_READ_AND_WRITE
    IMPLICIT NONE
    INTEGER, PARAMETER :: READ_TUNNEL  = 21
    INTEGER, PARAMETER :: WRITE_TUNNEL = 22
    INTEGER :: END_OF_FILE = -1, CURRENT_LINE = 0
    INTEGER, PARAMETER :: CHARACTER_LENGTH = 10000
    INTEGER, PARAMETER :: CHARACTER_LENGTH2 = CHARACTER_LENGTH*2
    CHARACTER (LEN = CHARACTER_LENGTH2) :: THE_READ_LINE
    INTEGER :: INPUT_FILE_WIDETH=0, INPUT_FILE_LENGTH=0
    INTEGER :: HOW_LAST_PERIOD_CHARACTER = -1


    CONTAINS


    INTEGER FUNCTION THE_READ_UNIT()
               IMPLICIT NONE
               THE_READ_UNIT = READ_TUNNEL
               RETURN
    END FUNCTION THE_READ_UNIT


    INTEGER FUNCTION THE_WRITE_UNIT()
               IMPLICIT NONE
               THE_WRITE_UNIT = WRITE_TUNNEL
               RETURN
    END FUNCTION THE_WRITE_UNIT


    INTEGER FUNCTION THE_MAX_LENGTH_OF_A_LINE()
               IMPLICIT NONE
               THE_MAX_LENGTH_OF_A_LINE = CHARACTER_LENGTH
               RETURN
    END FUNCTION THE_MAX_LENGTH_OF_A_LINE


    SUBROUTINE READ_FILE_OPEN(A)
               IMPLICIT NONE
               CHARACTER (LEN = *) :: A
               OPEN(READ_TUNNEL, FILE=A)
               END_OF_FILE = -1
               CURRENT_LINE = 0
               RETURN
    END SUBROUTINE READ_FILE_OPEN


    SUBROUTINE WRITE_FILE_OPEN(A)
               IMPLICIT NONE
               CHARACTER (LEN = *) ::  A
               OPEN(WRITE_TUNNEL, FILE=A)
               RETURN
    END SUBROUTINE WRITE_FILE_OPEN


    SUBROUTINE READ_FILE_CLOSE()
               IMPLICIT NONE
               CLOSE(READ_TUNNEL)
               RETURN
    END SUBROUTINE READ_FILE_CLOSE


    SUBROUTINE WRITE_FILE_CLOSE()
               IMPLICIT NONE
               CLOSE(WRITE_TUNNEL)
               RETURN
    END SUBROUTINE WRITE_FILE_CLOSE


    SUBROUTINE CHECK_INPUT_FILE(A)
               IMPLICIT NONE
               INTEGER :: LT
               CHARACTER (LEN = *) :: A
               INPUT_FILE_WIDETH = 0
               CALL READ_FILE_OPEN(A)
               A_SIMPL_TEST_DO: DO
                      CALL READ_A_LINE()
                      IF(END_OF_FILE .EQ. 1) EXIT A_SIMPL_TEST_DO
                      LT = LEN_TRIM(THE_READ_LINE)
                      IF(LT .GE. CHARACTER_LENGTH) THEN
                         WRITE(*,*) 'The input file: '//TRIM(A)//' is very wide: ', LT
                         WRITE(*,*) 'Please double CHARACTER_LENGTH in the module FOR_READ_AND_WRITE,'
                         WRITE(*,*) '                                    recompile, and run me again.'
                         STOP
                      END IF
                      IF(INPUT_FILE_WIDETH .LT. LT) THEN
                         INPUT_FILE_WIDETH   =  LT
                         !WRITE(*,*) INPUT_FILE_WIDETH, LT, CURRENT_LINE
                      END IF
               END DO A_SIMPL_TEST_DO
               CALL READ_FILE_CLOSE()
               INPUT_FILE_LENGTH = CURRENT_LINE
               INPUT_FILE_WIDETH = INPUT_FILE_WIDETH + 10
               RETURN
    END SUBROUTINE CHECK_INPUT_FILE


    SUBROUTINE READ_A_LINE()
               IMPLICIT NONE
               10 FORMAT(A)
               END_OF_FILE = -1
               THE_READ_LINE = ' '
               READ(READ_TUNNEL, 10, END = 50) THE_READ_LINE
               CURRENT_LINE = CURRENT_LINE + 1
               RETURN
               50 END_OF_FILE = 1
               RETURN
    END SUBROUTINE READ_A_LINE


    SUBROUTINE WRITE_A_LINE()
               IMPLICIT NONE
               10 FORMAT(A)
               IF(LEN_TRIM(THE_READ_LINE) .LE. 0) THEN
                  WRITE(WRITE_TUNNEL, *)
               ELSE
                  WRITE(WRITE_TUNNEL, 10) TRIM(THE_READ_LINE)
               END IF
               RETURN
    END SUBROUTINE WRITE_A_LINE


    SUBROUTINE SIMPLY_READ_AND_WRITE_A_LINE()
               IMPLICIT NONE
               CALL READ_A_LINE()
               IF(END_OF_FILE .NE. 1) CALL WRITE_A_LINE()
               RETURN
    END SUBROUTINE SIMPLY_READ_AND_WRITE_A_LINE


    SUBROUTINE SIMPLY_READ_AND_WRITE_N_LINES(N)
               IMPLICIT NONE
               INTEGER :: I, N
               DO I = 1, N, 1
                  CALL SIMPLY_READ_AND_WRITE_A_LINE()
               END DO
               RETURN
    END SUBROUTINE SIMPLY_READ_AND_WRITE_N_LINES


    SUBROUTINE READ_AN_AREA(A, L, N)
               IMPLICIT NONE
               INTEGER             :: I, L, N
               CHARACTER (LEN = L) :: A(N)
               IF(L .LT. INPUT_FILE_WIDETH) THEN
                     WRITE(*,*) 'SORRY:          CHARACTER ARRAY A IS TOO NARROW IN READ_AN_AREA().'
                     WRITE(*,*) 'SORRY: ', CURRENT_LINE, L, INPUT_FILE_WIDETH, ' IN READ_AN_AREA().'
                     STOP
               END IF
               IF(END_OF_FILE .EQ. 1) THEN
                     WRITE(*,*) 'SORRY:                 FILE ALREADY ENDED IN READ_AN_AREA().'
                     WRITE(*,*) 'SORRY: ', N, CURRENT_LINE, END_OF_FILE, ' IN READ_AN_AREA().'
                     STOP
               END IF
               DO I = 1, N, 1
                  CALL READ_A_LINE()
                  IF(END_OF_FILE .EQ. 1) THEN
                     WRITE(*,*) 'SORRY: ', N, CURRENT_LINE, END_OF_FILE, ' IN READ_AN_AREA().'
                     STOP
                  END IF
                  A(I)(1:L) = THE_READ_LINE(1:L)
               END DO
               RETURN
    END SUBROUTINE READ_AN_AREA


    SUBROUTINE WRITE_AN_AREA(A, L, N)
               IMPLICIT NONE
               INTEGER             :: I, L, N
               CHARACTER (LEN = L) :: A(N)
               10 FORMAT(A)
               DO I = 1, N, 1
                  IF(LEN_TRIM(A(I)) .LE. 0) THEN
                     WRITE(WRITE_TUNNEL, *)
                  ELSE
                     WRITE(WRITE_TUNNEL, 10) TRIM(A(I))
                  END IF
               END DO
               RETURN
    END SUBROUTINE WRITE_AN_AREA


    SUBROUTINE CHRACTER_ARRAY_EXPANSION(A, L, O, N)
               IMPLICIT NONE
               INTEGER             :: L, O, N
               CHARACTER (LEN = L), ALLOCATABLE :: A(:), B(:)
               ALLOCATE (B(O))
               B = A
               DEALLOCATE (A)
               ALLOCATE (A(N))
               A = ' '
               A(1:O)(1:L) = B(1:O)(1:L)
               DEALLOCATE (B)
               RETURN
    END SUBROUTINE CHRACTER_ARRAY_EXPANSION


    SUBROUTINE READ_AND_COPY_A_LINE(LINES_TO_STORE, L, O, I)
               IMPLICIT NONE
               INTEGER :: L, O, N, I
               CHARACTER(LEN=L), ALLOCATABLE:: LINES_TO_STORE(:)
               IF(L .LT. INPUT_FILE_WIDETH) THEN
                     WRITE(*,*) 'SORRY:            CHARACTER ARRAY IS TOO NARROW IN READ_AND_COPY_A_LINE().'
                     WRITE(*,*) 'SORRY: ', CURRENT_LINE, L, INPUT_FILE_WIDETH, ' IN READ_AND_COPY_A_LINE().'
                     STOP
               END IF
               IF(END_OF_FILE .EQ. 1) THEN
                     WRITE(*,*) 'SORRY:                 FILE ALREADY ENDED IN READ_AND_COPY_A_LINE().'
                     WRITE(*,*) 'SORRY: ', N, CURRENT_LINE, END_OF_FILE, ' IN READ_AND_COPY_A_LINE().'
                     STOP
               END IF
               CALL READ_A_LINE()
               IF(END_OF_FILE .EQ. 1) THEN
                  WRITE(*,*) 'SORRY, NOT ENOUGH LINES TO BE READ: ', CURRENT_LINE, L, I, O, ' IN READ_AND_COPY_A_LINE().'
                  STOP
               END IF
               IF(O.LT.I) THEN
                  N = I + 10
                  CALL CHRACTER_ARRAY_EXPANSION(LINES_TO_STORE, L, O, N)
                  O = N
               END IF
               LINES_TO_STORE(I)(1:L) = THE_READ_LINE(1:L)
               RETURN
    END SUBROUTINE READ_AND_COPY_A_LINE


    SUBROUTINE PERIOD_ARRAY_EXPANSION(A, N)
               IMPLICIT NONE
               INTEGER          :: NOLD(3), NNEW(3), N
               DOUBLE PRECISION, ALLOCATABLE  :: A(:, :, :), B(:, :, :)
               NOLD = SHAPE (A)
               ALLOCATE (B (NOLD(1), NOLD(2), NOLD(3)))
               B = A
               DEALLOCATE (A)
               NNEW = NOLD
               NNEW(3) = N
               ALLOCATE (A (NNEW(1), NNEW(2), NNEW(3)))
               A = 0.0D0
               A(1:NOLD(1), 1:NOLD(2), 1:NOLD(3)) = B(1:NOLD(1), 1:NOLD(2), 1:NOLD(3))
               DEALLOCATE (B)
               RETURN
    END SUBROUTINE PERIOD_ARRAY_EXPANSION


    SUBROUTINE RR_ARRAY_EXPANSION(A, N)
               IMPLICIT NONE
               INTEGER          :: NOLD(2), NNEW(2), N
               DOUBLE PRECISION, ALLOCATABLE  :: A(:, :), B(:, :)
               NOLD = SHAPE (A)
               ALLOCATE (B (NOLD(1), NOLD(2)))
               B = A
               DEALLOCATE (A)
               NNEW = NOLD
               NNEW(2) = N
               ALLOCATE (A (NNEW(1), NNEW(2)))
               A = 0.0D0
               A(1:NOLD(1), 1:NOLD(2)) = B(1:NOLD(1), 1:NOLD(2))
               DEALLOCATE (B)
               RETURN
    END SUBROUTINE RR_ARRAY_EXPANSION


    SUBROUTINE II_ARRAY_EXPANSION(A, N)
               IMPLICIT NONE
               INTEGER          :: NOLD(2), NNEW(2), N
               INTEGER, ALLOCATABLE  :: A(:, :), B(:, :)
               NOLD = SHAPE (A)
               ALLOCATE (B (NOLD(1), NOLD(2)))
               B = A
               DEALLOCATE (A)
               NNEW = NOLD
               NNEW(2) = N
               ALLOCATE (A (NNEW(1), NNEW(2)))
               A = 0
               A(1:NOLD(1), 1:NOLD(2)) = B(1:NOLD(1), 1:NOLD(2))
               DEALLOCATE (B)
               RETURN
    END SUBROUTINE II_ARRAY_EXPANSION


    SUBROUTINE I_ARRAY_EXPANSION(A, N)
               IMPLICIT NONE
               INTEGER          :: NOLD, N
               INTEGER, ALLOCATABLE  :: A(:), B(:)
               NOLD = SIZE (A)
               ALLOCATE (B (NOLD))
               B = A
               DEALLOCATE (A)
               ALLOCATE (A (N))
               A = 0
               A(1:NOLD) = B(1:NOLD)
               DEALLOCATE (B)
               RETURN
    END SUBROUTINE I_ARRAY_EXPANSION


    SUBROUTINE A_VALUE_TO_AN_EXPANDABLE_ARRAY(A, SIZE, TO, VALUE)
               IMPLICIT NONE
               INTEGER :: SIZE, TO, VALUE
               INTEGER, ALLOCATABLE :: A(:)
               IF(TO .GT. SIZE) THEN
                  SIZE = TO + 10
                  CALL I_ARRAY_EXPANSION(A, SIZE)
               END IF
               A(TO) = VALUE
               RETURN
    END SUBROUTINE A_VALUE_TO_AN_EXPANDABLE_ARRAY


    SUBROUTINE LAST_DECIMAL_DIGIT(A, I, L)
               IMPLICIT NONE
               INTEGER :: I, J, K, L, M, N
               CHARACTER(LEN = *):: A
               J = LEN_TRIM(A)
               IF(J .LT. I) THEN
                  WRITE(*,*) "The line # ", CURRENT_LINE
                  WRITE(*,*) TRIM(A)
                  WRITE(*,*) 'of the input file has nothing after ', I, ' characters.'
                  STOP
               END IF

               CHECK_DO: DO M = J, I, -1
                    IF(A(M:M) .EQ. '.') THEN
                       IF(HOW_LAST_PERIOD_CHARACTER .LE. 2) THEN
                          TALK_DO: DO
                             WRITE(*,*) "The last effective character in the line # ", CURRENT_LINE
                             WRITE(*,*) 'of the input file is the period "." characters: '
                             WRITE(*,*) A(1:M)
                             WRITE(*,*) 'What do you want to do?'
                             WRITE(*,*) "     ignored  for now,  reply me with 1 please;"
                             WRITE(*,*) "     as digit for now,  reply me with 2 please;"
                             WRITE(*,*) "     always   ignored,  reply me with 3 please;"
                             WRITE(*,*) "     always  as digit,  reply me with 4 please;"
                             WRITE(*,*) "     STOP me, then you will change the routine "
                             WRITE(*,*) '     "LAST_DECIMAL_DIGIT()", reply me with 5 please.'
                             READ(*,*) N
                             IF(N .EQ. 5) STOP
                             IF((N .LE. 4) .AND. (N .GE. 1)) EXIT TALK_DO
                          END DO TALK_DO
                          IF(N .GE. 3) HOW_LAST_PERIOD_CHARACTER = N
                          IF(N .EQ. 2*(N/2)) THEN
                             L = M
                             RETURN
                          ELSE
                             CYCLE CHECK_DO
                          END IF
                       ELSE IF(HOW_LAST_PERIOD_CHARACTER .EQ. 3) THEN
                             CYCLE CHECK_DO
                       ELSE IF(HOW_LAST_PERIOD_CHARACTER .EQ. 4) THEN
                             L = M
                             RETURN
                       END IF
                    END IF

                    K = IACHAR(A(M:M))
                    IF((K .GE. 48) .AND. (K .LE. 57)) THEN
                        L = M
                        RETURN
                    END IF
               END DO CHECK_DO

               WRITE(*,*) "Sorry no decimal digit is found in the line # ", CURRENT_LINE
               WRITE(*,*) TRIM(A)
               WRITE(*,*) 'of the input file after ', I, ' characters.'
               STOP

               RETURN
    END SUBROUTINE LAST_DECIMAL_DIGIT


END MODULE FOR_READ_AND_WRITE




MODULE KEY_WORD_GROUPS
    IMPLICIT NONE
    INTEGER            :: CURRENT_GROUP
    INTEGER, PARAMETER :: NUMBER_OF_KEY_WORD_GROUPS = 3
    CHARACTER(LEN = *), PARAMETER:: IDENTIFIER_WORDS(NUMBER_OF_KEY_WORD_GROUPS)&
                    = (/ 'A_very_nice_molecule_to_be_preocessed             ', &
                         'G E O M E T R Y    I N    X - Y - Z    F O R M A T', &
                         'Coordinates in Geometry Cycle                     '  &
                                                                                /)
    LOGICAL, PARAMETER :: IDENTIFIER_LINE_KEPT(NUMBER_OF_KEY_WORD_GROUPS) = (/ &
                                              .FALSE.,                         &
                                              .TRUE.,                          &
                                              .TRUE.                           &
                                                                                /)
    INTEGER, PARAMETER :: NUM_OF_FOLLOWING_LINES_SKIPPED(NUMBER_OF_KEY_WORD_GROUPS) &
                                 = (/          0,                                   &
                                               2,                                   &
                                               1                                    &
                                                                                     /)
END MODULE KEY_WORD_GROUPS




SUBROUTINE DIRECTION()
    USE KEY_WORD_GROUPS
    IMPLICIT NONE
    SELECT CASE (CURRENT_GROUP)
    CASE (1)
              CALL NICE_MOLECULE_PROCESSING()
    CASE (2)
              CALL XYZ_FORMAT()
    CASE (3)
              CALL GEOMETRY()

    CASE DEFAULT
              WRITE(*,*) 'UNACCEPTABLE "CURRENT_GROUP": ', CURRENT_GROUP
              WRITE(*,*) 'STOPPED IN SUBROUTINE DIRECTION().'
              STOP
    END SELECT
    RETURN
END SUBROUTINE DIRECTION




SUBROUTINE NICE_MOLECULE_PROCESSING()
    USE FOR_READ_AND_WRITE
    USE KEY_WORD_GROUPS
    IMPLICIT NONE
    CHARACTER (LEN = * ), PARAMETER :: SUBGROUP_KEY_WORDS='ATOMS'
    CHARACTER (LEN = * ), PARAMETER :: END_SUBGROUP_KEY_WORDS='END'
    CHARACTER (LEN = * ), PARAMETER :: ADDITIONAL_KEY_WORDS='Lattice'
    INTEGER, PARAMETER :: NOT_FOUND = -1, JUST_FOUND = 0
    INTEGER, PARAMETER :: DIMENSIONS = 3, INITIAL_SIZE = 2
    CHARACTER(LEN=INPUT_FILE_WIDETH), ALLOCATABLE :: READ_LINES_OF_THE_REGION(:)
    INTEGER, ALLOCATABLE :: LINES_TO_BE_CONVERTED(:)
    LOGICAL :: SUBGROUP_FOUND
    INTEGER :: ADDITIONAL_FOUND, FIRST_ADDITIONAL_LINE, SIZE_OF_LINES_TO_BE_CONVERTED
    INTEGER :: I, J, K, REGION_SIZE, TOTAL_OF_LINES_TO_BE_CONVERTED, TOTAL_LINES_OF_THE_REGION
    DOUBLE PRECISION    :: PERIOD_VECTORS(DIMENSIONS, DIMENSIONS)
    DOUBLE PRECISION, ALLOCATABLE :: R(:,:), F(:,:), R_BACK(:,:)
    DOUBLE PRECISION :: RELATIVE_ERROR_ALLOWED = 0.1D-12

    REGION_SIZE = INITIAL_SIZE
    SIZE_OF_LINES_TO_BE_CONVERTED = INITIAL_SIZE
    !No worried about the sizes of the following two arrays. They will be expanded automatically, if needed.
    ALLOCATE (READ_LINES_OF_THE_REGION(REGION_SIZE))
    ALLOCATE (LINES_TO_BE_CONVERTED(SIZE_OF_LINES_TO_BE_CONVERTED))
    TOTAL_OF_LINES_TO_BE_CONVERTED = 0
    ADDITIONAL_FOUND = NOT_FOUND
    SUBGROUP_FOUND = .FALSE.

    READ_AND_ANALYZE_INPUT_FILE_DO: DO I = 1, INPUT_FILE_LENGTH, 1

        CALL READ_AND_COPY_A_LINE(READ_LINES_OF_THE_REGION, INPUT_FILE_WIDETH, REGION_SIZE, I)

        IF(INDEX(READ_LINES_OF_THE_REGION(I), TRIM(SUBGROUP_KEY_WORDS)) .GT. 0) THEN
           SUBGROUP_FOUND = .TRUE.
           CYCLE
        END IF
        IF(INDEX(READ_LINES_OF_THE_REGION(I), TRIM(END_SUBGROUP_KEY_WORDS)) .GT. 0) THEN
           SUBGROUP_FOUND = .FALSE.
           CYCLE
        END IF
        IF(INDEX(READ_LINES_OF_THE_REGION(I), TRIM(ADDITIONAL_KEY_WORDS)) .GT. 0) THEN
           ADDITIONAL_FOUND = JUST_FOUND
           CYCLE
        END IF

        IF(SUBGROUP_FOUND) THEN
           TOTAL_OF_LINES_TO_BE_CONVERTED = TOTAL_OF_LINES_TO_BE_CONVERTED + 1
           CALL A_VALUE_TO_AN_EXPANDABLE_ARRAY(LINES_TO_BE_CONVERTED, SIZE_OF_LINES_TO_BE_CONVERTED, &
                                                                     TOTAL_OF_LINES_TO_BE_CONVERTED, I)
           CYCLE
        END IF

        IF( ADDITIONAL_FOUND .EQ. JUST_FOUND) FIRST_ADDITIONAL_LINE = I
        IF( ADDITIONAL_FOUND .GE. JUST_FOUND) ADDITIONAL_FOUND = ADDITIONAL_FOUND + 1
        IF((ADDITIONAL_FOUND .EQ. DIMENSIONS)) THEN
            TOTAL_LINES_OF_THE_REGION = I
            EXIT READ_AND_ANALYZE_INPUT_FILE_DO
        END IF

    END DO READ_AND_ANALYZE_INPUT_FILE_DO


    ALLOCATE (R(DIMENSIONS,      TOTAL_OF_LINES_TO_BE_CONVERTED))
    ALLOCATE (F(DIMENSIONS,      TOTAL_OF_LINES_TO_BE_CONVERTED))
    ALLOCATE (R_BACK(DIMENSIONS, TOTAL_OF_LINES_TO_BE_CONVERTED))

    READ_R_DATA_DO: DO I = 1, TOTAL_OF_LINES_TO_BE_CONVERTED, 1
        J = LINES_TO_BE_CONVERTED(I)
        CALL LAST_DECIMAL_DIGIT(READ_LINES_OF_THE_REGION(J), 1, K)
        READ(READ_LINES_OF_THE_REGION(J)(1:K), *) R(1:DIMENSIONS, I)
        !WRITE(*, *) I, R(1:DIMENSIONS, I)
    END DO READ_R_DATA_DO

    READ_THE_PERIOD_DO: DO I = 1, DIMENSIONS, 1
        CALL LAST_DECIMAL_DIGIT(READ_LINES_OF_THE_REGION(FIRST_ADDITIONAL_LINE + I - 1), 1, K)
        READ(READ_LINES_OF_THE_REGION(FIRST_ADDITIONAL_LINE + I - 1)(1:K), *) PERIOD_VECTORS(1:DIMENSIONS, I)
        !WRITE(*, *) I, PERIOD_VECTORS(1:DIMENSIONS, I)
    END DO READ_THE_PERIOD_DO

    CALL CARTESIAN_TO_FRACTIONAL_CONVERT(PERIOD_VECTORS, R, F, TOTAL_OF_LINES_TO_BE_CONVERTED)
    CALL CARTESIAN_TO_FRACTIONAL_CHECK(  PERIOD_VECTORS, R, F, R_BACK, &
                TOTAL_OF_LINES_TO_BE_CONVERTED, RELATIVE_ERROR_ALLOWED)


    WRITE_R_DATA_DO: DO I = 1, TOTAL_OF_LINES_TO_BE_CONVERTED, 1
        J = LINES_TO_BE_CONVERTED(I)
        READ_LINES_OF_THE_REGION(J) = ' '
        WRITE(READ_LINES_OF_THE_REGION(J), '(3F13.8)') F(1:DIMENSIONS, I)
        !WRITE(*, *) I, F(1:DIMENSIONS, I)
    END DO WRITE_R_DATA_DO


    !WRITE(*,*) FIRST_ADDITIONAL_LINE, TOTAL_OF_LINES_TO_BE_CONVERTED, TOTAL_LINES_OF_THE_REGION
    !WRITE(*,*) LINES_TO_BE_CONVERTED
    CALL WRITE_AN_AREA(READ_LINES_OF_THE_REGION, INPUT_FILE_WIDETH, TOTAL_LINES_OF_THE_REGION)
    DEALLOCATE (READ_LINES_OF_THE_REGION, LINES_TO_BE_CONVERTED, R, F, R_BACK)

    RETURN
END SUBROUTINE NICE_MOLECULE_PROCESSING




SUBROUTINE XYZ_FORMAT()
    USE FOR_READ_AND_WRITE
    USE KEY_WORD_GROUPS
    IMPLICIT NONE
    CHARACTER (LEN = * ), PARAMETER :: ADDITIONAL_KEY_WORDS='VEC'
    CHARACTER (LEN = * ), PARAMETER :: ENDING_KEY_WORDS='Total nr. of atoms:'
    INTEGER, PARAMETER :: DIMENSIONS = 3, R_STARTING_POSITION = 5, P_STARTING_POSITION = 6
    CHARACTER(LEN=INPUT_FILE_WIDETH), ALLOCATABLE :: READ_LINES_OF_THE_REGION(:)
    INTEGER :: I, J, K, L, FIRST_ADDITIONAL_LINE
    INTEGER :: REGION_SIZE, TOTAL_OF_LINES_TO_BE_CONVERTED, TOTAL_LINES_OF_THE_REGION
    DOUBLE PRECISION    :: PERIOD_VECTORS(DIMENSIONS, DIMENSIONS)
    DOUBLE PRECISION, ALLOCATABLE :: R(:,:), F(:,:), R_BACK(:,:)
    DOUBLE PRECISION :: RELATIVE_ERROR_ALLOWED = 0.1D-12

    REGION_SIZE = 2
    !No worried about the size of the following array. It will be expanded automatically, if needed.
    ALLOCATE (READ_LINES_OF_THE_REGION(REGION_SIZE))

    FIRST_ADDITIONAL_LINE = -1
    READ_AND_ANALYZE_INPUT_FILE_DO: DO I = 1, INPUT_FILE_LENGTH, 1
        CALL READ_AND_COPY_A_LINE(READ_LINES_OF_THE_REGION, INPUT_FILE_WIDETH, REGION_SIZE, I)
        IF(INDEX(READ_LINES_OF_THE_REGION(I), TRIM(ADDITIONAL_KEY_WORDS)) .GT. 0) THEN
           IF(FIRST_ADDITIONAL_LINE .LE. 0) FIRST_ADDITIONAL_LINE = I
           CYCLE
        END IF
        IF(INDEX(READ_LINES_OF_THE_REGION(I), TRIM(ENDING_KEY_WORDS)) .GT. 0) THEN
           TOTAL_LINES_OF_THE_REGION = I
           EXIT READ_AND_ANALYZE_INPUT_FILE_DO
        END IF
    END DO READ_AND_ANALYZE_INPUT_FILE_DO

    TOTAL_OF_LINES_TO_BE_CONVERTED = FIRST_ADDITIONAL_LINE - 1
    ALLOCATE (R(DIMENSIONS,      TOTAL_OF_LINES_TO_BE_CONVERTED))
    ALLOCATE (F(DIMENSIONS,      TOTAL_OF_LINES_TO_BE_CONVERTED))
    ALLOCATE (R_BACK(DIMENSIONS, TOTAL_OF_LINES_TO_BE_CONVERTED))

    READ_R_DATA_DO: DO I = 1, TOTAL_OF_LINES_TO_BE_CONVERTED, 1
        CALL LAST_DECIMAL_DIGIT(READ_LINES_OF_THE_REGION(I), R_STARTING_POSITION, K)
        READ(READ_LINES_OF_THE_REGION(I)(R_STARTING_POSITION:K), *) R(1:DIMENSIONS, I)
        !WRITE(*, *) I, R(1:DIMENSIONS, I)
    END DO READ_R_DATA_DO

    READ_THE_PERIOD_DO: DO I = 1, DIMENSIONS, 1
        CALL LAST_DECIMAL_DIGIT(READ_LINES_OF_THE_REGION(FIRST_ADDITIONAL_LINE + I - 1), P_STARTING_POSITION, K)
        READ(READ_LINES_OF_THE_REGION(FIRST_ADDITIONAL_LINE + I - 1)(P_STARTING_POSITION:K), *) &
                                                                  PERIOD_VECTORS(1:DIMENSIONS, I)
        !WRITE(*, *) I, PERIOD_VECTORS(1:DIMENSIONS, I)
    END DO READ_THE_PERIOD_DO

    CALL CARTESIAN_TO_FRACTIONAL_CONVERT(PERIOD_VECTORS, R, F, TOTAL_OF_LINES_TO_BE_CONVERTED)
    CALL CARTESIAN_TO_FRACTIONAL_CHECK(  PERIOD_VECTORS, R, F, R_BACK, &
                TOTAL_OF_LINES_TO_BE_CONVERTED, RELATIVE_ERROR_ALLOWED)


    WRITE_R_DATA_DO: DO I = 1, TOTAL_OF_LINES_TO_BE_CONVERTED, 1
              READ_LINES_OF_THE_REGION(I)(R_STARTING_POSITION:) = ' '
        WRITE(READ_LINES_OF_THE_REGION(I)(R_STARTING_POSITION:), '(3F16.9)') F(1:DIMENSIONS, I)
        !WRITE(*, *) I, F(1:DIMENSIONS, I)
    END DO WRITE_R_DATA_DO


    !WRITE(*,*) FIRST_ADDITIONAL_LINE, TOTAL_OF_LINES_TO_BE_CONVERTED, TOTAL_LINES_OF_THE_REGION
    CALL WRITE_AN_AREA(READ_LINES_OF_THE_REGION, INPUT_FILE_WIDETH, TOTAL_LINES_OF_THE_REGION)
    DEALLOCATE (READ_LINES_OF_THE_REGION, R, F, R_BACK)

    RETURN
END SUBROUTINE XYZ_FORMAT




SUBROUTINE GEOMETRY()
    USE FOR_READ_AND_WRITE
    USE KEY_WORD_GROUPS
    IMPLICIT NONE
    CHARACTER (LEN = * ), PARAMETER :: ADDITIONAL_KEY_WORDS='Lattice Vectors'
    INTEGER, PARAMETER :: DIMENSIONS = 3, NOT_FOUND = -1, JUST_FOUND = 0
    INTEGER, PARAMETER :: R_STARTING_POSITION = 11, P_STARTING_POSITION = 1
    CHARACTER(LEN=INPUT_FILE_WIDETH), ALLOCATABLE :: READ_LINES_OF_THE_REGION(:)
    INTEGER :: I, J, K, L, ADDITIONAL_FOUND, FIRST_ADDITIONAL_LINE
    INTEGER :: REGION_SIZE, TOTAL_OF_LINES_TO_BE_CONVERTED, TOTAL_LINES_OF_THE_REGION
    DOUBLE PRECISION    :: PERIOD_VECTORS(DIMENSIONS, DIMENSIONS)
    DOUBLE PRECISION, ALLOCATABLE :: R(:,:), F(:,:), R_BACK(:,:)
    DOUBLE PRECISION :: RELATIVE_ERROR_ALLOWED = 0.1D-12

    REGION_SIZE = 2
    !No worried about the size of the following array. It will be expanded automatically, if needed.
    ALLOCATE (READ_LINES_OF_THE_REGION(REGION_SIZE))

    ADDITIONAL_FOUND = NOT_FOUND
    READ_AND_ANALYZE_INPUT_FILE_DO: DO I = 1, INPUT_FILE_LENGTH, 1
        CALL READ_AND_COPY_A_LINE(READ_LINES_OF_THE_REGION, INPUT_FILE_WIDETH, REGION_SIZE, I)
        IF(INDEX(READ_LINES_OF_THE_REGION(I), TRIM(ADDITIONAL_KEY_WORDS)) .GT. 0) THEN
           ADDITIONAL_FOUND = JUST_FOUND
           CYCLE
        END IF
        IF( ADDITIONAL_FOUND .EQ. JUST_FOUND) FIRST_ADDITIONAL_LINE = I
        IF( ADDITIONAL_FOUND .GE. JUST_FOUND) ADDITIONAL_FOUND = ADDITIONAL_FOUND + 1
        IF((ADDITIONAL_FOUND .EQ. DIMENSIONS)) THEN
            TOTAL_LINES_OF_THE_REGION = I
            EXIT READ_AND_ANALYZE_INPUT_FILE_DO
        END IF
    END DO READ_AND_ANALYZE_INPUT_FILE_DO

    TOTAL_OF_LINES_TO_BE_CONVERTED = FIRST_ADDITIONAL_LINE - 2
    ALLOCATE (R(DIMENSIONS,      TOTAL_OF_LINES_TO_BE_CONVERTED))
    ALLOCATE (F(DIMENSIONS,      TOTAL_OF_LINES_TO_BE_CONVERTED))
    ALLOCATE (R_BACK(DIMENSIONS, TOTAL_OF_LINES_TO_BE_CONVERTED))

    READ_R_DATA_DO: DO I = 1, TOTAL_OF_LINES_TO_BE_CONVERTED, 1
        CALL LAST_DECIMAL_DIGIT(READ_LINES_OF_THE_REGION(I), R_STARTING_POSITION, K)
        READ(READ_LINES_OF_THE_REGION(I)(R_STARTING_POSITION: K), *) R(1:DIMENSIONS, I)
        !WRITE(*, *) I, R(1:DIMENSIONS, I)
    END DO READ_R_DATA_DO

    READ_THE_PERIOD_DO: DO I = 1, DIMENSIONS, 1
        CALL LAST_DECIMAL_DIGIT(READ_LINES_OF_THE_REGION(FIRST_ADDITIONAL_LINE + I - 1), P_STARTING_POSITION, K)
        READ(READ_LINES_OF_THE_REGION(FIRST_ADDITIONAL_LINE + I - 1)(P_STARTING_POSITION:K), *) &
                                                                   PERIOD_VECTORS(1:DIMENSIONS, I)
        !WRITE(*, *) I, PERIOD_VECTORS(1:DIMENSIONS, I)
    END DO READ_THE_PERIOD_DO

    CALL CARTESIAN_TO_FRACTIONAL_CONVERT(PERIOD_VECTORS, R, F, TOTAL_OF_LINES_TO_BE_CONVERTED)
    CALL CARTESIAN_TO_FRACTIONAL_CHECK(  PERIOD_VECTORS, R, F, R_BACK, &
                TOTAL_OF_LINES_TO_BE_CONVERTED, RELATIVE_ERROR_ALLOWED)


    WRITE_R_DATA_DO: DO I = 1, TOTAL_OF_LINES_TO_BE_CONVERTED, 1
              READ_LINES_OF_THE_REGION(I)(R_STARTING_POSITION:) = ' '
        WRITE(READ_LINES_OF_THE_REGION(I)(R_STARTING_POSITION:), '(3F15.9)') F(1:DIMENSIONS, I)
        !WRITE(*, *) I, F(1:DIMENSIONS, I)
    END DO WRITE_R_DATA_DO


    !WRITE(*,*) FIRST_ADDITIONAL_LINE, TOTAL_OF_LINES_TO_BE_CONVERTED, TOTAL_LINES_OF_THE_REGION
    CALL WRITE_AN_AREA(READ_LINES_OF_THE_REGION, INPUT_FILE_WIDETH, TOTAL_LINES_OF_THE_REGION)
    DEALLOCATE (READ_LINES_OF_THE_REGION, R, F, R_BACK)

    RETURN
END SUBROUTINE GEOMETRY




PROGRAM DATA_SECTIONS_TO_BE_CONVERTED
    USE FOR_READ_AND_WRITE
    USE KEY_WORD_GROUPS
    IMPLICIT NONE
    INTEGER :: DATA_SECTIONS_TRIED_PER_GROUP(NUMBER_OF_KEY_WORD_GROUPS)

    CALL CHECK_INPUT_FILE('in.dat')
    WRITE(*,*) INPUT_FILE_WIDETH, INPUT_FILE_LENGTH

    CALL  READ_FILE_OPEN('in.dat')
    CALL WRITE_FILE_OPEN('out.dat')
    DATA_SECTIONS_TRIED_PER_GROUP = 0
    MAIN_LOOP: DO
               CALL READ_A_LINE()
               IF(END_OF_FILE .EQ. 1) EXIT MAIN_LOOP
               DO CURRENT_GROUP = 1, NUMBER_OF_KEY_WORD_GROUPS ,1
                  IF(INDEX(THE_READ_LINE, TRIM(IDENTIFIER_WORDS(CURRENT_GROUP))) .GT. 0) THEN
                     DATA_SECTIONS_TRIED_PER_GROUP(CURRENT_GROUP) = DATA_SECTIONS_TRIED_PER_GROUP(CURRENT_GROUP) + 1
                     IF(IDENTIFIER_LINE_KEPT(CURRENT_GROUP)) CALL WRITE_A_LINE()
                     CALL  SIMPLY_READ_AND_WRITE_N_LINES(NUM_OF_FOLLOWING_LINES_SKIPPED(CURRENT_GROUP))
                     CALL  DIRECTION()
                     CYCLE MAIN_LOOP
                  END IF
               END DO
               CALL WRITE_A_LINE()
    END DO MAIN_LOOP
    CALL  READ_FILE_CLOSE()
    CALL WRITE_FILE_CLOSE()
    WRITE(*,*) "Total ", CURRENT_LINE, " lines of the input file have been processed."
    WRITE(*,*) "Sections found to be or have been converted per group: "
    WRITE(*,*) "    ", DATA_SECTIONS_TRIED_PER_GROUP
    WRITE(*,*) "    then total: ", SUM(DATA_SECTIONS_TRIED_PER_GROUP)

    STOP
END PROGRAM DATA_SECTIONS_TO_BE_CONVERTED



