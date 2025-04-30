*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE TRI (N,NTAB,RTAB)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Auteur : D.Martin (Mars 1994)
*
* Procedure triant le tableau NTAB selon l'ordre croissant,
* les permutations etant aussi effectuees sur RTAB
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL    RTAB(*)
      INTEGER N,NTAB(*)
*
      INTEGER I,J,AUX1
      REAL    AUX2
      LOGICAL PERMUT
*
      I=1
1     PERMUT=.FALSE.
      I=I+1
*
      DO 2 J=N,I,-1
        IF (NTAB(J).LT.NTAB(J-1)) THEN
          AUX1=NTAB(J-1)
          NTAB(J-1)=NTAB(J)
          NTAB(J)=AUX1
          AUX2=RTAB(J-1)
          RTAB(J-1)=RTAB(J)
          RTAB(J)=AUX2
          PERMUT=.TRUE.
        END IF
2     CONTINUE
*   
      IF (PERMUT) GOTO 1  
      END    
