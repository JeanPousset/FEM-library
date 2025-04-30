*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE CDESSE (NBLIGN,ADPRCL,NUMCOL,ADSUCL,MATRIS,SECMBR
     &                  ,NUDDIR,VALDIR,ADPRC0,NUMCO0,MATRI0,SECMB0) 
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Auteur : D.Martin (Mars 1995)
*
* Procedure de prise en compte des d.l. Dirichlet et de creation d'un 
* systeme lineaire apres rangement des coefficients d'une ligne par numeros
* de colonne croissant
*
*-- Arguments d'entree --
*------- Systeme lineaire non ordonne avant elimination
* NBLIGN nombre de lignes de la matrice
* ADPRCL adresses des 1er coeff. des lignes de la partie triang. inf. stricte
* NUMCOL numeros de colonnes des coeff. de la partie triang. inf. stricte
* ADSUCL adresse du coeff. suivant sur la meme ligne
* MATRIS tableau contenant la matrice A
* SECMBR second membre en entree
* NUDDIR indicateurs des d.l. de Dirichlet
* VALDIR valeurs de blocage des d.l. Dirichlet
*-- Arguments de sortie --
*--------- Systeme ordonne apres elimination des d.l. Dirichlet
* ADPRC0 adresses des 1er coeff. des lignes de la partie triang. inf. stricte
* NUMCO0 numeros de colonnes des coeff. de la partie triang. inf. stricte
* MATRI0 tableau contenant la matrice en sortie A0
* SECMBR second membre en sortie
*
*-- Remarques --
* Les tableaux ADPRCL et ADPRC0 peuvent coincider
* Les tableaux SECMBR et SECMB0 peuvent coincider
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER NBLIGN,ADPRCL(*),NUMCOL(*),ADSUCL(*),NUDDIR(*)
     &       ,ADPRC0(*),NUMCO0(*)   
      REAL    MATRIS(*),SECMBR(*),MATRI0(*),SECMB0(*),VALDIR(*) 
*
      INTEGER I,J,ADMATR,ADMAT0,ADMAT1,ADMATX
*----------------------------------------------------------------------
      DO 1 I = 1,NBLIGN
         SECMB0(I) = SECMBR(I)
         MATRI0(I) = MATRIS(I)
1     CONTINUE
*
      IF (NUDDIR(1).LT.0) THEN
         SECMB0(1) = MATRI0(1)*VALDIR(1)
      ELSEIF (NUDDIR(1).LE.0) THEN
         SECMB0(1)=0.
      ENDIF
*     ADMAT0 adresse dans la partie triang. inferieure stricte de A0
      ADMAT0 = 1
*     ADMATX adresse dans la partie triang. inferieure stricte de A
      ADMATX = ADPRCL(1)
*
      DO 100 I = 2,NBLIGN
*
*       Adresse du premier coefficient de la ligne
        ADMATR      = ADMATX
        ADMATX      = ADPRCL(I)
        ADPRC0(I-1) = ADMAT0
*
        IF (NUDDIR(I).LE.0) THEN
*
*          Cas des lignes Dirichlet
*
*          Bi = Aii * Ui au second membre
           SECMB0 (I) = 0.
           IF (NUDDIR(I).LT.0) SECMB0(I) = MATRIS(I)*VALDIR(I)
*
10         IF (ADMATR.NE.0) THEN
              J = NUMCOL(ADMATR)
*              
*             Report de la valeur de blocage au second membre
*             des lignes non Dirichlet; les coefficients sont
*             consideres sur les colonne i : Bj = Bj - Aji * Ui
              IF (NUDDIR(J).GT.0) 
     &        SECMB0(J) = SECMB0(J)-MATRIS(NBLIGN+ADMATR)*VALDIR(I)
              ADMATR = ADSUCL(ADMATR)
              GOTO 10
           ENDIF
*
        ELSE
*
*          Cas des lignes non Dirichlet
*
*          Adresse du premier coefficient de la ligne dans A0
           ADMAT1=ADMAT0
20         IF (ADMATR.NE.0) THEN
             J = NUMCOL(ADMATR)
             IF (NUDDIR(J).LT.0) THEN
*               Cas d'une colonne Dirichlet non homogene :
*               Report de la donnee au second membre et*               elimination du coefficient Bi = Bi - Aij * Uj
                 SECMB0(I) = SECMB0(I)-MATRIS(NBLIGN+ADMATR)*VALDIR(J)    
             ELSEIF (NUDDIR(J).GT.0) THEN
*               Cas des colonnes non Dirichlet : coefficient conserve
                MATRI0(NBLIGN+ADMAT0) = MATRIS(NBLIGN+ADMATR)
                NUMCO0(ADMAT0) = J 
                ADMAT0         = ADMAT0+1 
             ENDIF
             ADMATR = ADSUCL(ADMATR)
             GOTO 20
           ENDIF
*
*          Tri des colonnes de la ligne par ordre croissant
*
           IF (ADMAT0.GT.ADMAT1) 
     &     CALL TRI(ADMAT0-ADMAT1,NUMCO0(ADMAT1),MATRI0(NBLIGN+ADMAT1))
*   
         ENDIF  
100   CONTINUE 
      ADPRC0(NBLIGN) = ADMAT0
      END
