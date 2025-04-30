************************************************************************
      SUBROUTINE AFFSMO(NBLIGN,ADPRC0,NUMCO0,MATRI0,SECMB0)
************************************************************************
* Imprime la matrice et second membre assembles selon la S.M.O.
*
*** Arguments d'entree ***  
* NBLIGN : nombre de lignes de la matrice
* ADPRC0 : Tableau destine a recevoir les adresses du premier coefficient 
*          non nul de chaque ligne de la matrice
* NUMCO0 : Tableau destine a recevoir les numeros de colonne de la S.M.O.
* MATRI0 : Tableau destine a recevoir les coefficients diagonaux 
*          et les coefficients non logiquement nuls de la partie
*          triangulaire inferieure stricte
* SECMB0 : Second membre assemble de la S.M.O.
*** Arguments de sortie ***
************************************************************************
      INTEGER NBLIGN
* TABLEAUX :
      INTEGER ADPRC0(*),NUMCO0(*)
      REAL MATRI0(*),SECMB0(*)
* VARIABLES LOCALES :
      INTEGER KDEB,KFIN,IAD,IADNXT,J,K
*
*-----------------------------------------------------------------------     
*
      WRITE(*,*)' NOMBRE DE LIGNES :',NBLIGN
1     WRITE(*,*)' NUMEROS DE LA PREMIERE ET DE LA DERNIERE LIGNE ?'
      WRITE(*,*)' TAPER Q POUR TERMINER.'
      READ(*,*,ERR=30)KDEB,KFIN
      IF ((KDEB.LE.0).OR.(KFIN.LT.KDEB).OR.(KFIN.GT.NBLIGN)) GOTO 1
      WRITE(*,2000)
      IAD=1
      IF (KDEB.GT.1) IAD=ADPRC0(KDEB-1)
      DO 20 K=KDEB,KFIN
        IADNXT=ADPRC0(K)
        WRITE(*,1000)K,SECMB0(K),
     &               (MATRI0(J+NBLIGN),J=IAD,IADNXT-1),MATRI0(K)
        WRITE(*,1001)(NUMCO0(J),J=IAD,IADNXT-1),K
        IAD=IADNXT
20    CONTINUE
      GOTO 1
* FORMATS
*
1000  FORMAT(T2,I5,TR1,E12.4,(T20,8E12.4))
1001  FORMAT((T20,8I12))
2000  FORMAT(T2,'LIGNE',TR5,'SECMB0',TR6,
     &          'MATRI0(coefficients'/t31,'+ colonnes)')
30    END
