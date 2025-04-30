*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ASSMAT (I,J,X,ADPRCL,NUMCOL,ADSUCL,LMATRI,NEXTAD)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Auteur : D.Martin (Mars 1995)
*
* Assemblage de la partie triangulaire inferieure stricte de la matrice :
*
*   -creation de l'element A(i,j) (i > j) avec la valeur X s'il n'existe pas
*
*   -incrementation de l'element A(i,j) de la valeur X sinon
*
*-- Arguments --
* I,J    numeros de ligne et de colonnes
* X      valeur du coefficients A(i,j)
* ADPRCL adresses des 1er coeff. des lignes de la partie triang. inf. stricte
* NUMCOL numeros de colonnes des coeff. de la partie triang. inf. stricte
* ADSUCL adresse du coeff. suivant sur la meme ligne
* LMATRI tableau contenant la partie triangulaire inferieure stricte
* NEXTAD prochaine adresse libre dans la partie triangulaire inferieure stricte
*-- Remarques --
* Le depassement de capacite des tableaux NUMCOL,ADSUCL et LMATRI n'est pas
* teste dans cette procedure (NEXTAD > taille de LMATRI)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER   ADPRCL(*),NUMCOL(*),ADSUCL(*),NEXTAD
      REAL      LMATRI(*),X
*
      INTEGER I,J,IAD,NXT
*-----------------------------------------------------------------------
      IF (I.LE.J) THEN
         WRITE (*,*) '*Bug* La procedure AssMat n''est utilisee que pour
     & l''assemblage de la partie triangulaire stricte. Vous etes sur le
     & coefficient (',I,',',J,').'
         STOP
      ENDIF
*
      IAD=ADPRCL(I-1)
      IF (IAD.GT.0) THEN
*
*        La ligne I a deja ete rencontree : recherche de la colonne J
1        IF (NUMCOL(IAD).NE.J) THEN
*
*           Parcours des elements de la ligne I
            NXT=ADSUCL(IAD)
            IF (NXT.LE.0) GOTO 2
            IAD=NXT
            GOTO 1
         ELSE
*
*           L'element A(i,j) existe : on l'incremente et on s'casse
            LMATRI(IAD)=LMATRI(IAD)+X
            RETURN
         ENDIF
*
*        Creation de la colonne J
2        ADSUCL(IAD)=NEXTAD
      ELSE
*
*        Creation de la ligne I
         ADPRCL(I-1)=NEXTAD
      ENDIF
*
*     Creation du coefficient A(i,j) (nouvelle ligne ou colonne) 
      NUMCOL(NEXTAD)=J
      LMATRI(NEXTAD)=X
      ADSUCL(NEXTAD)=0
      NEXTAD=NEXTAD+1
      END
