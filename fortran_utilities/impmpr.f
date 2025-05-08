*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE IMPMPR (IMPFCH,RANG,PROFIL,D,L)
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
* Auteur : D.Martin (Juillet 91)
* Derniere modification (D.Martin) le 4 octobre 1994
*  
* Impression d'une matrice stockee Profil inferieur 
*  
*** Arguments d'entree *** 
* RANG   rang de la matrice
* PROFIL tableau des pointeurs dans la mat. du 1er coeff. non nul des lignes
* D      vecteur contenant la diagonale de la matrice
* L      vecteur contenant la partie triangulaire inferieure de la matrice
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER IMPFCH,RANG,PROFIL(*)
      REAL    D(*),L(*) 
*
      INTEGER IAD,IADNXT,LAD,KOL,LIG
*----------------------------------------------------------------------
      IAD=PROFIL(1)
      WRITE (IMPFCH,10000)
      WRITE (IMPFCH,10001) 1,1,D(1)
      DO 1 LIG=2,RANG
          IADNXT=PROFIL(LIG)
          KOL=LIG-IADNXT+IAD
          WRITE (IMPFCH,10001) LIG,KOL,(L(LAD),LAD=IAD,IADNXT-1),D(LIG)
          IAD=IADNXT
1     CONTINUE
*-----------------------------------------------------------------------
10000 FORMAT(T2,' Numero',T15,'  Numero',T30,'Coefficients'/
     &       T2,'de ligne',T15,'de colonne')
10001 FORMAT(I8,T15,I8,(T28,10E10.3))
                                                                    END !ImpMpr 
