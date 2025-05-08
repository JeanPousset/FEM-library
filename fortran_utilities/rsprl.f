*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE rsprl(rang,profil,d,l,b,x)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*cette procedure resoud un systeme lineaire triangulaire inferieur range
*sous forme PROFIL 
*
*designation des variables:
*  rang   = nombre de lignes du systeme
*  profil = tableau des pointeurs du profil inferieur
*  d      = diagonale de la matrice L
*  l      = profil inferieur strict de la matrice L
*  b      = vecteur second membre
*  x      = vecteur solution
*
*  epsil  = seuil de singularite de L
*
*parametres formels
      integer rang,profil(*)
      real    d(*),l(*),b(*),x(*),epsil
*
*variables locales
      integer iad,iadnxt,ligne,nbcoef,col,lad
      real    rcoef
*-----------------------------------------------------------------------
      epsil=1.E-10
      if (abs(d(1)).lt.epsil) goto 99991
      x(1)=b(1)/d(1)
      iad=profil(1)
      do 10 ligne=2,rang
        iadnxt=profil(ligne)
        nbcoef=iadnxt-iad
*=======test servant uniquement dans ltlpr lorsque le nombre de coeffs
*=======de la ligne globale de la matrice est > ou = au numero de la 
*=======ligne partielle de la sous matrice
        if (nbcoef.ge.ligne) nbcoef=ligne-1
*=======fin
        col=ligne-nbcoef
        rcoef=0.
        do 5 lad=iadnxt-nbcoef,iadnxt-1
          rcoef=rcoef+l(lad)*x(col)
          col=col+1
5       continue
        if (abs(d(ligne)).lt.epsil) goto 99991
        x(ligne)=( b(ligne)-rcoef )/d(ligne)
        iad=iadnxt
10    continue
      return
*-----------------------------------------------------------------------
99991 write(*,*) ' LA MATRICE L N''EST PAS INVERSIBLE.'
      stop
      end
