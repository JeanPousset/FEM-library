*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE rspru(rang,profil,d,l,b,x)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*cette procedure resoud un systeme lineaire triangulaire superieur range
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
      integer iad,iadnxt,ligne,xad,lad
      real    xconnu
*-----------------------------------------------------------------------
      epsil=1.e-10
      do 4 ligne=1,rang
4       x(ligne)=b(ligne)
*
      if (abs(d(rang)).lt.epsil) goto 99991
      x(rang)=x(rang)/d(rang)
      iadnxt=profil(rang)
      do 10 ligne=rang-1,1,-1
        iad=profil(ligne)
        xad=ligne
        xconnu=x(ligne+1)
        do 5 lad=iadnxt-1,iad,-1   
          x(xad)=x(xad)-l(lad)*xconnu
          xad=xad-1
5       continue
        if (abs(d(ligne)).lt.epsil) goto 99991
        x(ligne)=x(ligne)/d(ligne)   
        iadnxt=iad
10    continue
      return
*-----------------------------------------------------------------------
99991 write(*,*) ' LA MATRICE L N''EST PAS INVERSIBLE.'
      stop
      end
