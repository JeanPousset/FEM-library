*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE ltlpr(rang,profil,ad,al,eps,ld,ll)
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
* cette procedure effectue la factorisation A=LLt
* d'une matrice A symetrique definie positive
* stockee selon le stockage PROFIL
*
* procedure utilise
* rsprl resolution d'un systeme triangulaire inferieur
*       stockee sous forme PROFIL
*
* designation des variables:
*  rang   = nombre de lignes du systeme
*  profil = tableau des pointeurs du profil inferieur
*  ad     = diagonale de la matrice A
*  ld     = diagonale de la matrice L
*  al     = profil inferieur strict de la matrice A
*  ll     = profil inferieur strict de la matrice L
*  eps    = seuil de singularite de A
*
* parametres formels
      integer rang,profil(*)
      real    ad(*),al(*),ld(*),ll(*),eps
*
* variables locales
      integer iad,iadnxt,ligne,nbcoef,col,llad
      real    rcoef
*-----------------------------------------------------------------------
      if (ad(1).lt.eps) goto 99991
      ld(1)=sqrt(ad(1))
      iadnxt=1
      do 10 ligne=2,rang
        iad=iadnxt
        iadnxt=profil(ligne)
        nbcoef=iadnxt-iad
*=======resolution du sous systeme
        if (nbcoef.gt.0) then
          col=ligne-nbcoef 
          call rsprl(nbcoef,profil(col),ld(col),ll,al(iad),ll(iad))
        endif
*=======fin
        rcoef=ad(ligne)
        do 5 llad=iad,iadnxt-1
5         rcoef=rcoef-ll(llad)**2
        if (rcoef.lt.eps) goto 99991
        ld(ligne)=sqrt(rcoef)
10    continue
      return
*-----------------------------------------------------------------------
99991 write(*,*) ' LA MATRICE A N''EST PAS INVERSIBLE.'
      stop
*
      end
