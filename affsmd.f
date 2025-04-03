************************************************************************
      subroutine AFFSMD (nblign,adprcl,numcol,adsucl,matris,secmbr,
     &nuddir,valdir)
      integer nblign
      integer adprcl(*),numcol(*),adsucl(*), nuddir(*)
      real matris(*), secmbr(*), valdir(*)
* Variables locales
      integer i,j, Lmin, Lmax
* -------------------------------------------------------------------------
*  Affichage de lignes de la matrice du systeme lineaire stocke selon la
*  S.M.D. (structure morse desordonnee).
*  Les tableaux sont supposes etre de taille suffisante.
*
*  Arguments d'entree :
* nblign : nombre de lignes de la matrice
* adprcl : adresses (dans lmatri) des premiers coefficients des lignes de
*          la partie triangulaire inferieure stricte de la matrice
* numcol : numeros de colonnes des coeff. de la partie triang. inf. stricte
* adsucl : adresse (dans lmatri) du coefficient suivant sur la meme ligne
* matris : tableau contenant la diagonale de la matrice (dmatri) puis la
*          partie triangulaire inferieure stricte (lmatri) sous forme S.M.D.
* secmbr : tableau contenant le second membre
* nuddir : tableau de description des conditions de Dirichlet
* valdir : valeur de blocage aux noeuds Dirichlet non homogenes
*
*  Arguments de sortie :
* -------------------------------------------------------------------------
      write (*,*) '****** Affichage de la S.M.D. ******'
   10 write (*,*) 'Donner, dans l''intervalle [1,',nblign,'], les ',
     &'indices min et max des lignes a afficher. Taper q pour arreter.'
      read (*,*,err=60) Lmin, Lmax
      if (Lmin .gt. Lmax) goto 10
      if (Lmin.lt.1 .or. Lmax.gt.nblign) goto 10
*
      write (*,30)
      do i=Lmin, Lmax
       write (*,40) i, secmbr(i),nuddir(i),valdir(i)
      enddo
*
      write (*,*) '=== MATRICE ==========================='
      if (Lmin.eq.1) then
       write (*,*) '--- Ligne 1 ---'
       write (*,50) 1, matris(1)
       Lmin = 2
      endif
      do i=Lmin, Lmax
        write (*,*) '--- Ligne ', i, ' ---'
        j = adprcl(i-1)
   20   if (j.eq.0) goto 25
        write (*,50) numcol(j), matris(nblign+j)
        j = adsucl(j)
        goto 20
   25   write (*,50) i, matris(i)
      enddo
      goto 10
*
   30 format(/,t2,'Ligne',tr5,'SecMbr',tr4,'Nuddir',tr5,'Valdir')
   40 format(t2,i5,tr1,e12.4,tr3,i5,tr1,e12.4)
   50 format('Col ',i3,' : ',g14.6)
   60 end
