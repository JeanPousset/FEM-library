************************************************************************
      subroutine AFFSOL (nblign,coor,u,uex,impfch)
      integer nblign
      real coor(2,nblign),u(nblign),uex(nblign)
* Variables locales
      integer i, n1, n2
      real r, eps, nordq,noruq, nordm,norum, errqua, errmax
      logical ouvert
      parameter (eps = 1.e-5)
* -------------------------------------------------------------------------
*  Si impfch est positif, affichage dans le fichier associe a l'unite
*  logique impfch des erreurs relatives (quadratique et maximum) et de la
*  dimension caracteristique du maillage (pas h = distance entre les deux
*  premiers noeuds).
*  Nota : Ces valeurs sont accumulees dans le fichier lors d'appels
*         successifs a cette fonction au cours d'une meme execution,
*         mais aussi lors d'executions successives. Il faut renommer, vider
*         ou supprimer le fichier pour "repartir a zero".
*  Si impfch est negatif ou nul, affichage a l'ecran, pour chaque noeud, de
*  la solution calculee, de la solution exacte et de l'erreur relative sous
*  forme du quotient 2 * |u(i)-uex(i)| / |u(i)+uex(i)| (affichage de la
*  valeur -1 si ce quotient est indefini), suivi de l'affichage des erreurs
*  relatives (quadratique et maximum).
*
*  Arguments d'entree :
* nblign : nombre de lignes de la matrice
* coor   : coordonnees des noeuds
* u      : solution calculee
* uex    : solution exacte
* impfch : numero de l'unite logique du fichier qui regroupe les resultats
*
*  Arguments de sortie : neant
* -------------------------------------------------------------------------
*
*  Calcul de l'erreur quadratique relative et de l'erreur maximum relative
      nordq = 0.
      noruq = 0.
      nordm = 0.
      norum = 0.
      do i=1,nblign
        r = u(i)-uex(i)
        nordq = nordq + r*r
        noruq = noruq + uex(i)*uex(i)
        nordm = max(nordm,abs(r))
        norum = max(norum,abs(uex(i)))
      enddo
      errqua = -1.
      errmax = -1. 
      if (noruq .gt. eps) errqua = sqrt(nordq) / sqrt(noruq)
      if (norum .gt. eps) errmax = nordm / norum
*
*  Impressions
      if (impfch .gt. 0) then
       nordq = coor(1,2)-coor(1,1)
       noruq = coor(2,2)-coor(2,1)
       r = sqrt(nordq*nordq+noruq*noruq)
       inquire (unit=impfch, opened=ouvert)
       if (.not. ouvert) then
         open (unit=impfch,status='unknown',position='append')
       endif
       write (impfch,30) errqua, errmax, r
      else
       write (*,*) 'Affichage des resultats selon le format :'
       write (*,*) ' x y SolApprochee SolExacte ErrRelative'
       write (*,*) 'Nombre de noeuds : ', nblign
   10  write (*,*) 'Donner le premier et le dernier numeros de noeud.',
     & ' (q pour terminer)'
       read (*,*,err=20) n1, n2
       if (n1.lt.1 .or. n2.gt.nblign) goto 10
       do i=n1,n2
         r = abs(u(i)+uex(i))
         if (r .gt. eps) then
          r = 2  * abs(u(i)-uex(i)) / r
         else
          r = -1.
         endif
         write (*,30) coor(1,i), coor(2,i), u(i), uex(i), r
       enddo
       goto 10
   20  continue
*
       write (*,*) 'Erreur quadratique relative : ', errqua
       write (*,*) 'Erreur maximum relative     : ', errmax
      endif
*
   30 format(5(1x,E13.5))
      end
