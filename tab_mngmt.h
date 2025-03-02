#ifndef FEM_PROJECT_TAB_MNGMT_H
#define FEM_PROJECT_TAB_MNGMT_H

/*
--------------------------------------------------------------------------------
  Cette fonction alloue de la memoire pour stocker une matrice de dimensions
  dim1 x dim2 de type float. La fonction alloue un tableau de pointeurs (ptr)
  dont chacun des dim1 elements pointe vers une zone de dim2 elements.

  La fonction renvoie NULL en cas d'erreur lors de l'allocation.

  La liberation de la memoire allouee peut etre faite par appel a freetab.
  Remarque :
    Les elements utiles de la matrice sont ranges consecutivement en memoire,
    en dim1 blocs de dim2 elements chacun. Il s'ensuit que cet ensemble peut
    etre transmis a une procedure Fortran via la valeur de ptr[0], comme un seul
    tableau de dimensions (dim2, dim1). On a donc la correspondance d'adressage
    suivante (ptr[i][j] et tab(j,i) designent le meme element) :
      - langage C : ptr[i][j], i=0,dim1-1, j=0,dim2-1,
      - Fortran   : tab(j,i),  i=1,dim1,   j=1,dim2.
    En pratique, etant donnees le sous-programme Fortran :
      subroutine trucF (tab,nbli,nbco)
      integer nbli, nbco
      real tab(nbli, nbco)
    et la fonction C :
      void trucC(float **ptr)
    l'allocation de memoire et les appels se feront de la maniere suivante :
      ptr = matF_alloc(nbco, nbli);
      FORTRANNAME(trucF) (ptr[0], nbli, nbco);
      trucC (ptr);
    La reference aux elements de la matrice se font alors par la notation a
    deux indices indiquee ci-dessus.
--------------------------------------------------------------------------------
*/
float **matF_alloc(int dim1, int dim2);


/*
    Same but with a Integer tab
*/
int **matI_alloc(int dim1, int dim2);

/*
    Once again float but every element is set to 0.0 (calloc)
*/
float **matF_alloc0(int dim1, int dim2);



/*
--------------------------------------------------------------------------------
  Cette fonction libere la memoire allouee par matF_alloc or matI_alloc.
--------------------------------------------------------------------------------
*/
void free_mat(void *ptr);

#endif //FEM_PROJECT_TAB_MNGMT_H
