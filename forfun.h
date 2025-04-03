#ifndef FORFUN_H
#define FORFUN_H

/* ---------- Macro to handle name of Fortran subprogramms ----------- */

#ifdef C_NO_TRAILING_UNDERSCORE
     /* case of compilator xlf (IBM, Apple) */
     #define FORTRANNAME(x) x
#else
     /* Case of almost all other compilators */
     #define FORTRANNAME(x) x##_
#endif
/* -------------------------------------------------------------------------- */


/*
           Interface with Fortran procedures of the project
           ================================================
*/

void FORTRANNAME(affsmd)
     (const int *nblign, const int *adprcl, const int *numcol,
      const int *adsucl, const float *matris, const float *secmbr,
      const int *nuddir, const float *valdir);

void FORTRANNAME(affsmo)
     (const int *NBLIGN, const int *ADPRCO, const int *NUMCO0,
      const float *MATRI0, const float *SECMB0);

void FORTRANNAME(affsol)
     (const int *nblign, const float*coor, const float*u, const float*uex,
      const int *impfch);

void FORTRANNAME(assmat)
     (const int *I, const int *J, const float *X, int *ADPRCL, int *NUMCOL,
      int *ADSUCL,  float *LMATRI, int *NEXTAD);

void FORTRANNAME(cdesse)
     (const int *NBLIGN, const int *ADPRCL, const int *NUMCOL,
      const int *ADSUCL, const float *MATRIS, const float *SECMBR,
      const int *NUDDIR, const float *VALDIR,
      int *ADPRC0, int *NUMCO0, float *MATRI0, float *SECMB0);

void FORTRANNAME(impmpr)
     (const int*IMPFCH, const int *RANG, const int *PROFIL, const float *D,
      const float *L);

void FORTRANNAME(ltlpr)
     (const int*rang, const int *profil, const float *ad, const float *al,
      const float *eps, float *ld, float *ll);

void FORTRANNAME(rsprl)
     (const int*rang, const int *profil, const float *d, const float *l,
      const float *b, float *x);

void FORTRANNAME(rspru)
     (const int*rang, const int *profil, const float *d, const float *l,
      const float *b, float *x);

void FORTRANNAME(tri)
     (const int *N, int *NTAB, float *RTAB);
#endif
