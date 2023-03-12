#ifndef RKF45_H_INCLUDED
#define RKF45_H_INCLUDED

__MATHSUITE __JBURKARDT void   r8_fehl ( void ( ityp, ityp [], ityp b[]), const register dim_typ neqn,
  ityp [static neqn], const register ityp, const register ityp, ityp [static neqn], ityp [static neqn], ityp [static neqn], ityp [static neqn], ityp [static neqn], ityp [static neqn], ityp [static neqn] );
__MATHSUITE __JBURKARDT short   r8_rkf45 ( void ( ityp, ityp [], ityp [] ), int neqn,
  ityp [static neqn], ityp [static neqn], ityp *, const register ityp, ityp *,const register ityp, int );

#endif // RKF45_H_INCLUDED
