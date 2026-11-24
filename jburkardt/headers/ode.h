#ifndef WRAPPER_ODE_H_INCLUDED
#define WRAPPER_ODE_H_INCLUDED

#define ode_intrp(a,b,c,d,e,f,g,h,i) FUNCNAME_ODEINTRP(C_ITPITIT2PIT2DT2PIT(a,b,c,d,e,f,g,h,i))

__MATHSUITE __JBURKARDT void * FUNCNAME_ODEINTRP(void *);
__MATHSUITE __JBURKARDT  void   ode_step
(
  ityp *,
  ityp [],
  void f ( ityp t, ityp y[], ityp yp[] ),
  const register dim_typ neqn,
  ityp *,
  ityp *,
  ityp [static neqn],
  int *,
  ityp *,
  dim_typ *,
  int *,
  int *,
  ityp [static neqn<<4],
  ityp [static neqn],
  ityp [static neqn],
  ityp [static 12],
  ityp [static 12],
  ityp [static 12],
  ityp [static 13],
  ityp [static 12],
  ityp [static 12],
  ityp [static 13],
  int *,
  int *,
  int *
);

#endif
