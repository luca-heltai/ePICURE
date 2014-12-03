#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include <math.h>

static double max(double a, double b)
{
  if (a < b)
  {
    return b;
  }
  else 
  {
    return a;
  }
}

static double min(double a, double b)
{
  if (a < b)
  {
    return a;
  }
  else 
  {
    return b;
  }
}

static char _Bas_module__doc__[] = "_Bas module. Version 0.1\n\
\n\
This module implements low level NURBS functions.\n\
\n";

static double **matrix(int nrows, int ncols) 
{
  int row;
  double **mat;

  mat = (double**) malloc(nrows*sizeof(double*));
  for (row = 0; row < nrows; row++)
    mat[row] = (double*) malloc(ncols*sizeof(double));
  return mat;
}

static void freematrix(double **mat, int nrows) {
    int row;
    for (row=0; row<nrows; row++)
        free(mat[row]);
    free(mat);
}

// Compute logarithm of the gamma function
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg214.
static double _gammaln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6] = {76.18009172947146,-86.50532032291677,
                          24.01409824083091,-1.231739572450155,
                          0.12086650973866179e-2, -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

// computes ln(n!)
// Numerical Recipes in C
// Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
static double _factln(int n)
{
  static int ntop = 0;
  static double a[101];
  
  if (n <= 1) return 0.0;
  while (n > ntop)
  {
    ++ntop;
    a[ntop] = _gammaln(ntop+1.0);
  }
  return a[n];
}

static char bincoeff__doc__[] =
"Computes the binomial coefficient.\n\
\n\
 ( n )      n!\n\
 (   ) = --------\n\
 ( k )   k!(n-k)!\n\
\n\
 Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.\n";

static double _bincoeff(int n, int k)
{
  return floor(0.5+exp(_factln(n)-_factln(k)-_factln(n-k)));
}

static PyObject * _Bas_bincoeff(PyObject *self, PyObject *args)
{
    int n, k;
    double ret;
    if(!PyArg_ParseTuple(args, "ii", &n, &k))
        return NULL;
    ret = _bincoeff(n, k);
    return Py_BuildValue("d",ret);
}

// Find the knot span of the parametric point u. 
//
// INPUT:
//
//   n - number of control points - 1
//   p - spline degree       
//   u - parametric point    
//   U - knot sequence
//
// RETURN:
//
//   s - knot span
//
// Algorithm A2.1 from 'The NURBS BOOK' pg68.

static char findspan__doc__[] =
"// Find the knot span of the parametric point u. \n\
//\n\
// INPUT:\n\
//\n\
//   n - number of control points - 1\n\
//   p - spline degree       \n\
//   u - parametric point    \n\
//   U - knot sequence\n\
//\n\
// RETURN:\n\
//\n\
//   s - knot span\n\
//\n\
// Algorithm A2.1 from 'The NURBS BOOK' pg68.\n\
 \n";

static int _findspan(int n, int p, double u, double *U)
{
  int low, high, mid;
   
  //printf ("n=%d, p=%d, u=%g\n", n, p ,u);
  // special case
  if (u == U[n+1]) return(n);
    
  // do binary search
  low = p;
  high = n + 1;
  mid = (low + high) / 2;
  while (u < U[mid] || u >= U[mid+1])
  {
    if (u < U[mid])
      high = mid;
    else
      low = mid;
    mid = (low + high) / 2; // 
  }  

  return(mid);
}

static PyObject * _Bas_findspan(PyObject *self, PyObject *args)
{
    int n, p, ret;
    double u;
    PyObject *input_U;
    PyArrayObject  *U;
    if(!PyArg_ParseTuple(args, "iidO", &n, &p, &u, &input_U))
        return NULL;

    U = (PyArrayObject *) PyArray_ContiguousFromObject(input_U, NPY_DOUBLE, 1, 1);
    if(U == NULL)
        return NULL;

    ret = _findspan(n, p, u, (double *)PyArray_DATA(U));
    return Py_BuildValue("i", ret);
}

// Basis Function. 
//
// INPUT:
//
//   i - knot span  ( from FindSpan() )
//   p - spline degree
//   u - parametric point
//   U - knot sequence
//
// OUTPUT:
//
//   N - Basis functions vector[p+1]
//
// Algorithm A2.2 from 'The NURBS BOOK' pg70.
//
static char basisfuns__doc__[] =
"Compute basis functions which are not zero for given u.\n\
\n\
 INPUT:\n\
   i - knot span  ( from FindSpan() )\n\
   p - spline degree\n\
   u - parametric point\n\
   U - knot sequence\n\
\n\
 OUTPUT:\n\
\n\
   N - Basis functions vector[p+1]\n\
\n\
 Algorithm A2.2 from 'The NURBS BOOK' pg70.\n\
 \n";

static void _basisfuns(int i, double u, int p, double *U, double *N)
{
  int j,r;
  double saved, temp;

  // work space
  double *left  = (double*) malloc((p+1)*sizeof(double));
  double *right = (double*) malloc((p+1)*sizeof(double));
  
  // printf("i=%d\n",i);
  N[0] = 1.0;
  for (j = 1; j <= p; j++)
  {
    left[j]  = u - U[i+1-j];
    right[j] = U[i+j] - u;
    saved = 0.0;
    
    for (r = 0; r < j; r++)
    {
      temp = N[r] / (right[r+1] + left[j-r]);
      N[r] = saved + right[r+1] * temp;
      saved = left[j-r] * temp;
    } 

    N[j] = saved;
  }
  
//   printf("-- in N[0]=%g\n",N[0]);
//   printf("-- in N[1]=%g\n",N[1]);
//   printf("-- in N[2]=%g\n",N[2]);
  free(left);
  free(right);
}

static PyObject * _Bas_basisfuns(PyObject *self, PyObject *args)
{
    int i, p;
    double u;
    npy_intp dim[2];
    PyObject *input_U;
    PyArrayObject *N, *U;

    if(!PyArg_ParseTuple(args, "iidO", &i, &p, &u, &input_U))
        return NULL;

    U = (PyArrayObject *) PyArray_ContiguousFromObject(input_U, NPY_DOUBLE, 1, 1);
    if(U == NULL)
        return NULL;

    dim[0] = p+1;
    dim[1] = 0; 
    //  N array is p+1 in size
    N = (PyArrayObject *) PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    _basisfuns(i, u, p, (double *)PyArray_DATA(U), (double *)PyArray_DATA(N));
    // printf("*N->data[0]=%g \n", (double *)N->data[0]);
    // printf("*N->data[1]=%g \n", (double *)N->data[1]);
    //return (PyObject *)N;
    // return Py_BuildValue("O", (PyObject *)N );
    return PyArray_Return(N);
}

static char bspeval__doc__[] =
"Evaluation of univariate B-Spline. \n\
\n\
INPUT:\n\
\n\
 d - spline degree       integer\n\
 c - control points      double  matrix(mc,nc)\n\
 k - knot sequence       double  vector(nk)\n\
 u - parametric points   double  vector(nu)\n\
\n\
OUTPUT:\n\
\n\
   p - evaluated points    double  matrix(mc,nu)\n\
\n\
Modified version of Algorithm A3.1 from 'The NURBS BOOK' pg82.\n\
\n";

static void _bspeval(int d, double *ctrl, int mc, int nc, double *k, int nk, double *u,
            int nu, double *pnt)
{
  int i, s, tmp1, row, col;
  double tmp2;

  // space for the basis functions
  double *N = (double*) malloc((d+1)*sizeof(double));

  // for each parametric point u[col]
  for (col = 0; col < nu; col++)
  {
    // find the span of u[col]
    s = _findspan(nc-1, d, u[col], k);
    _basisfuns(s, u[col], d, k, N);
    
    tmp1 = s - d;
    for (row = 0; row < mc; row++)
    {
      tmp2 = 0.0;   
      for (i = 0; i <= d; i++)
        tmp2 += N[i] * ctrl[row*nc+tmp1+i];
      pnt[row*nu+col] = tmp2;
    }
  }
  free(N);
}

static PyObject * _Bas_bspeval(PyObject *self, PyObject *args)
{
    int d;
    npy_intp dim[2], mc, nc, nu;
    double *ctrldat, *pntdat, *kdat, *udat;
    PyObject *input_ctrl, *input_k, *input_u;
    PyArrayObject *ctrl, *k, *u, *pnt;

    if(!PyArg_ParseTuple(args, "iOOO", &d, &input_ctrl, &input_k, &input_u))
        return NULL;

    ctrl = (PyArrayObject *) PyArray_ContiguousFromAny(input_ctrl, NPY_DOUBLE, 2, 2);
    if(ctrl == NULL)
        return NULL;
    k = (PyArrayObject *) PyArray_ContiguousFromAny(input_k, NPY_DOUBLE, 1, 1);
    if(k == NULL)
        return NULL;
    u = (PyArrayObject *) PyArray_ContiguousFromAny(input_u, NPY_DOUBLE, 1, 1);
    if(u == NULL)
        return NULL;

    nu = PyArray_DIM(u, 0);
    mc = PyArray_DIM(ctrl, 0);
    nc = PyArray_DIM(ctrl, 1);
    dim[0] = mc;
    dim[1] = nu;
    pnt = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);

    ctrldat = (double *)PyArray_DATA(ctrl);
    pntdat  = (double *)PyArray_DATA(pnt);
    kdat    = (double *)PyArray_DATA(k);
    udat    = (double *)PyArray_DATA(u);
    _bspeval(d, ctrldat, mc, nc, kdat, PyArray_DIM(k, 0), udat, nu, pntdat);

    Py_DECREF(ctrl);
    Py_DECREF(k);
    Py_DECREF(u);
    return PyArray_Return(pnt);
}

static char dersbasisfuns__doc__[] =
" Compute Non-zero basis functions and their derivatives.\n\n\
 INPUT:\n\n\
   d  - spline degree          integer\n\
   k  - knot sequence          double  vector(nk)\n\
   u  - parametric point       double\n\
   s  - knot span              integer\n\
   n  - number of derivatives  integer\n\n\
 OUTPUT:\n\n\
   dN - Basis functions        double  matrix(n+1,d+1)\n\
        and derivatives up\n\
        to the nth derivative\n\
        (n < d)\n\n\
 Algorithm A2.3 from 'The NURBS BOOK' pg72.\n";

static void _dersbasisfuns(int d, double *k, int nk, double u, int s, int n, double *ders)
{
    // ders dimensions: (n+1, d+1)

    int i,j,r,s1,s2,rk,pk,j1,j2,dp1;
    double temp, saved, der;
    double **ndu, **a, *left, *right;

    dp1 = d+1;
    ndu = matrix(d+1, d+1);
    a = matrix(2, d+1);
    left = (double *) malloc((d+1)*sizeof(double));
    right = (double *) malloc((d+1)*sizeof(double));

    ndu[0][0] = 1.0;

    for (j=1; j<=d; j++) {
        left[j] = u - k[s+1-j];
        right[j] = k[s+j]-u;
        saved = 0.0;
        for (r=0; r<j; r++) {
            ndu[j][r] = right[r+1] + left[j-r];
            temp = ndu[r][j-1]/ndu[j][r];

            ndu[r][j] = saved + right[r+1]*temp;
            saved = left[j-r]*temp;
        }
        ndu[j][j] = saved;
    }

    for (j=0; j<=d; j++)
        ders[0*dp1+j] = ndu[j][d];

    for (r=0; r<=d; r++) {
        s1 = 0; s2 = 1;
        a[0][0] = 1.0;

        for (i = 1; i <= n; i++) {
            der = 0.0;
            rk = r-i;  pk = d-i;

            if (r >= i) {
                a[s2][0] = a[s1][0] / ndu[pk+1][rk];
                der = a[s2][0] * ndu[rk][pk];
            }  
            if (rk >= -1)
                j1 = 1;
            else
                j1 = -rk;
            if (r-1 <= pk)
                j2 = i-1;
            else
                j2 = d-r;

            for (j = j1; j <= j2; j++) {
                a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
                der += a[s2][j] * ndu[rk+j][pk];
            }
            if (r <= pk) {
                a[s2][i] = -a[s1][i-1] / ndu[pk+1][r];
                der += a[s2][i] * ndu[r][pk];
            }
            ders[i*dp1+r] = der;
            j = s1; s1 = s2; s2 = j;
        }        
    }
  
    r = d;
    for (i = 1; i <= n; i++) {
        for (j = 0; j <= d; j++)
            ders[i*dp1+j] *= r;
        r *= d-i;
    }

    freematrix(ndu, d+1);
    freematrix(a, 2);
    free(left);
    free(right);
}

static PyObject * _Bas_dersbasisfuns(PyObject *self, PyObject *args)
{
    int d, s, nk, n;
    double u;
    double *udat, *dNdat;
    npy_intp dim[2];
    PyObject *input_U;
    PyArrayObject *dN, *U;

    if (!PyArg_ParseTuple(args, "iOdii", &d, &input_U, &u, &s, &n))
        return NULL;

    U = (PyArrayObject *) PyArray_ContiguousFromAny(input_U, NPY_DOUBLE, 1, 1);
    if (U == NULL)
        return NULL;
    udat = (double *)PyArray_DATA(U);

    nk = PyArray_DIM(U, 0);
    dim[0] = n+1; 
    dim[1] = d+1;
    dN = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    dNdat = (double *)PyArray_DATA(dN);

    _dersbasisfuns(d, udat, nk, u, s, n, dNdat);
    return PyArray_Return(dN);
}

static char bspdeval__doc__[] =
"Evaluate a B-Spline derivative curve.\n\
\n\
INPUT:\n\
\n\
 d - spline degree       integer\n\
 c - control points      double  matrix(mc,nc)\n\
 k - knot sequence       double  vector(nk)\n\
 u - parametric point    double\n\
 n - nth derivative      integer\n\
\n\
OUTPUT:\n\
\n\
 p - evaluated points    double  matrix(mc, n+1)\n\
\n\
Modified version of Algorithm A3.2 from 'The NURBS BOOK' pg93.\n\
\n";

static void _bspdeval(int d, double *c, int mc, int nc, double *k, int nk, 
             double u, int n, double *p)
{
    int row, col, j, s, dp1, np1;
    int du = min(d,n);
    double *dNdat;
    PyArrayObject *dN;
    npy_intp dim[2];

    np1 = n+1;
    dp1 = d+1;
    dim[0] = np1;
    dim[1] = dp1;

    dN = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    dNdat = (double *)PyArray_DATA(dN);

    for (col = d+1; col < np1; col++)
        for (row = 0; row < mc; row++)
            p[row*np1 + col] = 0.0;

    s = _findspan(nc-1, d, u, k);
    _dersbasisfuns(d, k, nk, u, s, n, dNdat);

    for (col = 0; col <= du; col++) {
        for (row = 0; row < mc; row++) {
            p[row*np1 + col] = 0.0;
            for (j = 0; j <= d; j++)
                p[row*np1 + col] += dNdat[col*dp1+j] * c[row*nc+s-d+j];
        }
    }
    PyArray_free(dN);
}

static PyObject * _Bas_bspdeval(PyObject *self, PyObject *args)
{
    int d, n;
    npy_intp dim[2], mc, nc;
    double u, *ctrldat, *pntdat, *kdat;
    PyObject *input_ctrl, *input_k;
    PyArrayObject *ctrl, *k, *pnt;
    if(!PyArg_ParseTuple(args, "iOOdi", &d, &input_ctrl, &input_k, &u, &n))
        return NULL;
    ctrl = (PyArrayObject *) PyArray_ContiguousFromAny(input_ctrl, NPY_DOUBLE, 2, 2);
    if(ctrl == NULL)
        return NULL;
    k = (PyArrayObject *) PyArray_ContiguousFromAny(input_k, NPY_DOUBLE, 1, 1);
    if(k == NULL)
        return NULL;
    mc = PyArray_DIM(ctrl, 0);
    nc = PyArray_DIM(ctrl, 1);
    dim[0] = mc;
    dim[1] = n + 1;
    pnt = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    ctrldat = (double *)PyArray_DATA(ctrl);
    pntdat  = (double *)PyArray_DATA(pnt);
    kdat    = (double *)PyArray_DATA(k);
    _bspdeval(d, ctrldat, mc, nc, kdat, PyArray_DIM(k, 0), u, n, pntdat);
    Py_DECREF(ctrl);
    Py_DECREF(k);
    return PyArray_Return(pnt);
}

static char bspkntins__doc__[] =
"Insert Knot into a B-Spline.\n\
\n\
INPUT:\n\
\n\
 d - spline degree       integer\n\
 c - control points      double  matrix(mc,nc)\n\
 k - knot sequence       double  vector(nk)\n\
 u - new knots           double  vector(nu)\n\
\n\
OUTPUT:\n\
\n\
 ic - new control points double  matrix(mc,nc+nu)\n\
 ik - new knot sequence  double  vector(nk+nu)\n\
\n\
Modified version of Algorithm A5.4 from 'The NURBS BOOK' pg164.\n\
\n";

static void _bspkntins(int d, double *ctrl, int mc, int nc, double *k, int nk, 
              double *u, int nu, double *ictrl, double *ik)
{
    // ctrl dimensions:  (mc, nc)
    // ictrl dimensions: (mc, nci)

    int a, b, r, l, i, j, m, n, s, q, ind;
    const int nci = nc+nu;
    double alfa;

    n = nc - 1;
    r = nu - 1;

    m = n + d + 1;
    a = _findspan(n, d, u[0], k);
    b = _findspan(n, d, u[r], k);
    ++b;

    for (q = 0; q < mc; q++) {
        for (j = 0; j <= a-d; j++) ictrl[q*nci+j] = ctrl[q*nc+j];
        for (j = b-1; j <= n; j++) ictrl[q*nci+j+r+1] = ctrl[q*nc+j];
    }
    for (j = 0; j <= a; j++)   ik[j] = k[j];
    for (j = b+d; j <= m; j++) ik[j+r+1] = k[j];

    i = b + d - 1;
    s = b + d + r;
    for (j = r; j >= 0; j--) {
        while (u[j] <= k[i] && i > a) {
            for (q = 0; q < mc; q++)
                ictrl[q*nci+s-d-1] = ctrl[q*nc+i-d-1];
            ik[s] = k[i];
            --s;
            --i;
        }
        for (q = 0; q < mc; q++)
           ictrl[q*nci+s-d-1] = ictrl[q*nci+s-d];
        for (l = 1; l <= d; l++) {
            ind = s - d + l;
            alfa = ik[s+l] - u[j];
            if (fabs(alfa) == 0.0)
                for (q = 0; q < mc; q++)
                    ictrl[q*nci+ind-1] = ictrl[q*nci+ind];
            else {
                alfa /= (ik[s+l] - k[i-d+l]);
                for (q = 0; q < mc; q++)
                    ictrl[q*nci+ind-1] = alfa*ictrl[q*nci+ind-1] + (1.0-alfa)*ictrl[q*nci+ind];
            }
        }

        ik[s] = u[j];
        --s;
    }
}

static PyObject * _Bas_bspkntins(PyObject *self, PyObject *args)
{
    int d;
    npy_intp mc, nc, nk, nu, dim[2];
    //double **ctrlmat, **icmat;
    double *ctrldat, *icdat, *kdat, *ikdat, *udat;
    PyObject *input_ctrl, *input_k, *input_u;
    PyArrayObject *ctrl, *k, *u, *ic, *ik;
    if(!PyArg_ParseTuple(args, "iOOO", &d, &input_ctrl, &input_k, &input_u))
        return NULL;
    ctrl = (PyArrayObject *) PyArray_ContiguousFromAny(input_ctrl, NPY_DOUBLE, 2, 2);
    if(ctrl == NULL)
        return NULL;
    k = (PyArrayObject *) PyArray_ContiguousFromAny(input_k, NPY_DOUBLE, 1, 1);
    if(k == NULL)
        return NULL;
    u = (PyArrayObject *) PyArray_ContiguousFromAny(input_u, NPY_DOUBLE, 1, 1);
    if(u == NULL)
        return NULL;
    mc = PyArray_DIM(ctrl, 0);
    nc = PyArray_DIM(ctrl, 1);
    nk = PyArray_DIM(k, 0);
    nu = PyArray_DIM(u, 0);
    dim[0] = mc;
    dim[1] = nc + nu;
    ic = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    dim[0] = nk + nu;
    ik = (PyArrayObject *) PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    //ctrlmat = vec2mat(ctrl->data, mc, nc);
    //icmat = vec2mat(ic->data, mc, nc + nu);
    ctrldat = (double *)PyArray_DATA(ctrl);
    icdat   = (double *)PyArray_DATA(ic);
    kdat    = (double *)PyArray_DATA(k);
    ikdat   = (double *)PyArray_DATA(ik);
    udat    = (double *)PyArray_DATA(u);
    _bspkntins(d, ctrldat, mc, nc, kdat, nk, udat, nu, icdat, ikdat);
    //free(icmat);
    //free(ctrlmat);
    Py_DECREF(ctrl);
    Py_DECREF(k);
    Py_DECREF(u);
    return Py_BuildValue("(OO)", (PyObject *)ic, (PyObject *)ik);
}

static char bspdegelev__doc__[] =
"Degree elevate a B-Spline t times.\n\
\n\
INPUT:\n\
\n\
 n,p,U,Pw,t\n\
\n\
OUTPUT:\n\
\n\
 nh,Uh,Qw\n\
\n\
Modified version of Algorithm A5.9 from 'The NURBS BOOK' pg206.\n\
\n";

static void _bspdegelev(int d, double *ctrl, int mc, int nc, double *k, int nk, 
               int t, int *nh, double *ictrl, double *ik)
{
    // ctrl dimensions:  (mc, nc)
    // ictrl dimensions: (mc, nci)
  int i, j, q, s, m, ph, ph2, mpi, mh, r, a, b, cind, oldr, mul;
  int n, lbz, rbz, save, tr, kj, first, kind, last, bet, ii;
  const int nci = nc*(t+1);
  double inv, ua, ub, numer, den, alf, gam;
  double **bezalfs, **bpts, **ebpts, **Nextbpts, *alfs; 

  n = nc - 1;

  bezalfs = matrix(d+1,d+t+1);
  bpts = matrix(mc,d+1);
  ebpts = matrix(mc,d+t+1);
  Nextbpts = matrix(mc,d);
  alfs = (double *) malloc(d*sizeof(double));

  m = n + d + 1;
  ph = d + t;
  ph2 = ph / 2;

  // compute bezier degree elevation coefficeients  
  bezalfs[0][0] = bezalfs[d][ph] = 1.0;

  for (i = 1; i <= ph2; i++)
  {
    inv = 1.0 / _bincoeff(ph,i);
    mpi = min(d,i);
    
    for (j = max(0,i-t); j <= mpi; j++)
      bezalfs[j][i] = inv * _bincoeff(d,j) * _bincoeff(t,i-j);
  }    
  
  for (i = ph2+1; i <= ph-1; i++)
  {
    mpi = min(d, i);
    for (j = max(0,i-t); j <= mpi; j++)
      bezalfs[j][i] = bezalfs[d-j][ph-i];
  }       

  mh = ph;
  kind = ph+1;
  r = -1;
  a = d;
  b = d+1;
  cind = 1;
  ua = k[0];
  for (ii = 0; ii < mc; ii++)
    ictrl[ii*nci+0] = ctrl[ii*nc+0];
  
  for (i = 0; i <= ph; i++)
    ik[i] = ua;
    
  // initialise first bezier seg
  for (i = 0; i <= d; i++)
    for (ii = 0; ii < mc; ii++)
      bpts[ii][i] = ctrl[ii*nc+i];

  // big loop thru knot vector
  while (b < m)
  {
    i = b;
    while (b < m && k[b] == k[b+1])
      b++;

    mul = b - i + 1;
    mh += mul + t;
    ub = k[b];
    oldr = r;
    r = d - mul;
    
    // insert knot u(b) r times
    if (oldr > 0)
      lbz = (oldr+2) / 2;
    else
      lbz = 1;

    if (r > 0)
      rbz = ph - (r+1)/2;
    else
      rbz = ph;  

    if (r > 0)
    {
      // insert knot to get bezier segment
      numer = ub - ua;
      for (q = d; q > mul; q--)
        alfs[q-mul-1] = numer / (k[a+q]-ua);
      for (j = 1; j <= r; j++)  
      {
        save = r - j;
        s = mul + j;            

        for (q = d; q >= s; q--)
          for (ii = 0; ii < mc; ii++)
            bpts[ii][q] = alfs[q-s]*bpts[ii][q]+(1.0-alfs[q-s])*bpts[ii][q-1];

        for (ii = 0; ii < mc; ii++)
          Nextbpts[ii][save] = bpts[ii][d];
      }  
    }
    // end of insert knot

    // degree elevate bezier
    for (i = lbz; i <= ph; i++)
    {
      for (ii = 0; ii < mc; ii++)
        ebpts[ii][i] = 0.0;
      mpi = min(d, i);
      for (j = max(0,i-t); j <= mpi; j++)
        for (ii = 0; ii < mc; ii++)
          ebpts[ii][i] = ebpts[ii][i] + bezalfs[j][i]*bpts[ii][j];
    }
    // end of degree elevating bezier

    if (oldr > 1)
    {
      // must remove knot u=k[a] oldr times
      first = kind - 2;
      last = kind;
      den = ub - ua;
      bet = (ub-ik[kind-1]) / den;
      
      // knot removal loop
      for (tr = 1; tr < oldr; tr++)
      {        
        i = first;
        j = last;
        kj = j - kind + 1;
        while (j - i > tr)
        {
          // loop and compute the new control points
          // for one removal step
          if (i < cind)
          {
            alf = (ub-ik[i])/(ua-ik[i]);
            for (ii = 0; ii < mc; ii++)
              ictrl[ii*nci+i] = alf*ictrl[ii*nci+i] + (1.0-alf)*ictrl[ii*nci+i-1];
          }
          if (j >= lbz)
          {
            if (j-tr <= kind-ph+oldr)
            {  
              gam = (ub-ik[j-tr]) / den;
              for (ii = 0; ii < mc; ii++)
                ebpts[ii][kj] = gam*ebpts[ii][kj] + (1.0-gam)*ebpts[ii][kj+1];
            }
            else
            {
              for (ii = 0; ii < mc; ii++)
                ebpts[ii][kj] = bet*ebpts[ii][kj] + (1.0-bet)*ebpts[ii][kj+1];
            }
          }
          i++;
          j--;
          kj--;
        }      
        
        first--;
        last++;
      }                    
    }
    // end of removing knot n=k[a]
                  
    // load the knot ua
    if (a != d)
      for (i = 0; i < ph-oldr; i++)
      {
        ik[kind] = ua;
        kind++;
      }

    // load ctrl pts into ic
    for (j = lbz; j <= rbz; j++)
    {
      for (ii = 0; ii < mc; ii++)
        ictrl[ii*nci+cind] = ebpts[ii][j];
      cind++;
    }
    
    if (b < m)
    {
      // setup for next pass thru loop
      for (j = 0; j < r; j++)
        for (ii = 0; ii < mc; ii++)
          bpts[ii][j] = Nextbpts[ii][j];
      for (j = r; j <= d; j++)
        for (ii = 0; ii < mc; ii++)
          bpts[ii][j] = ctrl[ii*nc+b-d+j];
      a = b;
      b++;
      ua = ub;
    }
    else
      // end knot
      for (i = 0; i <= ph; i++)
        ik[kind+i] = ub;
  }                  
  // end while loop   
  
  *nh = mh - ph - 1;

  freematrix(bezalfs, d+1);
  freematrix(bpts, mc);
  freematrix(ebpts, mc);
  freematrix(Nextbpts, mc);
  free(alfs);
}

static PyObject * _Bas_bspdegelev(PyObject *self, PyObject *args)
{
    int d, t, nh;
    npy_intp mc, nc, nk, dim[2];
    double *ctrldat, *icdat, *kdat, *ikdat;
    PyObject *input_ctrl, *input_k;
    PyArrayObject *ctrl, *k, *ic, *ik;
    if(!PyArg_ParseTuple(args, "iOOi", &d, &input_ctrl, &input_k, &t))
        return NULL;
    ctrl = (PyArrayObject *) PyArray_ContiguousFromAny(input_ctrl, NPY_DOUBLE, 2, 2);
    if(ctrl == NULL)
        return NULL;
    k = (PyArrayObject *) PyArray_ContiguousFromAny(input_k, NPY_DOUBLE, 1, 1);
    if(k == NULL)
        return NULL;
    mc = PyArray_DIM(ctrl, 0);
    nc = PyArray_DIM(ctrl, 1);
    nk = PyArray_DIM(k, 0);
    dim[0] = mc;
    dim[1] = nc*(t + 1);
    ic = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    dim[0] = (t + 1)*nk;
    ik = (PyArrayObject *) PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    ctrldat = (double *)PyArray_DATA(ctrl);
    icdat   = (double *)PyArray_DATA(ic);
    kdat    = (double *)PyArray_DATA(k);
    ikdat   = (double *)PyArray_DATA(ik);
    _bspdegelev(d, ctrldat, mc, nc, kdat, nk, t, &nh, icdat, ikdat);
    Py_DECREF(ctrl);
    Py_DECREF(k);
    return Py_BuildValue("(OOi)", (PyObject *)ic, (PyObject *)ik, nh);
}

static char bspbezdecom__doc__[] =
"Decompose a B-Spline to Bezier segments.\n\
\n\
INPUT:\n\
\n\
 n,p,U,Pw\n\
\n\
OUTPUT:\n\
\n\
 Qw\n\
\n\
Modified version of Algorithm A5.6 from 'The NURBS BOOK' pg173.\n\
\n";

static void _bspbezdecom(int d, double *ctrl, int mc, int nc, double *k, int nk, 
               double *ictrl, int nci)
{
    // ctrl dimensions:  (mc, nc)
    // ictrl dimensions: (mc, nci)
  int i, j, s, m, r, a, b, mul, n, nb, ii, save, q;
  double ua, ub, numer;
  double *alfs; 

  n = nc - 1;

  alfs = (double *) malloc(d*sizeof(double));

  m = n + d + 1;
  a = d;
  b = d+1;
  ua = k[0];
  nb = 0;
  
  // initialise first bezier seg
  for (i = 0; i <= d; i++)
    for (ii = 0; ii < mc; ii++)
      ictrl[ii*nci+i] = ctrl[ii*nc+i];  

  // big loop thru knot vector
  while (b < m)
  {
    i = b;
    while (b < m && k[b] == k[b+1])
      b++;

    mul = b - i + 1;
    ub = k[b];
    r = d - mul;
    
    // insert knot u(b) r times
    if (r > 0)
    {
      // insert knot to get bezier segment
      numer = ub - ua;
      for (q = d; q > mul; q--)
        alfs[q-mul-1] = numer / (k[a+q]-ua);
      for (j = 1; j <= r; j++)  
      {
        save = r - j;
        s = mul + j;            

        for (q = d; q >= s; q--)
          for (ii = 0; ii < mc; ii++)
            ictrl[ii*nci+q+nb] = alfs[q-s]*ictrl[ii*nci+q+nb]+(1.0-alfs[q-s])*ictrl[ii*nci+q-1+nb];

        for (ii = 0; ii < mc; ii++)
          ictrl[ii*nci+save+nb+d+1] = ictrl[ii*nci+d]; 
      }  
    }
    // end of insert knot
    nb += d;
    if (b < m)
    {
      // setup for next pass thru loop
      for (j = r; j <= d; j++)
        for (ii = 0; ii < mc; ii++)
          ictrl[ii*nci+j+nb] = ctrl[ii*nc+b-d+j];
      a = b;
      b++;
      ua = ub;
    }
  }                 
  // end while loop   
  
  free(alfs);
}

static PyObject * _Bas_bspbezdecom(PyObject *self, PyObject *args)
{
    int i, b, c, d, m;
    npy_intp mc, nc, nk, dim[2];
    double *ctrldat, *icdat, *kdat;
    PyObject *input_ctrl, *input_k;
    PyArrayObject *ctrl, *k, *ic;
    if(!PyArg_ParseTuple(args, "iOO", &d, &input_ctrl, &input_k))
        return NULL;
    ctrl = (PyArrayObject *) PyArray_ContiguousFromAny(input_ctrl, NPY_DOUBLE, 2, 2);
    if(ctrl == NULL)
        return NULL;
    k = (PyArrayObject *) PyArray_ContiguousFromAny(input_k, NPY_DOUBLE, 1, 1);
    if(k == NULL)
        return NULL;
    mc = PyArray_DIM(ctrl, 0);
    nc = PyArray_DIM(ctrl, 1);
    nk = PyArray_DIM(k, 0);
    
    kdat    = (double *)PyArray_DATA(k);

    i = d + 1;
    c = 0;
    m = nk - d - 1;
    while (i < m) {
        b = 1;
        while (i < m && kdat[i] == kdat[i+1]) {
            b++;
            i++;
        }
        if(b < d) 
            c = c + (d - b); 
        i++;
    }
    dim[0] = mc;
    dim[1] = nc+c;
    ic = (PyArrayObject *) PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    ctrldat = (double *)PyArray_DATA(ctrl);
    icdat   = (double *)PyArray_DATA(ic);
    _bspbezdecom(d, ctrldat, mc, nc, kdat, nk, icdat, dim[1]);
    Py_DECREF(ctrl);
    Py_DECREF(k); 
    return Py_BuildValue("O", ic);
}

static PyMethodDef _Bas_methods[] =
{
    {"bincoeff", _Bas_bincoeff, METH_VARARGS, bincoeff__doc__},
    {"bspeval", _Bas_bspeval, METH_VARARGS, bspeval__doc__},
    {"bspdeval", _Bas_bspdeval, METH_VARARGS, bspdeval__doc__},
    {"bspkntins", _Bas_bspkntins, METH_VARARGS, bspkntins__doc__},
    {"bspdegelev", _Bas_bspdegelev, METH_VARARGS, bspdegelev__doc__},
    {"bspbezdecom", _Bas_bspbezdecom, METH_VARARGS, bspbezdecom__doc__},
    // se:
    {"basisfuns", _Bas_basisfuns, METH_VARARGS, basisfuns__doc__},
    {"dersbasisfuns", _Bas_dersbasisfuns, METH_VARARGS, dersbasisfuns__doc__},
    {"findspan", _Bas_findspan, METH_VARARGS, findspan__doc__},
    {NULL, NULL}
};

PyMODINIT_FUNC init_Bas(void)
{
    PyObject *m;
    m = Py_InitModule3("_Bas", _Bas_methods, _Bas_module__doc__);
    import_array();
}
    
