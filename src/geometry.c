// geometry.c 01/09/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.xm
*/

#include "dutils.h" // DA RENDERE VISIBILE SIA AL COMPILATORE CHE AL LINKER

#define PREC (5e-16f)
#define ISZERO(x) (fabs(x)<PREC)

__MSUTIL_ inline void toupper_s(char *string)
{
    const size_t len = strlen(string);
    #pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = toupper(string[i]);
    return;
}

__MSUTIL_ inline void tolower_s(char *string)
{
    const size_t len = strlen(string);
    #pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = tolower(string[i]);
    return;
}

__MSNATIVE_ void strundsc(const char *string, char name[])
{
    char *fname;
    char str[MIN_STRING];
    strcpy(str, string);
    tolower_s(str);
    fname = replace(str, BLANK_STRING, FILENAMES_SEPERATOR);
    strcpy(name, fname);
    free(fname);
    return;
}

__MSNATIVE_ void strboolize(const char *string, char name[])
{
	char str[INFO_STRING] = NULL_CHAR;
	strcpy(str, string);
	char *fname = strtok(str, BLANK_STRING);
	fname[0] = tolower(fname[0]);
	// strcpy(name, fname);
	dim_typ i=0;
	do
	{
		if(i++)
			fname[0] = toupper(fname[0]);
		strcat(name, fname);
	}
	while((fname=strtok(NULL,BLANK_STRING)));
	fname = NULL;
	return;
}

__MSNATIVE_ inline void strfnm(const char *string, char file_name[static MAX_PATH_LENGTH])
{
    strundsc(string, file_name);
    strcat(file_name, "."DEFAULT_HELP_FILE_EXTENSION);
    return;
}

/*  Bit counter by Ratko Tomic */
__MSUTIL_ inline int __export countbits(long i)
{
      i = ((i & 0xAAAAAAAAL) >>  1) + (i & 0x55555555L);
      i = ((i & 0xCCCCCCCCL) >>  2) + (i & 0x33333333L);
      i = ((i & 0xF0F0F0F0L) >>  4) + (i & 0x0F0F0F0FL);
      i = ((i & 0xFF00FF00L) >>  8) + (i & 0x00FF00FFL);
      i = ((i & 0xFFFF0000L) >> 16) + (i & 0x0000FFFFL);
      return (int)i;
}

__MSUTIL_ inline int __export ucountbits(unsigned long num)
{
    int count;
    for(count = 0; num; ++count, num &= num - 1);
    return count;
}

__MSUTIL_ char __export *strrev(char *str)
{
    char *p1, *p2;

    if (! str || ! *str)
        return str;

    for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
        *p1 ^= *p2 ^= *p1 ^= *p2;

    return str;
}

__MSUTIL_ char __export *replace(char const * const original, char const * const pattern, char const * const replacement)
{
  size_t const replen = strlen(replacement);
  size_t const patlen = strlen(pattern);
  size_t const orilen = strlen(original);

  size_t patcnt = 0;
  const char * oriptr;
  const char * patloc;

  // find how many times the pattern occurs in the original string
  for (oriptr = original; (patloc = strstr(oriptr, pattern)); oriptr = patloc + patlen)
    ++patcnt;

  {
    // allocate memory for the new string
    size_t const retlen = orilen + patcnt * (replen - patlen);
    char * const returned = (char *) malloc( sizeof(char) * (retlen + 1) );

    if (returned != NULL)
    {
      // copy the original string,
      // replacing all the instances of the pattern
      char * retptr = returned;
      for (oriptr = original; (patloc = strstr(oriptr, pattern)); oriptr = patloc + patlen)
      {
        size_t const skplen = patloc - oriptr;
        // copy the section until the occurence of the pattern
        strncpy(retptr, oriptr, skplen);
        retptr += skplen;
        // copy the replacement
        strncpy(retptr, replacement, replen);
        retptr += replen;
      }
      // copy the rest of the string.
      strcpy(retptr, oriptr);
    }
    return returned;
  }
}

__MSUTIL_ bool __export __system file_exists(const char * filename)
{

    FILE *file = NULL;

    if((file = fopen(filename, "r")))
    {
        fclose(file);
        return true;
    }
    return false;
}

__MSNATIVE_ bool __system readFile(const char path[static MAX_PATH_LENGTH])
{
    char c;
    FILE *fp = NULL;

    if(!(fp = fopen(path, "r")))
        return false;

    SHOWPAUSEMESSAGE();

    printf("\n%s:\n", path);
    PRINTL();

    while((c = fgetc(fp)) != EOF)
    {
        putchar(c);
        if(catchPause()) return true;
    }

    fclose(fp);
    PRINT2N();
    PRINTL();

    return true;
}

__MSNATIVE_ bool __system printFile(const char path[static MAX_PATH_LENGTH])
{
    if(!file_exists(path))
    {
        printErr(2, "Non-existent File:\n%s.\n\n", path);
        return false;
    }

    char str[MAX_PATH_LENGTH+INFO_STRING];

    #if WINOS
        printf2(COLOR_CREDITS, "Enter Device Name on which you want to print the File:\n%s.\n", path);
        printf2(COLOR_CREDITS, "Enter %c for Back.\n\n", SCANFEXIT_CHAR);

        char dname[MAX_STRING];

        while(scanf("%s", dname) != 1 || dname[0] == SCANFEXIT_CHAR)
        {
            if(dname[0] == SCANFEXIT_CHAR) return false;
            printErr(1, "Inserted string refers to an Invalid Device name.");
        }
        sprintf(str, "print /D:%s %s", dname, path);
    #else
        sprintf(str, "lpr -#1 -h -sP "DEFAULT_LINUX_SPOOLFOLDER" %s", path);
    #endif

    system(str);

    return true;
}

__MSNATIVE_ inline bool __system writeFile(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp;
    if((fp = checkForFHErrors(path, "w")) == NULL)
        return false;

    fclose(fp);
    return true;
}

__MSNATIVE_ inline FILE * __system checkForFHErrors(const char path[static MAX_PATH_LENGTH], char mode[static 1])
{
    FILE *fp;

    if((fp = fopen(path, mode)) == NULL)
        printErr(2, "Not possible to %s File:\n%s", mode[0] == 'r' ? "read":"write/create", path);

    return fp;
}

__MSNATIVE_ inline bool __system frename(const char name[static MAX_PATH_LENGTH], const char newname[static MAX_PATH_LENGTH])
{
    int err;

    if((err = rename(name, newname)))
    {
        printErr(err, "An error occurred during File Renaming Process: \n%s", name);
        return false;
    }

    return true;
}

// Conversion of matrix to upper triangular
// by Bibek Subedi original, adapted by me
__MSUTIL_ bool matrixUTConv(ityp **mat, dim_typ dimq)
{
    dim_typ i, j, k;

    ityp ratio;

    for(i = 0; i < dimq; ++i)
    {
        if( ISZERO(mat[i][i]) ) return false;
        for(j = 0; j < dimq; ++j)
            if(j>i)
            {
                ratio = mat[j][i]/mat[i][i];
                for(k = 0; k < dimq; ++k)
                    mat[j][k] -= ratio * mat[i][k];
            }
    }

    return true;
}

__MSNATIVE_ inline const ityp sarrus(ityp **mat)
{
    return ((mat[0][0]*mat[1][1]*mat[2][2])+(mat[0][1]*mat[1][2]*mat[2][0])+(mat[0][2]*mat[1][0]*mat[2][1])) -
            ((mat[0][2]*mat[1][1]*mat[2][0])+(mat[0][0]*mat[1][2]*mat[2][1])+(mat[0][1]*mat[1][0]*mat[2][2]));
}

#define CARLUCCI_DIM 4

// It calculates the determinant of a nXn matrix
// by up-triangularizing the square matrix passed
__MSNATIVE_ __MSUTIL_ ityp det(ityp **mat, dim_typ dimq, bool *flag)
{

    if(dimq > 0 && dimq < CARLUCCI_DIM) return checkStdMat(mat, dimq);

    // Conversion of matrix to upper triangular
    ityp D = 1; // storage for determinant

    if(matrixUTConv(mat, dimq))
    {
        dim_typ i;

        for(i = 0; i < dimq; ++i)
            D *= mat[i][i];
    }
    else
    {
        (*flag) = true;
        ityp * S = NULL;
        ityp ** V = NULL;

        S = malloc(sizeof(ityp)*dimq);
        errMem(S, MAX_VAL);

        if(!matrixAlloc(&V, (dim_typ2){dimq, dimq}))
        {
            free(S);
            matrixFree(&V, dimq);
            D = MAX_VAL;
        }
        else
        {
            dsvd(mat, (dim_typ2){dimq, dimq}, S, mat);
            D = fabs(productory(dimq, false, S));
            free(S);
        }

    }

    return D;
}

__MSNATIVE_ ityp trace(ityp **mat, dim_typ dimq)
{
    dim_typ i;
    ityp res;

    res = 0;

    for(i=0; i<dimq; ++i)
        res += mat[i][i];

    return res;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ ityp carlucci(ityp **mat)
{
 ityp det = 0.00;
 for(int i = 0; i <= CARLUCCI_DIM+1; ++i)
  if(i < CARLUCCI_DIM)
   det += i < CARLUCCI_DIM ? (pow(-1,i)*mat[CARLUCCI_DIM-2][i]*mat[CARLUCCI_DIM-1][(i+1)%CARLUCCI_DIM] + pow(-1,i+1)*mat[CARLUCCI_DIM-2][(i+1)%CARLUCCI_DIM]*mat[CARLUCCI_DIM-1][i])*(mat[0][(CARLUCCI_DIM-2+i)%CARLUCCI_DIM]*mat[1][(CARLUCCI_DIM-1+i)%CARLUCCI_DIM] - mat[0][(CARLUCCI_DIM-1+i)%CARLUCCI_DIM]*mat[1][(CARLUCCI_DIM-2+i)%CARLUCCI_DIM]) :
       (mat[CARLUCCI_DIM-2][i-CARLUCCI_DIM+2]*mat[CARLUCCI_DIM-1][i-CARLUCCI_DIM] - mat[CARLUCCI_DIM-2][i-CARLUCCI_DIM]*mat[CARLUCCI_DIM-1][i-CARLUCCI_DIM+2])*(mat[0][1-(i-CARLUCCI_DIM)]*mat[1][(CARLUCCI_DIM-1)-(i-CARLUCCI_DIM)] - mat[0][(CARLUCCI_DIM-1)-(i-CARLUCCI_DIM)]*mat[1][1-(i-CARLUCCI_DIM)]);

 return det;

}

__MSSHELL_WRAPPER_ __MSNATIVE_ ityp checkStdMat(ityp **a, dim_typ n)
{
    ityp det = 0.00;
    switch(n)
    {

        case 1:
            det = a[0][0];
            break;
        case 2:
            det = ((a[0][0]*a[1][1]) - (a[0][1]*a[1][0]));
            break;
        case 3:
            det = sarrus(a); // apply SARRUS algorithm
            break;
        case 4:
            det = carlucci(a);
            break;
    }

    return det;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ bool randomMatrix(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])  
{
    dim_typ range;
    printf2(COLOR_CREDITS, "\n\nEnter a non-negative integer to set pseudo-random numbers Range.\n");

    if(PARSING_SYSTEM_ALLOWED)
        PRINTHOWTOBACKMESSAGE();

    ityp tmp;

    while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
        (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != (range = (dim_typ)tmp) || range < 1 || range > SHRT_MAX)
    {
        CLEARBUFFER();
        if(access(exitHandle) == 1) continue;
        if(exitHandleCheck) return false;
        printErr(33, "Invalid inserted Value");
    }

    CLEARBUFFER();

    // RANDOMIZING THE MATRIX
    dim_typ i, j;

	#pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
            matrix[i][j] = random(range);

    printf2(COLOR_SYSTEM, "\n\n[%hu X %hu] Randomized Matrix with Range: %hu is:\n", dim[RAWS], dim[COLUMNS], range);
    PRINTL();

    printMatrix(stdout, matrix, dim);

    return true;
}

__MSNATIVE_ void transpose(ityp **matrix, ityp **matrix2, const register dim_typ dim[static MAX_DIMENSIONS])  
{

    dim_typ i, j;
    #pragma omp parallel for
    for(i=0; i<dim[COLUMNS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[RAWS]; ++j)
            matrix2[i][j] = matrix[j][i];

    return;
}


/// thanks to: http://www.di.unipi.it/~bozzo/fino/appunti/3/lr.c
/// for this implementation of LU Decomposition.
/*
   Calcola la fattorizzazione LU della matrice passata.
   In pratica esegue la prima fase dell'algoritmo di eliminazione
   di Gauss. La matrice "a" viene via via trasformata in forma
   triangolare superiore (diventa U) e, memorizzando tutti i
   valori di "m" si ottiene la parte inferiore di L.
   Restituisce 0 se non riesce a fattorizzare.
*/

__MSUTIL_ bool __export FattLU(dim_typ n, ityp **c, ityp **l, ityp **a)
{
    ityp m, pivot;
    dim_typ i, j, k;

    // Copia di C su A.
    if(!equalMatrix(&a, c, (dim_typ2){n, n}))
        return false;

    for (k=0; k<n; ++k)
    {
        pivot = a[k][k];
        if ( ISZERO(pivot) ) return false;
        for (i=k+1; i<n; ++i)
        {
            l[i][k] = m = a[i][k] / pivot ;
            for (j=k; j<n; ++j)
                a[i][j] -= m*a[k][j] ;
        }
    }

  /* Adesso "a" contiene U, e "l" la parte inferiore di L */

	#pragma omp parallel for
    for (i=0; i<n; ++i) /* Completa L con zeri ed uni */
    {
        l[i][i]=1.00;
        #pragma omp parallel for
        for (j=i+1; j<n; ++j)
            l[i][j]=0.00;
    }

  return true;
}


/*
thanks to: Bibek Subedi:
http://programming-technique.blogspot.it/2011/09/numerical-methods-inverse-of-nxn-matrix.html
for this part of code, which I renamed, modified and adapted to this program
*/

__MSUTIL_ bool __export invertMatrix(ityp **matrix, dim_typ n)
{

    dim_typ i, j;
    ityp a;
	
	#pragma omp parallel for
    for(i = 0; i < n; ++i)
   		#pragma omp parallel for
        for(j = n; j < 2*n; ++j)
            matrix[i][j] = i == (j-n);

    if(!matrixUTConv(matrix, n))
        return false;

	#pragma omp parallel for
    for(i = 0; i < n; ++i)
    {
        const ityp a = matrix[i][i];
        for(j = 0; j < 2*n; ++j)
            matrix[i][j] /= a;
    }

    return true;
}

__MSUTIL_ static inline _MS__private __export ityp PYTHAG(ityp a, ityp b)
{
    const ityp at = fabs(a), bt = fabs(b);
    ityp ct;
    return (at > bt ? at * sqrt(1.00 + ((ct=bt/at) * ct)) : (bt > 0.00 ? bt * sqrt(1.00 + ((ct=at/bt) * ct)) : 0.00));
}

// Jacobi Singular Value Decomposition (SVD)
__MSUTIL_ bool __export dsvd(ityp **a, const register dim_typ dim[static MAX_DIMENSIONS], ityp *w, ityp **v)
{
    // checking whether m < n and correct it
    // by transposing the matrix into n,m matrix
    // causes its correctly dsvd decomposition...

    const register dim_typ m = dim[RAWS], n = dim[COLUMNS];

    if(m < n)
    {
        ityp ** tmp = NULL;
        ityp res;
        if(!matrixAlloc(&tmp, (dim_typ2){n, m}))
            return false;
        transpose(a, tmp, dim); // m, n);
        res = dsvd(tmp, (dim_typ2){n, m}, w, v);
        matrixFree(&tmp, n);
        return res;
    }

    bool flag;
    int i, its, j, jj, k, l, nm;
    ityp c, f, h, s, x, y, z;
    ityp anorm = 0.00, g = 0.00, scale = 0.00;
    ityp *rv1;

    rv1 = malloc(sizeof(ityp)*n);
    errMem(rv1, false);

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; ++i)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.00;
        if (i < m)
        {
            for (k = i; k < m; ++k)
                scale += fabs(a[k][i]);
            if (scale)
            {
                for (k = i; k < m; ++k)
                {
                    a[k][i] = (a[k][i]/scale);
                    s += (a[k][i] * a[k][i]);
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                if (i != n - 1)
                {
                    for (j = l; j < n; ++j)
                    {
                        for (s = 0.00, k = i; k < m; ++k)
                            s += (a[k][i] * a[k][j]);
                        f = s / h;
                        for (k = i; k < m; ++k)
                            a[k][j] += (f * a[k][i]);
                    }
                }
                for (k = i; k < m; ++k)
                    a[k][i] = a[k][i]*scale;
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.00;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; ++k)
                scale += fabs(a[i][k]);
            if (scale)
            {
                for (k = l; k < n; ++k)
                {
                    a[i][k] = a[i][k]/scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k < n; ++k)
                    rv1[k] = a[i][k] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; ++j)
                    {
                        for (s = 0.00, k = l; k < n; ++k)
                            s += a[j][k] * a[i][k];
                        for (k = l; k < n; ++k)
                            a[j][k] += s * rv1[k];
                    }
                }
                for (k = l; k < n; ++k)
                    a[i][k] = a[i][k]*scale;
            }
        }
        register ityp reg;
        if((reg = (fabs(w[i]) + fabs(rv1[i]))) < anorm)
            anorm = reg;
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; --i)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; ++j)
                    v[j][i] = ((a[i][j] / a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.00, k = l; k < n; ++k)
                        s += a[i][k] * v[k][j];
                    for (k = l; k < n; ++k)
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j < n; ++j)
                v[i][j] = v[j][i] = 0.00;
        }
        v[i][i] = 1.00;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; --i)
    {
        l = i + 1;
        g = w[i];
        if (i < n - 1)
            for (j = l; j < n; ++j)
                a[i][j] = 0.00;
        if (g)
        {
            g = 1.00 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; ++j)
                {
                    for (s = 0.00, k = l; k < m; ++k)
                        s += a[k][i] * a[k][j];
                    f = (s / a[i][i]) * g;
                    for (k = i; k < m; ++k)
                        a[k][j] += f * a[k][i];
                }
            }
            for (j = i; j < m; ++j)
                a[j][i] = a[j][i]*g;
        }
        else
        {
            for (j = i; j < m; ++j)
                a[j][i] = 0.00;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; --k)
    {                             /* loop over singular values */
        for (its = 0; its < 30; ++its)
        {                         /* loop over allowed iterations */
            flag = true;
            for (l = k; l >= 0; --l)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = false;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.00;
                s = 1.00;
                for (i = l; i <= k; ++i)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h;
                        h = 1.00 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; ++j)
                        {
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = (y * c + z * s);
                            a[j][i] = (z * c - y * s);
                        }
                    }
                }
            }
            z = w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.00)
                {              /* make singular value nonnegative */
                    w[k] = (-z);
                    for (j = 0; j < n; ++j)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= (access(curLayout)->max_dsvd_iterations/1000))
            {
                printErr(33, "No convergence after %d! iterations", access(curLayout)->max_dsvd_iterations);
                free((void*) rv1);
                return false;
            }

            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.00 * h * y);
            g = PYTHAG(f, 1.00);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.00;
            for (j = l; j <= nm; ++j)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; ++jj)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = (x * c + z * s);
                    v[jj][i] = (z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z)
                {
                    z = 1.00 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; ++jj)
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = (y * c + z * s);
                    a[jj][i] = (z * c - y * s);
                }
            }
            rv1[l] = 0.00;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free((void*) rv1);
    return true;
}

// RANK Calculator of a nXm Matrix using
// Jacobi Singular Value Decomposition METHOD (SVD)
__MSNATIVE_ dim_typ __export rank(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    ityp * S = NULL;
    ityp ** V = NULL;

    const register ityp maxv = dim[dim[RAWS] >= dim[COLUMNS]];

    S = malloc(sizeof(ityp)*maxv);
    errMem(S, USHRT_MAX);

    if(!matrixAlloc(&V, (dim_typ2){maxv, maxv}))
    {
        free(S);
        return USHRT_MAX;
    }

    dsvd(matrix, dim, S, V);
    matrixFree(&V, maxv);

    register dim_typ rnk = maxv;

    for(dim_typ i=0; i<maxv; ++i)
        if(ISZERO(S[i])) -- rnk;

    free(S);

    return rnk;
}

__MSNATIVE_ void __export matrixToVector(ityp **mdim_array, const register dim_typ dim[static MAX_DIMENSIONS], ityp vector[static dim[RAWS]*dim[COLUMNS]], bool direction)
{
    dim_typ i, j;
    int k;

    for(i=k=0; i<dim[RAWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            if(!direction) vector[k++] = mdim_array[i][j];
            else mdim_array[i][j] = vector[k++];
    return;
}

__MSNATIVE_ void __export vectorToMatrix(const register dim_typ dim[static MAX_DIMENSIONS], ityp vector[static dim[RAWS]*dim[COLUMNS]], ityp **mdim_array)
{
    matrixToVector(mdim_array, dim, vector, VECTOR_TO_MATRIX);
    return;
}

__MSNATIVE_ inline void __system _flushLogBuf(logObj * const which_log)
{
    strcpy(which_log->buffer, NULL_CHAR);
    return;
}

__MSNATIVE_ __WINCALL void __system _editLog(const char path[static MAX_PATH_LENGTH])
{
    #if (WINOS)
    if(isnSett(BOOLS_SYSLOGSECURITYCHECK))
    {
        char str[PATHCONTAINS_STRING] = NULL_CHAR;
        sprintf(str, "notepad.txt %s", path);
        system(str);
    }
    else
    #endif
    {
        printf2(COLOR_SYSTEM, "You entered Sequential Log Editing:\n%s.\nPress CTRL + C to exit.\n\n", path);

        FILE *fp = NULL;

        if((fp = checkForFHErrors(path, "a")) == NULL)
            return;

        fputs(NULL_CHAR, fp);

        char str[MAX_BUFSIZ] = NULL_CHAR;
        for( ; ; )
        {
            signal(SIGINT, (__p_sig_fn_t) sigexit);
            gets(str);
            if(access(sigresult))
            {
                fclose(fp);
                printf2(COLOR_SYSTEM, "\nYou exited Sequential Log Editing:\n%s.\n\n", path);
                access(sigresult) = false;
                return;
            }
            fputs(str, fp);
            /// writelog functions
        }
        fclose(fp); /// in order to provide unaligned strongly dangerous goto.
    }
    return;
}

__MSNATIVE_ void __system logCheck(logObj * const which_log, const char *format, const char path[static MAX_PATH_LENGTH])
{
    if(strlen(which_log->buffer) < which_log->buflen)
        strcat(which_log->buffer, format);
    else
    {
        FILE *fp;

        if((fp = checkForFHErrors(path, "a")) == NULL)
            return;

        logWrite(fp, which_log);
        fclose(fp);
    }

    return;
}

__MSNATIVE_ inline void __system prependTimeToString(char *string, const bool condition)
{
    if(access(curLayout) && (isSett(BOOLS_PRINTTIME) || condition))
    {
        char str_time[MAX_BUFSIZ+INFO_STRING];
        time_t ora;
        ora = time(NULL);
        sprintf(str_time, "[%s] - %s", asctime(localtime(&ora)), string);
        strcpy(string, str_time);
    }
}

__MSNATIVE_ void __system printErr(const int err, const char *format, ...)
{
    CLEARBUFFER();
    PRINTL();

    _sprint(format, COLOR_ERROR);

    errno = err;
    perror(".\nERRORE");

    PRINTL();

    return;
}

__MSNATIVE_ void __system sprint(const char *format, ...)
{
    _sprint(format, COLOR_SYSTEM);
    return;
}

__MSNATIVE_ void __system fprintf2(FILE *fp, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);

    char str[MAX_BUFSIZ];
    vsprintf(str, format, ap);

	#if WINOS
	    const bool cond = fp == stdout;
	    if(cond)
	        SetColor(COLOR_USER);
	#endif

    fprintf(fp, str);

    #if WINOS
        if(cond)
            SetDefaultColor();
    #endif


    if(getItemsListNo(LOGS) != STARTING_LOGSNO)
    {
        nodelist * const cur_log = listNo(access(lists)[LOGS].cur_item, LOGS);
        logCheck(((logObj *)(cur_log->data)), str, cur_log->path);
    }

    va_end(ap);

    return;
}


__MSNATIVE_ void __system printf2(const sel_typ col, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);

    char str[MAX_BUFSIZ];
    static sel_typ col_cache = MAX_COLORS;
    const bool cond = getItemsListNo(LOGS) != STARTING_LOGSNO;

	#if WINOS
   		SetColor(col);	
   	#endif
    vsprintf(str, format, ap);
    if(col != col_cache)
        prependTimeToString(str, cond);
    printf(str);
    SetDefaultColor();

    if(cond)
    {
        nodelist * const cur_log = listNo(access(lists)[LOGS].cur_item, LOGS);
        logCheck(((logObj *)(cur_log->data)), str, cur_log->path);
    }

    va_end(ap);
    col_cache = col;
    return;
}

__MSNATIVE_ bool __system scanf2(sel_typ count, const char *format, ...)
{
    int scanner;

    scanner = 0;

    va_list ap;
    va_start(ap, format);

    access(exitHandle) = INVALID_EXITHANDLE;
    signal(SIGINT, (__p_sig_fn_t) sigproc);
    scanner = vscanf(format, ap);

    va_end(ap);

    if(access(sigresult))
        access(exitHandle) = EXITHANDLE_GETCMD;

    if(!scanner)
        access(exitHandle) = EXITHANDLE_EXIT;

    return scanner == count;
}

__MSNATIVE_ _MS__private void __system printMatrix(FILE *fp, ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])  
{

    const bool assert = fp == stdout;

    if(assert)
        PRINTN();

    dim_typ  i, j;

    for(i=0; i<dim[RAWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
            fprintf2(fp, OUTPUT_CONVERSION_FORMAT, matrix[i][j]);
            fprintf2(fp, "; ");
            if(j >= dim[COLUMNS]-1)
               fputc('\n', fp);
        }

    if(assert)
    {
        PRINT2N();
        PRINTN();
        
        if(access(lmpMatrix) && !access(lmpMatrix)->matrix)
        	matrixFree(&(access(lmpMatrix)->matrix), access(lmpMatrix)->dim[RAWS]);
        
        if(access(lmpMatrix) && !matrixAlloc(&(access(lmpMatrix)->matrix), dim))
        	resetLmpMatrix();

        if(matrix && access(lmpMatrix)->matrix && equalMatrix(&(access(lmpMatrix)->matrix), matrix, dim))
        {
        	/// printf("\n\nROMUALDO\n\n")
            access(lmpMatrix)->dim[RAWS] = dim[RAWS];
            access(lmpMatrix)->dim[COLUMNS] = dim[COLUMNS];
        }
    }

    return;
}

__MSNATIVE_ bool __system __export parse(char expr[], ityp *res)
{

    volatile int err;
    exprObj * exp = INIT_OBJLIST;

    err = exprCreate(&exp, access(func_list), access(exprVars)->var_list, access(const_list), NULL, 0);

    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Expr Creation Error.\n");
        exprFree(exp);
        return false;
    }

    if(expr[strlen(expr)-1] != TERMINATING_CHAR)
        strcat(expr, TERMINATING_STRING);

    err = exprParse(exp, expr);
    if(err != EXPR_ERROR_NOERROR)
    {
        int start, end;
        exprGetErrorPosition(exp, &start, &end);
        sprint("Parse Error (%d,%d).\n", start, end);
        exprFree(exp);
        return false;
    }

    ityp val;
    err = exprEval(exp, &val);

    if(err != EXPR_ERROR_NOERROR)
    {
        sprint("Eval Error: %d.\n", err);
        exprFree(exp);
        return false;
    }

    (*res) = val;
    exprFree(exp);

    return true;
}

/*
// STATIC MATRIX EXTRACTION SYSTEM
// this system is faster than the actual one,
// but is slightly deprecated.
// This went promptly discarded.
bool extractMat(ityp ***matrix, dim_typ *righe, dim_typ *colonne, char item_path[static MAX_PATH_LENGTH])
{

    FILE *fp;

    fp = NULL;


    if((fp = checkForFHErrors(item_path, "r")) == NULL)
       return false;

    fflush(fp);

    char str[MAX_BUFSIZ] = NULL_CHAR; // char str[MAX_BUFSIZ] = NULL_CHAR;


    exprObj *e = INIT_OBJLIST;
    int start, end;
    time_t t1, t2;
    double diff;
    diff = 0;
    volatile int err;

    size_t len;

    (*righe) = (*colonne) = 0;

    if(fscanf(fp, "%hu", righe) != 1) // fgets(str, MAX_FILE_LINES, fp) == NULL)
    {
        fclose(fp);
        return false;
    }

    // len = strlen(str)-1;
    // str[len] = '\0';

    // (*righe) = strtod(str, NULL);

    if(fscanf(fp, "%hu", colonne) != 1) // fgets(str, MAX_FILE_LINES, fp) == NULL)
    {
        fclose(fp);
        return false;
    }

    // len = strlen(str)-1;
    // str[len] = '\0';

    // (*colonne) = strtod(str, NULL);

    if(!matrixAlloc(matrix, (*righe), (*colonne)))
    {
        fclose(fp);
        return false;
    }

    dim_typ i, j;

    for(i=0; i<(*righe); i++)
    {
        for(j=0; j<(*colonne); j++)
        {
            // printf("\nCURRENT RAWCOLUMN: [%hu, %hu]\n\n", i, j);
            if(fgets(str, MAX_FILE_LINES, fp) == NULL)
            {
                fclose(fp);
                return false;
            }

            err = exprCreate(&e, suite.exprVars.func_list, suite.exprVars.var_list, suite.exprVars.const_list, NULL, 0);

            if(err != EXPR_ERROR_NOERROR)
            {
                printf("Expr Creation Error.\n");
                exprFree(e);
                continue;
            }


            err = exprParse(e, str);
            if(err != EXPR_ERROR_NOERROR)
            {
                exprGetErrorPosition(e, &start, &end);
                printf("Parse Error (%d,%d).\n", start, end);
                exprFree(e);
                continue;
            }

            if(suite.bools[BOOLS_SHOWDIFFTIME].state)
                t1 = time(NULL);

            ityp val;
            err = exprEval(e, &val);

            if(err != EXPR_ERROR_NOERROR)
                printf("Eval Error: %d.\n", err);

            if(suite.bools[BOOLS_SHOWDIFFTIME].state)
            {
                t2 = time(NULL);
                diff += difftime(t2, t1);
            }

            (*matrix)[i][j] = val;

            exprFree(e);

        }

    }

    fclose(fp);

    if(suite.bools[BOOLS_SHOWDIFFTIME].state)
    {
        PRINTL();
        printf("Average Time: %.*f;\n", suite.precision, (EXPRTYPE)diff);
        PRINTL();
    }

    return true;
}
*/

#define MINMAX_BUFFER_LEN MAX_STRING

__MSNATIVE_ bool __system __export extractMat(dim_typ which_mat)
{

    FILE *fp;
    time_t t1;
    ityp diff = 0.00;
    char str[MINMAX_BUFFER_LEN] = NULL_CHAR;

    // char *ij_element = NULL;

    dim_typ i;

    matrixObj * const tmp = malloc(sizeof(matrixObj));
    errMem(tmp, false);

    tmp->matrix = NULL;

    fp = NULL;

    t1 = time(NULL);
    if((fp = checkForFHErrors(listNo(which_mat, MATRICES)->path, "r")) == NULL)
        return false;

    fflush(fp);


    char *ij_element = NULL;


    for(tmp->dim[RAWS]=tmp->dim[COLUMNS]=INIT_DIM; fgets(str, sizeof(str), fp) != NULL; ++ tmp->dim[RAWS])  // fscanf(fp, "%[^\n]s", str)) // fgets(str, sizeof(str), fp) != NULL)
    {
        if(!(tmp->dim[RAWS]) && !matrixAlloc(&(tmp->matrix), (dim_typ2){1, lazy_exec ? 1 : access(curLayout)->stabilizer_factor}))
            return false;
        else
        {
            const dim_typ analog_raws = (tmp->dim[RAWS])+1;

            if(lazy_exec)
                tmp->matrix = realloc(tmp->matrix, sizeof(ityp *)*analog_raws);
            else if(!(analog_raws % access(curLayout)->stabilizer_factor))
                tmp->matrix = realloc(tmp->matrix, sizeof(ityp *) * (tmp->dim[RAWS]+access(curLayout)->stabilizer_factor));

            errMem(tmp->matrix, false);

            (tmp->matrix)[tmp->dim[RAWS]] = malloc((tmp->dim[COLUMNS])*sizeof(ityp));
            errMem((tmp->matrix)[tmp->dim[RAWS]], false);
        }


        for(i = 0, ij_element = strtok(str, MATRIXES_SEPERATOR_STRING); ij_element != NULL; ++i, ij_element = strtok(NULL, MATRIXES_SEPERATOR_STRING))
        {
            if(strrchr(ij_element, '\n') != NULL)
            {
                ij_element = NULL;
                break;
            }

            if(!(tmp->dim[RAWS]))
            {
                const dim_typ analog_columns = (tmp->dim[COLUMNS])+1;
                if(lazy_exec)
                    (tmp->matrix)[0] = realloc((tmp->matrix)[0], sizeof(ityp)*analog_columns);
                else if(!(analog_columns % access(curLayout)->stabilizer_factor))
                    (tmp->matrix)[0] = realloc((tmp->matrix)[0], sizeof(ityp)*((tmp->dim[COLUMNS]) + access(curLayout)->stabilizer_factor));

                errMem((tmp->matrix)[0], false);
            }

            if(isSett(BOOLS_MATRIXPARSING) && !parse(ij_element, &((tmp->matrix)[tmp->dim[RAWS]][i])))
                continue;
            else
                (tmp->matrix)[tmp->dim[RAWS]][i] = strtod(ij_element, NULL);

            if(!(tmp->dim[RAWS]))
               ++(tmp->dim[COLUMNS]);
        }
    }

    fclose(fp);
    diff = difftime(time(NULL), t1);

    if(!(lazy_exec))
    {
        (tmp->matrix) = realloc(tmp->matrix, sizeof(ityp*)*(tmp->dim[RAWS]));
        errMem(tmp->matrix, false);
    }

    if(isSett(BOOLS_SHOWDIFFTIME))
    {
        PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", DEFAULT_PRECISION, (EXPRTYPE)diff);
        PRINTL();
    }

    // REDIRECTING WHICH_MAT LISTBOX ITEM
    listNo(which_mat, MATRICES)->data = tmp;

    return true;
}

__MSNATIVE_ bool __system __export matrixToken(const char string[], ityp ***matrix, dim_typ *righe, dim_typ *colonne)
{

    char target[MAX_BUFSIZ];

    strcpy(target, string);

    char *line;
    char *token;
    char buf[MAX_STRING];

    /// char str2[MAX_BUFSIZ];
    dim_typ i;

    line = token = NULL;

    for((*righe) = (*colonne) = INIT_DIM, line = strtok(target, TERMINATING_STRING); line != NULL; ++ (*righe), line = strtok (line + strlen (line) + 1, TERMINATING_STRING))
    {
        /* String to scan is in buf string..token */
        strncpy(buf, line, sizeof(buf));

        if((!(*righe)) && !matrixAlloc(matrix, (dim_typ2){1, lazy_exec ? 1 : access(curLayout)->stabilizer_factor}))
            return false;
        else
        {

            const dim_typ analog_raws = ((*righe))+1;

            if(lazy_exec)
                (*matrix) = realloc((*matrix), sizeof(ityp *)*analog_raws);
            else if(!(analog_raws % access(curLayout)->stabilizer_factor))
                (*matrix) = realloc((*matrix), sizeof(ityp *) * (((*righe))+access(curLayout)->stabilizer_factor));

            errMem((*matrix), false);

            (*matrix)[(*righe)] = malloc(((*colonne))*sizeof(ityp));
            errMem((*matrix)[(*righe)], false);
        }

        for(i=0, token=strtok(buf,","); token != NULL; ++ i, token = strtok (token + strlen (token) + 1, ","))
        {

            if(!((*righe)))
            {
                const dim_typ analog_columns = ((*colonne))+1;
                if(lazy_exec)
                    (*matrix)[0] = realloc((*matrix)[0], sizeof(ityp)*analog_columns);
                else if(!(analog_columns % access(curLayout)->stabilizer_factor))
                    (*matrix)[0] = realloc((*matrix)[0], sizeof(ityp)*(((*colonne))+ access(curLayout)->stabilizer_factor));

                errMem((*matrix)[0], false);
            }

            char token2[strlen(token)+1];
            sprintf(token2, "%s;", token);

            if(isSett(BOOLS_MATRIXPARSING) && !parse(token2, &((*matrix)[(*righe)][i])))
                continue;
            else
                (*matrix)[(*righe)][i] = strtod(token2, NULL);

            if(!((*righe)))
                ++ (*colonne);

        }
    }

    return true;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ const sprog * const __system searchProgram(const char cmdname[static SIGN_STRING])
{
    dim_typ i;

    for(i=0; i<MAX_PROGRAMMI; ++i)
        if(!strcmp(cmdname, main_menu[i].cmdname))
            return &main_menu[i];

    for(i=0; i<MAX_ADVCALC_PROGS; ++i)
        if(!strcmp(cmdname, adv_calc[i].cmdname))
            return &adv_calc[i];

    for(i=0; i<MAX_ALGEBRA_OPERATIONS; ++i)
        if(!strcmp(cmdname, alg_operations[i].cmdname))
            return &alg_operations[i];

    for(i=0; i<MAX_ENVSMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, envs_manager[i].cmdname))
            return &envs_manager[i];

    for(i=0; i<MAX_MATMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, mat_manager[i].cmdname))
            return &mat_manager[i];

    for(i=0; i<MAX_LOGSMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, logs_manager[i].cmdname))
            return &logs_manager[i];

    for(i=0; i<MAX_SYSLOGMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, syslog_manager[i].cmdname))
            return &syslog_manager[i];

    for(i=0; i<MAX_LAYOUTSMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, layouts_manager[i].cmdname))
            return &layouts_manager[i];

    for(i=0; i<MAX_SETTINGS; ++i)
        if(!strcmp(cmdname, change_settings[i].cmdname))
            return &change_settings[i];

    for(i=0; i<MAX_COLSMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, cols_manager[i].cmdname))
            return &cols_manager[i];

    for(i=0; i<MAX_LFSMANAGER_PROGS; ++i)
        if(!strcmp(cmdname, lfs_manager[i].cmdname))
            return &lfs_manager[i];

    return NULL;

}

__MSUTIL_ inline int __export cmpfunc(const void * a, const void * b)
{
   return ( *(ityp*)a - *(ityp*)b );
}

__MSNATIVE_ inline void __system _changeAlgebraDims(const dim_typ dim)
{
    access(curLayout)->algebra = dim;
    access(matrixSumFunc) = suite_c.matrixSumFuncs[dim];
    access(matrixProdFunc) = suite_c.matrixProdFuncs[dim];
    access(matrixKProdFunc) = suite_c.matrixKProdFuncs[dim];
    return;
}

__MSNATIVE_ inline void __system _showUsage(const sprog * const prog)
{
    printf2(COLOR_ERROR, "\nUSAGE: ");
    printf2(COLOR_SYSTEM, "%s %s;\n", prog->cmdname, prog->usage);
    return;
}

__MSNATIVE_ void __system printUsage(const sprog * const prog)
{
    _showUsage(prog);
    prog->program_function(0, NULL);
    return;
}

__MSNATIVE_ void __system prepareToExit(void)
{
	flushAllMemoizersBuffers();
	
    #if WINOS
        _backupColFile();
    #endif

    if(isSett(BOOLS_ITEMSAUTOSAVING))
        for(dim_typ i=0; i<MAX_LISTS; ++i)
            updAll(i);
    
    printf2(COLOR_CREDITS, EXIT_MESSAGE);
    PRINTL();
    return;
 }

__MSNATIVE_ inline void __system safeExit(const volatile int exval)
{
    prepareToExit();
    exit(exval);
    return;
}

__MSNATIVE_ void __system _handleCmdLine(const register sel_typ argc, char ** argv)
{
    // catch _MSS_CMD exception
    if(!strcmp(argv[0], _MSS_CMD))
    {
        access(mss) = true;
        printf2(COLOR_USER, "\nScripting Mode has been enabled.\n\n");
        return;
    }

    // catch EXIT_CMD exception
    if(argc == MAX_DIMENSIONS && !strcmp(argv[0], EXIT_CMD))
    {
        ityp tmp;
        if(PARSING_SYSTEM_ALLOWED)
        {

            if(!parse(argv[1], &tmp))
            {
                printErr(1, "Parse Error on "EXIT_CMD" command.");
                return;
            }
        }
        else
            tmp = strtod(argv[1], NULL);

        if(tmp < INT_MIN || tmp > INT_MAX)
        {
            printErr(33, EXIT_CMD" accepts only integers between %d and %d", INT_MIN, INT_MAX);
            return;
        }

        safeExit(tmp);
    }

    const sprog * const prog = searchProgram(argv[0]);

    if(prog)
        prog->program_function(argc-1, &argv[1]);
    else
        printErr(1, "SubProgram hasn't been found in Program Macro-List");

    return;
}

__MSNATIVE_ bool __system _execScriptFiles(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp = NULL;

    if((fp = checkForFHErrors(path, "r")))
    {
        dim_typ i;
        char str[MAX_FILE_LINES];

        for( ; fgets(str, sizeof(str), fp) ; )
        {
            char *cmdtab[MAX_ARGS];
            const size_t len = strlen(str)-1;

            if(str[len] == '\n')
                str[len] = '\0';

            if(access(mss))
            {
                char mss_apnt[MAX_FILE_LINES+SIGN_STRING];
                sprintf(mss_apnt, CMD_BASECALC" %s", str);
                strcpy(str, mss_apnt);
            }

            if(!str[0])
                continue;



            for(cmdtab[0]=strtok(str,BLANK_STRING),i=0; cmdtab[i] != NULL; cmdtab[++i] = strtok(NULL, BLANK_STRING));

            _handleCmdLine(i, cmdtab);
        }

        fclose(fp);
        return true;
    }

    return false;
}

__MSNATIVE_ bool __system _lfLoader(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp = NULL;

    if((fp = checkForFHErrors(path, "r")))
    {
        sel_typ mode;
        char str[MAX_PATH_LENGTH] = NULL_CHAR;

        for( ; fgets(str,sizeof(str),fp) ; )
        {
            const size_t len = strlen(str)-1;

            if(str[len] == '\n')
                str[len] = '\0';

            if((mode = checkItemTypeByExtension(strrchr(str, '.')+1)) != MAX_LISTS)
				listInsertProc(str, mode);
        }

        fclose(fp);
        return true;
    }

    return false;
}

__MSNATIVE_ bool __system _lfCreate(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp;

    if((fp = checkForFHErrors(path, "w")))
    {

        dim_typ i, j;
        dim_typ itemsno;

        for(i=0; i<MAX_LISTS; ++i)
            if((itemsno = getItemsListNo(i)) != STARTING_ITEMSNO)
                for(j=0; j<itemsno; ++j)
                    fputs(listNo(j, i)->path, fp);

        fclose(fp);
        return true;
    }

    return true;
}

__MSUTIL_ static inline void ftoa(char *string, const float value, const fsel_typ prec)
{
	sprintf(string, "%.*f", prec, value);
	return;
}
#ifdef XMLCALL

	__MSUTIL_ XMLCALL inline xmlDoc * __system __export xmlInit(const char file_name[static XML_FILENAMES_LENGTH], xmlXPathContext ** xpathCtx)
	{
		xmlInitParser();
		LIBXML_TEST_VERSION;
		xmlDoc * doc = xmlParseFile( file_name );
		(*xpathCtx) = xmlXPathNewContext( doc );
		return doc;
	}
	
	__MSUTIL_ XMLCALL inline void __system __export xmlExit(const char file_name[static XML_FILENAMES_LENGTH], xmlDoc ** doc, xmlXPathObject ** xpathObj, xmlXPathContext ** xpathCtx)
	{
		xmlSaveFileEnc(file_name, (*doc), XML_ENCODING);
		xmlFreeDoc((*doc));
    	xmlCleanupParser();
    	(*doc) = (xmlDoc*) NULL;
    	(*xpathObj) = (xmlXPathObject*) NULL;
    	(*xpathCtx) = (xmlXPathContext*) NULL;
    	return;
	}
	
	__MSUTIL_ XMLCALL inline void __system __export xmlWriteInt(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const int value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		char tmp[SIGN_STRING] = NULL_CHAR;
		itoa(value, tmp, sizeof(tmp));
		xmlNodeSetContent(node, BAD_CAST tmp);
		xmlXPathFreeObject((*xpathObj));
		return;
	}
	
	__MSUTIL_ XMLCALL inline void __system __export xmlWriteBool(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const bool value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlNodeSetContent(node, BAD_CAST suite_c.bools_identifiers[value]);
		xmlXPathFreeObject((*xpathObj));
		return;
	}
	
	__MSUTIL_ XMLCALL inline void __system __export xmlWriteFloat(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const float value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		char tmp[SIGN_STRING] = NULL_CHAR;
		ftoa(tmp, value, DEFAULT_PRECISION);
		xmlNodeSetContent(node, BAD_CAST tmp);
		xmlXPathFreeObject((*xpathObj));
		return;
	}
	
	__MSUTIL_ XMLCALL inline void __system __export xmlWriteString(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const char string[static MAX_XML_FIELDSTRINGS])
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlNodeSetContent(node, BAD_CAST string);
		xmlXPathFreeObject((*xpathObj));
		return;
	}
	
	__MSUTIL_ XMLCALL inline int __system __export xmlGetInt(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		return atoi(xmlNodeGetContent(node));
	}
	
	__MSUTIL_ XMLCALL inline bool __system __export xmlGetBool(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		return strcmp(xmlNodeGetContent(node), suite_c.bools_identifiers[false]);
	}
	
	__MSUTIL_ XMLCALL inline float __system __export xmlGetFloat(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		return atof(xmlNodeGetContent(node));
	}
	
	__MSUTIL_ XMLCALL inline void __system __export xmlGetString(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, char string[static MAX_XML_FIELDSTRINGS])
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		strcpy(string, xmlNodeGetContent(node));
		xmlXPathFreeObject((*xpathObj));
		return;
	}
	
	__MSNATIVE_ XMLCALL void __system getProgramSettings(dim_typ which_layout)
	{
		dim_typ i;
	    nodelist * const tmp = listNo(which_layout, LAYOUTS);
	    layoutObj * const cur_layout = (layoutObj *)(tmp->data);
	    
	    xmlXPathContext * xpathCtx = (xmlXPathContext*) NULL;
		xmlXPathObject * xpathObj = (xmlXPathObject*) NULL;
		xmlDoc * doc = xmlInit(tmp->path, &xpathCtx);
		char ex_char[MAX_XML_FIELDSTRINGS];
		sprintf(ex_char, "%c", cur_layout->exit_char);
		xmlWriteString(&xpathObj, xpathCtx, "/settings/generalSettings/exitChar", ex_char);
		
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/programPrecision", cur_layout->precision);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/stabilizerFactor", cur_layout->stabilizer_factor);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/minStirlingRequiresNumber", cur_layout->min_stirlingrequires_number);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/algebra", cur_layout->algebra);
		xmlWriteFloat(&xpathObj, xpathCtx, "/settings/generalSettings/outlierConstant", cur_layout->outlier_constant);
	
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxRaws", cur_layout->matrix_max_raws);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxColumns", cur_layout->matrix_max_columns);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxDSVDIterations", cur_layout->max_dsvd_iterations);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxSimplexIterations", cur_layout->max_simplex_iterations);
	
	
		#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
		for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		{
			char str[2*INFO_STRING] = NULL_CHAR;
			char strboolized[INFO_STRING] = NULL_CHAR;
			strboolize(suite_c.memoizers_names[i], strboolized);
			strboolized[0] = toupper(strboolized[0]);
			sprintf(str, "/settings/memoizerOptions/max%sMemoizableIndex", strboolized);
			xmlWriteInt(&xpathObj, xpathCtx, str, cur_layout->max_memoizable_indeces[i]);
		}
		
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/baseConversions/minBase", cur_layout->basecalc_minbase);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxBase", cur_layout->basecalc_maxbase);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxChangebaseBinaryConvnum", cur_layout->max_changebase_binary_convnum);
	
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/minDimension", cur_layout->min_newton_difftables_dim);
		xmlWriteInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/maxDimension", cur_layout->max_newton_difftables_dim);
	
	    xmlWriteInt(&xpathObj, xpathCtx, "/settings/romanNumbers/minProcessableNumber", cur_layout->min_roman_number);
	    xmlWriteInt(&xpathObj, xpathCtx, "/settings/romanNumbers/maxProcessableNumber", cur_layout->max_roman_number);
	    
	    xmlWriteInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/minRaws", cur_layout->pascal_triangle_min_raws);
	    xmlWriteInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/maxRaws", cur_layout->pascal_triangle_max_raws);
	    
	    /// write some stuffs...
	
		#pragma omp parallel for
	    for(i=0; i<MAX_BOOL_SETTINGS; ++i)
	    {   
	    	char name[4*MIN_STRING] = NULL_CHAR;
			char strboolized[2*MIN_STRING] = NULL_CHAR;
			strboolize(suite_c.bools_names[i], strboolized);
	    	sprintf(name, "/settings/booleanKeys/%s", strboolized);
	        xmlWriteBool(&xpathObj, xpathCtx, name, (cur_layout->bools & suite_c.bools[i].bmask) == suite_c.bools[i].bmask);
	    }
	    
	    xmlExit(tmp->path, &doc, &xpathObj, &xpathCtx);
	    return;
	}
	
	__MSNATIVE_ XMLCALL void __system resetProgramSettings(layoutObj * const tmp, const char path[static MAX_PATH_LENGTH])
	{
		dim_typ i;
		xmlXPathContext * xpathCtx = (xmlXPathContext*) NULL;
	    xmlXPathObject * xpathObj = (xmlXPathObject*) NULL;
	    xmlDoc * doc = xmlInit(path, &xpathCtx);
	
		char ex_char[1];
		xmlGetString(&xpathObj, xpathCtx, "/settings/generalSettings/exitChar", ex_char);
		tmp->exit_char = ex_char[0];
		tmp->precision = xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/programPrecision");
		tmp->stabilizer_factor = xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/stabilizerFactor");
		tmp->min_stirlingrequires_number = xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/minStirlingRequiresNumber");
		tmp->algebra = xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/algebra");
		tmp->outlier_constant = xmlGetFloat(&xpathObj, xpathCtx, "/settings/generalSettings/outlierConstant");
			
		tmp->matrix_max_raws = xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxRaws");
		tmp->matrix_max_columns = xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxColumns");
		tmp->max_dsvd_iterations = xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxDSVDIterations");
		tmp->max_simplex_iterations = xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxSimplexIterations");
	
		
		#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
		for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		{
			char str[2*INFO_STRING] = NULL_CHAR;
			char strboolized[INFO_STRING] = NULL_CHAR;
			strboolize(suite_c.memoizers_names[i], strboolized);
			strboolized[0] = toupper(strboolized[0]);
			sprintf(str, "/settings/memoizerOptions/max%sMemoizableIndex", strboolized);
			tmp->max_memoizable_indeces[i] = xmlGetInt(&xpathObj, xpathCtx, str);
		}
	
	
		tmp->basecalc_minbase = xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/minBase");
		tmp->basecalc_maxbase = xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxBase");
		tmp->max_changebase_binary_convnum = xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxChangebaseBinaryConvnum");
		
		tmp->min_newton_difftables_dim = xmlGetInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/minDimension");
		tmp->max_newton_difftables_dim = xmlGetInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/maxDimension");
		
		tmp->min_roman_number = xmlGetInt(&xpathObj, xpathCtx, "/settings/romanNumbers/minProcessableNumber");
		tmp->max_roman_number = xmlGetInt(&xpathObj, xpathCtx, "/settings/romanNumbers/maxProcessableNumber");
		
		tmp->pascal_triangle_min_raws = xmlGetInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/minRaws");
		tmp->pascal_triangle_max_raws = xmlGetInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/maxRaws");
		
		/// printf("\n\nMAX PASCALTRIANGLE RAWS: %hu\n\n", tmp->pascal_triangle_max_raws);
	
	    volatile bool tmp_bool = false;
		
		#pragma omp parallel for
	    for(i=0; i<MAX_BOOL_SETTINGS; ++i)
	    {
			char name[4*MIN_STRING] = NULL_CHAR;
			char strboolized[2*MIN_STRING] = NULL_CHAR;
	        strboolize(suite_c.bools_names[i], strboolized);
	        sprintf(name, "/settings/booleanKeys/%s", strboolized);
	        tmp_bool = xmlGetBool(&xpathObj, xpathCtx, name);
	        if(tmp_bool)
	            tmp->bools |= suite_c.bools[i].bmask;
	    }
	    
	    xmlExit(path, &doc, &xpathObj, &xpathCtx);
	
	    return;
	}
	
	
#endif


#if WINOS

	#ifdef XMLCALL
		__MSNATIVE_ XMLCALL void __system _colFileLoader(const char path[static MAX_PATH_LENGTH])
	    {
	        static bool once_executed = false;
	
	        if(once_executed)
			{
				if(isSett(BOOLS_ITEMSAUTOSAVING))
	            	_backupColFile();
			}
			else
	        	once_executed = true;
	        
	        xmlXPathContext * xpathCtx = (xmlXPathContext*) NULL;
	        xmlXPathObject * xpathObj = (xmlXPathObject*) NULL;
	        xmlDoc * doc = xmlInit(access(colors_path), &xpathCtx);
			        
	        COLOR_DEFAULT = xmlGetInt(&xpathObj, xpathCtx, "/colors/defaultColor");
	        COLOR_ERROR = xmlGetInt(&xpathObj, xpathCtx, "/colors/errorsColor");
	        COLOR_CREDITS = xmlGetInt(&xpathObj, xpathCtx, "/colors/creditsColor");
	        COLOR_USER = xmlGetInt(&xpathObj, xpathCtx, "/colors/userColor");
	        COLOR_SYSTEM = xmlGetInt(&xpathObj, xpathCtx, "/colors/systemColor");
	        COLOR_AUTHOR = xmlGetInt(&xpathObj, xpathCtx, "/colors/authorColor");
	        
			xmlExit(access(colors_path), &doc, &xpathObj, &xpathCtx);
	        return;
	    }
	    
	    __MSNATIVE_ XMLCALL void __system _backupColFile(void)
		{
			xmlXPathContext * xpathCtx = (xmlXPathContext*) NULL;
			xmlXPathObject * xpathObj = (xmlXPathObject*) NULL;
			xmlDoc * doc = xmlInit(access(colors_path), &xpathCtx);
			
			xmlWriteInt(&xpathObj, xpathCtx, "/colors/defaultColor", COLOR_DEFAULT);
			xmlWriteInt(&xpathObj, xpathCtx, "/colors/errorsColor", COLOR_ERROR);
			xmlWriteInt(&xpathObj, xpathCtx, "/colors/creditsColor", COLOR_CREDITS);
			xmlWriteInt(&xpathObj, xpathCtx, "/colors/userColor", COLOR_USER);
			xmlWriteInt(&xpathObj, xpathCtx, "/colors/systemColor", COLOR_SYSTEM);
			xmlWriteInt(&xpathObj, xpathCtx, "/colors/authorColor", COLOR_AUTHOR);
			
			xmlExit(access(colors_path), &doc, &xpathObj, &xpathCtx);
		    return;
		}
	#endif
    
    __MSUTIL_ __WINCALL inline BOOL WINAPI __system __export SetExitButtonState(const bool state)
	{
		return DeleteMenu(GetSystemMenu(GetConsoleWindowNT(),state),6,MF_BYPOSITION);
	}
	
    __MSUTIL_ inline const char * const __system __export getFilename(const char path[static MAX_PATH_LENGTH])
    {
        const char * const stkr = strrchr(path, '\\');
        return stkr ? stkr+1 : path;
    }

    __MSUTIL_ __WINCALL HWND WINAPI __system __export GetConsoleWindowNT()
    {

        // declare function pointer type

        typedef HWND WINAPI (*GetConsoleWindowT)(void);

        // declare one such function pointer

        GetConsoleWindowT GetConsoleWindow;

        // get a handle on kernel32.dll

        const HMODULE hK32Lib = GetModuleHandle(TEXT("KERNEL32.DLL"));

        // assign procedure address to function pointer

        GetConsoleWindow = (GetConsoleWindowT)GetProcAddress(hK32Lib,TEXT("GetConsoleWindow"));

        // check if the function pointer is valid

        // since the function is undocumented

        if ( GetConsoleWindow == NULL )
            return NULL;
        // call the undocumented function

        return GetConsoleWindow();

    }

    __MSUTIL_ __WINCALL bool __system __export windowsFileHandler(char *name, const char *ext_pattern, const char default_ext[static MAX_EXTENSION_LENGTH], bool mode)
    {
        bool result;

        result = false;
            // HWND hwnd;
        OPENFILENAME ofn;
        char szFileName[MAX_PATH] = NULL_CHAR;
        // const bool assert = __pmode__ == ENVS_OPEN;

        ZeroMemory(&ofn, sizeof(ofn));

        ofn.lStructSize = sizeof(ofn); // GUARDATE LA NOTA SOTTO

        ofn.lpstrFilter = ext_pattern; // SPERIMENTALE

        // ofn.hwndOwner = hwnd;

        ofn.lpstrFile = szFileName;
        ofn.nMaxFile = MAX_PATH;
        ofn.Flags = mode ? (OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY) :
                    (OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT);
        ofn.lpstrDefExt = default_ext;

        if((result = mode ? GetOpenFileName(&ofn) : GetSaveFileName(&ofn)))
            strcpy(name, ofn.lpstrFile);

        return result;
    }
    
#endif

__MSUTIL_ inline void __system updInfo(void)
{

    char title[MAXX_STRING];

    sprintf(title, PROG__NAME" - [ %s ] - [ %s ] - [ %s ] - (%s) - [ %s ]", getItemsListNo(ENVS) ? getFilename(listNo(access(lists)[ENVS].cur_item, ENVS)->path) : NULL_CHAR, getItemsListNo(MATRICES) ? getFilename(listNo(access(lists)[MATRICES].cur_item, MATRICES)->path) : NULL_CHAR,
                    getItemsListNo(LOGS) ? getFilename(listNo(access(lists)[LOGS].cur_item, LOGS)->path) : NULL_CHAR, access(sysLogPath), getItemsListNo(LAYOUTS) ? getFilename(listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path) : NULL_CHAR);

	#if WINOS
	    if(!SetConsoleTitle(title))
	        printErr(22, "SetConsoleTitle failed with error: %lu", GetLastError());
	#else
		printf("%c]0;%s%c", '\033', title, '\007');
	#endif

    return;
}

__MSUTIL_ void __system __export SetColor(const sel_typ ForgC)
{
	#if WINOS
	    const HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	    CONSOLE_SCREEN_BUFFER_INFO csbi;
	
	    if(GetConsoleScreenBufferInfo(hStdOut, &csbi))
	        SetConsoleTextAttribute(hStdOut, (csbi.wAttributes & 0xF0) + (ForgC & 0x0F));
	#endif
	
	// Linux Part WIP
	// actually this results in a NOP

    return;
}

__MSUTIL_ inline bool __export min_cmpfunc(const register ityp a, const register ityp b)
{
    return (a<b);
}

__MSUTIL_ inline bool __export max_cmpfunc(const register ityp a, const register ityp b)
{
    return (a>b);
}

__MSNATIVE_ ityp __export MINMAX(const register dim_typ dim, const ityp vector[static dim], const bool mode, dim_typ *idx)
{
    ityp tmp;
    dim_typ i;
    bool (* const cfunc)(const register ityp, const register ityp) = mode ? max_cmpfunc : min_cmpfunc;
    if(idx) (*idx) = 0;
    for(i=1,tmp=vector[0]; i<dim; ++i)
        if(cfunc(vector[i], tmp))
        {
            tmp = vector[i];
            if(idx) (*idx) = i;
        }


    return tmp;
}

__MSNATIVE_ inline bool __system __export isDomainForbidden(ityp val, bool mode)
{
    if(TYPE_DOMAIN(val))
    {
        printf2(COLOR_ERROR, "\n%sPUT OVERFLOW ERROR.\n\n", mode ? "OUT":"IN");
        return true;
    }
    return false;
}

__MSUTIL_ inline void _MS__private __system __export free_foreach(ityp ***matrix, const dim_typ algebra_units, const dim_typ dimraws, bool mode)
{
    if(mode)
        return;
    #pragma omp parallel for num_threads(algebra_units)
    for(dim_typ i=0; i<algebra_units; ++i)
        matrixFree(&matrix[i], dimraws);
    return;
}

__MSUTIL_ inline void _MS__private __system __export free_foreach2(ityp ***matrix, const dim_typ dim)
{
	#pragma omp parallel for
    for(dim_typ i=0; i<dim; ++i)
        matrixFree(&matrix[i], dim);
    return;
}

__MSUTIL_ inline void _MS__private __system __export free_foreach3(ityp ***matrix, const dim_typ algebra_units)
{
	#pragma omp parallel for num_threads(algebra_units)
    for(dim_typ i=0; i<algebra_units; ++i)
        free(matrix[i]);
    return;
}


__MSUTIL_ inline bool __system __export checkErrMem(const void * pntr)
{
    if(!pntr)
    {
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        printErr(12, "An error occurred during Heap Dynamic Memory Allocation");
        sprint("\n(Sub)Program Terminating...\n\n");
        system("PAUSE");
        return true;
    }
    return false;
}

__MSNATIVE_ bool __system __export matrixAlloc(ityp ***matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    if(!(*matrix))
        (*matrix) = NULL;

    (*matrix) = calloc(sizeof(ityp *), dim[RAWS]);
    errMem((*matrix), false);

    dim_typ i;

    for(i=0; i<dim[RAWS]; ++i)
    {
        (*matrix)[i] = calloc(sizeof(ityp), dim[COLUMNS]);
        errMem((*matrix)[i], false);
    }

    return true;
}

__MSNATIVE_ void __system __export _matrixFree(ityp ***matrix, dim_typ raws, bool mode)
{
    if(mode || !(*matrix))
        return;

    dim_typ i;
    #pragma omp parallel for
    for(i=0; i<raws; ++i)
        free((*matrix)[i]);
    free((*matrix));
    (*matrix) = NULL; // to avoid dangling references,
    // even if all this pointer passed to this function
    // are generally locally allocated into static void functions.
    return;
}

__MSNATIVE_ bool __system __export equalMatrix(ityp ***matrix1, ityp **matrix2, const register dim_typ dim[static MAX_DIMENSIONS])  
{
    (*matrix1) = realloc((*matrix1), sizeof(ityp*)*dim[RAWS]);
    errMem((*matrix1), false);

    dim_typ i;

    for(i=0; i<dim[RAWS]; ++i)
    {
        (*matrix1)[i] = realloc((*matrix1)[i], sizeof(ityp)*dim[COLUMNS]);
        errMem((*matrix1)[i], false);
    }

    dim_typ j;

    // Phisically equalling matrix1 values to matrix2 ones
    #pragma omp parallel for
    for(i=0; i<dim[RAWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
            (*matrix1)[i][j] = matrix2[i][j];

    return true;
}

__MSNATIVE_ _MS__private inline void __system resetLmpMatrix(void)
{
	access(lmpMatrix)->matrix = NULL;
    access(lmpMatrix)->dim[RAWS] = STARTING_LMP_RAWS;
    access(lmpMatrix)->dim[COLUMNS] = STARTING_LMP_COLUMNS;
    return;
}

__MSNATIVE_ inline void __system __export _flushMemoizersBuffers(sel_typ mode)
{

    free(access(sysMem)[mode].memoizer);
    access(sysMem)[mode].memoizer = NULL;
    access(sysMem)[mode].current_max_index = 0;
    sprint("%s Function Memoizer has been properly flushed.\n\n", suite_c.memoizers_names[mode]);
    return;
}

__MSNATIVE_ inline void __system __export flushAllMemoizersBuffers(void)
{
	for(dim_typ i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		_flushMemoizersBuffers(i);
		
	return;
}

__MSNATIVE_ dim_typ __system __export selectListItem(dim_typ dim, bool mode, const char *string, const char list[static dim][MIN_STRING])
{
    dim_typ item;
    dim_typ i;

    PRINTN();
    printf2(COLOR_CREDITS, string);
    printf2(COLOR_CREDITS, ":\n");

    PRINTL();

    if(mode)
    {
        for(i=0; i<dim; ++i)
            printf("- %hu: %s;\n", i, list[i]); // ext_math.functions[i].name);
        printf("- %hu: Back to previous SubProgram.\n", i);

        PRINTL();

        ityp item2;

        while((PARSING_SYSTEM_ALLOWED ? (isNullVal((item2 = requires(NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS)))) :
            (!scanf2(1, INPUT_CONVERSION_FORMAT, &item2))) || item2 != (item = (dim_typ)item2) || item < 0 || item > dim)
            printErr(1, "Invalid inserted value");

        CLEARBUFFER();
    }
    else
    {
        for(i=0; i<dim; ++i)
            printf("- %c: %s;\n", i+'A', list[i]);
        printf("- %c: Back to previous SubProgram.\n", i+'A');

        PRINTL();

        sel_typ item2;

        do item2 = toupper(getch());
        while(item2 < 'A' || item2 > dim+'A');

        item = item2-'A';
    }

    return item;
}

__MSNATIVE_ void __system viewProgramSettings(dim_typ which_layout)
{
    printf2(COLOR_SYSTEM, "\nCurrent Settings Layout:\n\n");

	dim_typ i;
    nodelist * const tmp = listNo(which_layout, LAYOUTS);
    layoutObj * const cur_layout = ((layoutObj *)(tmp->data));

    PRINTL();

	sprint("- Exit CHAR: %c;\n", cur_layout->exit_char);
    sprint("- PROGRAM PRECISION: %hu;\n", cur_layout->precision);
    sprint("- STABILIZER FACTOR: %hu;\n", cur_layout->stabilizer_factor);
    sprint("- Min Stirling-Requires NUMBER: %hu;\n", cur_layout->min_stirlingrequires_number);
    sprint("- ALGEBRA: %s;\n", suite_c.algebra_elements_names[cur_layout->algebra]);
    sprint("- Outlier Constant: %.*f;\n", DEFAULT_PRECISION, cur_layout->outlier_constant);

    sprint("- Matrices MAX RAWS: %hu;\n", cur_layout->matrix_max_raws);
    sprint("- Matrices MAX COLUMNS: %hu;\n", cur_layout->matrix_max_columns);
    sprint("- Max DSVD Iterations: %hu;\n", cur_layout->max_dsvd_iterations);
    sprint("- Max SIMPLEX METHOD Iterations: %hu;\n", cur_layout->max_simplex_iterations);


	for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
	{
		char str[SIGN_STRING] = NULL_CHAR;
		strcpy(str, suite_c.memoizers_names[i]);
		toupper_s(str);
		sprint("- Max %s Memoizable Index: %hu;\n", str, cur_layout->max_memoizable_indeces[i]);
	}

    sprint("- MIN Processable BASE %hu;\n", cur_layout->basecalc_minbase);
    sprint("- MAX Processable BASE: %hu;\n", cur_layout->basecalc_maxbase);
    sprint("- MAX PROCESSABLE Decimal2Bin NUM: %d;\n", cur_layout->max_changebase_binary_convnum);

    sprint("- NEWTON DIFFTABLES Minimum Dimension: %hu;\n", cur_layout->min_newton_difftables_dim);
    sprint("- NEWTON DIFFTABLES Maximum Dimension: %hu;\n", cur_layout->max_newton_difftables_dim);

    sprint("- MIN Processable ROMAN NUMBER: %hu;\n", cur_layout->min_roman_number);
    sprint("- MAX Processable ROMAN  NUMBER: %hu;\n", cur_layout->max_roman_number);

    sprint("- Pascal's Triangle MIN RAWS: %hu;\n", cur_layout->pascal_triangle_min_raws);
    sprint("- Pascal's Triangle MAX RAWS: %hu.\n", cur_layout->pascal_triangle_max_raws);

    printf2(COLOR_USER, "\n*** Boolean Settings ***\n\n");

    PRINTL();

    for(i=0; i<MAX_BOOL_SETTINGS; ++i)
        sprint("- %s: %s;\n", suite_c.bools_names[i], (cur_layout->bools & suite_c.bools[i].bmask) == suite_c.bools[i].bmask ? "Enabled":"Disabled");

    PRINTL();

    return;
}

__MSNATIVE_ void __system setProgramSettings(dim_typ which_layout)
{
    nodelist * item_data = listNo(which_layout, LAYOUTS);
    layoutObj * const tmp = malloc(sizeof(layoutObj));

    errMem(tmp, VSPACE);

    resetProgramSettings(tmp, item_data->path);

    item_data->data = tmp;

    return;
}

__MSNATIVE_ inline bool __system __export catchPause()
{
    signal(SIGINT, (__p_sig_fn_t) sigexit);
    if(access(sigresult))
    {
        printf("\nPress any key to continue Program Action.\n");
        printf("or %c to stop Execution.\n", access(curLayout)->exit_char);
        access(sigresult) = false;
        if(getch() == access(curLayout)->exit_char) return true;
    }

    return false;
}

__MSNATIVE_ inline void __system logPrint(logObj * const which_log)
{
    printf("\nSelected Log Content:\n\n");
    SHOWPAUSEMESSAGE();

    size_t i;

    for(i=0; i<which_log->buflen; ++i)
    {
        putchar(which_log->buffer[i]);
        if(catchPause()) return;
    }

    PRINTL();
    PRINTN();

    return;
}

__MSNATIVE_ inline void __system logWrite(FILE *fpnt, logObj * const which_log)
{

    if(fpnt == stdout)
        logPrint(which_log);
    else
    {
        size_t i;
        for(i=0; i<which_log->buflen; ++i)
            putc(which_log->buffer[i], fpnt);
        _flushLogBuf(which_log);
    }

    return;

}

// void getVarList(FILE *fpnt);
__MSNATIVE_ void __system getVarListEx(FILE *fpnt, dim_typ which_env)
{
    void *cookie = NULL;
    EXPRTYPE vval;
    char *vname;

    if(fpnt == stdout)
        printf("\nVariable list items:\n\n");

    exprValList * const tmp = ((exprValList *)((exprType*)(listNo(which_env, ENVS)->data))->var_list);

    cookie = exprValListGetNext(tmp, &vname, &vval, NULL, NULL);
    // cookie = exprValListGetNext(vlist, &vname, &vval, NULL, NULL);

    do
    {

        if(fpnt != stdout && (!strcmp(vname, DEFAULT_ENVS_ANSVALNAME)||!strcmp(vname, "global")))
            break;

        // fputs("\n", fpnt);
        char str[MIN_STRING];
        // sprintf(str,"%s = %.*f;\n", vname, suite.precision, vval);
        sprintf(str,"%s = %.*f;\n", vname, access(curLayout)->precision, vval);
        fputs(str, fpnt);
        cookie = exprValListGetNext(tmp/*suite.exprVars->var_list*/, &vname, &vval, NULL, cookie);
    }
    while(cookie);

    return;
}

__MSNATIVE_ ityp __system __export requires(const char *cmd_string, const char *string, const char *res_string, const unsigned options)
{

    char buf[MAX_BUFSIZ];
    exprObj *e = INIT_OBJLIST;
    volatile int err;
    ityp val;
    ityp diff;
    time_t t1;
    jmp_buf jumper;
    int start, end;


    /* Set error buffer */
    err = setjmp(jumper);
    if(err && e)
        exprFree(e);


    /* Gets an expression */
    if(cmd_string)
        strcpy(buf, cmd_string);
    else
    {
        printf2(COLOR_CREDITS, string);
        PRINT2N();
        gets(buf);
    }

    access(exitHandle) = INVALID_EXITHANDLE;
    signal(SIGINT, (__p_sig_fn_t) sigproc); // RIPRISTINO HANDLER ALLA SIG_DFL

    // To check whether User wants to exit from requiringINPUT or Not.
    if(access(sigresult) || (!strcmp(buf, MATRIXGET_COMMAND)) || (!strcmp(buf, MATRIXSET_COMMAND)) || (!strcmp(buf, MATRIXBACK_COMMAND)) || !buf[0])
    {

        if(access(sigresult) || !strcmp(buf, MATRIXGET_COMMAND))
        {
            if(!strcmp(buf, MATRIXGET_COMMAND))
                sigproc();
            access(exitHandle) = EXITHANDLE_GETCMD;
        }

        else if(!strcmp(buf, MATRIXSET_COMMAND))
            access(exitHandle) = EXITHANDLE_SETCMD;

        else if(!strcmp(buf, MATRIXBACK_COMMAND))
            access(exitHandle) = EXITHANDLE_BACKCMD;

        if(!buf[0])
            access(exitHandle) = EXITHANDLE_EXIT;

        access(sigresult) = false;
        // signal(SIGINT, (__p_sig_fn_t) sigproc); // REIMPOSTO HANDLER DEL SEGNALE SIGINT ALLA MIA FUNZIONE
        return NULL_VAL;
    }

    if(buf[strlen(buf)-1] != TERMINATING_CHAR)
        strcat(buf, TERMINATING_STRING);

    PRINTN();

    // Creates expression
    err = exprCreate(&e, access(func_list), access(exprVars)->var_list, access(const_list), NULL, 0);
    if(err != EXPR_ERROR_NOERROR)
        {
        sprint("Expr Creation Error.\n");
        longjmp(jumper, err);
        }

    /* Parse expr */
    err = exprParse(e, buf);
    if(err != EXPR_ERROR_NOERROR)
        {
        exprGetErrorPosition(e, &start, &end);
        sprint("Parse Error (%d,%d).\n", start, end);

        longjmp(jumper, err);
        }

    /* Eval expression */
    const bool assert = (isSett(BOOLS_SHOWDIFFTIME) && (options & PARSER_SHOWDIFFTIME) == PARSER_SHOWDIFFTIME);

    if(assert)
        t1 = time(NULL);

    err = exprEval(e, &val);

    if(err != EXPR_ERROR_NOERROR)
        sprint("Eval Error: %d.\n", err);

    if((options & PARSER_SHOWRESULT) == PARSER_SHOWRESULT)
    {

        PRINTN();
        printf2(COLOR_USER, res_string);
        printf2(COLOR_USER, ": ");
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, val);
        printf2(COLOR_USER, ".\n\n");

    }

    if(assert)
    {
        PRINTL();
        diff = difftime(time(NULL), t1);
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", DEFAULT_PRECISION, (EXPRTYPE)diff);
        PRINTL();
    }

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS) && (options & PARSER_SAVERESULT) == PARSER_SAVERESULT)
        *(access(exprVars)->e_ANS) = val;



    /* Enumerate the items in the variable list */

    /// const bool ols_assert = getItemsListNo(ENVS) != STARTING_ENVSNO;

    // does this has to be called each time we parse and evaluate a new expression???
    // this depends on BOOLS_RESETLISTS boolean suite structure-capsulated bit-fielded macrovar
    if(isSett(BOOLS_RESETLISTS)) // here must be written these piece of code
        refreshExprEvalLists();
    // freeExprEvalLists();

    if(getItemsListNo(ENVS) != STARTING_ENVSNO)
    {
        if(isSett(BOOLS_SHOWVARLIST) && (options & PARSER_SHOWVARLIST) == PARSER_SHOWVARLIST)
            getVarList(stdout);
        if(lazy_exec && isSett(BOOLS_ITEMSAUTOSAVING))
            saveItem(access(lists)[ENVS].cur_item, ENVS);
    }

    exprFree(e);


    /* We are done */
    // longjmp(jumper, -1);

    // signal(SIGINT, (__p_sig_fn_t) sigproc); // REIMPOSTO HANDLER DEL SEGNALE SIGINT ALLA MIA FUNZIONE
    return val;
}

__MSNATIVE_ bool __system insertDims(dim_typ *righe, dim_typ *colonne)
{
    const dim_typ old_dims[MAX_DIMENSIONS] =
    {
        (*righe),
        (*colonne)
    };

    printf2(COLOR_CREDITS, "\nEnter RAWS and COLUMNS as expected:\n\
[RAWS%sCOLUMNS].\n", PARSING_SYSTEM_ALLOWED ? "]\n[" : " ");

    if(PARSING_SYSTEM_ALLOWED)
        PRINTHOWTOBACKMESSAGE();

    ityp tmp;
    ityp tmp2;

    tmp = tmp2 = 0.00;

    // uint64_t tmp3[MAX_DIMENSIONS];

    while(PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted RAWS", PARSER_SHOWRESULT))) || tmp != ((*righe) = (dim_typ)tmp) || (*righe) < 1 || (*righe) > access(curLayout)->matrix_max_raws) ||
        (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted COLUMNS", PARSER_SHOWRESULT))) || tmp != ((*colonne) = (dim_typ)tmp) || (*colonne) < 1 || (*colonne) > access(curLayout)->matrix_max_columns) :
        (!scanf2(2, "%lf %lf", &tmp, &tmp2)) || tmp != ((*righe) = (dim_typ)tmp) || tmp2 != ((*colonne) = (dim_typ)tmp2) || (*righe) < 1 || (*colonne) < 1 || (*righe) > access(curLayout)->matrix_max_raws || (*colonne) > access(curLayout)->matrix_max_columns)
    {
        CLEARBUFFER();
        if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
        if(exitHandleCheck) // if(tmp3[RAWS] == NULL_VAL || tmp3[COLUMNS] == NULL_VAL)
        {
            (*righe) = old_dims[RAWS];
            (*colonne) = old_dims[COLUMNS];
            return false;
        }
        printErr(33, "Invalid [RAWS COLUMNS] format.\nYou have to insert non-negative RAWS and COLUMNS,\n\
and must be respectively less than: %hu and %hu", access(curLayout)->matrix_max_raws, access(curLayout)->matrix_max_columns);
    }
    return true;
}

__MSNATIVE_ bool __system insertDim(dim_typ *dim, bool mode)
{
    const dim_typ old_dim = (*dim);
    dim_typ max_dim;

    PRINTN();

    if(mode > COLUMNS)
    {
        max_dim = MAX_RIGHE_PER_COLONNE;
        printf2(COLOR_CREDITS, "Enter Quad Matrix DIMENSION.");
    }
    else
    {
        max_dim = mode ? access(curLayout)->matrix_max_columns : access(curLayout)->matrix_max_raws;
        printf2(COLOR_CREDITS, "Enter Matrix %s.", mode ? "COLUMNS":"RAWS");
    }

    if(PARSING_SYSTEM_ALLOWED)
        PRINTHOWTOBACKMESSAGE();

    ityp tmp;

    // uint64_t tmp2;

    while((PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted DIMENSION is", PARSER_SHOWRESULT)))) :
          (!scanf2(1, INPUT_CONVERSION_FORMAT, &tmp))) || tmp != ((*dim) = (dim_typ)tmp) || (*dim) < 1 || (*dim) > max_dim)
    {
        CLEARBUFFER();
        if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
        if(exitHandleCheck)
        {
            (*dim) = old_dim;
            return false;
        }
        printErr(33, "Invalid inserted Value.\nMust be a non-negative integer less than: %hu", max_dim);
    }

    // (*dim) = tmp2;



    return true;
}

// La seguente funzione sarebbe stata la funzione Handler
// del segnale SIGINT. Il problema e' che dovrebbe essere
// sempre e continuamente richiamata la funzione signal
// per far sì che funzioni correttamente. Meglio evitare.


__MSNATIVE_ void __system sigproc(void)
{
    sigexit();
    /// signal(SIGINT, (__p_sig_fn_t) sigproc);
    if(getItemsListNo(MATRICES) != STARTING_MATNO && access(curMatrix)->matrix && access(lmpMatrix)->matrix)
    {
        if(isSett(BOOLS_ITEMSAUTOSAVING) && !saveItem(access(lists)[MATRICES].cur_item, MATRICES))
            return;
        access(curMatrix)->matrix = access(lmpMatrix)->matrix;
        access(curMatrix)->dim[RAWS] = access(lmpMatrix)->dim[RAWS];
        access(curMatrix)->dim[COLUMNS] = access(lmpMatrix)->dim[COLUMNS];
        // equalMatrix(&suite.current_matrix, suite.last_matrix_printed, suite.lmp_raws, suite.lmp_columns);
        sprint("\nCurrent Matrix has been equalled\nto Last Matrix Printed.\n\n");
    }
    return;
}

__MSNATIVE_ inline void __system sigexit(void)
{
    access(sigresult) = true;
    return;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ void showNewtonDifferenceTable(dim_typ n, ityp x[static n], ityp y[access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool mode)
{

    printf2(COLOR_SYSTEM, "\n***********%s Difference Table ***********\n", mode ? "FORWARD" : "BACKWARD");

    //display Difference Table

    dim_typ i, j;

    for(i=0;i<n;++i)

    {

        PRINTT();
        printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, x[i]);
        printf2(COLOR_USER, "\t%.*f", access(curLayout)->precision, x[i]);

        for(j=0; (mode ? (j<(n-i)):(j<=i)) ;++j)
        {
            PRINTT();
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, y[i][j]);
        }

        PRINTN();

    }

    return;
}

__MSNATIVE_ bool __system insertNMMatrix(ityp ***matrix, const register dim_typ dim[static MAX_DIMENSIONS])  
{
    if(!matrixAlloc(matrix, dim))
    {
        #if WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        return false;
    }

    volatile char tmp=0;
    dim_typ i, j;
    dim_typ start_col_index;


    // we must seek for backtracking...
    // BY Inserting somewhere a 'return false' statement.
    // (FINAL BACKTRACKING RESPONSE)
    for(i=start_col_index=0; tmp != -1 && i<dim[RAWS]; ++i)
        for(j=start_col_index; tmp != -1 && j<dim[COLUMNS]; ++j)
        {
            while((tmp = insertElement((*matrix), (dim_typ2){i, j}, false)) != 1 && tmp != -2 && !(char_insert))
                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                    if(access(curMatrix)->dim[RAWS] != dim[RAWS] || access(curMatrix)->dim[COLUMNS] != dim[COLUMNS])
                        printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Raws and %hu Columns", dim[RAWS], dim[COLUMNS]);
                    else
                    {
                        if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                        {
                            #if WINOS
                                SetExitButtonState(DISABLED);
                            #endif // WINOS
                            return false;
                        }
                        sprint("\nYou're correctly using Current Matrix.\n\n");
                        break;
                    }
                else
                {
                    if(!i)
                    {
                        matrixFree(matrix, dim[RAWS]);
                        #if WINOS
                            SetExitButtonState(ENABLED);
                        #endif // WINOS
                        return false;
                    }
                    printErr(1, "You cannot enter characters");
                }

            if(tmp == -2)
            {
                if(j > 0)
                    j -= 2;
                else if(i > 0)
                {
                    i -= 2;
                    start_col_index = dim[COLUMNS]-1;
                    break;
                }
                else
                {
                    matrixFree(matrix, dim[RAWS]);
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    return false;
                }
            }

        }

    return true;
}

__MSNATIVE_ volatile char __system insertElement(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS], bool square)
{

    if(square && (dim[RAWS] == MAX_RIGHE_PER_COLONNE || dim[COLUMNS] == MAX_RIGHE_PER_COLONNE))
    {
        printErr(33, "MAX RAWS per COLUMNS Reached");
        return 0;
    }

    if(dim[RAWS] == access(curLayout)->matrix_max_raws)
    {
        printErr(33, "MAX RAWS Reached");
        return 0;
    }

    if(dim[COLUMNS] == access(curLayout)->matrix_max_columns)
    {
        printErr(33, "MAX COLUMNS Reached");
        return 0;
    }

    // PRINTN();
    printf2(COLOR_CREDITS, "\nEnter [%hu,%hu] Matrix Element.\n", dim[RAWS]+1, dim[COLUMNS]+1);

    sel_typ tmp;
    tmp = 1;
    // tmp = true;

    do
    {
        CLEARBUFFER();

        char str[MIN_STRING];

        sprintf(str, "[%hu,%hu] Matrix Element correctly inserted", dim[RAWS]+1, dim[COLUMNS]+1);


        // tmp = scanf(INPUT_CONVERSION_FORMAT, &matrix[riga][colonna]);
        // if((matrix[riga][colonna] = requires(NULL, str, restr, true, false, false, true, false)) == NULL_VAL) tmp=false;
        // if(tmp == NULL_VAL) tmp = 0;

        // FARE ATTENZIONE PER IL GETSET SYSTEM CARATTERIZZATO DAI VALORI DI RITORNO NULL_VALN della funzione requires(NULL, ...)


        tmp = PARSING_SYSTEM_ALLOWED ? (!((matrix[dim[RAWS]][dim[COLUMNS]] = requires(NULL, NULL_CHAR, str, PARSER_SHOWRESULT)) == NULL_VAL)) :
                scanf2(1, INPUT_CONVERSION_FORMAT, &matrix[dim[RAWS]][dim[COLUMNS]]);// , printf("%s: ", str), printf(OUTPUT_CONVERSION_FORMAT, matrix[riga][colonna]), printf(".\n\n");




        CLEARBUFFER();

        if(access(exitHandle) == EXITHANDLE_SETCMD)
            return -1;

        if(access(exitHandle) == EXITHANDLE_BACKCMD)
            return -2;

        if(access(exitHandle) == EXITHANDLE_EXIT && (!dim[RAWS]) && (!dim[COLUMNS]))
            return -3;

        // tmp = PARSING_SYSTEM_ALLOWED ? (!(isNullVal((matrix[riga][colonna] = requires(NULL, NULL_CHAR, str, true, false, false, true, false))))) :
            // scanf(INPUT_CONVERSION_FORMAT, &matrix[riga][colonna]);// , printf2("%s: ", str), printf2(OUTPUT_CONVERSION_FORMAT, matrix[riga][colonna]), printf(".\n\n");

        if(isnSett(BOOLS_MATRIXPARSING))
        {
            printf2(COLOR_USER, "%s: ", str);
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, matrix[dim[RAWS]][dim[COLUMNS]]);
            printf2(COLOR_USER, ".\n\n");
        }

        if(!(matrix[dim[RAWS]][dim[COLUMNS]]) && dcheck && INVERSE_OPS && ((__pmode__ > 5 && __pmode__ < 11)||__pmode__ == ALGEBRA_SCALARDIVISIONMATRIX))
           printErr(33, "You cannot enter a 0 because program is performing a Division somewhere");

    }
    while((tmp != 1 && !(dim[RAWS]) && !(dim[COLUMNS])) || (!(matrix[dim[RAWS]][dim[COLUMNS]]) && dcheck && INVERSE_OPS && ((__pmode__ > 5 && __pmode__ < 11)
    ||__pmode__ == ALGEBRA_SCALARDIVISIONMATRIX))||isDomainForbidden(matrix[dim[RAWS]][dim[COLUMNS]], INPUT));

    CLEARBUFFER();

    return tmp;
}

__MSNATIVE_ inline volatile bool __system checkBackTracking(volatile char tmp, dim_typ *colonna)
{
    if(tmp < -1)
    {
        if((*colonna) > 0)
            (*colonna) -= 2;
        else
            return false;
    }
    return true;
}

__MSNATIVE_ inline volatile __system sel_typ checkBackTracking2(volatile char tmp, dim_typ *riga, dim_typ *colonna, dim_typ *start_col_index, dim_typ colonne)
{
    if(tmp == -2)
    {
        if((*colonna) > 0)
        {
            (*colonna) -= 2;
            return 0;
        }
        else if((*riga) > 1)
        {

            (*riga) -= 2;
            (*start_col_index) = colonne-1;
            return 1;
        }
        return 2;
    }

    return 3;
}

__MSNATIVE_ bool __system __export enterMatrix(ityp ***matrix, dim_typ *righe, dim_typ *colonne, bool square, bool view)
{

    printf2(COLOR_CREDITS, "\n\nEnter%s Matrix", square ? " Quad" : "");
    printf2(COLOR_CREDITS, " by following procedure you'll see.\n");

    volatile char tmp;
    dim_typ start_col_index;

    tmp = 1;
    start_col_index = 0;

    #if WINOS
        SetExitButtonState(DISABLED);
    #endif // WINOS


    if(isSett(BOOLS_INSERTMODE))
    {
        printf2(COLOR_CREDITS, "And when you reach desired raws %s columns dimensions, %s.\n\n",
               square ? "=":"and", isSett(BOOLS_MATRIXPARSING) ? "press ENTER":"insert an alphanumeric value");

        (*matrix) = malloc(2*sizeof(ityp *));
        errMem((*matrix), false);


        (*matrix)[0] = malloc(sizeof(ityp)*(lazy_exec ? 1 : access(curLayout)->stabilizer_factor));

        (*righe) = 1;

        for(*colonne = 0; tmp; ++ (*colonne))
        {

            if(tmp < -1)
            {
                // ONE STEP FORWARD AND TWO STEP BACK
                // BACKTRACKING FORMULA!!!
                if((*colonne) > 1)
                    (*colonne) -= 2;
                else
                {
                    matrixFree(matrix, (*righe));
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    return false;
                }
            }

            const dim_typ analog_columns = (*colonne)+1;
            if(lazy_exec)
                (*matrix)[0] = realloc((*matrix)[0], sizeof(ityp)*analog_columns);
            else if(!(analog_columns % access(curLayout)->stabilizer_factor))
                (*matrix)[0] = realloc((*matrix)[0], sizeof(ityp)*((*colonne) + access(curLayout)->stabilizer_factor));

            errMem((*matrix)[0], false);

            if((tmp = insertElement((*matrix), (dim_typ2){0, (*colonne)}, square)) == -1 && getItemsListNo(MATRICES) != STARTING_MATNO)
            {
                if(square && access(curMatrix)->dim[RAWS] != access(curMatrix)->dim[COLUMNS])
                {
                    printErr(1, "You cannot use Current Matrix because\nit isn't a Quad one");
                    (*colonne) --;
                }
                else
                {
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                    {
                        matrixFree(matrix, (*righe));
                        return false;
                    }
                    sprint("\nYou are correctly using Current Matrix.\n\n");
                    (*righe) = access(curMatrix)->dim[RAWS];
                    (*colonne) = access(curMatrix)->dim[COLUMNS];
                    return true;
                }
            }
        }

        tmp = 1;
        (*colonne) --;

        (*matrix)[0] = realloc((*matrix)[0], sizeof(ityp)*(*colonne));
        errMem((*matrix)[0], false);


        for(*righe = 1; (square ? *righe < *colonne : (tmp && __pmode__ != ALGEBRA_SCALARPRODUCT)); ++(*righe))
        {

            dim_typ i;
            const dim_typ analog_raws = (*righe)+1;

            if(lazy_exec)
                (*matrix) = realloc((*matrix), sizeof(ityp *)*analog_raws);
            else if(!(analog_raws % access(curLayout)->stabilizer_factor))
                (*matrix) = realloc((*matrix), sizeof(ityp *) * ((*righe)+access(curLayout)->stabilizer_factor));

            errMem((*matrix), false);

            (*matrix)[*righe] = malloc((*colonne)*sizeof(ityp));
            errMem((*matrix)[*righe], false);

            for(i=start_col_index; tmp && i<*colonne; ++i)
            {
                /// if(tmp == -2) return false; // simple backtracking on second or major raws...
                // BACKTRACKING ON SECOND or MAJOR RAWS EXPERIMENTAL FORMULA

                while((tmp = insertElement((*matrix), (dim_typ2){*righe, i}, square)) != 1 && tmp != -2 && !(char_insert) && i)
                    if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                        if(square && access(curMatrix)->dim[RAWS] != access(curMatrix)->dim[COLUMNS])
                            printErr(1, "You cannot use Current Matrix because\nit isn't a Quad one");
                        else
                        {
                            #if WINOS
                                SetExitButtonState(ENABLED);
                            #endif
                            if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                            {
                                matrixFree(matrix, (*righe));
                                return false;
                            }
                            sprint("\nYou're correctly using Current Matrix.\n\n");
                            (*righe) = access(curMatrix)->dim[RAWS];
                            (*colonne) = access(curMatrix)->dim[COLUMNS];
                            return true;
                        }
                    else
                        printErr(1, "You cannot enter an alphanumeric value\nif you haven't entered all elements of last raw");

                if(tmp == -2)
                {
                    if(i > 0)
                        i -= 2;
                    else if((*righe) > 1)
                    {

                        (*righe) -= 2;
                        start_col_index = (*colonne)-1;
                        break;
                    }
                    else
                    {
                        matrixFree(matrix, (*righe));
                        #if WINOS
                            SetExitButtonState(ENABLED);
                        #endif // WINOS
                        return false;
                    }
                }
            }

        }

        if(!(square) && __pmode__ != ALGEBRA_SCALARPRODUCT)
        {
           (*righe) --;
            free((*matrix)[*righe]);
        }

        if(!(lazy_exec))
        {
            (*matrix) = realloc((*matrix), sizeof(ityp*)*((*righe)+1));
            errMem((*matrix), false);
        }

    }
    else
    {


        if(square)
        {
            if(!insertDim(righe, MAX_DIMENSIONS))
            {
                matrixFree(matrix, (*righe));
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif
                return false;
            }
            (*colonne) = (*righe);
        }

        else if(!insertDims(righe, colonne)) return false;

        CLEARBUFFER();

        if(!matrixAlloc(matrix, (dim_typ2){*righe, *colonne}))
        {
            matrixFree(matrix, (*righe));
            #if WINOS
                SetExitButtonState(ENABLED);
            #endif
            return false;
        }

        dim_typ i;

        for(i=0; i<(*colonne); ++i)
        {
            while((tmp = insertElement((*matrix), (dim_typ2){0, i}, square)) != 1 && tmp != -2 && !(char_insert))
                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                    if(access(curMatrix)->dim[RAWS] != *righe || access(curMatrix)->dim[COLUMNS] != *colonne)
                        printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Raws and %hu Columns", righe, colonne);
                    else
                    {
                        #if WINOS
                            SetExitButtonState(ENABLED);
                        #endif
                        if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                        {
                            matrixFree(matrix, (*righe));
                            return false;
                        }
                        sprint("\nYou're correctly using Current Matrix.\n\n");
                        return true;
                    }
                else
                    printErr(1, "You cannot insert characters");

            if(!checkBackTracking(tmp, &i))
            {
                matrixFree(matrix, (*righe));
                #if WINOS
                    SetExitButtonState(ENABLED);
                #endif
                return false;
            }
        }

        dim_typ j;
        volatile sel_typ back_tracking;

        for(i = 1; i<(*righe) && __pmode__ != ALGEBRA_SCALARPRODUCT; ++i)
            for(j=start_col_index; j<(*colonne); ++j)
            {
                while((tmp = insertElement((*matrix), (dim_typ2){i, j}, square)) != 1 && tmp != -2 && !(char_insert))
                    if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                        if(access(curMatrix)->dim[RAWS] != *righe || access(curMatrix)->dim[COLUMNS] != *colonne)
                            printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Raws and %hu Columns", righe, colonne);
                        else
                        {
                            #if WINOS
                                SetExitButtonState(ENABLED);
                            #endif
                            if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                            {
                                matrixFree(matrix, (*righe));
                                return false;
                            }
                            sprint("\nYou're correctly using Current Matrix.\n\n");
                            return true;
                        }
                    else
                        printErr(1, "You cannot insert characters");

                if((back_tracking = checkBackTracking2(tmp, &i, &j, &start_col_index, (*colonne))) == 1)
                    break;

                else if(back_tracking == 2)
                {
                    matrixFree(matrix, (*righe));
                    #if WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    return false;
                }



                if(tmp == -2)
                {
                    if(j > 0)
                        j -= 2;
                    else if(i > 1)
                    {

                        i -= 2;
                        start_col_index = (*colonne)-1;
                        break;
                    }
                    else
                    {
                        matrixFree(matrix, (*righe));
                        #if WINOS
                            SetExitButtonState(ENABLED);
                        #endif
                        return false;
                    }
                }

            }

    }


    if(view)
    {
        printf2(COLOR_SYSTEM, "\nInserted [%hu X %hu] Matrix is:\n\n", *righe, *colonne);
        printMatrix(stdout, (*matrix), (dim_typ2){*righe, *colonne});
    }

    return true;
}
