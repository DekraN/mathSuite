// geometry.c 04/10/2014 Marco Chiarelli aka DekraN
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

    #ifdef WINOS
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

    (void) system(str);

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
__MSUTIL_ bool matrixUTConv(ityp *restrict mat, dim_typ dimq)
{
    dim_typ i, j, k;

    ityp ratio;

    for(i = 0; i < dimq; ++i)
    {
        if( ISZERO(*(mat + dimq*i + i)) ) return false;
        for(j = 0; j < dimq; ++j)
            if(j>i)
            {
                ratio = *(mat + dimq*j + i)/ *(mat + dimq*i + i);
                for(k = 0; k < dimq; ++k)
                    *(mat + dimq*j + k) -= ratio * *(mat + dimq*i + k);
            }
    }

    return true;
}

enum
{
	SARRUS_DIM = 3,
	CARLUCCI_DIM
};

__MSNATIVE_ inline const ityp sarrus(ityp *restrict mat)
{
	return (((*(mat) * *(mat + 4) * *(mat + 8))+(*(mat + 1) * *(mat + 5) * *(mat + 6))+(*(mat + 2) * *(mat + 3) * *(mat + 7))) -
	((*(mat + 2) * *(mat + 4) * *(mat + 6))+(*(mat) * *(mat + 5) * *(mat + 7))+(*(mat + 1) * *(mat + 3) * *(mat + 8))));
}

// It calculates the determinant of a nXn matrix
// by up-triangularizing the square matrix passed
__MSNATIVE_ __MSUTIL_ ityp det(ityp *restrict mat, dim_typ dimq, bool *flag)
{

    if(dimq > 0 && dimq < CARLUCCI_DIM) return checkStdMat(mat, dimq);

    // Conversion of matrix to upper triangular
    ityp D = 1; // storage for determinant

    if(matrixUTConv(mat, dimq))
    {
        dim_typ i;

        for(i = 0; i < dimq; ++i)
            D *= *(mat + dimq*i + i);
    }
    else
    {
    	if(flag)
        	(*flag) = true;
        ityp * S = NULL;
        ityp * V = NULL; // ityp ** V = NULL;

        S = malloc(sizeof(ityp)*dimq);
        errMem(S, MAX_VAL);

        if(!matrixAlloc(&V, (dim_typ2){dimq, dimq}))
        {
            free(S);
        	matrixFree(&V); // matrixFree(&V, dimq);
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

__MSNATIVE_ ityp _matrixTrace(ityp *restrict mat, dim_typ dimq)
{
    dim_typ i;
    ityp res = 0.00;

    for(i=0; i<dimq; ++i)
        res += *(mat + dimq*i + i);

    return res;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ ityp carlucci(ityp *restrict mat)
{
 ityp det = 0.00;
 for(int i = 0; i < CARLUCCI_DIM; ++i)
   det += i < CARLUCCI_DIM ? (pow(-1,i) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + i) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) + (i+1)%CARLUCCI_DIM) + pow(-1,i+1) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + (i+1)%CARLUCCI_DIM) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) +i))*(*(mat + (CARLUCCI_DIM-2+i)%CARLUCCI_DIM) * *(mat + CARLUCCI_DIM + (CARLUCCI_DIM-1+i)%CARLUCCI_DIM) - *(mat + (CARLUCCI_DIM-1+i)%CARLUCCI_DIM) * *(mat + CARLUCCI_DIM + (CARLUCCI_DIM-2+i)%CARLUCCI_DIM)) :
       (*(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + i-CARLUCCI_DIM+2) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) + i-CARLUCCI_DIM) - *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-2) + i-CARLUCCI_DIM) * *(mat + CARLUCCI_DIM*(CARLUCCI_DIM-1) + i-CARLUCCI_DIM+2))*(*(mat + 1-(i-CARLUCCI_DIM)) * *(mat + CARLUCCI_DIM + (CARLUCCI_DIM-1)-(i-CARLUCCI_DIM)) - *(mat + (CARLUCCI_DIM-1)-(i-CARLUCCI_DIM)) * *(mat + CARLUCCI_DIM + 1-(i-CARLUCCI_DIM)));

 return det;

}

__MSSHELL_WRAPPER_ __MSNATIVE_ ityp checkStdMat(ityp *restrict a, dim_typ n)
{
    ityp det = 0.00;
    switch(n)
    {

        case 1:
            det = *a;
            break;
        case 2:
            det = ((*a * *(a + n + 1)) - (*(a + 1) * *(a + n)));
            break;
        case 3:
            det = sarrus(a); // apply SARRUS algorithm
            break;
        case CARLUCCI_DIM:
            det = carlucci(a);
            break;
    }

    return det;
}

__MSSHELL_WRAPPER_ __MSNATIVE_ bool randomMatrix(ityp *restrict matrix, const register dim_typ dim[static MAX_DIMENSIONS])
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
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
            *(matrix + dim[COLUMNS]*i + j) = random(range);

    printf2(COLOR_SYSTEM, "\n\n[%hu X %hu] Randomized Matrix with Range: %hu is:\n", dim[ROWS], dim[COLUMNS], range);
    PRINTL();

    printMatrix(stdout, matrix, dim);

    return true;
}

__MSNATIVE_ void transpose(ityp *restrict matrix, ityp *restrict matrix2, const register dim_typ dim[static MAX_DIMENSIONS])
{

    dim_typ i, j;
    #pragma omp parallel for
    for(i=0; i<dim[COLUMNS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[ROWS]; ++j)
        	*(matrix2 + dim[ROWS]*i + j) = *(matrix + dim[COLUMNS]*j + i);

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

__MSUTIL_ bool __export FattLU(dim_typ n, ityp *restrict c, ityp *restrict l, ityp * a)
{
    ityp m, pivot;
    dim_typ i, j, k;

    // Copia di C su A.
    if(!equalMatrix(&a, c, (dim_typ2){n, n}))
        return false;

    for (k=0; k<n; ++k)
    {
        pivot = *(a + n*k + k);
        if ( ISZERO(pivot) ) return false;
        for (i=k+1; i<n; ++i)
        {
            *(l + n*i + k) = m = *(a + n*i + k) / pivot ;
            for (j=k; j<n; ++j)
                *(a + n*i + j) -= m * *(a + n*k + j);
        }
    }

  /* Adesso "a" contiene U, e "l" la parte inferiore di L */

	#pragma omp parallel for
    for (i=0; i<n; ++i) /* Completa L con zeri ed uni */
    {
        *(l + n*i + i)=1.00;
        #pragma omp parallel for
        for (j=i+1; j<n; ++j)
            *(l + n*i + j)=0.00;
    }

  return true;
}


/*
thanks to: Bibek Subedi:
http://programming-technique.blogspot.it/2011/09/numerical-methods-inverse-of-nxn-matrix.html
for this part of code, which I renamed, modified and adapted to this program
*/

__MSUTIL_ bool __export invertMatrix(ityp *restrict matrix, dim_typ n)
{

    dim_typ i, j;
    ityp a;

    const register dim_typ n2 = n<<1;

	#pragma omp parallel for
    for(i = 0; i < n; ++i)
   		#pragma omp parallel for
        for(j = n; j < n2; ++j)
            *(matrix + n*i + j) = i == (j-n);

    if(!matrixUTConv(matrix, n))
        return false;

	#pragma omp parallel for
    for(i = 0; i < n; ++i)
    {
        const ityp a = *(matrix + n*i + i);
        for(j = 0; j < n2; ++j)
            *(matrix + n*i + j) /= a;
    }

    return true;
}


/*
   Find the cofactor matrix of a square matrix
*/
__MSUTIL_ bool __export CoFactor(ityp *restrict a, ityp *restrict b, dim_typ n)
{
	dim_typ i,j,ii,jj,i1,j1;
	const register dim_typ nminus1 = n-1; 
	ityp _det = 0.00;
	ityp *c = NULL;
	
	if(!matrixAlloc(&c, (dim_typ2){nminus1,nminus1}))
		return false;
		
	for (j=0; j<n; j++) 
  		for (i=0; i<n; i++)
  		{

     		/* Form the adjoint a_ij */
     		i1 = 0;
   			for (ii=0; ii<n; ++ii)
	 		{
    			if(ii == i)
           			continue;
        		j1 = 0;
        		for(jj=0; jj<n; ++jj)
				{
	           		if(jj == j)
	              		continue;
	           		*(c + nminus1*i1 + j1) = *(a + n*ii + jj);
	           		++ j1;
	        	}
        		++ i1;
     		}
		     /* Calculate the determinant */
		     _det = det(c, n-1, NULL);
	     	/* Fill in the elements of the cofactor */
	     	*(b + n*i + j) = pow(-1.0,i+j+2.0) * _det;
  		}

	matrixFree(&c);
	return true;
}

__MSUTIL_ inline bool _MS__private __export adjoint(ityp *restrict a, ityp *restrict b, dim_typ n)
{
	if(CoFactor(a,b,n))
	{
		transpose(b,a,(dim_typ2){n,n});
		return true;
	}
	return false;
}

__MSUTIL_ static inline _MS__private __export ityp PYTHAG(ityp a, ityp b)
{
    const ityp at = fabs(a), bt = fabs(b);
    ityp ct;
    return (at > bt ? at * sqrt(1.00 + ((ct=bt/at) * ct)) : (bt > 0.00 ? bt * sqrt(1.00 + ((ct=at/bt) * ct)) : 0.00));
}

// Jacobi Singular Value Decomposition (SVD)
__MSUTIL_ bool __export dsvd(ityp *restrict a, const register dim_typ dim[static MAX_DIMENSIONS], ityp *w, ityp *v)
{
    // checking whether m < n and correct it
    // by transposing the matrix into n,m matrix
    // causes its correctly dsvd decomposition...

    const register dim_typ m = dim[ROWS], n = dim[COLUMNS];

    if(m < n)
    {
        ityp * tmp = NULL;
        ityp res;
        if(!matrixAlloc(&tmp, (dim_typ2){n, m}))
            return false;
        transpose(a, tmp, dim); // m, n);
        res = dsvd(tmp, (dim_typ2){n, m}, w, v);
        matrixFree(&tmp);
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
                scale += fabs(*(a + n*k + i));
            if (scale)
            {
                for (k = i; k < m; ++k)
                {
                    *(a + n*k + i) /= scale;
                    s += (*(a + n*k + i) * *(a + n*k + i));
                }
                f = *(a + n*i + i);
				g = -SIGN(sqrt(s), f);
                h = f * g - s;
                *(a + n*i + i) = f - g;
                if (i != n - 1)
                {
                    for (j = l; j < n; ++j)
                    {
                        for (s = 0.00, k = i; k < m; ++k)
                            s += (*(a + n*k + i) * *(a + n*k + j));
                        f = s / h;
                        for (k = i; k < m; ++k)
                            *(a + n*k + j) += (f * *(a + n*k + i));
                    }
                }
                for (k = i; k < m; ++k)
                    *(a + n*k + i) *= scale;
            }
        }
        w[i] = scale * g;

        /* right-hand reduction */
        g = s = scale = 0.00;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; ++k)
                scale += fabs(*(a + n*i + k));
            if (scale)
            {
                for (k = l; k < n; ++k)
                {
                    *(a + n*i + k) /= scale;
                    s += *(a + n*i + k) * *(a + n*i + k);
                }
                f = *(a + n*i + l);
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                *(a + n*i + l) = f - g;
                for (k = l; k < n; ++k)
                    rv1[k] = *(a + n*i + k) / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; ++j)
                    {
                        for (s = 0.00, k = l; k < n; ++k)
                            s += (*(a + n*j +k) * *(a + n*i + k));
                        for (k = l; k < n; ++k)
                            *(a + n*j + k) += s * rv1[k];
                    }
                }
                for (k = l; k < n; ++k)
                    *(a + n*i + k) *= scale;
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
                    *(v + n*j + i) = (*(a + n*i + j) / *(a + n*i + l) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.00, k = l; k < n; ++k)
                        s += *(a + n*i + k) * *(v + n*k + j);
                    for (k = l; k < n; ++k)
                        *(v + n*k + j) += s * *(v + n*k + i);
                }
            }
            for (j = l; j < n; ++j)
                *(v + n*i + j) = *(v + n*j + i) = 0.00;
        }
        *(v + n*i + i) = 1.00;
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
                *(a + n*i + j) = 0.00;
        if (g)
        {
            g = 1.00 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; ++j)
                {
                    for (s = 0.00, k = l; k < m; ++k)
                        s += *(a + n*k + i) * *(a + n*k + j);
                    f = (s / *(a + n*i + i)) * g;
                    for (k = i; k < m; ++k)
                        *(a + n*k + j) += f * *(a + n*k + i);
                }
            }
            for (j = i; j < m; ++j)
                *(a + n*j + i) *= g;
        }
        else
        {
            for (j = i; j < m; ++j)
                *(a + n*j + i) = 0.00;
        }
        ++ *(a + n*i + i);
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
                            y = *(a + n*j + nm);
                            z = *(a + n*j + i);
                            *(a + n*j + nm) = (y * c + z * s);
                            *(a + n*j + i) = (z * c - y * s);
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
                        *(v + n*j + k) = - *(v + n*j + k);
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
                    x = *(v + n*jj + j);
                    z = *(v + n*jj + i);
                    *(v + n*jj + j) = (x * c + z * s);
                    *(v + n*jj + i) = (z * c - x * s);
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
                    y = *(a + n*jj + j);
                    z = *(a + n*jj + i);
                    *(a + n*jj + j) = (y * c + z * s);
                    *(a + n*jj + i) = (z * c - y * s);
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
__MSNATIVE_ dim_typ __export rank(ityp *restrict matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    ityp * S = NULL;
    ityp * V = NULL;

    const register ityp maxv = dim[dim[ROWS] >= dim[COLUMNS]];

    S = malloc(sizeof(ityp)*maxv);
    errMem(S, USHRT_MAX);

    if(!matrixAlloc(&V, (dim_typ2){maxv, maxv}))
    {
        free(S);
        return USHRT_MAX;
    }

    dsvd(matrix, dim, S, V);
    matrixFree(&V);

    register dim_typ rnk = maxv;

    for(dim_typ i=0; i<maxv; ++i)
        if(ISZERO(S[i])) -- rnk;

    free(S);

    return rnk;
}

__MSNATIVE_ inline void __system _flushLogBuf(logObj * const which_log)
{
    strcpy(which_log->buffer, NULL_CHAR);
    return;
}

__MSNATIVE_ __WINCALL void __system _editLog(const char path[static MAX_PATH_LENGTH])
{
    #ifdef WINOS
    if(isnSett(BOOLS_SYSLOGSECURITYCHECK))
    {
        char str[PATHCONTAINS_STRING] = NULL_CHAR;
        sprintf(str, "notepad.txt %s", path);
        (void) system(str);
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
        	#ifdef WINOS
            	signal(SIGINT, (__p_sig_fn_t) sigexit);
            #else
            	signal(SIGINT, (__sighandler_t) sigexit);
            #endif
            (void) gets(str);
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
    
    SetColor(COLOR_ERROR);
    
    errno = err;
    perror("\nERRORE");

    _sprint(format, COLOR_ERROR);

    SetDefaultColor();

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

    const bool cond = fp == stdout;
    if(cond)
        SetColor(COLOR_USER);

    fprintf(fp, str);

    if(cond)
        SetDefaultColor();


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

    SetColor(col);
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
    #ifdef WINOS
    	signal(SIGINT, (__p_sig_fn_t) sigproc);
    #else
    	signal(SIGINT, (__sighandler_t) sigproc);
    #endif
    scanner = vscanf(format, ap);

    va_end(ap);

    if(access(sigresult))
        access(exitHandle) = EXITHANDLE_GETCMD;

    if(!scanner)
        access(exitHandle) = EXITHANDLE_EXIT;

    return scanner == count;
}

__MSNATIVE_ _MS__private void __system printMatrix(FILE *fp, ityp *matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{

    const bool assert = fp == stdout;

    if(assert)
        PRINTN();

    dim_typ  i, j;

    for(i=0; i<dim[ROWS]; ++i)
	{
    	if(assert && isSett(BOOLS_PRINTROWSLABELS))
        	fprintf2(fp, "R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
            fprintf2(fp, OUTPUT_CONVERSION_FORMAT, *(matrix + dim[COLUMNS]*i + j));
            fprintf2(fp, "; ");
            if(j >= dim[COLUMNS]-1)
            	fputc('\n', fp);

        }
	}

    if(assert)
    {
        PRINT2N();
        PRINTN();

        if(access(lmpMatrix) && access(lmpMatrix)->matrix)
        	matrixFree(&(access(lmpMatrix)->matrix));

        if(access(lmpMatrix) && !matrixAlloc(&(access(lmpMatrix)->matrix), dim))
        	resetLmpMatrix();

        if(access(lmpMatrix)->matrix && equalMatrix(&(access(lmpMatrix)->matrix), matrix, dim))
        {
            access(lmpMatrix)->dim[ROWS] = dim[ROWS];
            access(lmpMatrix)->dim[COLUMNS] = dim[COLUMNS];
        }
    }

    return;
}

__MSNATIVE_ bool __system __export parse(char expr[], ityp *res)
{

    int err;
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

#define MINMAX_BUFFER_LEN MAX_STRING

__MSNATIVE_ bool __system __export extractMat(dim_typ which_mat)
{

    FILE *fp;
    struct timeval tvBegin;
    char str[MINMAX_BUFFER_LEN] = NULL_CHAR;

    // char *ij_element = NULL;

    dim_typ i;

    matrixObj * const tmp = malloc(sizeof(matrixObj));
    errMem(tmp, false);

    tmp->matrix = NULL;

    fp = NULL;

    gettimeofday(&tvBegin, NULL);
    if((fp = checkForFHErrors(listNo(which_mat, MATRICES)->path, "r")) == NULL)
        return false;

    fflush(fp);


	dim_typ analog_rows, analog_columns = 1;
    char *ij_element = NULL;


    for(tmp->dim[ROWS]=tmp->dim[COLUMNS]=INIT_DIM; fgets(str, sizeof(str), fp) != NULL; ++ tmp->dim[ROWS])  // fscanf(fp, "%[^\n]s", str)) // fgets(str, sizeof(str), fp) != NULL)
    {
        if(!(tmp->dim[ROWS]) && !matrixAlloc(&(tmp->matrix), (dim_typ2){1, 1}))
            return false;
        else
        {
            analog_rows = (tmp->dim[ROWS])+1;

            tmp->matrix = realloc(tmp->matrix, sizeof(ityp)*analog_rows*analog_columns);
            errMem(tmp->matrix, false);
        }


        for(i = 0, ij_element = strtok(str, MATRIXES_SEPERATOR_STRING); ij_element != NULL; ++i, ij_element = strtok(NULL, MATRIXES_SEPERATOR_STRING))
        {
            if(strrchr(ij_element, '\n') != NULL)
            {
                ij_element = NULL;
                break;
            }

            if(!(tmp->dim[ROWS]))
        	{
                tmp->matrix = realloc(tmp->matrix, sizeof(ityp)*analog_rows*analog_columns);
                errMem(tmp->matrix, false);
        	}

            if(isSett(BOOLS_MATRIXPARSING) && !parse(ij_element, (tmp->matrix) + tmp->dim[COLUMNS]*tmp->dim[ROWS] + i))
                continue;
            else
                *((tmp->matrix) + tmp->dim[COLUMNS]*tmp->dim[ROWS] + i) = strtod(ij_element, NULL);

            if(!(tmp->dim[ROWS]))
               analog_columns = ++(tmp->dim[COLUMNS]) +1;
        }
    }

    fclose(fp);

    if(isSett(BOOLS_SHOWDIFFTIME))
    {
        PRINTL();
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
    }

    // REDIRECTING WHICH_MAT LISTBOX ITEM
    listNo(which_mat, MATRICES)->data = tmp;

    return true;
}

__MSNATIVE_ bool __system __export matrixToken(const char string[], ityp **matrix, dim_typ *righe, dim_typ *colonne)
{

    char target[MAX_BUFSIZ];

    strcpy(target, string);

    char *line;
    char *token;
    char buf[MAX_STRING];

    /// char str2[MAX_BUFSIZ];
    dim_typ i;
	dim_typ analog_rows, analog_columns = 1;

    line = token = NULL;

    for((*righe) = (*colonne) = INIT_DIM, line = strtok(target, TERMINATING_STRING); line != NULL; ++ (*righe), line = strtok (line + strlen (line) + 1, TERMINATING_STRING))
    {
        /* String to scan is in buf string..token */
        strncpy(buf, line, sizeof(buf));

        if((!(*righe)) && !matrixAlloc(matrix, (dim_typ2){1, 1}))
            return false;
        else
        {

            analog_rows = ((*righe))+1;

            (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*analog_columns);
            errMem((*matrix), false);
        }

        for(i=0, token=strtok(buf,","); token != NULL; ++ i, token = strtok (token + strlen (token) + 1, ","))
        {

            if(!((*righe)))
            {
                (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*analog_columns);
                errMem((*matrix), false);
            }

            char token2[strlen(token)+1];
            sprintf(token2, "%s;", token);

            if(isSett(BOOLS_MATRIXPARSING) && !parse(token2, (*matrix) + (*colonne)*(*righe) + i))
                continue;
            else
            	*((*matrix) + (*colonne)*(*righe) + i) = strtod(token2, NULL);

            if(!((*righe)))
                analog_columns = ++ (*colonne) +1;

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
            
    for(i=0; i<MAX_MSSMANAGER_PROGS; ++i)
    	if(!strcmp(cmdname, mss_manager[i].cmdname))
    		return &mss_manager[i];

    return NULL;

}

__MSUTIL_ inline int __export cmpfunc(const void * a, const void * b)
{
   return ( *(ityp*)a - *(ityp*)b );
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
    _backupColFile();

    if(isSett(BOOLS_ITEMSAUTOSAVING))
        for(dim_typ i=0; i<MAX_LISTS; ++i)
            updAll(i);

    printf2(COLOR_CREDITS, EXIT_MESSAGE);
    PRINTL();
    return;
 }

__MSNATIVE_ inline void __system safeExit(const int exval)
{
    prepareToExit();
    exit(exval);
    return;
}

__MSNATIVE_ void __system _handleCmdLine(const sel_typ argc, char ** argv)
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
                sprintf(mss_apnt, CMD_BCALC" %s", str);
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

#ifdef XMLCALL

	#ifndef WINOS
		__MSNATIVE_ __MSUTIL_ static inline void itoa(const int val, char string[])
		{
			sprintf(string, "%d", val);
			return;
		}
	#endif
	
	__MSUTIL_ static inline void ftoa(char *string, const float value, const fsel_typ prec)
	{
		sprintf(string, "%.*f", prec, value);
		return;
	}

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

	__MSUTIL_ XMLCALL inline bool __system __export xmlWriteInt(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const int value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			char tmp[SIGN_STRING] = NULL_CHAR;
			#ifdef WINOS
				itoa(value, tmp, sizeof(tmp));
			#else
				itoa(value, tmp);
			#endif
			xmlNodeSetContent(node, BAD_CAST tmp);
			return true;
		}	
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlWriteBool(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const bool value)
	{
		static const char bools_identifiers[MAX_DIMENSIONS][SIGN_STRING] =
		{
			"false",
			"true"
		};
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			xmlNodeSetContent(node, BAD_CAST bools_identifiers[value]);
			return true;
		}
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlWriteFloat(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const float value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			char tmp[SIGN_STRING] = NULL_CHAR;
			ftoa(tmp, value, DEFAULT_PRECISION);
			xmlNodeSetContent(node, BAD_CAST tmp);
			return true;
		}
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlWriteString(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const char string[static MAX_XML_FIELDSTRINGS])
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			xmlNodeSetContent(node, BAD_CAST string);
			return true;
		}
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlGetInt(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, int * value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			(*value) = atoi(xmlNodeGetContent(node));
			return true;
		}
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlGetBool(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, bool * value)
	{
		static const char bools_identifiers[MAX_DIMENSIONS][SIGN_STRING] =
		{
			"false",
			"true"
		};
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			(*value) = strcmp(xmlNodeGetContent(node), bools_identifiers[false]);
			return true;
		}
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlGetFloat(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, float * value)
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			(*value) = atof(xmlNodeGetContent(node));
			return true;
		}
		return false;
	}

	__MSUTIL_ XMLCALL inline bool __system __export xmlGetString(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, char string[static MAX_XML_FIELDSTRINGS])
	{
		(*xpathObj) = xmlXPathEvalExpression( BAD_CAST nodeAddress, xpathCtx );
		xmlNode * node = (*xpathObj)->nodesetval->nodeTab[0];
		xmlXPathFreeObject((*xpathObj));
		if(node)
		{
			strcpy(string, xmlNodeGetContent(node));
			return true;
		}
		return false;
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
		if(!xmlWriteString(&xpathObj, xpathCtx, "/settings/generalSettings/exitChar", ex_char))
			printErr(5, "In Writing /settings/generalSettings/exitChar field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/programPrecision", cur_layout->precision))
			printErr(5, "In Writing /settings/generalSettings/programPrecision field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/stabilizerFactor", cur_layout->stabilizer_factor))
			printErr(5, "In Writing /settings/generalSettings/stabilizerFactor field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/minStirlingNumber", cur_layout->min_stirling_number))
			printErr(5, "In Writing /settings/generalSettings/minStirlingNumber field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/generalSettings/algebra", cur_layout->algebra))
			printErr(5, "In Writing /settings/generalSettings/algebra field");
		if(!xmlWriteFloat(&xpathObj, xpathCtx, "/settings/generalSettings/outlierConstant", cur_layout->outlier_constant))
			printErr(5, "In Writing /settings/generalSettings/outlierConstant field");

		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxRows", cur_layout->matrix_max_rows))
			printErr(5, "In Writing /settings/matricesOptions/maxRows field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxColumns", cur_layout->matrix_max_columns))
			printErr(5, "In Writing /settings/matricesOptions/maxColumns field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/blockSize", cur_layout->block_size))
			printErr(5, "In Writing /settings/matricesOptions/blockSize field");

		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/minOSMMDim", cur_layout->min_osmm_dim))
			printErr(5, "In Writing /settings/matricesOptions/minOSMMDim field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/minStrassenDim", cur_layout->min_strassen_dim))
			printErr(5, "In Writing /settings/matricesOptions/minStrassenDim field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxEigenValuesIterations", cur_layout->max_eigvalues_iterations))
			printErr(5, "In Writing /settings/matricesOptions/maxEigenValuesIterations field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxDSVDIterations", cur_layout->max_dsvd_iterations))
			printErr(5, "In Writing /settings/matricesOptions/maxDSVDIterations field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxSimplexIterations", cur_layout->max_simplex_iterations))
			printErr(5, "In Writing /settings/matricesOptions/maxSimplexIterations field");


		#ifdef WINOS
			#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
		#endif
		for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		{
			char str[DINFO_STRING] = NULL_CHAR;
			char strboolized[INFO_STRING] = NULL_CHAR;
			strboolize(suite_c.memoizers_names[i], strboolized);
			strboolized[0] = toupper(strboolized[0]);
			sprintf(str, "/settings/memoizerOptions/max%sMemoizableIndex", strboolized);
			if(!xmlWriteInt(&xpathObj, xpathCtx, str, cur_layout->max_memoizable_indices[i]))
				printErr(5, "In Writing %s field", str);
		}

		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/baseConversions/minBase", cur_layout->basecalc_minbase))
			printErr(5, "In Writing /settings/baseConversions/minBase field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxBase", cur_layout->basecalc_maxbase))
			printErr(5, "In Writing /settings/baseConversions/maxBase field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxChangebaseBinaryConvnum", cur_layout->max_changebase_binary_convnum))
			printErr(5, "In Writing /settings/baseConversions/maxChangebaseBinaryConvnum field");

		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/minDimension", cur_layout->min_newton_difftables_dim))
			printErr(5, "In Writing /settings/newtonDifferenceTables/minDimension field");
		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/maxDimension", cur_layout->max_newton_difftables_dim))
			printErr(5, "In Writing /settings/newtonDifferenceTables/maxDimension field");

	    if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/romanNumbers/minProcessableNumber", cur_layout->min_roman_number))
	    	printErr(5, "In Writing /settings/romanNumbers/minProcessableNumber field");
	    if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/romanNumbers/maxProcessableNumber", cur_layout->max_roman_number))
	    	printErr(5, "In Writing /settings/romanNumbers/maxProcessableNumber field");

	    if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/minRows", cur_layout->pascal_triangle_min_rows))
	    	printErr(5, "In Writing /settings/pascalsTriangle/minRows field");
	    if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/maxRows", cur_layout->pascal_triangle_max_rows))
	    	printErr(5, "In Writing /settings/pascalsTriangle/maxRows field");

	    /// write some stuffs...

		#ifdef WINOS
			#pragma omp parallel for
		#endif
	    for(i=0; i<MAX_BOOL_SETTINGS; ++i)
	    {
	    	char name[MIN_STRING<<MAX_DIMENSIONS] = NULL_CHAR;
			char strboolized[MIN_STRING<<1] = NULL_CHAR;
			strboolize(suite_c.bools_names[i], strboolized);
	    	sprintf(name, "/settings/booleanKeys/%s", strboolized);
	        if(!xmlWriteBool(&xpathObj, xpathCtx, name, (cur_layout->bools & suite_c.bools[i].bmask) == suite_c.bools[i].bmask))
	        	printErr(5, "In Writing %s field", name);
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
		if(xmlGetString(&xpathObj, xpathCtx, "/settings/generalSettings/exitChar", ex_char))
			tmp->exit_char = ex_char[0];
		else
			printErr(5, "In Parsing /settings/generalSettings/exitChar field");
		int tmp_int;
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/programPrecision", &tmp_int))
			tmp->precision = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/generalSettings/programPrecision field");
			tmp->precision = DEFAULT_PRECISION;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/stabilizerFactor", &tmp_int))
			tmp->stabilizer_factor = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/generalSettings/stabilizerFactor field");
			tmp->stabilizer_factor = DEFAULT_STABILIZER_FACTOR;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/minStirlingNumber", &tmp_int))
			tmp->min_stirling_number = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/generalSettings/minStirlingNumber field");
			tmp->min_stirling_number = DEFAULT_MIN_STIRLING_NUMBER;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/generalSettings/algebra", &tmp_int))
			tmp->algebra = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/generalSettings/algebra field");
			tmp->algebra = DEFAULT_ALGEBRA;
		}
		float tmp_float;
		if(xmlGetFloat(&xpathObj, xpathCtx, "/settings/generalSettings/outlierConstant", &tmp_float))
			tmp->outlier_constant = tmp_float;
		else
		{
			printErr(5, "In Parsing /settings/generalSettings/outlierConstant field");
			tmp->outlier_constant = DEFAULT_OUTLIER_CONSTANT;
		}

	
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxRows", &tmp_int))
			tmp->matrix_max_rows = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxRows field");
			tmp->matrix_max_rows = MAX_RIGHE;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxColumns", &tmp_int))
			tmp->matrix_max_columns = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxColumns field");
			tmp->matrix_max_columns = MAX_COLONNE;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/blockSize", &tmp_int))
			tmp->block_size = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/blockSize field");
			tmp->block_size = DEFAULT_BLOCKSIZE;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/minOSMMDim", &tmp_int))
			tmp->min_osmm_dim = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/minOSMMDim field");
			tmp->min_osmm_dim = DEFAULT_MINOSMMDIM;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/minStrassenDim", &tmp_int))
			tmp->min_strassen_dim = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/minStrassenDim field");
			tmp->min_strassen_dim = DEFAULT_MINSTRASSENDIM;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxEigenValuesIterations", &tmp_int))
			tmp->max_eigvalues_iterations = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxEigenValuesIterations field");
			tmp->max_eigvalues_iterations = DEFAULT_MAX_EIGVALUES_ITERATIONS;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxDSVDIterations", &tmp_int))
			tmp->max_dsvd_iterations = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxDSVDIterations field");
			tmp->max_dsvd_iterations = DEFAULT_MAX_DSVD_ITERATIONS;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxSimplexIterations", &tmp_int))
			tmp->max_simplex_iterations = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxSimplexIterations field");
			tmp->max_simplex_iterations = DEFAULT_MAX_SIMPLEXMETHOD_ITERATIONS;
		}


		#ifdef WINOS
			#pragma omp parallel for num_threads(MAX_MEMOIZABLE_FUNCTIONS)
		#endif
		for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		{
			char str[DINFO_STRING] = NULL_CHAR;
			char strboolized[INFO_STRING] = NULL_CHAR;
			strboolize(suite_c.memoizers_names[i], strboolized);
			strboolized[0] = toupper(strboolized[0]);
			sprintf(str, "/settings/memoizerOptions/max%sMemoizableIndex", strboolized);
			if(xmlGetInt(&xpathObj, xpathCtx, str, &tmp_int))
				tmp->max_memoizable_indices[i] = tmp_int;
			else
			{
				printErr(5, "In Parsing %s field", str);
				tmp->max_memoizable_indices[i] = suite_c.max_memoizable_indices[i];
			}
		}


		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/minBase", &tmp_int))
			tmp->basecalc_minbase = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/baseConversions/minBase field");
			tmp->basecalc_minbase = BCALC_CAMBIAMENTODIBASE_MINBASE;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxBase", &tmp_int))
			tmp->basecalc_maxbase = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/baseConversions/maxBase field");
			tmp->basecalc_maxbase = BCALC_CAMBIAMENTODIBASE_MAXBASE;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxChangebaseBinaryConvnum", &tmp_int))
			tmp->max_changebase_binary_convnum = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/baseConversions/maxChangebaseBinaryConvnum field");
			tmp->max_changebase_binary_convnum = MAX_MINBASE_CONVERTIBLE_NUM;
		}

		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/minDimension", &tmp_int))
			tmp->min_newton_difftables_dim = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/newtonDifferenceTables/minDimension field");
			tmp->min_newton_difftables_dim = MIN_NEWTON_DIFFTABLES_DIM;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/newtonDifferenceTables/maxDimension", &tmp_int))
			tmp->max_newton_difftables_dim = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/newtonDifferenceTables/maxDimension field");
			tmp->max_newton_difftables_dim = MAX_NEWTON_DIFFTABLES_DIM;
		}

		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/romanNumbers/minProcessableNumber", &tmp_int))
			tmp->min_roman_number = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/romanNumbers/minProcessableNumber field");
			tmp->min_roman_number = MIN_PROCESSABLE_ROMAN_NUMBER;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/romanNumbers/maxProcessableNumber", &tmp_int))
			tmp->max_roman_number = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/romanNumbers/maxProcessableNumber field");
			tmp->max_roman_number = MAX_PROCESSABLE_ROMAN_NUMBER;
		}

		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/minRows", &tmp_int))
			tmp->pascal_triangle_min_rows = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/pascalsTriangle/minRows field");
			tmp->pascal_triangle_min_rows = MIN_PASCALTRIANGLE_ROWS;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/pascalsTriangle/maxRows", &tmp_int))
			tmp->pascal_triangle_max_rows = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/pascalsTriangle/maxRows field");
			tmp->pascal_triangle_max_rows = MAX_PASCALTRIANGLE_ROWS;
		}
		
		bool tmp_bool = false;                           

		#ifdef WINOS
			#pragma omp parallel for
		#endif
	    for(i=tmp->bools=0; i<MAX_BOOL_SETTINGS; ++i)
	    {
			char name[MIN_STRING<<MAX_DIMENSIONS] = NULL_CHAR;
			char strboolized[DINFO_STRING] = NULL_CHAR;
	        strboolize(suite_c.bools_names[i], strboolized);
	        sprintf(name, "/settings/booleanKeys/%s", strboolized);
	        if(!xmlGetBool(&xpathObj, xpathCtx, name, &tmp_bool))
	        {
	        	printErr(5, "In Parsing %s field", name);
	        	tmp_bool = suite_c.bools[i].default_val;
	    	}
	        if(tmp_bool)
	            tmp->bools |= suite_c.bools[i].bmask;
	    }

	    xmlExit(path, &doc, &xpathObj, &xpathCtx);

	    return;
	}

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

			int tmp_int;
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/defaultColor", &tmp_int))
	        	COLOR_DEFAULT = tmp_int;
	        else
	        {
	        	printf2(_COLOR_ERROR, "ERROR: In Parsing /colors/defaultColor field");
	        	COLOR_DEFAULT = DEFAULT_COLOR;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/errorsColor", &tmp_int))
	        	COLOR_ERROR = tmp_int;
	        else
	        {
	        	printf2(_COLOR_ERROR, "ERROR: In Parsing /colors/errorsColor field");
	        	COLOR_ERROR = _COLOR_ERROR;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/creditsColor", &tmp_int))
	        	COLOR_CREDITS = tmp_int;
	        else
	        {
	        	printf2(_COLOR_ERROR, "ERROR: In Parsing /colors/creditsColor field");
	        	COLOR_CREDITS = _COLOR_CREDITS;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/userColor", &tmp_int))
	        	COLOR_USER = tmp_int;
	        else
	        {
	        	printf2(_COLOR_ERROR, "ERROR: In Parsing /colors/userColor field");
	        	COLOR_USER = _COLOR_USER;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/systemColor", &tmp_int))
	        	COLOR_SYSTEM = tmp_int;
	        else
	        {
	        	printf2(_COLOR_ERROR, "ERROR: In Parsing /colors/systemColor field");
	        	COLOR_SYSTEM = _COLOR_SYSTEM;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/authorColor", &tmp_int))
	        	COLOR_AUTHOR = tmp_int;
	        else
	        {
	        	printf2(_COLOR_ERROR, "ERROR: In Parsing /colors/authorColor field");
	        	COLOR_AUTHOR = _COLOR_AUTHOR;
	        }

			xmlExit(access(colors_path), &doc, &xpathObj, &xpathCtx);
	        return;
	    }

	    __MSNATIVE_ XMLCALL void __system _backupColFile(void)
		{
			xmlXPathContext * xpathCtx = (xmlXPathContext*) NULL;
			xmlXPathObject * xpathObj = (xmlXPathObject*) NULL;
			xmlDoc * doc = xmlInit(access(colors_path), &xpathCtx);

			if(!xmlWriteInt(&xpathObj, xpathCtx, "/colors/defaultColor", COLOR_DEFAULT))
				printErr(5, "In Writing /colors/defaultColor field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/colors/errorsColor", COLOR_ERROR))
				printErr(5, "In Writing /colors/errorsColor field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/colors/creditsColor", COLOR_CREDITS))
				printErr(5, "In Writing /colors/creditsColor field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/colors/userColor", COLOR_USER))
				printErr(5, "In Writing /colors/userColor field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/colors/systemColor", COLOR_SYSTEM))
				printErr(5, "In Writing /colors/systemColor field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/colors/authorColor", COLOR_AUTHOR))
				printErr(5, "In Writing /colors/authorColor field");

			xmlExit(access(colors_path), &doc, &xpathObj, &xpathCtx);
		    return;
		}

#endif


#ifdef WINOS

    __MSUTIL_ __WINCALL inline BOOL WINAPI __system __export SetExitButtonState(const bool state)
	{
		return DeleteMenu(GetSystemMenu(GetConsoleWindowNT(),state),6,MF_BYPOSITION);
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

#else
	__MSUTIL_ __MSNATIVE_ __system __export int getch( )
    {
        struct termios oldt,
             newt;
        int            ch;
        tcgetattr( STDIN_FILENO, &oldt );
        newt = oldt;
        newt.c_lflag &= ~( ICANON | ECHO );
        tcsetattr( STDIN_FILENO, TCSANOW, &newt );
        ch = getchar();
        tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
        return ch;
    }
#endif

__MSUTIL_ inline const char * const __system __export getFilename(const char path[static MAX_PATH_LENGTH])
{
    const char * const stkr = strrchr(path, '\\');
    return stkr ? stkr+1 : path;
}

__MSUTIL_ inline void __system updInfo(void)
{

    char title[MAXX_STRING];

    sprintf(title, PROG__NAME" - [ %s ] - [ %s ] - [ %s ] - (%s) - [ %s ]", getItemsListNo(ENVS) ? getFilename(listNo(access(lists)[ENVS].cur_item, ENVS)->path) : NULL_CHAR, getItemsListNo(MATRICES) ? getFilename(listNo(access(lists)[MATRICES].cur_item, MATRICES)->path) : NULL_CHAR,
                    getItemsListNo(LOGS) ? getFilename(listNo(access(lists)[LOGS].cur_item, LOGS)->path) : NULL_CHAR, access(sysLogPath), getItemsListNo(LAYOUTS) ? getFilename(listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path) : NULL_CHAR);

	#ifdef WINOS
	    if(!SetConsoleTitle(title))
	        printErr(22, "SetConsoleTitle failed with error: %lu", GetLastError());
	#else
		printf("%c]0;%s%c", '\033', title, '\007');
	#endif

    return;
}

__MSUTIL_ void __system __export SetColor(const sel_typ ForgC)
{
	#ifdef WINOS
	    const HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	    CONSOLE_SCREEN_BUFFER_INFO csbi;

	    if(GetConsoleScreenBufferInfo(hStdOut, &csbi))
	        SetConsoleTextAttribute(hStdOut, (csbi.wAttributes & 0xF0) + (ForgC & 0x0F));
	#else
		static const char xcolors[MAX_COLORS][SIGN_STRING] =
        {
	        "\e[0;30m",
	        "\e[0;34m",
	        "\e[0;32m",
	        "\e[0;36m",
	        "\e[0;31m",
	        "\e[0;35m",
	        "\e[0;33m",
	        "\e[0;37m",
	        "\e[0;90m",
	        "\e[0;94m",
	        "\e[0;92m",
	        "\e[0;96m",
	        "\e[0;91m",
	        "\e[0;95m",
	        "\e[0;93m",
	        "\e[0;97m"
        };
        printf(xcolors[ForgC]);
	#endif

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

__MSNATIVE_ __MSUTIL_ ityp __system __export getDiffTime(struct timeval * t1)
{
	struct timeval result, t2;
	gettimeofday(&t2, NULL);
	long int diff =
	(t2.tv_usec + 1000000 * t2.tv_sec) -
	(t1->tv_usec + 1000000 * t1->tv_sec);
	result.tv_sec = diff / 1000000;
	result.tv_usec = diff % 1000000;
	char str[MINMIN_STRING];
	sprintf(str, "%ld.%06ld\n", (long int) result.tv_sec, (long int) result.tv_usec);
	return atof(str);
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

__MSUTIL_ inline void _MS__private __system __export free_foreach(ityp **matrix, const dim_typ algebra_units, bool mode)
{
    if(mode)
        return;
    #pragma omp parallel for num_threads(algebra_units)
    for(dim_typ i=0; i<algebra_units; ++i)
        matrixFree(&matrix[i]);
    return;
}

__MSUTIL_ inline void _MS__private __system __export free_foreach2(ityp **matrix, const dim_typ algebra_units)
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
        #ifdef WINOS
            SetExitButtonState(ENABLED);
        #endif // WINOS
        printErr(12, "An error occurred during Heap Dynamic Memory Allocation");
        sprint("\n(Sub)Program Terminating...\n\n");
        #ifdef WINOS
            (void) system("PAUSE");
        #endif
        return true;
    }
    return false;
}

__MSNATIVE_ bool __system __export matrixAlloc(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    if(!(*matrix))
        (*matrix) = NULL;

    (*matrix) = calloc(dim[ROWS]*dim[COLUMNS], sizeof(ityp));
    errMem((*matrix), false);

    return true;
}

__MSNATIVE_ void __system __export _matrixFree(ityp **matrix, bool mode)
{
    if(mode || !(*matrix))
        return;

    free((*matrix));
    (*matrix) = NULL; // to avoid dangling references,
    // even if all this pointer passed to this function
    // are generally locally allocated into static void functions.
    return;
}

__MSNATIVE_ bool __system __export equalMatrix(ityp **matrix1, ityp *matrix2, const register dim_typ dim[static MAX_DIMENSIONS])
{
    (*matrix1) = realloc((*matrix1), sizeof(ityp)*dim[ROWS]*dim[COLUMNS]);
    errMem((*matrix1), false);

    dim_typ i, j;

    // Phisically equalling matrix1 values to matrix2 ones
    for(i=0; i<dim[ROWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            *((*matrix1) + dim[COLUMNS]*i + j) = *(matrix2 + dim[COLUMNS]*i + j);

    return true;
}

__MSNATIVE_ _MS__private inline void __system resetLmpMatrix(void)
{
	access(lmpMatrix)->matrix = NULL;
    access(lmpMatrix)->dim[ROWS] = STARTING_LMP_ROWS;
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
		if(access(sysMem)[i].current_max_index)
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
    sprint("- Min Stirling NUMBER: %hu;\n", cur_layout->min_stirling_number);
    sprint("- ALGEBRA: %s;\n", suite_c.algebra_elements_names[cur_layout->algebra]);
    sprint("- Outlier Constant: %.*f;\n", DEFAULT_PRECISION, cur_layout->outlier_constant);

    sprint("- Matrices MAX ROWS: %hu;\n", cur_layout->matrix_max_rows);
    sprint("- Matrices MAX COLUMNS: %hu;\n", cur_layout->matrix_max_columns);
    sprint("- Matrices BLOCK SIZE: %hu;\n", cur_layout->block_size);
    sprint("- Matrices Min OSMM Dimension: %hu;\n", cur_layout->min_osmm_dim);
    sprint("- Matrices Min Strassen Dimension: %hu;\n", cur_layout->min_strassen_dim);
    sprint("- Max EIGENVALUES FINDER Iterations: %hu;\n", cur_layout->max_eigvalues_iterations);
    sprint("- Max DSVD Iterations: %hu;\n", cur_layout->max_dsvd_iterations);
    sprint("- Max SIMPLEX METHOD Iterations: %hu;\n", cur_layout->max_simplex_iterations);


	for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
	{
		char str[SIGN_STRING] = NULL_CHAR;
		strcpy(str, suite_c.memoizers_names[i]);
		toupper_s(str);
		sprint("- Max %s Memoizable Index: %hu;\n", str, cur_layout->max_memoizable_indices[i]);
	}

    sprint("- MIN Processable BASE %hu;\n", cur_layout->basecalc_minbase);
    sprint("- MAX Processable BASE: %hu;\n", cur_layout->basecalc_maxbase);
    sprint("- MAX PROCESSABLE Decimal2Bin NUM: %d;\n", cur_layout->max_changebase_binary_convnum);

    sprint("- NEWTON DIFFTABLES Minimum Dimension: %hu;\n", cur_layout->min_newton_difftables_dim);
    sprint("- NEWTON DIFFTABLES Maximum Dimension: %hu;\n", cur_layout->max_newton_difftables_dim);

    sprint("- MIN Processable ROMAN NUMBER: %hu;\n", cur_layout->min_roman_number);
    sprint("- MAX Processable ROMAN  NUMBER: %hu;\n", cur_layout->max_roman_number);

    sprint("- Pascal's Triangle MIN ROWS: %hu;\n", cur_layout->pascal_triangle_min_rows);
    sprint("- Pascal's Triangle MAX ROWS: %hu.\n", cur_layout->pascal_triangle_max_rows);

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
	#ifdef WINOS
    	signal(SIGINT, (__p_sig_fn_t) sigexit);
    #else
    	signal(SIGINT, (__sighandler_t) sigexit);
    #endif
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
    int err;
    ityp val;
	struct timeval tvBegin;
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
        (void) gets(buf);
    }

    access(exitHandle) = INVALID_EXITHANDLE;
    #ifdef WINOS
    	signal(SIGINT, (__p_sig_fn_t) sigproc);
    #else
    	signal(SIGINT, (__sighandler_t) sigproc);
    #endif

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
    	gettimeofday(&tvBegin, NULL);

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
        printf2(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
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
    return val;
}

__MSNATIVE_ bool __system insertDims(dim_typ *righe, dim_typ *colonne)
{
    const dim_typ old_dims[MAX_DIMENSIONS] =
    {
        (*righe),
        (*colonne)
    };

    printf2(COLOR_CREDITS, "\nEnter ROWS and COLUMNS as expected:\n\
[ROWS%sCOLUMNS].\n", PARSING_SYSTEM_ALLOWED ? "]\n[" : " ");

    if(PARSING_SYSTEM_ALLOWED)
        PRINTHOWTOBACKMESSAGE();

    ityp tmp;
    ityp tmp2;

    tmp = tmp2 = 0.00;

    // uint64_t tmp3[MAX_DIMENSIONS];

    while(PARSING_SYSTEM_ALLOWED ? (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted ROWS", PARSER_SHOWRESULT))) || tmp != ((*righe) = (dim_typ)tmp) || (*righe) < 1 || (*righe) > access(curLayout)->matrix_max_rows) ||
        (isNullVal((tmp = requires(NULL, NULL_CHAR, "Inserted COLUMNS", PARSER_SHOWRESULT))) || tmp != ((*colonne) = (dim_typ)tmp) || (*colonne) < 1 || (*colonne) > access(curLayout)->matrix_max_columns) :
        (!scanf2(2, "%lf %lf", &tmp, &tmp2)) || tmp != ((*righe) = (dim_typ)tmp) || tmp2 != ((*colonne) = (dim_typ)tmp2) || (*righe) < 1 || (*colonne) < 1 || (*righe) > access(curLayout)->matrix_max_rows || (*colonne) > access(curLayout)->matrix_max_columns)
    {
        CLEARBUFFER();
        if(access(exitHandle) == EXITHANDLE_GETCMD) continue;
        if(exitHandleCheck) // if(tmp3[ROWS] == NULL_VAL || tmp3[COLUMNS] == NULL_VAL)
        {
            (*righe) = old_dims[ROWS];
            (*colonne) = old_dims[COLUMNS];
            return false;
        }
        printErr(33, "Invalid [ROWS COLUMNS] format.\nYou have to insert non-negative ROWS and COLUMNS,\n\
and must be respectively less than: %hu and %hu", access(curLayout)->matrix_max_rows, access(curLayout)->matrix_max_columns);
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
        max_dim = mode ? access(curLayout)->matrix_max_columns : access(curLayout)->matrix_max_rows;
        printf2(COLOR_CREDITS, "Enter Matrix %s.", mode ? "COLUMNS":"ROWS");
    }

	PRINTN();
	
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
    if(getItemsListNo(MATRICES) != STARTING_MATNO && access(curMatrix)->matrix && access(lmpMatrix)->matrix)
    {
        if(isSett(BOOLS_ITEMSAUTOSAVING) && !saveItem(access(lists)[MATRICES].cur_item, MATRICES))
            return;
        access(curMatrix)->matrix = access(lmpMatrix)->matrix;
        access(curMatrix)->dim[ROWS] = access(lmpMatrix)->dim[ROWS];
        access(curMatrix)->dim[COLUMNS] = access(lmpMatrix)->dim[COLUMNS];
        // equalMatrix(&suite.current_matrix, suite.last_matrix_printed, suite.lmp_rows, suite.lmp_columns);
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

__MSNATIVE_ bool __system insertNMMatrix(ityp **matrix, const register dim_typ dim[static MAX_DIMENSIONS])
{
    if(!matrixAlloc(matrix, dim))
    {
        #ifdef WINOS
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
    for(i=start_col_index=0; tmp != -1 && i<dim[ROWS]; ++i)
        for(j=start_col_index; tmp != -1 && j<dim[COLUMNS]; ++j)
        {
            while((tmp = insertElement((*matrix), (dim_typ2){i, j}, dim[COLUMNS], false)) != 1 && tmp != -2 && !(char_insert))
                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                    if(access(curMatrix)->dim[ROWS] != dim[ROWS] || access(curMatrix)->dim[COLUMNS] != dim[COLUMNS])
                        printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Rows and %hu Columns", dim[ROWS], dim[COLUMNS]);
                    else
                    {
                        if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                        {
                            #ifdef WINOS
                                SetExitButtonState(ENABLED);
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
                        matrixFree(matrix);
                        #ifdef WINOS
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
                    matrixFree(matrix);
                    #ifdef WINOS
                        SetExitButtonState(ENABLED);
                    #endif // WINOS
                    return false;
                }
            }

        }

    return true;
}

__MSNATIVE_ volatile char __system insertElement(ityp *restrict matrix, const register dim_typ dim[static MAX_DIMENSIONS], const register dim_typ columns, bool square)
{

    if(square && (dim[ROWS] == MAX_RIGHE_PER_COLONNE || dim[COLUMNS] == MAX_RIGHE_PER_COLONNE))
    {
        printErr(33, "MAX ROWS per COLUMNS Reached");
        return 0;
    }

    if(dim[ROWS] == access(curLayout)->matrix_max_rows)
    {
        printErr(33, "MAX ROWS Reached");
        return 0;
    }

    if(dim[COLUMNS] == access(curLayout)->matrix_max_columns)
    {
        printErr(33, "MAX COLUMNS Reached");
        return 0;
    }

    // PRINTN();
    printf2(COLOR_CREDITS, "\nEnter [%hu,%hu] Matrix Element.\n", dim[ROWS]+1, dim[COLUMNS]+1);

    sel_typ tmp;
    tmp = 1;
    // tmp = true;

    do
    {
        CLEARBUFFER();

        char str[MIN_STRING];

        sprintf(str, "[%hu,%hu] Matrix Element correctly inserted", dim[ROWS]+1, dim[COLUMNS]+1);


        // tmp = scanf(INPUT_CONVERSION_FORMAT, &matrix[riga][colonna]);
        // if((matrix[riga][colonna] = requires(NULL, str, restr, true, false, false, true, false)) == NULL_VAL) tmp=false;
        // if(tmp == NULL_VAL) tmp = 0;

        // FARE ATTENZIONE PER IL GETSET SYSTEM CARATTERIZZATO DAI VALORI DI RITORNO NULL_VALN della funzione requires(NULL, ...)


        tmp = PARSING_SYSTEM_ALLOWED ? (!((*(matrix + columns*dim[ROWS] + dim[COLUMNS]) = requires(NULL, NULL_CHAR, str, PARSER_SHOWRESULT)) == NULL_VAL)) :
                scanf2(1, INPUT_CONVERSION_FORMAT, matrix + columns*dim[ROWS] + dim[COLUMNS]);




        CLEARBUFFER();

        if(access(exitHandle) == EXITHANDLE_SETCMD)
            return -1;

        if(access(exitHandle) == EXITHANDLE_BACKCMD)
            return -2;

        if(access(exitHandle) == EXITHANDLE_EXIT && (!dim[ROWS]) && (!dim[COLUMNS]))
            return -3;

        // tmp = PARSING_SYSTEM_ALLOWED ? (!(isNullVal((matrix[riga][colonna] = requires(NULL, NULL_CHAR, str, true, false, false, true, false))))) :
            // scanf(INPUT_CONVERSION_FORMAT, &matrix[riga][colonna]);// , printf2("%s: ", str), printf2(OUTPUT_CONVERSION_FORMAT, matrix[riga][colonna]), printf(".\n\n");

        if(isnSett(BOOLS_MATRIXPARSING))
        {
            printf2(COLOR_USER, "%s: ", str);
            printf2(COLOR_USER, OUTPUT_CONVERSION_FORMAT, *(matrix + columns*dim[ROWS] + dim[COLUMNS]));
            printf2(COLOR_USER, ".\n\n");
        }

        if(!(*(matrix + columns*dim[ROWS] + dim[COLUMNS])) && dcheck && INVERSE_OPS && ((__pmode__ >= ALGOPS_MATRIXMULTIPLICATION && __pmode__ <= ALGOPS_DOTPRODUCT )|| __pmode__ == ALGOPS_SCALARDIVISIONMATRIX))
           printErr(33, "You cannot enter a 0 because program is performing a Division somewhere");

    }
    while((tmp != 1 && !(dim[ROWS]) && !(dim[COLUMNS])) || (!(*(matrix + columns*dim[ROWS] + dim[COLUMNS])) && dcheck && INVERSE_OPS && ((__pmode__ >= ALGOPS_DOTPRODUCT && __pmode__ <= ALGOPS_SCALARDIVISIONMATRIX)
    ||__pmode__ == ALGOPS_SCALARDIVISIONMATRIX))||isDomainForbidden(*(matrix + columns*dim[ROWS] + dim[COLUMNS]), INPUT));

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

__MSNATIVE_ bool __system __export enterMatrix(ityp **matrix, dim_typ *righe, dim_typ *colonne, bool square, bool view)
{

    printf2(COLOR_CREDITS, "\n\nEnter%s Matrix", square ? " Quad" : "");
    printf2(COLOR_CREDITS, " by following procedure you'll see.\n");

    volatile char tmp;
    dim_typ start_col_index;

    tmp = 1;
    start_col_index = 0;

    #ifdef WINOS
        SetExitButtonState(DISABLED);
    #endif // WINOS


    if(isSett(BOOLS_INSERTMODE))
    {
        printf2(COLOR_CREDITS, "And when you reach desired rows %s columns dimensions, %s.\n\n",
               square ? "=":"and", isSett(BOOLS_MATRIXPARSING) ? "press ENTER":"insert an alphanumeric value");

        (*matrix) = malloc(MAX_DIMENSIONS*sizeof(ityp));
        errMem((*matrix), false);

        (*righe) = 1;

        dim_typ analog_rows = MAX_DIMENSIONS;
        dim_typ analog_columns;

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
                    matrixFree(matrix);
                    #ifdef WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    return false;
                }
            }

            analog_columns = (*colonne)+1;

            (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*analog_columns);

            errMem((*matrix), false);

            if((tmp = insertElement((*matrix), (dim_typ2){0, (*colonne)}, *colonne, square)) == -1 && getItemsListNo(MATRICES) != STARTING_MATNO)
            {
                if(square && access(curMatrix)->dim[ROWS] != access(curMatrix)->dim[COLUMNS])
                {
                    printErr(1, "You cannot use Current Matrix because\nit isn't a Quad one");
                    (*colonne) --;
                }
                else
                {
                    #ifdef WINOS
                        SetExitButtonState(ENABLED);
                    #endif
                    if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                    {
                        matrixFree(matrix);
                        return false;
                    }
                    sprint("\nYou are correctly using Current Matrix.\n\n");
                    (*righe) = access(curMatrix)->dim[ROWS];
                    (*colonne) = access(curMatrix)->dim[COLUMNS];
                    return true;
                }
            }
        }

        tmp = 1;
        (*colonne) --;

        (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*(*colonne));
        errMem((*matrix), false);


        for(*righe = 1; (square ? *righe < *colonne : (tmp && __pmode__ != ALGOPS_DOTPRODUCT)); ++(*righe))
        {
            dim_typ i;
            analog_rows = (*righe)+1;

            (*matrix) = realloc((*matrix), sizeof(ityp)*analog_rows*(*colonne));
            errMem((*matrix), false);

            for(i=start_col_index; tmp && i<*colonne; ++i)
            {
                /// if(tmp == -2) return false; // simple backtracking on second or major rows...
                // BACKTRACKING ON SECOND or MAJOR ROWS EXPERIMENTAL FORMULA

                while((tmp = insertElement((*matrix), (dim_typ2){*righe, i}, *colonne, square)) != 1 && tmp != -2 && !(char_insert) && i)
                    if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                        if(square && access(curMatrix)->dim[ROWS] != access(curMatrix)->dim[COLUMNS])
                            printErr(1, "You cannot use Current Matrix because\nit isn't a Quad one");
                        else
                        {
                            #ifdef WINOS
                                SetExitButtonState(ENABLED);
                            #endif
                            if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                            {
                                matrixFree(matrix);
                                return false;
                            }
                            sprint("\nYou're correctly using Current Matrix.\n\n");
                            (*righe) = access(curMatrix)->dim[ROWS];
                            (*colonne) = access(curMatrix)->dim[COLUMNS];
                            return true;
                        }
                    else
                        printErr(1, "You cannot enter an alphanumeric value\nif you haven't entered all elements of last row");

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
                        matrixFree(matrix);
                        #ifdef WINOS
                            SetExitButtonState(ENABLED);
                        #endif // WINOS
                        return false;
                    }
                }
            }

        }

        if(!(square) && __pmode__ != ALGOPS_DOTPRODUCT)
           (*righe) --;

    }
    else
    {


        if(square)
        {
            if(!insertDim(righe, MAX_DIMENSIONS))
            {
                matrixFree(matrix);
                #ifdef WINOS
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
            matrixFree(matrix);
            #ifdef WINOS
                SetExitButtonState(ENABLED);
            #endif
            return false;
        }

        dim_typ i;

        for(i=0; i<(*colonne); ++i)
        {
            while((tmp = insertElement((*matrix), (dim_typ2){0, i}, *colonne, square)) != 1 && tmp != -2 && !(char_insert))
                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                    if(access(curMatrix)->dim[ROWS] != *righe || access(curMatrix)->dim[COLUMNS] != *colonne)
                        printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Rows and %hu Columns", righe, colonne);
                    else
                    {
                        #ifdef WINOS
                            SetExitButtonState(ENABLED);
                        #endif
                        if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                        {
                            matrixFree(matrix);
                            return false;
                        }
                        sprint("\nYou're correctly using Current Matrix.\n\n");
                        return true;
                    }
                else
                    printErr(1, "You cannot insert characters");

            if(!checkBackTracking(tmp, &i))
            {
                matrixFree(matrix);
                #ifdef WINOS
                    SetExitButtonState(ENABLED);
                #endif
                return false;
            }
        }

        dim_typ j;
        volatile sel_typ back_tracking;

        for(i = 1; i<(*righe) && __pmode__ != ALGOPS_DOTPRODUCT; ++i)
            for(j=start_col_index; j<(*colonne); ++j)
            {
                while((tmp = insertElement((*matrix), (dim_typ2){i, j}, *colonne, square)) != 1 && tmp != -2 && !(char_insert))
                    if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                        if(access(curMatrix)->dim[ROWS] != *righe || access(curMatrix)->dim[COLUMNS] != *colonne)
                            printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Rows and %hu Columns", righe, colonne);
                        else
                        {
                            #ifdef WINOS
                                SetExitButtonState(ENABLED);
                            #endif
                            if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim))
                            {
                                matrixFree(matrix);
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
                    matrixFree(matrix);
                    #ifdef WINOS
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
                        matrixFree(matrix);
                        #ifdef WINOS
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

