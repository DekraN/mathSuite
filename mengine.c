// mengine.c 04/10/2014 Marco Chiarelli aka DekraN
/*
WARNING!!! This program is intended to be used, so linked at the compilation,
exclusively with main.c of my suite program! I do not assume any responsibilities
about the use with any other code-scripts.xm
*/

#include "dutils.h" // DA RENDERE VISIBILE SIA AL COMPILATORE CHE AL LINKER

 inline void toupper_s(char *string)
{
    const register size_t len = strlen(string);
    #pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = toupper(string[i]);
    return;
}

 inline void tolower_s(char *string)
{
    const register size_t len = strlen(string);
    #pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = tolower(string[i]);
    return;
}

inline void asteriskize(char *string)
{
	const register size_t len = strlen(string);
	#pragma omp parallel for
    for(dim_typ i=0; i<len; ++i)
        string[i] = ASTERISK_CHAR;
    return; 
}
__MATHSUITE void strundsc(const char *string, char name[])
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

__MATHSUITE void strboolize(const char *string, char name[])
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

__MATHSUITE inline void strfnm(const char *string, char file_name[static MAX_PATH_LENGTH])
{
    strundsc(string, file_name);
    strcat(file_name, "."DEFAULT_HELP_FILE_EXTENSION);
    return;
}

/*  Bit counter by Ratko Tomic */
 inline int  countbits(long i)
{
      i = ((i & 0xAAAAAAAAL) >>  1) + (i & 0x55555555L);
      i = ((i & 0xCCCCCCCCL) >>  2) + (i & 0x33333333L);
      i = ((i & 0xF0F0F0F0L) >>  4) + (i & 0x0F0F0F0FL);
      i = ((i & 0xFF00FF00L) >>  8) + (i & 0x00FF00FFL);
      i = ((i & 0xFFFF0000L) >> 16) + (i & 0x0000FFFFL);
      return (int)i;
}

 inline int  ucountbits(unsigned long num)
{
    int count;
    for(count = 0; num; ++count, num &= num - 1);
    return count;
}

 char  *strrev(char *str)
{
    char *p1, *p2;

    if (! str || ! *str)
        return str;

    for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
        *p1 ^= *p2 ^= *p1 ^= *p2;

    return str;
}

 char  *replace(char const * const original, char const * const pattern, char const * const replacement)
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

 bool   file_exists(const char * filename)
{

    FILE *file = NULL;

    if((file = fopen(filename, "r")))
    {
        fclose(file);
        return true;
    }
    return false;
}

__MATHSUITE bool  readFile(const char path[static MAX_PATH_LENGTH])
{
    char c;
    FILE *fp = NULL;

    if(!(fp = fopen(path, "r")))
        return false;

    SHOWPAUSEMESSAGE();
    msyprintf(COLOR_SYSTEM, "\n%s:\n", path);
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

__MATHSUITE bool  printFile(const char path[static MAX_PATH_LENGTH])
{
    if(!file_exists(path))
    {
        printErr(2, "Non-existing File:\n%s", path);
        return false;
    }

    char str[MAX_PATH_LENGTH+INFO_STRING];

    #ifdef WINOS
        msprintf(COLOR_CREDITS, "Enter Device Name on which you want to print the File:\n%s.\n", path);
        msprintf(COLOR_CREDITS, "Enter %c for Back.\n\n", SCANFEXIT_CHAR);

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

__MATHSUITE inline bool  writeFile(const char path[static MAX_PATH_LENGTH])
{
    FILE *fp;
    if((fp = checkForFHErrors(path, "w")) == NULL)
        return false;

    fclose(fp);
    return true;
}

__MATHSUITE inline FILE *  checkForFHErrors(const char path[static MAX_PATH_LENGTH], char mode[static 1])
{
    FILE *fp;

    if((fp = fopen(path, mode)) == NULL)
        printErr(2, "Not possible to %s File:\n%s", mode[0] == 'r' ? "read":"write/create", path);

    return fp;
}

__MATHSUITE inline bool  frename(const char name[static MAX_PATH_LENGTH], const char newname[static MAX_PATH_LENGTH])
{
    int err;

    if((err = rename(name, newname)))
    {
        printErr(err, "An error occurred during File Renaming Process: \n%s", name);
        return false;
    }

    return true;
}

// Conversion of matrix to uppe triangular
// by Bibek Subedi original, adapted by me
 bool matrixUTConv(mpfr_t *restrict mat, dim_typ dimq, sel_typ *perms)
{
    dim_typ i, j, k, m, n;
    mpfr_t ratio; 

	mpfr_init(ratio); 
    for(i = 0; i < dimq; ++i)
    {
    	if( mpfr_zero_p(*(mat + dimq*i + i)) )
		{
			for(j=i+1,m=i-1; j<dimq && m<dimq; ++j, --m)
			{
				// checking for useful permutations
				// and eventually swapping really the two vectors
				// the ratio variable is used here as a tmpbuf
				if(mpfr_zero_p(*(mat + dimq*j + i)) == 0 && mpfr_zero_p(*(mat + dimq*i + j)) == 0)
				{
					for(k=0; k<dimq; ++k)
					{
						mpfr_set(ratio, *(mat + dimq*i + k), MPFR_RNDN);
						mpfr_set(*(mat + dimq*i + k), *(mat + dimq*j + k), MPFR_RNDN);
						mpfr_set(*(mat + dimq*j + k), ratio, MPFR_RNDN);
					}
					if(perms)
						++ *perms;
					break;
				}
				if(mpfr_zero_p(*(mat + dimq*m + i)) == 0 && mpfr_zero_p(*(mat + dimq*i + m)) == 0)
				{
					for(n=0; n<dimq; ++n)
					{
						mpfr_set(ratio, *(mat + dimq*i + n), MPFR_RNDN);
						mpfr_set(*(mat + dimq*i + n), *(mat + dimq*m + n), MPFR_RNDN);
						mpfr_set(*(mat + dimq*m + n), ratio, MPFR_RNDN);
					}
					if(perms)
						++ *perms;
					break;
				}
			}
		}
		//return false;
        for(j = 0; j < dimq; ++j)
            if(j>i)
            {
                mpfr_div(ratio, *(mat + dimq*j + i), *(mat + dimq*i + i), MPFR_RNDN);
                for(k = 0; k < dimq; ++k)
                {
                	mpfr_mul(ratio, ratio, *(mat + dimq*i + k), MPFR_RNDN);
                	mpfr_sub(*(mat + dimq*j + k), *(mat + dimq*j + k), ratio, MPFR_RNDN);
				}
            }
    }

	mpfr_clear(ratio);
    return true;
}

__MATHSUITE void det_3d(mpfr_t rop, mpfr_t *restrict a)
{
	mpfr_t tmp, tmp2, tmp3; 
	mpfr_set_ui(rop, 0, MPFR_RNDN);
	
	#pragma omp parallel sections num_threads(3) private(tmp, tmp2, tmp3)
	{
		#pragma omp section
		{
			mpfr_init_set(tmp, a[1+1*3], MPFR_RNDN);
			mpfr_mul(tmp, tmp, a[2+2*3], MPFR_RNDN);
			mpfr_init_set(tmp2, a[1+2*3], MPFR_RNDN); 
			mpfr_mul(tmp2, tmp2, a[2+1*3], MPFR_RNDN);
			mpfr_init(tmp3);
			mpfr_sub(tmp3, tmp, tmp2, MPFR_RNDN);
			mpfr_mul(tmp3, tmp3, a[0+0*3], MPFR_RNDN);
			#pragma omp critical
				mpfr_add(rop, rop, tmp3, MPFR_RNDN);	
			mpfr_clears(tmp, tmp2, tmp3, NULL);
		}
		
		#pragma omp section
		{
			mpfr_init_set(tmp, a[1+2*3], MPFR_RNDN);
			mpfr_mul(tmp, tmp, a[2+0*3], MPFR_RNDN);
			mpfr_init_set(tmp2, a[1+0*3], MPFR_RNDN); 
			mpfr_mul(tmp2, tmp2, a[2+2*3], MPFR_RNDN);
			mpfr_init(tmp3);
			mpfr_sub(tmp3, tmp, tmp2, MPFR_RNDN);
			mpfr_mul(tmp3, tmp3, a[0+1*3], MPFR_RNDN);
			#pragma omp critical
				mpfr_add(rop, rop, tmp3, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, NULL);	
		}
		
		#pragma omp section
		{
			mpfr_init_set(tmp, a[1+0*3], MPFR_RNDN);
			mpfr_mul(tmp, tmp, a[2+1*3], MPFR_RNDN);
			mpfr_init_set(tmp2, a[1+1*3], MPFR_RNDN); 
			mpfr_mul(tmp2, tmp2, a[2+0*3], MPFR_RNDN);
			mpfr_init(tmp3);
			mpfr_sub(tmp3, tmp, tmp2, MPFR_RNDN);
			mpfr_mul(tmp3, tmp3, a[0+2*3], MPFR_RNDN);
			#pragma omp critical
				mpfr_add(rop, rop, tmp3, MPFR_RNDN);
			mpfr_clears(tmp, tmp2, tmp3, NULL);	
		}
		
	}
	
	return;
}

__MATHSUITE void det_4d(mpfr_t rop, mpfr_t *restrict a)
{
	mpfr_t tmp11, tmp12, tmp13,
		   tmp21, tmp22, tmp23,
		   tmp31, tmp32, tmp33,
		   tmp41, tmp42, tmp43;
	mpfr_t tmp, rop_tmp; 
		   
	mpfr_set_ui(rop, 0, MPFR_RNDN);
	
	#pragma omp parallel sections num_threads(4) private(rop_tmp, tmp)
	{
		#pragma omp section
		{
			mpfr_init_set(tmp, a[0+0*4], MPFR_RNDN); 
			mpfr_init_set_ui(rop_tmp, 0, MPFR_RNDN);
			#pragma omp parallel sections num_threads(3) private(tmp11, tmp12, tmp13)
			{
				#pragma omp section
				{
					mpfr_init_set(tmp11, a[2+2*4], MPFR_RNDN);
					mpfr_mul(tmp11, tmp11, a[3+3*4], MPFR_RNDN);
					mpfr_init_set(tmp12, a[2+3*4], MPFR_RNDN); 
					mpfr_mul(tmp12, tmp12, a[3+2*4], MPFR_RNDN);
					mpfr_init(tmp13);
					mpfr_sub(tmp13, tmp11, tmp12, MPFR_RNDN);
					mpfr_mul(tmp13, tmp13, a[1+1*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp13, MPFR_RNDN);	
					mpfr_clears(tmp11, tmp12, tmp13, NULL);
				}
				
				#pragma omp section
				{
					mpfr_init_set(tmp11, a[2+1*4], MPFR_RNDN);
					mpfr_mul(tmp11, tmp11, a[3+3*4], MPFR_RNDN);
					mpfr_init_set(tmp12, a[2+3*4], MPFR_RNDN); 
					mpfr_mul(tmp12, tmp12, a[3+1*4], MPFR_RNDN);
					mpfr_init(tmp13);
					mpfr_sub(tmp13, tmp11, tmp12, MPFR_RNDN);
					mpfr_mul(tmp13, tmp13, a[1+2*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_sub(rop_tmp, rop_tmp, tmp13, MPFR_RNDN);
					mpfr_clears(tmp11, tmp12, tmp13, NULL);	
				}
				
				#pragma omp section
				{
					// mpfr_t tmp, tmp2, tmp3; 
					mpfr_init_set(tmp11, a[2+1*4], MPFR_RNDN);
					mpfr_mul(tmp11, tmp11, a[3+2*4], MPFR_RNDN);
					mpfr_init_set(tmp12, a[2+2*4], MPFR_RNDN); 
					mpfr_mul(tmp12, tmp12, a[3+1*4], MPFR_RNDN);
					mpfr_init(tmp13);
					mpfr_sub(tmp13, tmp11, tmp12, MPFR_RNDN);
					mpfr_mul(tmp13, tmp13, a[1+3*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp13, MPFR_RNDN);
					mpfr_clears(tmp11, tmp12, tmp13, NULL);	
				}
				
			}
			mpfr_mul(tmp, tmp, rop_tmp, MPFR_RNDN);
			#pragma omp critical
				mpfr_add(rop, rop, tmp, MPFR_RNDN);	
			mpfr_clears(tmp, rop_tmp, NULL);
		}
		
		#pragma omp section
		{
			mpfr_init_set(tmp, a[0+1*4], MPFR_RNDN); 
			mpfr_init_set_ui(rop_tmp, 0, MPFR_RNDN);
			
			#pragma omp parallel sections num_threads(3) private(tmp21, tmp22, tmp23)
			{
				#pragma omp section
				{
					mpfr_init_set(tmp21, a[2+2*4], MPFR_RNDN);
					mpfr_mul(tmp21, tmp21, a[3+3*4], MPFR_RNDN);
					mpfr_init_set(tmp22, a[2+3*4], MPFR_RNDN); 
					mpfr_mul(tmp22, tmp22, a[3+2*4], MPFR_RNDN);
					mpfr_init(tmp23);
					mpfr_sub(tmp23, tmp21, tmp22, MPFR_RNDN);
					mpfr_mul(tmp23, tmp23, a[1+0*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp23, MPFR_RNDN);	
					mpfr_clears(tmp21, tmp22, tmp23, NULL);
				}
				
				#pragma omp section
				{
					mpfr_init_set(tmp21, a[2+0*4], MPFR_RNDN);
					mpfr_mul(tmp21, tmp21, a[3+3*4], MPFR_RNDN);
					mpfr_init_set(tmp22, a[2+3*4], MPFR_RNDN); 
					mpfr_mul(tmp22, tmp22, a[3+0*4], MPFR_RNDN);
					mpfr_init(tmp23);
					mpfr_sub(tmp23, tmp21, tmp22, MPFR_RNDN);
					mpfr_mul(tmp23, tmp23, a[1+2*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_sub(rop_tmp, rop_tmp, tmp23, MPFR_RNDN);
					mpfr_clears(tmp21, tmp22, tmp23, NULL);	
				}
				
				#pragma omp section
				{
					mpfr_init_set(tmp21, a[2+0*4], MPFR_RNDN);
					mpfr_mul(tmp21, tmp21, a[3+2*4], MPFR_RNDN);
					mpfr_init_set(tmp22, a[2+2*4], MPFR_RNDN); 
					mpfr_mul(tmp22, tmp22, a[3+0*4], MPFR_RNDN);
					mpfr_init(tmp23);
					mpfr_sub(tmp23, tmp21, tmp22, MPFR_RNDN);
					mpfr_mul(tmp23, tmp23, a[1+3*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp23, MPFR_RNDN);
					mpfr_clears(tmp21, tmp22, tmp23, NULL);	
				}
				
			}
			mpfr_mul(tmp, tmp, rop_tmp, MPFR_RNDN);
			#pragma omp critical
				mpfr_sub(rop, rop, tmp, MPFR_RNDN);	
			mpfr_clears(tmp, rop_tmp, NULL);
		}
		
		#pragma omp section
		{
			mpfr_init_set(tmp, a[0+2*4], MPFR_RNDN); 
			mpfr_init_set_ui(rop_tmp, 0, MPFR_RNDN);
			
			#pragma omp parallel sections num_threads(3) private(tmp31, tmp32, tmp33)
			{
				#pragma omp section
				{
					mpfr_init_set(tmp31, a[2+1*4], MPFR_RNDN);
					mpfr_mul(tmp31, tmp31, a[3+3*4], MPFR_RNDN);
					mpfr_init_set(tmp32, a[2+3*4], MPFR_RNDN); 
					mpfr_mul(tmp32, tmp32, a[3+1*4], MPFR_RNDN);
					mpfr_init(tmp33);
					mpfr_sub(tmp33, tmp31, tmp32, MPFR_RNDN);
					mpfr_mul(tmp33, tmp33, a[1+0*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp33, MPFR_RNDN);	
					mpfr_clears(tmp31, tmp32, tmp33, NULL);
				}
				
				#pragma omp section
				{
					mpfr_init_set(tmp31, a[2+0*4], MPFR_RNDN);
					mpfr_mul(tmp31, tmp31, a[3+3*4], MPFR_RNDN);
					mpfr_init_set(tmp32, a[2+3*4], MPFR_RNDN); 
					mpfr_mul(tmp32, tmp32, a[3+0*4], MPFR_RNDN);
					mpfr_init(tmp33);
					mpfr_sub(tmp33, tmp31, tmp32, MPFR_RNDN);
					mpfr_mul(tmp33, tmp33, a[1+1*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_sub(rop_tmp, rop_tmp, tmp33, MPFR_RNDN);
					mpfr_clears(tmp31, tmp32, tmp33, NULL);	
				}
				
				#pragma omp section
				{
					mpfr_init_set(tmp31, a[2+0*4], MPFR_RNDN);
					mpfr_mul(tmp31, tmp31, a[3+1*4], MPFR_RNDN);
					mpfr_init_set(tmp32, a[2+1*4], MPFR_RNDN); 
					mpfr_mul(tmp32, tmp32, a[3+0*4], MPFR_RNDN);
					mpfr_init(tmp33);
					mpfr_sub(tmp33, tmp31, tmp32, MPFR_RNDN);
					mpfr_mul(tmp33, tmp33, a[1+3*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp33, MPFR_RNDN);
					mpfr_clears(tmp31, tmp32, tmp33, NULL);	
				}
				
			}
			mpfr_mul(tmp, tmp, rop_tmp, MPFR_RNDN);
			#pragma omp critical
				mpfr_add(rop, rop, tmp, MPFR_RNDN);	
			mpfr_clears(tmp, rop_tmp, NULL);
		}
		
		#pragma omp section
		{
			mpfr_init_set(tmp, a[0+3*4], MPFR_RNDN); 
			mpfr_init_set_ui(rop_tmp, 0, MPFR_RNDN);
			
			#pragma omp parallel sections num_threads(3) private(tmp41, tmp42, tmp43)
			{
				#pragma omp section
				{
					mpfr_init_set(tmp41, a[2+1*4], MPFR_RNDN);
					mpfr_mul(tmp41, tmp41, a[3+2*4], MPFR_RNDN);
					mpfr_init_set(tmp42, a[2+2*4], MPFR_RNDN); 
					mpfr_mul(tmp42, tmp42, a[3+1*4], MPFR_RNDN);
					mpfr_init(tmp43);
					mpfr_sub(tmp43, tmp41, tmp42, MPFR_RNDN);
					mpfr_mul(tmp43, tmp43, a[1+0*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp43, MPFR_RNDN);	
					mpfr_clears(tmp41, tmp42, tmp43, NULL);
				}
				
				#pragma omp section
				{ 
					mpfr_init_set(tmp41, a[2+0*4], MPFR_RNDN);
					mpfr_mul(tmp41, tmp41, a[3+2*4], MPFR_RNDN);
					mpfr_init_set(tmp42, a[2+2*4], MPFR_RNDN); 
					mpfr_mul(tmp42, tmp42, a[3+0*4], MPFR_RNDN);
					mpfr_init(tmp43);
					mpfr_sub(tmp43, tmp41, tmp42, MPFR_RNDN);
					mpfr_mul(tmp43, tmp43, a[1+1*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_sub(rop_tmp, rop_tmp, tmp43, MPFR_RNDN);
					mpfr_clears(tmp41, tmp42, tmp43, NULL);	
				}
				
				#pragma omp section
				{
					mpfr_init_set(tmp41, a[2+0*4], MPFR_RNDN);
					mpfr_mul(tmp41, tmp41, a[3+1*4], MPFR_RNDN);
					mpfr_init_set(tmp42, a[2+1*4], MPFR_RNDN); 
					mpfr_mul(tmp42, tmp42, a[3+0*4], MPFR_RNDN);
					mpfr_init(tmp43);
					mpfr_sub(tmp43, tmp41, tmp42, MPFR_RNDN);
					mpfr_mul(tmp43, tmp43, a[1+2*4], MPFR_RNDN);
					#pragma omp critical
						mpfr_add(rop_tmp, rop_tmp, tmp43, MPFR_RNDN);
					mpfr_clears(tmp41, tmp42, tmp43, NULL);	
				}
				
			}
			mpfr_mul(tmp, tmp, rop_tmp, MPFR_RNDN);
			#pragma omp critical
				mpfr_sub(rop, rop, tmp, MPFR_RNDN);	
			mpfr_clears(tmp, rop_tmp, NULL);
		}
	}

	return;
}

__MATHSUITE void det_5d(mpfr_t rop, mpfr_t *restrict a)
{
	mpfr_t *b;
	mpfr_t tmp;
	dim_typ i, j, k;
	dim_typ2 dim =
	{
		4,
		4
	};

	mpfr_set_ui(rop, 0, MPFR_RNDN);
	
	#pragma omp parallel for num_threads(5) private(i, j, tmp, b)
    for ( k = 0; k < 5; ++k)
    {
    	(void) matrixAlloc(&b, dim); // FAST
    	mpfr_init(tmp);
    	#pragma omp parallel for num_threads(4)
        for ( i = 0; i < 4; ++i )
        	#pragma omp parallel for num_threads(4)
            for ( j = 0; j < 4; ++j )
            	mpfr_set(b[i+(j<<2)], a[i+1+(j+(j>=k))*5], MPFR_RNDN);

		det_4d(tmp, b);
		mpfr_mul(tmp, tmp, a[0+k*5], MPFR_RNDN);
		#pragma omp critical
		{
			if(k%2)
				mpfr_sub(rop, rop, tmp, MPFR_RNDN);
			else
				mpfr_add(rop, rop, tmp, MPFR_RNDN);
		}
		matrixFree(&b, dim);
		mpfr_clear(tmp);
    }
    
	return;
}

// It calculates the determinant of a nXn matrix
// by upper-triangularizing the square matrix passed
__MATHSUITE  void det(mpfr_t rop, mpfr_t *restrict mat, dim_typ dimq, bool *flag)
{
	switch(dimq)
    {
        case 1:
            mpfr_set(rop, *mat, MPFR_RNDN); 
            return;
        case 2:
        {
			mpfr_t tmp;
			mpfr_init(tmp);
			mpfr_mul(rop, *mat, *(mat + dimq + 1), MPFR_RNDN);
			mpfr_mul(tmp, *(mat + 1), *(mat + dimq), MPFR_RNDN);
			mpfr_sub(rop, rop, tmp, MPFR_RNDN); 
			mpfr_clear(tmp); 
            return;
        }
        case 3:
    	{
            det_3d(rop, mat);
            return;
    	}
        case 4:
        {
            det_4d(rop, mat);
            return;
        }
        case 5:
        {
        	det_5d(rop, mat);
			return;
        }
    }
	
	mpfr_t D;
	sel_typ perms = 0;
	
    // Conversion of matrix to upper triangular
	mpfr_set_ui(D, 1, MPFR_RNDN); // storage for determinant

    if(matrixUTConv(mat, dimq, &perms))
    {
    	dim_typ i;
		mpfr_set_si(rop, powi(-1, perms), MPFR_RNDN);
        for(i=0; i < dimq; ++i)
        	mpfr_mul(rop, rop, *(mat + dimq*i + i), MPFR_RNDN);
    }
    else
    {
    	mpfr_t * S;
    	mpfr_t * V;

    	dim_typ2 dimqd =
    	{
    		dimq,
    		dimq
    	};

    	matrixAlloc(&S, (dim_typ2){1,dimq});
    	mpfr_t tmp;
    	mpfr_init(tmp);

    	if(matrixAlloc(&V, dimqd) && !_matrixEigenValues(mat, S, V, dimq))
		{
			mpfr_set_ui(tmp, dimq, MPFR_RNDN);
			mpfr_product(tmp, tmp, false, S);
    		mpfr_set(rop, tmp, MPFR_RNDN);
    	}
    	/*
    	{
    		matrixFree(&S, ((dim_typ2){1,dimq}));
    		matrixFree(&V, dimqd);
    		if(flag)
    			(*flag) = true;

			if(!matrixAlloc(&S, (dim_typ2){dimq,1}))
			{
	        	mpfr_set_d(rop, MAX_VAL, MPFR_RNDN);
	        	return;
	        }

	        if(!matrixAlloc(&V, dimqd))
	        {
	        	matrixFree(&V, dimqd);
	        	mpfr_set_d(rop, MAX_VAL, MPFR_RNDN);
	        }
	        else
	        {
	            dsvd(mat, dimqd, S, V);
	            mpfr_t tmp;
	            mpfr_set_ui(tmp, dimq, MPFR_RNDN);
	            mpfr_product(rop, tmp, false, S);
	            mpfr_abs(rop, rop, MPFR_RNDN);
	            mpfr_clear(tmp);
	        }
			matrixFree(&S, ((dim_typ2){dimq,1}));
		}
		*/
    	else
    	{
    		mpfr_t tmp, tmp2;
			mpfr_t *svd_vec;
			mpfr_t *Lm;
			mpfr_t *Rm;
			
    		matrixFree(&S, ((dim_typ2){1,dimq}));
    		matrixFree(&V, dimqd);
    		
    		if(flag)
    			(*flag) = true;
    			
    		if(!matrixAlloc(&Lm, dimqd))
       		{
	        	mpfr_set_d(rop, MAX_VAL, MPFR_RNDN);
	        	return;
	        }
	        
	        if(!matrixAlloc(&Rm, dimqd))
		    {
		    	mpfr_set_d(rop, MAX_VAL, MPFR_RNDN);
		    	matrixFree(&Lm, dimqd);
		        return;
		    }
		    
		    if(!matrixAlloc(&svd_vec, (dim_typ2){dimq,1}))
		    {
		    	mpfr_set_d(rop, MAX_VAL, MPFR_RNDN);
		    	matrixFree(&Lm, dimqd);
		    	matrixFree(&Rm, dimqd);
		        return;
		    }
	        else
	        {
			    mpfr_init_set_d(tmp, SVD_EPSILON, MPFR_RNDN);
				mpfr_init_set_d(tmp2, SVD_TOLERANCE, MPFR_RNDN);
	        	svd(dimqd, SVD_WITHU, SVD_WITHV, tmp, tmp2, mat, svd_vec, Lm, Rm);
	            // dsvd(mat, dimqd, S, V);
	            mpfr_set_ui(tmp, dimq, MPFR_RNDN);
	            mpfr_product(rop, tmp, false, svd_vec);
	            mpfr_abs(rop, rop, MPFR_RNDN);
	            mpfr_clears(tmp, tmp2, NULL);
	            matrixFree(&Lm, dimqd);
		    	matrixFree(&Rm, dimqd);
		    	matrixFree(&svd_vec, ((dim_typ2){dimq,1}));
	        }
			// matrixFree(&S, ((dim_typ2){dimq,1}));
		}
    }

    return;
}

__MATHSUITE void _matrixTrace(mpfr_t rop, mpfr_t *restrict mat, dim_typ dimq)
{
	mpfr_set_ui(rop, 0, MPFR_RNDN);

    for(dim_typ i=0; i<dimq; ++i)
    	mpfr_add(rop, rop, *(mat + dimq*i + i), MPFR_RNDN);

    return;
}

__MSSHELL_WRAPPER_ __MATHSUITE bool randomMatrix(mpfr_t * matrix, const register dim_typ dim[static 2])
{
	mpfr_t tmp;
    dim_typ range;
    
    msprintf(COLOR_CREDITS, "\n\nEnter a non-negative integer to set pseudo-random numbers Range.\n");
    PRINTHOWTOBACKMESSAGE();

    while(requires(tmp, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(tmp) || mpfr_cmp_ui(tmp, (range = mpfr_get_ui(tmp, MPFR_RNDN))) || range < 1 || range > SHRT_MAX)
    {
		mpfr_clear(tmp);
        CLEARBUFFER();
        if(access(exitHandle) == 1) continue;
        if(exitHandleCheck)
			return false;
        printErr(33, "Invalid inserted Value");
    }
    
    mpfr_clear(tmp);
    CLEARBUFFER();

    // RANDOMIZING THE MATRIX
    dim_typ i, j;

	#pragma omp parallel for
    for(i=0; i<dim[ROWS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[COLUMNS]; ++j)
            mpfr_set_d(*(matrix + dim[COLUMNS]*i + j), random(range), MPFR_RNDN);

    msprintf(COLOR_SYSTEM, "\n\n[%hu X %hu] randomized Matrix with Range: %hu is:\n", dim[ROWS], dim[COLUMNS], range);
    PRINTL();
    printMatrix(stdout, &matrix, dim);
    return true;
}

__MATHSUITE void transpose(mpfr_t *restrict matrix, mpfr_t *restrict matrix2, const register dim_typ dim[static 2])
{

    dim_typ i, j;
    #pragma omp parallel for
    for(i=0; i<dim[COLUMNS]; ++i)
    	#pragma omp parallel for
        for(j=0; j<dim[ROWS]; ++j)
        	mpfr_set(*(matrix2 + dim[ROWS]*i + j), *(matrix + dim[COLUMNS]*j + i), MPFR_RNDN); 

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

 __MATHSUITE bool  FattLU(dim_typ n, mpfr_t *restrict c, mpfr_t *restrict l, mpfr_t * a)
{
    mpfr_t m, pivot;
    dim_typ i, j, k;

    // Copia di C su A.
    (void) equalMatrix(&a, c, (dim_typ2){n, n}, EQUALMATRIX_NOREALLOC);

	mpfr_inits(m, pivot, NULL); 
    for (k=0; k<n; ++k)
    {
        mpfr_set(pivot, *(a + n*k + k), MPFR_RNDN); 
        if ( mpfr_zero_p(pivot) ) return false;
        for (i=k+1; i<n; ++i)
        {
        	mpfr_div(m, *(a + n*i + k), pivot, MPFR_RNDN);
        	mpfr_set(*(l + n*i + k), m, MPFR_RNDN);
            for (j=k; j<n; ++j)
			{
            	mpfr_mul(m, m, *(a + n*k + j), MPFR_RNDN);
            	mpfr_sub(*(a + n*i + j), *(a + n*i + j), m, MPFR_RNDN); 
        	}
        }
    }

  /* Adesso "a" contiene U, e "l" la parte inferiore di L */

	#pragma omp parallel for
    for (i=0; i<n; ++i) /* Completa L con zeri ed uni */
    {
        mpfr_set_ui(*(l + n*i + i), 1, MPFR_RNDN);
        #pragma omp parallel for
        for (j=i+1; j<n; ++j)
            mpfr_set_ui(*(l + n*i + j), 0, MPFR_RNDN);
    }
    
  mpfr_clears(m, pivot, NULL); 
  return true;
}

 __MATHSUITE sel_typ  invertMatrix(mpfr_t *restrict matrix, mpfr_t *restrict matrix2, dim_typ n)
{
    mpfr_t dt;
    mpfr_init(dt);
    det(dt, matrix, n, NULL);
    
    if(mpfr_cmp_ui(dt, 0) == 0)
    	return INVERTMATRIX_SINGULAR;
    if(!adjoint(matrix, matrix2, n))
    	return INVERTMATRIX_ALLOCERROR;
    
    dim_typ i, j;
    
    #pragma omp parallel for
    for(i=0; i<n; ++i)	
    	#pragma omp parallel for
    	for(j=0; j<n; ++j)	
    		mpfr_div(*(matrix + n*i + j), *(matrix + n*i + j), dt, MPFR_RNDN);
    	
    	
    mpfr_clear(dt); 
	return INVERTMATRIX_SUCCESS;
}

 __MATHSUITE sel_typ  invertMatrixFast(mpfr_t *restrict matrix, mpfr_t * matrix2, dim_typ n)
 {

	dim_typ i, j, k;
	mpfr_t * matrix3;
	mpfr_t temp, tmp_prod;
	dim_typ2 dims =
	{
		n,
		n
	};
	
	for(i=0;i<n;++i)							
		for(j=0;j<n;++j)	
		{
			// WORST CASE: Use Backup Naive function
			if(i == j && mpfr_zero_p(*(matrix + n*i + i)))
			{
				matrixFree(&matrix2, dims);
				if(!matrixAlloc(&matrix2, dims))
					return INVERTMATRIX_ALLOCERROR;
					
				if(!equalMatrix(&matrix3, matrix, dims, EQUALMATRIX_REALLOC))
					return INVERTMATRIX_ALLOCERROR;
					
				const register sel_typ exitVal = invertMatrix(matrix3, matrix2, n);
					
				if(exitVal != INVERTMATRIX_SUCCESS)
					return exitVal;
					
				// matrixFree(&matrix2, dims);
				// matrix2 = matrix3;
				(void) equalMatrix(&matrix2, matrix3, dims, EQUALMATRIX_NOREALLOC);
				matrixFree(&matrix3, dims);
				return INVERTMATRIX_SUCCESS;
			}
			mpfr_set_ui(*(matrix2 + n*i + j), i == j, MPFR_RNDN); 
		}
	
	mpfr_inits(temp, tmp_prod, NULL);
	
	for(k=0;k<n;++k)	
	{							
		mpfr_set(temp, *(matrix + n*k + k), MPFR_RNDN);	
										
		for(j=0;j<n;++j)				
		{
			mpfr_div(*(matrix + n*k + j), *(matrix + n*k + j), temp, MPFR_RNDN);
			mpfr_div(*(matrix2 + n*k + j), *(matrix2 + n*k + j), temp, MPFR_RNDN);
		}	
									
		for(i=0;i<n;++i)
		{
			mpfr_set(temp, *(matrix + n*i + k), MPFR_RNDN);	
			for(j=0;j<n;++j)
			{
				if(i==k)
					break;
				mpfr_mul(tmp_prod, temp, *(matrix + k*n + j), MPFR_RNDN);
				mpfr_sub(*(matrix + n*i + j), *(matrix + n*i + j), tmp_prod, MPFR_RNDN);
				mpfr_mul(tmp_prod, temp, *(matrix2 + k*n + j), MPFR_RNDN);
				mpfr_sub(*(matrix2 + n*i + j), *(matrix2 + n*i + j), tmp_prod, MPFR_RNDN);
			}
		}	
	}
	
	mpfr_clears(temp, tmp_prod, NULL);
	return INVERTMATRIX_SUCCESS;
}

/*
   Find the cofactor matrix of a square matrix
*/
 bool  CoFactor(mpfr_t *restrict a, mpfr_t *restrict b, dim_typ n)
{
	dim_typ i,j,ii,jj,i1,j1;
	const register dim_typ nminus1 = n-1;
	mpfr_t _det;
	mpfr_t *c = NULL;

	if(!matrixAlloc(&c, (dim_typ2){nminus1,nminus1}))
		return false;
		
	mpfr_init(_det);

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
	           		mpfr_set(*(c + nminus1*i1 + j1), *(a + n*ii + jj), MPFR_RNDN);
	           		++ j1;
	        	}
        		++ i1;
     		}
		     /* Calculate the determinant */
		     det(_det, c, n-1, NULL);
	     	/* Fill in the elements of the cofactor */
	     	mpfr_mul_si(*(b + n*i + j), _det, pow(-1.0,i+j+2.0), MPFR_RNDN);
  		}
	
	mpfr_clear(_det);
	matrixFree(&c, ((dim_typ2){nminus1,nminus1}));
	return true;
}

 inline bool   adjoint(mpfr_t *restrict a, mpfr_t *restrict b, dim_typ n)
{
	if(CoFactor(a,b,n))
	{
		transpose(b,a,(dim_typ2){n,n});
		return true;
	}
	return false;
}

 static inline   void SIGN(mpfr_t res, mpfr_t a, mpfr_t b)
{
	mpfr_abs(res, a, MPFR_RNDN);
	if(mpfLZero(b))
		mpfr_neg(res, res, MPFR_RNDN);
	return;
}

 static inline   void PYTHAG(mpfr_t res, mpfr_t a, mpfr_t b)
{
	mpfr_t at, bt;
	mpfr_inits(at, bt, NULL);
	mpfr_abs(at, a, MPFR_RNDN);
	mpfr_abs(bt, b, MPFR_RNDN);
	if(mpfr_greater_p(at, bt))
	{
		mpfr_div(res, bt, at, MPFR_RNDN);
		mpfr_mul(res, res, res, MPFR_RNDN);
		mpfr_add_ui(res, res, 1, MPFR_RNDN);
		mpfr_sqrt(res, res, MPFR_RNDN);
		mpfr_mul(res, at, res, MPFR_RNDN);
	}
	else if(mpfGZero(bt))
	{
		mpfr_div(res, at, bt, MPFR_RNDN);
		mpfr_mul(res, res, res, MPFR_RNDN);
		mpfr_add_ui(res, res, 1, MPFR_RNDN);
		mpfr_sqrt(res, res, MPFR_RNDN);
		mpfr_mul(res, bt, res, MPFR_RNDN);
	}
	return;
}

/*  svd.c -- Singular value decomposition. Translated to 'C' from the
 *           original ALGOL code in "Handbook for Automatic Computation,
 *           vol. II, Linear Algebra", Springer-Verlag.
 *
 *  (C) 2000, C. Bond. All rights reserved.
 *
 *  This is almost an exact translation from the original, except that
 *  an iteration counter is added to prevent stalls. This corresponds
 *  to similar changes in other translations.
 *
 *  Returns an error code = 0, if no errors and 'k' if a failure to
 *  converge at the 'kth' singular value.
 * 
 */
 // Rewritten in MPFR by me.
 dim_typ  svd(const register dim_typ dim[static 2],const bool withu,const bool withv,mpfr_t eps,mpfr_t tol,
	mpfr_t *restrict a,mpfr_t *restrict q,mpfr_t *restrict u, mpfr_t *restrict v)
{
	
	const register dim_typ n = dim[ROWS];
	const register dim_typ m = dim[COLUMNS];
	
	if(n < m)
    {
        mpfr_t * transp;
        if(!matrixAlloc(&transp, (dim_typ2){m, n}))
            return false;
        transpose(a, transp, dim); 
        dim_typ ret = svd((dim_typ2){m, n}, withu, withv, eps, tol, transp, q, u, v);
        matrixFree(&transp, ((dim_typ2){m, n}));
        return ret;
    }
 
	dim_typ iter,retval=0;
	mpfr_t tmp, tmp2, tmp3, tmp4, i,j,k,l,l1,c,f,g,h,s,x,y,z;
	mpfr_inits(tmp, tmp2, tmp3, tmp4, i,j,k,l,l1,c,f,g,h,s,x,y,z, NULL);
	mpfr_t *e;
	(void) matrixAlloc(&e, (dim_typ2){n,1}); // it shouldn't be done such a thing, but this function is assumed to be extremely FAST
	
	/* Copy 'a' to 'u' */    
	for (mpfr_set_ui(i, 0, MPFR_RNDN); mpfr_cmp_ui(i,m) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
		for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
			mpfr_set(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], a[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], MPFR_RNDN);
			
	/* Householder's reduction to bidiagonal form. */
	mpfr_set_ui(g, 0, MPFR_RNDN);
	mpfr_set_ui(x, 0, MPFR_RNDN);   
	for (mpfr_set_ui(i, 0, MPFR_RNDN); mpfr_cmp_ui(i,n) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
	{
		mpfr_set(e[mpfr_get_ui(i, MPFR_RNDN)], g, MPFR_RNDN);
		mpfr_set_ui(s, 0, MPFR_RNDN);
		mpfr_add_ui(l, i, 1, MPFR_RNDN);
		for (mpfr_set(j, i, MPFR_RNDN); mpfr_cmp_ui(j,m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
		{
			mpfr_mul(tmp, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN); 
			mpfr_add(s, s, tmp, MPFR_RNDN);
		}
		if (mpfr_less_p(s, tol))
			mpfr_set_ui(g, 0, MPFR_RNDN);
		else 
		{
			mpfr_set(f, u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
			mpfr_sqrt(g, s, MPFR_RNDN);
			if(mpfGeZero(f))
				mpfr_neg(g, g, MPFR_RNDN);
			mpfr_mul(h, f, g, MPFR_RNDN);
			mpfr_sub(h, h, s, MPFR_RNDN);
			mpfr_sub(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], f, g, MPFR_RNDN); 
			for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
			{
				mpfr_set_ui(s, 0, MPFR_RNDN);
				for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k, m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
				{
					mpfr_mul(tmp, u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], MPFR_RNDN);
					mpfr_add(s, s, tmp, MPFR_RNDN);
				}
				mpfr_div(f, s, h, MPFR_RNDN);
				for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k, m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
				{
					mpfr_mul(tmp, f, u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
					mpfr_add(u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], tmp, MPFR_RNDN);
				}
			} /* end j */
		} /* end s */
		mpfr_set(q[mpfr_get_ui(i, MPFR_RNDN)], g, MPFR_RNDN);
		mpfr_set_ui(s, 0, MPFR_RNDN);
		for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
		{
			mpfr_mul(tmp, u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], MPFR_RNDN);
			mpfr_add(s, s, tmp, MPFR_RNDN);
		}
		if (mpfr_less_p(s, tol))
			mpfr_set_ui(g, 0, MPFR_RNDN);
		else
		{
			mpfr_set(f, u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)+1], MPFR_RNDN);
			mpfr_sqrt(g, s, MPFR_RNDN);
			if(mpfGeZero(f))
				mpfr_neg(g, g, MPFR_RNDN);
			mpfr_mul(h, f, g, MPFR_RNDN);
			mpfr_sub(h, h, s, MPFR_RNDN);
			mpfr_sub(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)+1], f, g, MPFR_RNDN);
			for (mpfr_set(j,l,MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
				mpfr_div(e[mpfr_get_ui(j, MPFR_RNDN)], u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], h, MPFR_RNDN);
			for (mpfr_set(j,l,MPFR_RNDN); mpfr_cmp_ui(j, m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
			{
				mpfr_set_ui(s, 0, MPFR_RNDN);
				for (mpfr_set(k,l,MPFR_RNDN); mpfr_cmp_ui(k,n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
				{
					mpfr_mul(tmp, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(k, MPFR_RNDN)], u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
					mpfr_add(s, s, tmp, MPFR_RNDN);
				}
				for (mpfr_set(k,l,MPFR_RNDN); mpfr_cmp_ui(k,n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
				{
					mpfr_mul(tmp, s, e[mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
					mpfr_add(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(k, MPFR_RNDN)], u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(k, MPFR_RNDN)], tmp, MPFR_RNDN);
				}
			} /* end j */
		} /* end s */
		mpfr_abs(tmp, e[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
		mpfr_abs(y, q[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
		mpfr_add(y, y, tmp, MPFR_RNDN);
		if(mpfr_greater_p(y, x))
			mpfr_set(x, y, MPFR_RNDN);
	} /* end i */
	
	/* accumulation of right-hand transformations */
	if (withv)
	{
		for (mpfr_set_ui(i, n-1, MPFR_RNDN); mpfGeZero(i); mpfr_sub_ui(i, i, 1, MPFR_RNDN))
		{
			if (!mpfr_zero_p(g)) 
			{
				mpfr_mul(h, u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)+1], g, MPFR_RNDN);
				for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
					mpfr_div(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], h, MPFR_RNDN);
				for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
				{
					mpfr_set_ui(s, 0, MPFR_RNDN);
					for (mpfr_set(k,l,MPFR_RNDN); mpfr_cmp_ui(k,n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
					{
						mpfr_mul(tmp, u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(k, MPFR_RNDN)], v[mpfr_get_ui(k, MPFR_RNDN)*m + mpfr_get_ui(j, MPFR_RNDN)], MPFR_RNDN);
						mpfr_add(s, s, tmp, MPFR_RNDN);
					}

					for (mpfr_set(k,l,MPFR_RNDN); mpfr_cmp_ui(k,n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
					{
						mpfr_mul(tmp, s, v[mpfr_get_ui(k, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
						mpfr_add(v[mpfr_get_ui(k, MPFR_RNDN)*m + mpfr_get_ui(j, MPFR_RNDN)], v[mpfr_get_ui(k, MPFR_RNDN)*m + mpfr_get_ui(j, MPFR_RNDN)], tmp, MPFR_RNDN);
					}	
				} /* end j */
			} /* end g */
			for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
			{
				mpfr_set_ui(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], 0, MPFR_RNDN);
				mpfr_set_ui(v[mpfr_get_ui(i, MPFR_RNDN)*m + mpfr_get_ui(j, MPFR_RNDN)], 0, MPFR_RNDN);
			}
			
			mpfr_set_ui(v[mpfr_get_ui(i, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], 1, MPFR_RNDN);
			mpfr_set(g, e[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
			mpfr_set(l, i, MPFR_RNDN);
		} /* end i */
	} /* end withv, parens added for clarity */
	
	/* accumulation of left-hand transformations */
	if (withu) 
	{
		for (mpfr_set_ui(i,n,MPFR_RNDN); mpfr_cmp_ui(i, m) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
		{
			for (mpfr_set_ui(j, n, MPFR_RNDN); mpfr_cmp_ui(j, m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
				mpfr_set_ui(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], 0, MPFR_RNDN);
			mpfr_set_ui(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], 1, MPFR_RNDN);
		}
		for (mpfr_set_ui(i, n-1, MPFR_RNDN); mpfGeZero(i); mpfr_sub_ui(i, i, 1, MPFR_RNDN))
		{
			mpfr_add_ui(l, i, 1, MPFR_RNDN);
			mpfr_set(g, q[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
			for (mpfr_set(j,l,MPFR_RNDN); mpfr_cmp_ui(j, m) < 0;mpfr_add_ui(j, j, 1, MPFR_RNDN))  /* upper limit was 'n' */
				mpfr_set_ui(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], 0, MPFR_RNDN);
			if (!mpfr_zero_p(g)) 
			{
				mpfr_mul(h, u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], g, MPFR_RNDN);
				for (mpfr_set(j,l,MPFR_RNDN); mpfr_cmp_ui(j, m) < 0;mpfr_add_ui(j, j, 1, MPFR_RNDN))
				{ /* upper limit was 'n' */
					mpfr_set_ui(s, 0, MPFR_RNDN);
					for (mpfr_set(k,l,MPFR_RNDN); mpfr_cmp_ui(k, m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
					{
						mpfr_mul(tmp, u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], MPFR_RNDN);
						mpfr_add(s, s, tmp, MPFR_RNDN);
					}
					mpfr_div(f, s, h, MPFR_RNDN);
					for (mpfr_set(k,i,MPFR_RNDN); mpfr_cmp_ui(k, m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
					{
						mpfr_mul(tmp, f, u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
						mpfr_add(u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], u[mpfr_get_ui(k, MPFR_RNDN)*n + mpfr_get_ui(j, MPFR_RNDN)], tmp, MPFR_RNDN);
					}
				} /* end j */
				for (mpfr_set(j,i,MPFR_RNDN); mpfr_cmp_ui(j, m) < 0;mpfr_add_ui(j, j, 1, MPFR_RNDN))
					mpfr_div(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], g, MPFR_RNDN);
			} /* end g */
			else 
				for (mpfr_set(j,i,MPFR_RNDN); mpfr_cmp_ui(j, m) < 0;mpfr_add_ui(j, j, 1, MPFR_RNDN))
					mpfr_set_ui(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], 0, MPFR_RNDN);
			mpfr_add_ui(u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], u[mpfr_get_ui(i, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], 1, MPFR_RNDN);
		} /* end i*/
	} /* end withu, parens added for clarity */
	
	/* diagonalization of the bidiagonal form */
	mpfr_mul(eps, eps, x, MPFR_RNDN);
	for (mpfr_set_ui(k, n-1, MPFR_RNDN); mpfGeZero(k); mpfr_sub_ui(k, k, 1, MPFR_RNDN))
	{
		iter = 0;
		test_f_splitting:
			for (mpfr_set(l,k, MPFR_RNDN); mpfGeZero(l); mpfr_sub_ui(l, l, 1, MPFR_RNDN))
			{
				mpfr_abs(tmp, e[mpfr_get_ui(l, MPFR_RNDN)], MPFR_RNDN);
				if(mpfr_lessequal_p(tmp, eps)) goto test_f_convergence;
				mpfr_abs(tmp, q[mpfr_get_ui(l, MPFR_RNDN)-1], MPFR_RNDN);
				if(mpfr_lessequal_p(tmp, eps)) goto cancellation;
			} /* end l */
		
		/* cancellation of e[l] if l > 0 */
		cancellation:
			mpfr_set_ui(c, 0, MPFR_RNDN);
			mpfr_set_ui(s, 1, MPFR_RNDN);
			mpfr_sub_ui(l1, l, 1, MPFR_RNDN);
			for (mpfr_set(i, l, MPFR_RNDN); mpfr_lessequal_p(i, k); mpfr_add_ui(i, i, 1, MPFR_RNDN))
			{
				mpfr_mul(f, s, e[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
				mpfr_mul(e[mpfr_get_ui(i, MPFR_RNDN)], e[mpfr_get_ui(i, MPFR_RNDN)], c, MPFR_RNDN);
				mpfr_abs(tmp, f, MPFR_RNDN);
				if(mpfr_lessequal_p(tmp, eps)) goto test_f_convergence;
				mpfr_set(g, q[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
				mpfr_mul(tmp, g, g, MPFR_RNDN);
				mpfr_mul(h, f, f, MPFR_RNDN);
				mpfr_add(h, h, tmp, MPFR_RNDN);
				mpfr_sqrt(h, h, MPFR_RNDN); 
				mpfr_set(q[mpfr_get_ui(i, MPFR_RNDN)], h, MPFR_RNDN);
				mpfr_div(s, f, h, MPFR_RNDN);
				mpfr_neg(s, s, MPFR_RNDN);
				mpfr_div(c, g, h, MPFR_RNDN);
				if (withu) 
				{
					for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j, m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
					{
						mpfr_set(y, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(l1, MPFR_RNDN)], MPFR_RNDN);
						mpfr_set(z, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
						mpfr_mul(tmp, z, s, MPFR_RNDN);
						mpfr_mul(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(l1, MPFR_RNDN)], y, c, MPFR_RNDN);
						mpfr_add(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(l1, MPFR_RNDN)], u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(l1, MPFR_RNDN)], tmp, MPFR_RNDN);
						mpfr_mul(tmp, z, c, MPFR_RNDN);
						mpfr_mul(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], y, s, MPFR_RNDN);
						mpfr_sub(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], tmp, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
					} /* end j */
				} /* end withu, parens added for clarity */
			} /* end i */
		test_f_convergence:
			mpfr_set(z, q[mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
			if (mpfr_equal_p(l, k)) goto convergence;

			/* shift from bottom 2x2 minor */
			if (++iter > 31) 
			{
				retval = mpfr_get_ui(k, MPFR_RNDN);
				mpfr_clears(tmp, tmp2, tmp3, tmp4, i,j,k,l,l1,c,f,g,h,s,x,y,z, NULL);
				matrixFree(&e, ((dim_typ2){n,1}));
				return retval;
			}
			mpfr_set(x, q[mpfr_get_ui(l, MPFR_RNDN)], MPFR_RNDN);
			mpfr_set(y, q[mpfr_get_ui(k, MPFR_RNDN)-1], MPFR_RNDN);
			mpfr_set(g, e[mpfr_get_ui(k, MPFR_RNDN)-1], MPFR_RNDN);
			mpfr_set(h, e[mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
			mpfr_mul_ui(tmp, h, 2, MPFR_RNDN);
			mpfr_mul(tmp, tmp, y, MPFR_RNDN);
			mpfr_add(tmp2, g, h, MPFR_RNDN);
			mpfr_sub(tmp3, g, h, MPFR_RNDN);
			mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);
			mpfr_add(tmp3, y, z, MPFR_RNDN);
			mpfr_sub(tmp4, y, z, MPFR_RNDN);
			mpfr_mul(tmp3, tmp3, tmp4, MPFR_RNDN);
			mpfr_add(tmp2, tmp2, tmp3, MPFR_RNDN);
			mpfr_div(f, tmp2, tmp, MPFR_RNDN);
			mpfr_mul(g, f, f, MPFR_RNDN);
			mpfr_add_ui(g, g, 1, MPFR_RNDN);
			mpfr_sqrt(g, g, MPFR_RNDN);
			if(f < 0)
				mpfr_sub(tmp, f, g, MPFR_RNDN); 
			else
				mpfr_add(tmp, f, g, MPFR_RNDN); 
			mpfr_div(tmp, y, tmp, MPFR_RNDN);
			mpfr_sub(tmp, tmp, h, MPFR_RNDN);
			mpfr_mul(tmp, tmp, h, MPFR_RNDN);
			mpfr_add(tmp2, x, z, MPFR_RNDN);
			mpfr_sub(tmp3, x, z, MPFR_RNDN);
			mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);
			mpfr_add(tmp, tmp, tmp2, MPFR_RNDN);
			mpfr_div(f, tmp, x, MPFR_RNDN);
			mpfr_set_ui(s, 1, MPFR_RNDN);
			
			/* next QR transformation */
			mpfr_set_ui(c, 1, MPFR_RNDN);
			
			for (mpfr_add_ui(i, l, 1, MPFR_RNDN); mpfr_lessequal_p(i,k); mpfr_add_ui(i, i, 1, MPFR_RNDN))
			{
				mpfr_set(g, e[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
				mpfr_set(y, q[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
				mpfr_mul(h, s, g, MPFR_RNDN);
				mpfr_mul(g, g, c, MPFR_RNDN);
				mpfr_mul(tmp, h, h, MPFR_RNDN);
				mpfr_mul(z, f, f, MPFR_RNDN);
				mpfr_add(z, z, tmp, MPFR_RNDN);
				mpfr_sqrt(z, z, MPFR_RNDN);
				mpfr_set(e[mpfr_get_ui(i, MPFR_RNDN)-1], z, MPFR_RNDN);
				mpfr_div(c, f, z, MPFR_RNDN);
				mpfr_div(s, h, z, MPFR_RNDN);
				mpfr_mul(tmp, g, s, MPFR_RNDN);
				mpfr_mul(f, x, c, MPFR_RNDN);
				mpfr_add(f, f, tmp, MPFR_RNDN);
				mpfr_mul(tmp, g, c, MPFR_RNDN);
				mpfr_mul(g, x, s, MPFR_RNDN);
				mpfr_neg(g, g, MPFR_RNDN);
				mpfr_add(g, g, tmp, MPFR_RNDN);
				mpfr_mul(h, y, s, MPFR_RNDN);
				mpfr_mul(y, y, c, MPFR_RNDN);
				
				if (withv) 
				{
					for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
					{
						mpfr_set(x, v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)-1], MPFR_RNDN);
						mpfr_set(z, v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
						mpfr_mul(tmp, z, s, MPFR_RNDN);
						mpfr_mul(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)-1], x, c, MPFR_RNDN);
						mpfr_add(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)-1], v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)-1], tmp, MPFR_RNDN);
						mpfr_mul(tmp, z, c, MPFR_RNDN);
						mpfr_mul(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], x, s, MPFR_RNDN);
						mpfr_sub(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], tmp, v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
					} /* end j */
				} /* end withv, parens added for clarity */
				mpfr_mul(tmp, h, h, MPFR_RNDN);
				mpfr_mul(z, f, f, MPFR_RNDN);
				mpfr_add(z, z, tmp, MPFR_RNDN);
				mpfr_sqrt(z, z, MPFR_RNDN);
				mpfr_set(q[mpfr_get_ui(i, MPFR_RNDN)-1], z, MPFR_RNDN);
				mpfr_div(c, f, z, MPFR_RNDN);
				mpfr_div(s, h, z, MPFR_RNDN);
				mpfr_mul(tmp, s, y, MPFR_RNDN);
				mpfr_mul(f, c, g, MPFR_RNDN);
				mpfr_add(f, f, tmp, MPFR_RNDN);
				mpfr_mul(tmp, c, y, MPFR_RNDN);
				mpfr_mul(x, s, g, MPFR_RNDN);
				mpfr_sub(x, tmp, x, MPFR_RNDN);
				if (withu) 
				{
					for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j, m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
					{
						mpfr_set(y, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)-1], MPFR_RNDN);
						mpfr_set(z, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
						mpfr_mul(tmp, z, s, MPFR_RNDN);
						mpfr_mul(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)-1], y, c, MPFR_RNDN);
						mpfr_add(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)-1], u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)-1], tmp, MPFR_RNDN);
						mpfr_mul(tmp, z, c, MPFR_RNDN);
						mpfr_mul(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], y, s, MPFR_RNDN);
						mpfr_sub(u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], tmp, u[mpfr_get_ui(j, MPFR_RNDN)*n + mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
					} /* end j */
				} /* end withu, parens added for clarity */
			} /* end i */
			mpfr_set_ui(e[mpfr_get_ui(l, MPFR_RNDN)], 0, MPFR_RNDN);
			mpfr_set(e[mpfr_get_ui(k, MPFR_RNDN)], f, MPFR_RNDN);
			mpfr_set(q[mpfr_get_ui(k, MPFR_RNDN)], x, MPFR_RNDN);
			goto test_f_splitting;
		convergence:
			if (mpfLZero(z)) 
			{
				/* q[k] is made non-negative */
				mpfr_neg(q[mpfr_get_ui(k,MPFR_RNDN)], z, MPFR_RNDN);
				if (withv) 
				{
					for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
						mpfr_neg(v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(k, MPFR_RNDN)], v[mpfr_get_ui(j, MPFR_RNDN)*m + mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
				} /* end withv, parens added for clarity */
			} /* end z */
	} /* end k */

	mpfr_clears(tmp, tmp2, tmp3, tmp4, i,j,k,l,l1,c,f,g,h,s,x,y,z, NULL);
	matrixFree(&e, ((dim_typ2){n,1}));
	return retval;
}

// Jacobi Singular Value Decomposition (SVD)
// written by me
 bool  dsvd(mpfr_t *a, const register dim_typ dim[static 2], mpfr_t *w, mpfr_t *v)
{
    // checking whether m < n and correct it
    // by transposing the matrix into n,m matrix
    // causes its correctly dsvd decomposition...
	const register dim_typ m = dim[ROWS];
    const register dim_typ n = dim[COLUMNS];

    if(m < n)
    {
        mpfr_t * transp;
        if(!matrixAlloc(&transp, (dim_typ2){n, m}))
            return false;
        transpose(a, transp, dim); 
        bool ret = dsvd(transp, (dim_typ2){n, m}, w, v);
        matrixFree(&transp, ((dim_typ2){n, m}));
        return ret;
    }
    
    bool flag;
    mpfr_t *rv1;
    mpfr_t tmp, tmp2, tmp3, tmp4, anorm, g, scale, i, its, j, jj, k, l, nm, c, f, h, s, x, y, z;
    
    if(!matrixAlloc(&rv1, (dim_typ2){n,1}))
    	return false;

    mpfr_inits(tmp, tmp2, tmp3, tmp4, anorm, g, scale, its, j, jj, k, l, nm, c, f, h, s, x, y, z, NULL);
    
    mpfr_set_ui(scale, 0, MPFR_RNDN); 
    mpfr_set_ui(s, 0, MPFR_RNDN);
	mpfr_set_ui(g, 0, MPFR_RNDN);
	
/* Householder reduction to bidiagonal form */
    for (mpfr_init_set_ui(i, 0, MPFR_RNDN); mpfr_cmp_ui(i, n) < 0; mpfr_add_ui(i, i, 1, MPFR_RNDN))
    {
        /* left-hand reduction */
        mpfr_add_ui(l, i, 1, MPFR_RNDN);
		mpfr_mul(rv1[mpfr_get_ui(i, MPFR_RNDN)], scale, g, MPFR_RNDN);
		mpfr_set_ui(scale, 0, MPFR_RNDN);
		mpfr_set_ui(s, 0, MPFR_RNDN);
		mpfr_set_ui(g, 0, MPFR_RNDN);
        if (mpfr_cmp_ui(i,m) < 0)
        {
            for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
            {
				mpfr_abs(tmp, *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
            	mpfr_add(scale, scale, tmp, MPFR_RNDN);
            }
            if (!mpfr_zero_p(scale))
            {
                for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                {
                    mpfr_div(*(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), scale, MPFR_RNDN);
                    mpfr_mul(tmp, *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
                    mpfr_add(s, s, tmp, MPFR_RNDN);
                }
                mpfr_set(f, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
				mpfr_sqrt(tmp, s, MPFR_RNDN);
				SIGN(g, tmp, f);
				mpfr_neg(g, g, MPFR_RNDN);
				mpfr_mul(tmp, f, g, MPFR_RNDN);
				mpfr_sub(h, tmp, s, MPFR_RNDN);
				mpfr_sub(*(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), f, g, MPFR_RNDN);
				mpfr_set_ui(tmp, n-1, MPFR_RNDN);
                if (mpfr_cmp(i, tmp))
                {
                	for(mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j, n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                    {
                        for (mpfr_set_ui(s, 0, MPFR_RNDN), mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                        {
                        	mpfr_mul(tmp, *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), MPFR_RNDN);
                        	mpfr_add(s, s, tmp, MPFR_RNDN);
                        }
                        mpfr_div(f, s, h, MPFR_RNDN);
                        for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
						{
							mpfr_mul(tmp, f, *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
							mpfr_add(*(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), tmp, MPFR_RNDN);
						}
                    }
                }
                for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                	mpfr_mul(*(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), scale, MPFR_RNDN);
            }
        }
        mpfr_mul(w[mpfr_get_ui(i, MPFR_RNDN)], scale, g, MPFR_RNDN);

        /* right-hand reduction */
        mpfr_set_ui(g, 0, MPFR_RNDN);
        mpfr_set_ui(s, 0, MPFR_RNDN);
        mpfr_set_ui(scale, 0, MPFR_RNDN);
        mpfr_set_ui(tmp, n-1, MPFR_RNDN);

        if (mpfr_cmp_ui(i,m) < 0 && mpfr_cmp(i, tmp))
        {
            for (mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k, n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
            {
            	mpfr_abs(tmp, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), MPFR_RNDN);
            	mpfr_add(scale, scale, tmp, MPFR_RNDN);
            }
            
            if (!mpfr_zero_p(scale))
            {
                for (mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k, n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                {
                	mpfr_div(*(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), scale, MPFR_RNDN);
                	mpfr_mul(tmp, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), MPFR_RNDN);
                	mpfr_add(s, s, tmp, MPFR_RNDN);
                }
                mpfr_set(f, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(l, MPFR_RNDN)), MPFR_RNDN);
                mpfr_sqrt(tmp, s, MPFR_RNDN);
                SIGN(g, tmp, f);
                mpfr_neg(g, g, MPFR_RNDN);
                mpfr_mul(h, f, g, MPFR_RNDN);
                mpfr_sub(h, h, s, MPFR_RNDN);
                mpfr_sub(*(a + n*mpfr_get_ui(i,MPFR_RNDN) + mpfr_get_ui(l, MPFR_RNDN)), f, g, MPFR_RNDN);
                for (mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k,n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                	mpfr_div(rv1[mpfr_get_ui(k, MPFR_RNDN)], *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), h, MPFR_RNDN);
                mpfr_set_ui(tmp, m-1, MPFR_RNDN);
                if (mpfr_cmp(i, tmp))
                {
                    for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                    {
                        for (mpfr_set_ui(s, 0, MPFR_RNDN), mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k, n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                        {
                        	mpfr_mul(tmp, *(a + n*mpfr_get_ui(j, MPFR_RNDN) +mpfr_get_ui(k, MPFR_RNDN)), *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), MPFR_RNDN);
                        	mpfr_add(s, s, tmp, MPFR_RNDN);
                        }
                        for (mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k, n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                    	{
                    		mpfr_mul(tmp, s, rv1[mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
                    		mpfr_add(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), *(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), tmp, MPFR_RNDN);
                    	}
                    }
                }
                for (mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k, n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                	mpfr_mul(*(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), scale, MPFR_RNDN);
            }
        }
        mpfr_t reg;
        mpfr_init(reg);
        mpfr_abs(tmp, rv1[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
        mpfr_abs(reg, w[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
        mpfr_add(reg, reg, tmp, MPFR_RNDN);
        if(mpfr_less_p(reg, anorm))
            mpfr_set(anorm, reg, MPFR_RNDN);
    	mpfr_clear(reg);
    }

    /* accumulate the right-hand transformation */
    for (mpfr_set_ui(i, n-1, MPFR_RNDN); mpfGeZero(i); mpfr_sub_ui(i, i, 1, MPFR_RNDN))
    {
    	mpfr_set_ui(tmp, n-1, MPFR_RNDN);
        if (mpfr_less_p(i, tmp))
        {
            if (!mpfr_zero_p(g))
            {
                for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j, n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                {
                	mpfr_div(tmp, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(l, MPFR_RNDN)), g, MPFR_RNDN);
                	mpfr_div(*(v + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), tmp, MPFR_RNDN);
                }
                    /* double division to avoid underflow */
                for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j, n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                {
                    for (mpfr_set_ui(s, 0, MPFR_RNDN), mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k,n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                    {
                    	mpfr_mul(tmp, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), *(v + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), MPFR_RNDN);
                    	mpfr_add(s, s, tmp, MPFR_RNDN);
                    }
                    for (mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k, n) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                    {
                    	mpfr_mul(tmp, s, *(v + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
                    	mpfr_add(*(v + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), *(v + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), tmp, MPFR_RNDN); 	
                    }
                }
            }
            for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
            {
            	mpfr_set_ui(*(v + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), 0, MPFR_RNDN);
            	mpfr_set_ui(*(v + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), 0, MPFR_RNDN);
            }
        }
        mpfr_set_ui(*(v + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), 1, MPFR_RNDN);
        mpfr_set(g, rv1[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
        mpfr_set(l, i, MPFR_RNDN);
    }

    /* accumulate the left-hand transformation */
    for (mpfr_set_ui(i, n-1, MPFR_RNDN); mpfGeZero(i); mpfr_sub_ui(i, i, 1, MPFR_RNDN))
    {
    	mpfr_add_ui(l, i, 1, MPFR_RNDN);
    	mpfr_set(g, w[mpfr_get_ui(i,MPFR_RNDN)], MPFR_RNDN);
    	mpfr_set_ui(tmp, n-1, MPFR_RNDN);
        if (mpfr_less_p(i, tmp))
            for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
            	mpfr_set_ui(*(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), 0, MPFR_RNDN);
        if (!mpfr_zero_p(g))
        {
        	mpfr_pow_si(g, g, -1, MPFR_RNDN);
			mpfr_set_ui(tmp, n-1, MPFR_RNDN);
            if (mpfr_cmp(i, tmp))
            {
                for (mpfr_set(j, l, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                {
                    for (mpfr_set_ui(s, 0, MPFR_RNDN), mpfr_set(k, l, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                	{
                		mpfr_mul(tmp, *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), MPFR_RNDN);
                		mpfr_add(s, s, tmp, MPFR_RNDN);
                	}
					mpfr_div(f, s, *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
					mpfr_mul(f, f, g, MPFR_RNDN);
                    for (mpfr_set(k, i, MPFR_RNDN); mpfr_cmp_ui(k,m) < 0; mpfr_add_ui(k, k, 1, MPFR_RNDN))
                    {
                    	mpfr_mul(tmp, f, *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
                    	mpfr_add(*(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), *(a + n*mpfr_get_ui(k, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), tmp, MPFR_RNDN);
                    }
                }
            }
            for (mpfr_set(j, i, MPFR_RNDN); mpfr_cmp_ui(j,m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
            	mpfr_mul(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), g, MPFR_RNDN);
        }
        else
        {
            for(mpfr_set(j, i, MPFR_RNDN); mpfr_cmp_ui(j,m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
            	mpfr_set_ui(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), 0, MPFR_RNDN);
        }
        mpfr_add_ui(*(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(i, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), 1, MPFR_RNDN);
    }

    /* diagonalize the bidiagonal form */
    for (mpfr_set_ui(k, n-1, MPFR_RNDN); mpfGeZero(k); mpfr_sub_ui(k, k, 1, MPFR_RNDN))
    {                             /* loop over singular values */
        for (mpfr_set_ui(its, 0, MPFR_RNDN); mpfr_cmp_ui(its, 30) < 0; mpfr_add_ui(its, its, 1, MPFR_RNDN))
        {                         /* loop over allowed iterations */
            flag = true;
            for (mpfr_set(l, k, MPFR_RNDN); mpfGeZero(l); mpfr_sub_ui(l, l, 1, MPFR_RNDN))
            {                     /* test for splitting */
            	mpfr_sub_ui(nm, l, 1, MPFR_RNDN);
            	mpfr_abs(tmp, rv1[mpfr_get_ui(l, MPFR_RNDN)], MPFR_RNDN);
            	mpfr_add(tmp, tmp, anorm, MPFR_RNDN);
                if (mpfr_cmp(tmp, anorm) == 0)
                {
                    flag = false;
                    break;
                }
                mpfr_abs(tmp, w[mpfr_get_ui(nm, MPFR_RNDN)], MPFR_RNDN);
                mpfr_add(tmp, tmp, anorm, MPFR_RNDN);
                if(mpfr_cmp(tmp, anorm) == 0)
                    break;
            }
            if (flag)
            {
            	mpfr_set_ui(c, 0, MPFR_RNDN);
            	mpfr_set_ui(s, 1, MPFR_RNDN);
                for (mpfr_set(i, l, MPFR_RNDN); mpfr_lessequal_p(i, k); mpfr_add_ui(i, i, 1, MPFR_RNDN))
                {
                	mpfr_mul(f, s, rv1[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
					mpfr_abs(tmp, f, MPFR_RNDN);
					mpfr_add(tmp, tmp, anorm, MPFR_RNDN);
					if(mpfr_cmp(tmp, anorm))
                    {
                    	mpfr_set(g, w[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
                    	PYTHAG(h, f, g);
                    	mpfr_set(w[mpfr_get_ui(i, MPFR_RNDN)], h, MPFR_RNDN);
                    	mpfr_pow_si(h, h, -1, MPFR_RNDN);
                    	mpfr_mul(c, g, h, MPFR_RNDN);
                    	mpfr_mul(s, f, h, MPFR_RNDN);
                    	mpfr_neg(s, s, MPFR_RNDN);
                        for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j,m) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                        {
                        	mpfr_set(y, *(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(nm, MPFR_RNDN)), MPFR_RNDN);
                        	mpfr_set(z, *(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
                        	mpfr_mul(tmp, z, s, MPFR_RNDN);
                        	mpfr_mul(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(nm, MPFR_RNDN)), y, c, MPFR_RNDN);
                        	mpfr_add(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(nm, MPFR_RNDN)), *(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(nm, MPFR_RNDN)), tmp, MPFR_RNDN);
                        	mpfr_mul(tmp, y, s, MPFR_RNDN);
                        	mpfr_mul(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), z, c, MPFR_RNDN);
                        	mpfr_sub(*(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), tmp, MPFR_RNDN);
                        }
                    }
                }
            }
            mpfr_set(z, w[mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
            if (mpfr_cmp(l, k) == 0)
            {                  /* convergence */
                if (mpfLZero(z))
                {              /* make singular value nonnegative */
            		mpfr_neg(w[mpfr_get_ui(k, MPFR_RNDN)], z, MPFR_RNDN);
                    for (mpfr_set_ui(j, 0, MPFR_RNDN); mpfr_cmp_ui(j,n) < 0; mpfr_add_ui(j, j, 1, MPFR_RNDN))
                    	mpfr_neg(*(v + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), *(v + n*mpfr_get_ui(j, MPFR_RNDN) + mpfr_get_ui(k, MPFR_RNDN)), MPFR_RNDN);
                }
                break;
            }
            if (mpfr_cmp_ui(its, (access(curLayout)->max_dsvd_iterations/1000)) >= 0)
            {
                printErr(33, "No convergence after %u! iterations", access(curLayout)->max_dsvd_iterations);
                mpfr_clears(tmp, tmp2, tmp3, tmp4, anorm, g, scale, i, its, j, jj, k, l, nm, c, f, h, s, x, y, z, NULL);
    			matrixFree(&rv1, ((dim_typ2){n,1}));
                free(rv1);
                return false;
            }

            /* shift from bottom 2 x 2 minor */
            mpfr_set(x, w[mpfr_get_ui(l, MPFR_RNDN)], MPFR_RNDN);
            mpfr_sub_ui(nm, k, 1, MPFR_RNDN);
            mpfr_set(y, w[mpfr_get_ui(nm, MPFR_RNDN)], MPFR_RNDN);
            mpfr_set(g, rv1[mpfr_get_ui(nm, MPFR_RNDN)], MPFR_RNDN);
            mpfr_set(h, rv1[mpfr_get_ui(k, MPFR_RNDN)], MPFR_RNDN);
            mpfr_mul_ui(tmp, h, 2, MPFR_RNDN);
            mpfr_mul(tmp, tmp, y, MPFR_RNDN);
            mpfr_add(tmp2, g, h, MPFR_RNDN);
            mpfr_sub(tmp3, g, h, MPFR_RNDN);
            mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);
            mpfr_add(tmp3, y, z, MPFR_RNDN);
            mpfr_sub(tmp4, y, z, MPFR_RNDN);
            mpfr_mul(tmp3, tmp3, tmp4, MPFR_RNDN);
            mpfr_add(tmp2, tmp2, tmp3, MPFR_RNDN);
			mpfr_div(f, tmp2, tmp, MPFR_RNDN);
			mpfr_set_ui(tmp, 1, MPFR_RNDN);
			PYTHAG(g, f, tmp);
			SIGN(tmp, g, f);
			mpfr_add(tmp, tmp, f, MPFR_RNDN);
			mpfr_div(tmp, y, tmp, MPFR_RNDN);
			mpfr_sub(tmp, tmp, h, MPFR_RNDN);
			mpfr_mul(tmp, tmp, h, MPFR_RNDN);
			mpfr_add(tmp2, x, z, MPFR_RNDN);
			mpfr_sub(tmp3, x, z, MPFR_RNDN); 
			mpfr_mul(tmp2, tmp2, tmp3, MPFR_RNDN);
			mpfr_add(tmp, tmp2, tmp, MPFR_RNDN);
			mpfr_div(f, tmp, x, MPFR_RNDN);

            /* next QR transformation */
            mpfr_set_ui(c, 1, MPFR_RNDN);
            mpfr_set_ui(s, 1, MPFR_RNDN);
            for (mpfr_set(j, l, MPFR_RNDN); mpfr_lessequal_p(j, nm); mpfr_add_ui(j, j, 1, MPFR_RNDN))
            {
            	mpfr_add_ui(i, j, 1, MPFR_RNDN);
            	mpfr_set(g, rv1[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
            	mpfr_set(y, w[mpfr_get_ui(i, MPFR_RNDN)], MPFR_RNDN);
            	mpfr_mul(h, s, g, MPFR_RNDN);
            	mpfr_mul(g, c, g, MPFR_RNDN);
            	PYTHAG(z, f, h);
            	mpfr_set(rv1[mpfr_get_ui(j, MPFR_RNDN)], z, MPFR_RNDN);
            	mpfr_div(c, f, z, MPFR_RNDN);
            	mpfr_div(s, h, z, MPFR_RNDN); 
            	mpfr_mul(tmp, g, s, MPFR_RNDN);
            	mpfr_mul(f, c, x, MPFR_RNDN);
            	mpfr_add(f, f, tmp, MPFR_RNDN);
            	mpfr_mul(tmp, x, s, MPFR_RNDN);
            	mpfr_mul(g, g, c, MPFR_RNDN);
            	mpfr_sub(g, g, tmp, MPFR_RNDN);
            	mpfr_mul(h, y, s, MPFR_RNDN);
            	mpfr_mul(y, y, c, MPFR_RNDN);
                for (mpfr_set_ui(jj, 0, MPFR_RNDN); mpfr_cmp_ui(jj, n) < 0; mpfr_add_ui(jj, jj, 1, MPFR_RNDN))
                {
                	mpfr_set(x, *(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), MPFR_RNDN);
                	mpfr_set(z, *(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
                	mpfr_mul(tmp, z, s, MPFR_RNDN);
                	mpfr_mul(*(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), x, c, MPFR_RNDN);
                	mpfr_add(*(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), *(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), tmp, MPFR_RNDN);
                	mpfr_mul(tmp, x, s, MPFR_RNDN);
                	mpfr_mul(*(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), z, c, MPFR_RNDN);
                	mpfr_sub(*(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(v + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), tmp, MPFR_RNDN);
                }
                PYTHAG(z, f, h);
                mpfr_set(w[mpfr_get_ui(j, MPFR_RNDN)], z, MPFR_RNDN);
                if (!mpfr_zero_p(z))
                {
                	mpfr_pow_si(z, z, -1, MPFR_RNDN);
                	mpfr_mul(c, f, z, MPFR_RNDN);
                	mpfr_mul(s, h, z, MPFR_RNDN);
                }
                mpfr_mul(tmp, s, y, MPFR_RNDN);
                mpfr_mul(f, c, g, MPFR_RNDN);
                mpfr_mul(f, f, tmp, MPFR_RNDN);
                mpfr_mul(tmp, s, g, MPFR_RNDN);
                mpfr_mul(x, c, y, MPFR_RNDN);
                mpfr_sub(x, x, tmp, MPFR_RNDN);
                for (mpfr_set_ui(jj, 0, MPFR_RNDN); mpfr_cmp_ui(jj,m) < 0; mpfr_add_ui(jj, jj, 1, MPFR_RNDN))
                {
                	mpfr_set(y, *(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), MPFR_RNDN);
                	mpfr_set(z, *(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), MPFR_RNDN);
                	mpfr_mul(tmp, z, s, MPFR_RNDN);
                	mpfr_mul(*(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), y, c, MPFR_RNDN);
                	mpfr_add(*(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), *(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(j, MPFR_RNDN)), tmp, MPFR_RNDN);
                	mpfr_mul(tmp, y, s, MPFR_RNDN);
                	mpfr_mul(*(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), z, c, MPFR_RNDN);
                	mpfr_sub(*(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), *(a + n*mpfr_get_ui(jj, MPFR_RNDN) + mpfr_get_ui(i, MPFR_RNDN)), tmp, MPFR_RNDN);
                }
            }
            mpfr_set_ui(rv1[mpfr_get_ui(l,MPFR_RNDN)], 0, MPFR_RNDN);
            mpfr_set(rv1[mpfr_get_ui(k,MPFR_RNDN)], f, MPFR_RNDN);
            mpfr_set(w[mpfr_get_ui(k,MPFR_RNDN)], x, MPFR_RNDN);
        }
    }
    mpfr_clears(tmp, tmp2, tmp3, tmp4, anorm, g, scale, i, its, j, jj, k, l, nm, c, f, h, s, x, y, z, NULL);
    matrixFree(&rv1, ((dim_typ2){n,1}));
    return true;
}

 dim_typ  absor(const register dim_typ dim[static 2], mpfr_t *restrict p, mpfr_t *restrict pp, mpfr_t ** R)
{
	mpfr_t tmp, tmp2;
	dim_typ i, j;
	dim_typ retval, k;
	mpfr_t *a;
	dim_typ2 dim2 =
	{
		dim[ABSOR_RDIM],
		dim[ABSOR_RDIM]
	};
	
	mpfr_init(tmp);
	if(!matrixAlloc(&a, dim2)) 
		return ABSOR_ERROR; 
	for(i=0; i<dim[ABSOR_DIM]; ++i)
		for(j=0; j<dim[ABSOR_RDIM]; ++j)
			for(k=0; k<dim[ABSOR_RDIM]; ++k)
			{
				mpfr_mul(tmp, p[i*dim[ABSOR_DIM]], pp[i*dim[ABSOR_DIM] + k], MPFR_RNDN);
				mpfr_add(a[j*dim[ABSOR_RDIM] + k], a[j*dim[ABSOR_RDIM] + k], tmp, MPFR_RNDN);
			}
				
	mpfr_t *Utransp;
	mpfr_t *svd_vec;
	mpfr_t *U;
	mpfr_t *S;
	mpfr_t *V;
	
	if(!matrixAlloc(&U, dim2))
		return ABSOR_ERROR;
	if(!matrixAlloc(&V, dim2))
	{
		matrixFree(&U, dim2);
		return ABSOR_ERROR;
	}
	if(!matrixAlloc(&Utransp, dim2))
	{
		matrixFree(&U, dim2);
		matrixFree(&V, dim2);
		return ABSOR_ERROR;
	}
	if(!matrixAlloc(R, dim2))
	{
		matrixFree(&U, dim2);
		matrixFree(&V, dim2);
		matrixFree(&Utransp, dim2);
		return ABSOR_ERROR;
	}	
	(void) matrixAlloc(&svd_vec, (dim_typ2){dim[ABSOR_RDIM],1});
	if(svd_vec == NULL)
	{
		matrixFree(&U, dim2);
		matrixFree(&V, dim2);
		matrixFree(&Utransp, dim2);
		matrixFree(R, dim2);
		return ABSOR_ERROR;
	}
	
	mpfr_init_set_d(tmp2, SVD_EPSILON, MPFR_RNDN);
	mpfr_init_set_d(tmp2, SVD_TOLERANCE, MPFR_RNDN);
	
	// const register dim_typ final_dim = dim[ABSOR_RDIM]*dim[ABSOR_RDIM];
	if(!(retval = svd(dim2, SVD_WITHU, SVD_WITHV, tmp, tmp2, a, svd_vec, U, V)))
	{

		transpose(U, Utransp, dim2);
		matrixFree(&U, dim2);
		_matrixMultiplication(&V, &Utransp, R, (dim_typ3){dim[ABSOR_RDIM],dim[ABSOR_RDIM],dim[ABSOR_RDIM]}, (dim_typ3){dim[ABSOR_RDIM],dim[ABSOR_RDIM],dim[ABSOR_RDIM]});
		matrixFree(&V, dim2);
		matrixFree(&Utransp, dim2);
	}
	mpfr_clears(tmp, tmp2, NULL);
	matrixFree(&a, dim2);
	return retval;
		
}

// RANK Calculator of a nXm Matrix using
// Jacobi Singular Value Decomposition METHOD (SVD)
__MATHSUITE dim_typ  rank(mpfr_t *restrict matrix, const register dim_typ dim[static 2])
{
    mpfr_t * S;
    mpfr_t * V;

    const register dim_typ maxv = dim[dim[ROWS] >= dim[COLUMNS]];

	if(!matrixAlloc(&S, (dim_typ2){maxv,1}))
		return USHRT_MAX;

    if(!matrixAlloc(&V, (dim_typ2){maxv, maxv}))
    {
        matrixFree(&S, ((dim_typ2){maxv,1}));
        return USHRT_MAX;
    }

    dsvd(matrix, dim, S, V);
    
    matrixFree(&V, ((dim_typ2){maxv, maxv}));
    register dim_typ rnk = maxv;
    mpfr_t tmp;
    
    mpfr_init(tmp);

    for(dim_typ i=0; i<maxv; ++i)
    {
    	mpfr_abs(tmp, S[i], MPFR_RNDN);
        if(mpfr_cmp_d(tmp, PREC) < 0) -- rnk;
   	}
        
    mpfr_clear(tmp);
	matrixFree(&S, ((dim_typ2){maxv,1}));
    return rnk;
}
/*
__MATHSUITE dim_typ  rank(mpfr_t *restrict matrix, const register dim_typ dim[static 2])
{
    mpfr_t tmp, tmp2;
	mpfr_t *svd_vec;
	mpfr_t *Lm;
	mpfr_t *Rm;
	
	const register dim_typ maxv = dim[dim[ROWS] >= dim[COLUMNS]];
	const register dim_typ secure_dim = MAX(dim[ROWS], dim[COLUMNS]);
    
    if(!matrixAlloc(&Lm, (dim_typ2){secure_dim, secure_dim}))
        return USHRT_MAX;
    
    if(!matrixAlloc(&Rm, (dim_typ2){secure_dim, secure_dim}))
    {
    	matrixFree(&Lm, ((dim_typ2){secure_dim, secure_dim}));
        return USHRT_MAX;
    }
    
    if(!matrixAlloc(&svd_vec, (dim_typ2){dim[ROWS],1}))
    {
    	matrixFree(&Lm, ((dim_typ2){secure_dim, secure_dim}));
    	matrixFree(&Rm, ((dim_typ2){secure_dim, secure_dim}));
        return USHRT_MAX;
    }
    
    mpfr_init_set_d(tmp, SVD_EPSILON, MPFR_RNDN);
	mpfr_init_set_d(tmp2, SVD_TOLERANCE, MPFR_RNDN);
	
	svd(dim, SVD_WITHU, SVD_WITHV, tmp, tmp2, matrix, svd_vec, Lm, Rm);
    // dsvd(matrix, dim, S, V);
    
    matrixFree(&Lm, ((dim_typ2){secure_dim, secure_dim}));
	matrixFree(&Rm, ((dim_typ2){secure_dim, secure_dim}));
    register dim_typ rnk = maxv;

    for(dim_typ i=0; i<maxv; ++i)
    {
    	mpfr_abs(tmp, svd_vec[i], MPFR_RNDN);
        if(mpfr_cmp_d(tmp, PREC) < 0) -- rnk;
   	}
        
    mpfr_clears(tmp, tmp2, NULL);
    matrixFree(&svd_vec, ((dim_typ2){dim[ROWS],1}));
    return rnk;
}
*/


__MATHSUITE inline void  _flushLogBuf(logObj * const which_log)
{
    strcpy(which_log->buffer, NULL_CHAR);
    return;
}

__MATHSUITE __WINCALL void  _editLog(const char path[static MAX_PATH_LENGTH])
{
#ifdef WINOS
    char str[PATHCONTAINS_STRING] = NULL_CHAR;
    sprintf(str, "notepad.txt %s", path);
    (void) system(str);
#else
    msprintf(COLOR_SYSTEM, "You entered Sequential Log Editing:\n%s.\nPress CTRL + C to exit.\n\n", path);

    FILE *fp = NULL;

    if((fp = checkForFHErrors(path, "a")) == NULL)
        return;

    fputs(NULL_CHAR, fp);

    char str[MAX_BUFSIZ] = NULL_CHAR;
    for( ; ; )
    {
        signal(SIGINT, sigexit);
        (void) gets(str);
        if(access(sigresult))
        {
            fclose(fp);
            msprintf(COLOR_SYSTEM, "\nYou exited Sequential Log Editing:\n%s.\n\n", path);
            access(sigresult) = false;
            return;
        }
        fputs(str, fp);
        /// writelog functions
    }
    fclose(fp); /// in order to provide unaligned strongly dangerous goto.
#endif
    return;
}

__MATHSUITE void  logCheck(logObj * const which_log, const char *format, const char path[static MAX_PATH_LENGTH])
{
    if(strlen(which_log->buffer)+strlen(format) >= which_log->buflen)
    {
        FILE *fp;

        if((fp = checkForFHErrors(path, "a")) == NULL)
            return;

        logWrite(fp, which_log);
        fclose(fp);
        _flushLogBuf(which_log);
        
    }
    
	strcat(which_log->buffer, format);
    return;
}

__MATHSUITE inline void  prependTimeToString(char *string)
{
    char str_time[MAX_BUFSIZ+INFO_STRING];
    time_t ora;
    ora = time(NULL);
    sprintf(str_time, "[%s] - %s", asctime(localtime(&ora)), string);
    strcpy(string, str_time);
    return;
}

__MATHSUITE void  printErr(const int err, const char *format, ...)
{
	va_list ap;
	va_start(ap, format);
	char str[MAX_BUFSIZ+INFO_STRING];         
	vsprintf(str, format, ap);
    errno = err;
    
    #ifndef __DISABLE_SERVER
    if(access(lastSessionSocket) == INVALID_SOCKET)
    {
    #endif
    	CLEARBUFFER();
    	PRINTL();
    	SetColor(COLOR_ERROR);
    	perror("\nERROR");
    	printf("%s.", str);
    	getch();
    	SetDefaultColor();
    #ifndef __DISABLE_SERVER
    }
    else if(access(server_mode))
	{
		if(strlen(access(serverBuffer))+strlen(str) >= access(serverBufferLength))
		{
			access(serverBufferLength) += access(serverBufferLength);
			access(serverBuffer) = realloc(access(serverBuffer), access(serverBufferLength));
			errMem(access(serverBuffer), VSPACE);
		}
    	strcat(access(serverBuffer), str); // send(access(lastSessionSocket), str, strlen(str), 0);
	}
	#endif

    // msyprintf(COLOR_ERROR, format, ap);
    
    PRINTL();
    va_end(ap);
    return;
}

__MATHSUITE void  msyprintf(const sel_typ col, const char *format, ...)
{
    va_list ap;                                                                  
    va_start(ap,format);                                                              
    char str[MAX_BUFSIZ+INFO_STRING];                                            
	                                                          	 
    vsprintf(str, format, ap);  
    
    #ifndef __DISABLE_SERVER
	if(access(lastSessionSocket) == INVALID_SOCKET)
	{
	#endif
		SetColor(col);   
    	printf(str);
    	SetDefaultColor();   
    #ifndef __DISABLE_SERVER
	}
    else if(access(server_mode))
	{
		if(strlen(access(serverBuffer))+strlen(str) >= access(serverBufferLength))
		{
			access(serverBufferLength) += access(serverBufferLength);
			access(serverBuffer) = realloc(access(serverBuffer), access(serverBufferLength));
			errMem(access(serverBuffer), VSPACE);
		}
    	strcat(access(serverBuffer), str); // send(access(lastSessionSocket), str, strlen(str), 0);
	}
	#endif                                          
	                                                 				 
    if(access(sysLog) != INIT_LOGOBJ) 											 
	{																			 
		prependTimeToString(str);                                                
		logCheck(access(sysLog), str, access(sysLogPath));               	 	 
	}																			 
                                                            
    va_end(ap);
    return;
}

__MATHSUITE void  fprintf2(FILE *fp, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);

    char str[MAX_BUFSIZ];
    vsprintf(str, format, ap);

    const bool cond = fp == stdout;
    
	if(cond)
	{	
		#ifndef __DISABLE_SERVER
		if(access(lastSessionSocket) == INVALID_SOCKET)
		{
		#endif
			SetColor(COLOR_USER);
    		printf(str);
    		SetDefaultColor();
    	#ifndef __DISABLE_SERVER
    	}
	    else if(access(server_mode))
		{
			if(strlen(access(serverBuffer))+strlen(str) >= access(serverBufferLength))
			{
				access(serverBufferLength) += access(serverBufferLength);
				access(serverBuffer) = realloc(access(serverBuffer), access(serverBufferLength));
				errMem(access(serverBuffer), VSPACE);
			}
	    	strcat(access(serverBuffer), str); // send(access(lastSessionSocket), str, strlen(str), 0);
		}
		#endif
	}
	else
    	fprintf(fp, str);
        
    if(getItemsListNo(LOGS) != STARTING_LOGSNO)
    {
        nodelist * const cur_log = listNo(access(lists)[LOGS].cur_item, LOGS);
        logCheck(((logObj *)(cur_log->data)), str, cur_log->path);
    }

    va_end(ap);
    return;
}


__MATHSUITE void  msprintf(const sel_typ col, const char *format, ...)
{
    va_list ap;
    va_start(ap, format);

    char str[MAX_BUFSIZ];

    vsprintf(str, format, ap);
    
    #ifndef __DISABLE_SERVER
    if(access(lastSessionSocket) == INVALID_SOCKET)
    {
    #endif
		SetColor(col);
    	printf(str);
    	SetDefaultColor();
    #ifndef __DISABLE_SERVER
    }
    else if(access(server_mode))
	{
		if(strlen(access(serverBuffer))+strlen(str) >= access(serverBufferLength))
		{
			access(serverBufferLength) += access(serverBufferLength);
			access(serverBuffer) = realloc(access(serverBuffer), access(serverBufferLength));
			errMem(access(serverBuffer), VSPACE);
		}
    	strcat(access(serverBuffer), str); // send(access(lastSessionSocket), str, strlen(str), 0);
	}
	#endif
    	
    /*
    else
    {
    	SetDefaultColor();
    	printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "Some Socket FD went corrupted");
    	va_end(ap);
    	return;
	}
	*/
    
    if(getItemsListNo(LOGS) != STARTING_LOGSNO)
    {
        nodelist * const cur_log = listNo(access(lists)[LOGS].cur_item, LOGS);
        logCheck(((logObj *)(cur_log->data)), str, cur_log->path);
    }
    
    va_end(ap);
    return;
}

__MATHSUITE bool  scanf2(sel_typ count, const char *format, ...)
{
	va_list ap;
    va_start(ap, format);

    int scanner = 0;

    access(exitHandle) = INVALID_EXITHANDLE;
    // signal(SIGINT, equalSpecialMatrix);
    scanner = vscanf(format, ap);

    va_end(ap);

    if(!scanner)
        access(exitHandle) = EXITHANDLE_EXIT;

    return scanner == count;
}

__MATHSUITE  void  _printMatrix(FILE *fp, mpfr_t **matrix, const register dim_typ dim[static 2], const char matname[])
{

    const bool assert = fp == stdout;

    if(assert)
        PRINTN();

	char buf[MAX_BUFSIZ];
    dim_typ  i, j;

	if(fp)
	{
		for(i=0; i<dim[ROWS]; ++i)
		{
			strcpy(buf, NULL_CHAR);
	    	if(assert && isSett(BOOLS_PRINTROWSLABELS))
	        	fprintf2(fp, "R%hu: ", i+1);
	        for(j=0; j<dim[COLUMNS]; ++j)
				mpfr_sprintf(buf, "%s"OUTPUT_CONVERSION_FORMAT"; ", buf, *((*matrix) + dim[COLUMNS]*i + j));
	        mpfr_sprintf(buf, "%s\n", buf);
	        fprintf2(fp, buf);
	        // fputc('\n', fp);
		}
	}
	#ifndef __DISABLE_DATABASE
		else
		{
			for(i=0; i<dim[ROWS]; ++i)
			{
				strcpy(buf, NULL_CHAR);
		        for(j=0; j<dim[COLUMNS]; ++j)
					mpfr_sprintf(buf, "%s"OUTPUT_CONVERSION_FORMAT",", buf, *((*matrix) + dim[COLUMNS]*i + j));
				buf[strlen(buf)-1] = '\0';
				writeMat(matname, buf);       
			}
		}
	#endif
    
    if(assert)
    {
        PRINT2N();
        PRINTN();
        
        matrixObj * const currentLmpMatrix = accessCurrentLmpMatrix;

		if(&(currentLmpMatrix->matrix) != matrix)
		{
			if(currentLmpMatrix && currentLmpMatrix->matrix)
        		matrixFree(&(currentLmpMatrix->matrix), currentLmpMatrix->dim);

	        if((!currentLmpMatrix->matrix) && equalMatrix(&(currentLmpMatrix->matrix), (*matrix), dim, EQUALMATRIX_REALLOC))
	        {
	            currentLmpMatrix->dim[ROWS] = dim[ROWS];
	            currentLmpMatrix->dim[COLUMNS] = dim[COLUMNS];
	        }
		}
		
        
    }

    return;
}

__MATHSUITE  void  printItypMatrix(ityp *matrix, const register dim_typ dim[static 2])
{

    PRINTN();
    dim_typ  i, j;

    for(i=0; i<dim[ROWS]; ++i)
	{
    	if(isSett(BOOLS_PRINTROWSLABELS))
        	msprintf(COLOR_USER, "R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
			msprintf(COLOR_USER, "%.*f", access(curLayout)->precision, *(matrix + dim[COLUMNS]*i + j));
            msprintf(COLOR_USER, "; ");
            if(j >= dim[COLUMNS]-1)
            	putchar('\n');

        }
	}

    return;
}

__MATHSUITE  void  printBoolMatrix(bool *matrix, const register dim_typ dim[static 2])
{

    PRINTN();
    dim_typ  i, j;

    for(i=0; i<dim[ROWS]; ++i)
	{
    	if(isSett(BOOLS_PRINTROWSLABELS))
        	msprintf(COLOR_USER, "R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
			msprintf(COLOR_USER, "%s", *(matrix + dim[COLUMNS]*i + j) ? "true":"false");
            msprintf(COLOR_USER, "; ");
            if(j >= dim[COLUMNS]-1)
            	putchar('\n');

        }
	}

    return;
}

__MATHSUITE  void  printUShortMatrix(dim_typ *matrix, const register dim_typ dim[static 2])
{

    PRINTN();
    dim_typ  i, j;

    for(i=0; i<dim[ROWS]; ++i)
	{
    	if(isSett(BOOLS_PRINTROWSLABELS))
        	msprintf(COLOR_USER, "R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
			msprintf(COLOR_USER, "%hu", *(matrix + dim[COLUMNS]*i + j));
            msprintf(COLOR_USER, "; ");
            if(j >= dim[COLUMNS]-1)
            	putchar('\n');

        }
	}

    return;
}

__MATHSUITE  void  printShortMatrix(short  *matrix, const register dim_typ dim[static 2])
{

    PRINTN();
    dim_typ  i, j;

    for(i=0; i<dim[ROWS]; ++i)
	{
    	if(isSett(BOOLS_PRINTROWSLABELS))
        	msprintf(COLOR_USER, "R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
			msprintf(COLOR_USER, "%h", *(matrix + dim[COLUMNS]*i + j));
            msprintf(COLOR_USER, "; ");
            if(j >= dim[COLUMNS]-1)
            	putchar('\n');

        }
	}

    return;
}

__MATHSUITE  void  printIntMatrix(int  *matrix, const register dim_typ dim[static 2])
{

    PRINTN();
    dim_typ  i, j;

    for(i=0; i<dim[ROWS]; ++i)
	{
    	if(isSett(BOOLS_PRINTROWSLABELS))
        	msprintf(COLOR_USER, "R%hu: ", i+1);
        for(j=0; j<dim[COLUMNS]; ++j)
        {   // remember to put comma and whitespace
            // in order to right-format matrix file-parsing system
			msprintf(COLOR_USER, "%d", *(matrix + dim[COLUMNS]*i + j));
            msprintf(COLOR_USER, "; ");
            if(j >= dim[COLUMNS]-1)
            	putchar('\n');

        }
	}

    return;
}

__MATHSUITE bool   _parse(char expr[], mpfr_t *res, exprValList * vlist)
{
	
    EXPRERRTYPE err;
    exprObj * exp = INIT_OBJLIST;

    err = exprCreate(&exp, access(func_list), vlist ? vlist : access(exprVars)->var_list, access(const_list), NULL, 0);

    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_SYSTEM, "Expr Creation Error.\n");
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
        msyprintf(COLOR_SYSTEM, "Parse Error (%d,%d).\n", start, end);
        exprFree(exp);
        return false;
    }

    mpfr_t val;
    err = exprEval(exp, val);

    if(err != EXPR_ERROR_NOERROR)
    {
        msyprintf(COLOR_ERROR, "Eval Error:\n%s.\n", suite_c.exprErrorString[err == EXPR_ERROR_UNKNOWN ? 0 : err+1  ]);
        exprFree(exp);
        return false;
    }

	if(res)
    	mpfr_init_set(*res, val, MPFR_RNDN);
    	
    mpfr_clear(val);
    exprFree(exp);
    return true;
}

#define MINMAX_BUFFER_LEN MAX_STRING

__MATHSUITE bool   extractMat(dim_typ which_mat)
{
    matrixObj * const tmp = malloc(sizeof(matrixObj));
    errMem(tmp, false);
    tmp->matrix = NULL;
        
    #ifndef __DISABLE_DATABASE
    if(!access(curLayout)->database.con)
    	dbEstablishConnection(0, NULL); // (char * [4]){NULL, NULL, NULL, NULL});
    	
    register bool fail = false;

	if(access(curLayout)->database.con)
	{
		char mat[MAX_BUFSIZ] = NULL_CHAR;
		char matname [MAX_PATH_LENGTH];
		strcpy(matname, listNo(which_mat, MATRICES)->path);
		(void) strtok(matname, EXTENSION_DOT);
		if((!readMat(matname, mat)) || !matrixToken(mat, &(tmp->matrix), &tmp->dim[ROWS], &tmp->dim[COLUMNS]))   
			fail = true;
	}
	
	if(fail || !access(curLayout)->database.con)
	{
	#endif
		dim_typ i;
		struct timeval tvBegin;
	    char str[MINMAX_BUFFER_LEN] = NULL_CHAR;
	    dim_typ analog_rows, analog_columns = 1;
	    char *ij_element = NULL;
		FILE *fp = NULL;
		
		if((fp = checkForFHErrors(listNo(which_mat, MATRICES)->path, "r")) == NULL)
       		return false;
		
	    fflush(fp);
	    gettimeofday(&tvBegin, NULL);
	
	    for(tmp->dim[ROWS]=tmp->dim[COLUMNS]=INIT_DIM; fgets(str, sizeof(str), fp) != NULL; ++ tmp->dim[ROWS])  // fscanf(fp, "%[^\n]s", str)) // fgets(str, sizeof(str), fp) != NULL)
	    {
	        if(!(tmp->dim[ROWS]) && !matrixAlloc(&(tmp->matrix), (dim_typ2){1, 1}))
	            return false;
	        else
	        {
	            analog_rows = (tmp->dim[ROWS])+1;
	
	            tmp->matrix = realloc(tmp->matrix, sizeof(mpfr_t)*analog_rows*analog_columns);
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
	                tmp->matrix = realloc(tmp->matrix, sizeof(mpfr_t)*analog_rows*analog_columns);
	                errMem(tmp->matrix, false);
	        	}
	
	            if(!parse(ij_element, (tmp->matrix) + tmp->dim[COLUMNS]*tmp->dim[ROWS] + i))
	                continue;
	            else
	                mpfr_set_d(*((tmp->matrix) + tmp->dim[COLUMNS]*tmp->dim[ROWS] + i), strtod(ij_element, NULL), MPFR_RNDN);
	
	            if(!(tmp->dim[ROWS]))
	               analog_columns = ++(tmp->dim[COLUMNS]) +1;
	        }
	    }
	
	    fclose(fp);
	    
	    if(isSett(BOOLS_SHOWDIFFTIME))
	    {
	        PRINTL();
	        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
	        PRINTL();
	    }
	#ifndef __DISABLE_DATABASE
	}
	#endif

    // REDIRECTING WHICH_MAT LISTBOX ITEM
    listNo(which_mat, MATRICES)->data = tmp;
    return true;
}


__MATHSUITE bool   matrixToken(const char string[], mpfr_t **matrix, dim_typ *righe, dim_typ *colonne)
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

            (*matrix) = realloc((*matrix), sizeof(mpfr_t)*analog_rows*analog_columns);
            errMem((*matrix), false);
        }

        for(i=0, token=strtok(buf,","); token != NULL; ++ i, token = strtok (token + strlen (token) + 1, ","))
        {

            if(!((*righe)))
            {
                (*matrix) = realloc((*matrix), sizeof(mpfr_t)*analog_rows*analog_columns);
                errMem((*matrix), false);
            }

            char token2[strlen(token)+1];
            sprintf(token2, "%s;", token);

            if(!parse(token2, (*matrix) + (*colonne)*(*righe) + i))
                continue;
            else
            	mpfr_set_d(*((*matrix) + (*colonne)*(*righe) + i), strtod(token2, NULL), MPFR_RNDN); 

            if(!((*righe)))
                analog_columns = ++ (*colonne) +1;

        }
    }

    return true;
}

__MATHSUITE inline void  delAll()
{
	for(uint64_t i=0; i<access(PRMSystem).currentIndex; ++i)
		if(accessContainer(i).pnt)
		{
			free(accessContainer(i).pnt);
			accessContainer(i).size = 0;
			accessContainer(i)._volatile = false;
		}
	access(PRMSystem).currentIndex = 0;
	return;
}

__MATHSUITE void  newFloat(const uint64_t siz, mpfr_t val, const bool _volatile, ityp * tmp)
{
	if(access(PRMSystem).currentIndex+1 == MAX_PRMCONTAINERSIZE)
    	if(isSett(BOOLS_ARRAYSAUTORESET))
    		delAll();
    	else
    	{
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    		return;
    	}
    else
    {
    	ityp * tmpval = tmp ? tmp : calloc(siz, sizeof(ityp));
    		
    	if(!tmpval)
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    	else
    	{
    		mpfr_set_ui(val, (access(PRMSystem).currentIndex) ++, MPFR_RNDN);
    		_accessContainer(val).pnt = tmpval;
    		_accessContainer(val).size = siz;
    		_accessContainer(val)._volatile = _volatile;
    	}
    }
	return;
}

__MATHSUITE void  newBool(const uint64_t siz, mpfr_t val, const bool _volatile, bool * tmp)
{
	if(access(PRMSystem).currentIndex+1 == MAX_PRMCONTAINERSIZE)
    	if(isSett(BOOLS_ARRAYSAUTORESET))
    		delAll();
    	else
    	{
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    		return;
    	}
    else
    {
    	ityp * tmpval = tmp ? tmp : calloc(siz, sizeof(bool));
    		
    	if(!tmpval)
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    	else
    	{
    		mpfr_set_ui(val, (access(PRMSystem).currentIndex) ++, MPFR_RNDN);
    		_accessContainer(val).pnt = tmpval;
    		_accessContainer(val).size = siz;
    		_accessContainer(val)._volatile = _volatile;
    	}
    }
	return;
}

__MATHSUITE void  newUShort(const uint64_t siz, mpfr_t val, const bool _volatile, dim_typ * tmp)
{
	if(access(PRMSystem).currentIndex+1 == MAX_PRMCONTAINERSIZE)
    	if(isSett(BOOLS_ARRAYSAUTORESET))
    		delAll();
    	else
    	{
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    		return;
    	}
    else
    {
    	ityp * tmpval = tmp ? tmp : calloc(siz, sizeof(dim_typ));
    		
    	if(!tmpval)
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    	else
    	{
    		mpfr_set_ui(val, (access(PRMSystem).currentIndex) ++, MPFR_RNDN);
    		_accessContainer(val).pnt = tmpval;
    		_accessContainer(val).size = siz;
    		_accessContainer(val)._volatile = _volatile;
    	}
    }
	return;
}

__MATHSUITE void  newShort(const uint64_t siz, mpfr_t val, const bool _volatile, short * tmp)
{
	if(access(PRMSystem).currentIndex+1 == MAX_PRMCONTAINERSIZE)
		if(isSett(BOOLS_ARRAYSAUTORESET))
    		delAll();
    	else
    	{
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    		return;
    	}
    else
    {
    	ityp * tmpval = tmp ? tmp : calloc(siz, sizeof(short));
    		
    	if(!tmpval)
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    	else
    	{
    		mpfr_set_ui(val, (access(PRMSystem).currentIndex) ++, MPFR_RNDN);
    		_accessContainer(val).pnt = tmpval;
    		_accessContainer(val).size = siz;
    		_accessContainer(val)._volatile = _volatile;
    	}
    }
	return;
}

__MATHSUITE void  newInt(const uint64_t siz, mpfr_t val, const bool _volatile, int * tmp)
{
	if(access(PRMSystem).currentIndex+1 == MAX_PRMCONTAINERSIZE)
		if(isSett(BOOLS_ARRAYSAUTORESET))
    		delAll();
    	else
    	{
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    		return;
    	}
    else
    {
    	ityp * tmpval = tmp ? tmp : calloc(siz, sizeof(int));
    		
    	if(!tmpval)
    		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
    	else
    	{
    		mpfr_set_ui(val, (access(PRMSystem).currentIndex) ++, MPFR_RNDN);
    		_accessContainer(val).pnt = tmpval;
    		_accessContainer(val).size = siz;
    		_accessContainer(val)._volatile = _volatile;
    	}
    }
	return;
}

__MATHSUITE void  del(const uint64_t i)
{
	free(accessContainer(i).pnt);
	accessContainer(i).size = 0;
	accessContainer(i)._volatile = false;
	if(isSett(BOOLS_ADJUSTARRAYSINDEX))
		for(uint64_t j=i; accessContainer(j).pnt == NULL && j<access(PRMSystem).currentIndex; ++j)
			if(j == access(PRMSystem).currentIndex-1)
				access(PRMSystem).currentIndex = i;
	return;
}

__MATHSUITE inline void  aCheckFloat(const uint64_t siz, mpfr_t val, ityp * tmp)
{
	if(tmp)
		newFloat(siz, val,false, tmp);
	else
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  aCheckBool(const uint64_t siz, mpfr_t val, bool * tmp)
{
	if(tmp)
		newBool(siz, val, false, tmp);
	else
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  aCheckUShort(const uint64_t siz, mpfr_t val, dim_typ * tmp)
{
	if(tmp)
		newUShort(siz, val, false, tmp);
	else
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  aCheckShort(const uint64_t siz, mpfr_t val, short * tmp)
{
	if(tmp)
		newShort(siz, val, false, tmp);
	else
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  aCheckInt(const uint64_t siz, mpfr_t val, int * tmp)
{
	if(tmp)
		newInt(siz, val, false, tmp);
	else
		mpfr_set_ui(val, MAX_PRMCONTAINERSIZE, MPFR_RNDN);
	return;
}

__MATHSUITE inline void  truncateVector(const register dim_typ dim, mpfr_t vec[static dim], ityp outputVec[static dim])
{
	for(dim_typ i=0; i<dim; ++i) 
		outputVec[i] = mpfr_get_d(vec[i], MPFR_RNDN);
	return;
}

__MATHSUITE inline ityp *  accessFloatContainerPointer(const mpfr_t a)
{
	return _accessContainer(a).pnt;
}

__MATHSUITE inline bool *  accessBoolContainerPointer(const mpfr_t a)
{ 
	return _accessContainer(a).pnt;
}

__MATHSUITE inline dim_typ *  accessUShortContainerPointer(const mpfr_t a)
{
	return _accessContainer(a).pnt;
}

__MATHSUITE inline short *  accessShortContainerPointer(const mpfr_t a)
{ 
	return _accessContainer(a).pnt;
}

__MATHSUITE inline int *  accessIntContainerPointer(const mpfr_t a)
{
	return _accessContainer(a).pnt;
}

static inline int stringcmp (const void *a, const void *b)
{
	return strcmp((const char*)a, (const char*)b);
}

__MSSHELL_WRAPPER_ __MATHSUITE sprog * const  searchProgram(char cmdname[static SIGN_STRING])
{
	static map_t cmdmap = NULL;
	
	if(cmdmap && cmdname[0] == SEARCHPROGRAM_FREECHAR)
	{ 
		dim_typ i;

		for(i=0; i<MAX_PROGRAMMI; ++i)
			hashmap_remove(cmdmap, main_menu[i].cmdname);
		
		#ifndef __DISABLE_MULTIUSER		
			for(i=0; i<MAX_MULTIUSER_PROGS; ++i)
				hashmap_remove(cmdmap, multiuser_prog[i].cmdname);
		#endif
					
	    for(i=0; i<MAX_ADVCALC_PROGS; ++i)
	    	hashmap_remove(cmdmap, adv_calc[i].cmdname);
	    	
	    #ifndef _DISABLE_CRYPTOGRAPHICHASH
		    for(i=0; i<MAX_CRYPTOGRAPHICHASH_PROGS; ++i)
		    	hashmap_remove(cmdmap, cryptographic_hash_prog[i].cmdname);
		#endif
		
		#ifndef _DISABLE_SPECIALMATRICES
		    for(i=0; i<MAX_SPECMAT_PROGS; ++i)
		    	hashmap_remove(cmdmap, specmat_prog[i].cmdname);
		#endif
	
	    for(i=0; i<MAX_ALGEBRA_OPERATIONS; ++i)
	    	hashmap_remove(cmdmap, alg_operations[i].cmdname);
	
		#ifndef __DISABLE_SYSTEM
		    for(i=0; i<MAX_SYSTEM_PROGS; ++i)
		    	hashmap_remove(cmdmap, system_prog[i].cmdname);
		#endif
	
		#ifndef __DISABLE_VARLISTMANAGER
		    for(i=0; i<MAX_ENVSMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, envs_manager[i].cmdname);
		#endif
	
		#ifndef __DISABLE_MATMANAGER
		    for(i=0; i<MAX_MATMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, mat_manager[i].cmdname);
		#endif
	            
	    #ifndef __DISABLE_LOGSMANAGER
		    for(i=0; i<MAX_LOGSMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, logs_manager[i].cmdname);
		#endif
	
		#ifndef __DISABLE_SYSLOGMANAGER
		    for(i=0; i<MAX_SYSLOGMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, syslog_manager[i].cmdname);
		#endif
	
		#ifndef __DISABLE_LAYOUTSMANAGER
		    for(i=0; i<MAX_LAYOUTSMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, layouts_manager[i].cmdname);
	    #endif
	
	    for(i=0; i<MAX_SETTINGS; ++i)
	    	hashmap_remove(cmdmap, change_settings[i].cmdname);
		
		#ifndef __DISABLE_COLSMANAGER
		    for(i=0; i<MAX_COLSMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, cols_manager[i].cmdname);
	    #endif
	
		#ifndef __DISABLE_LFSMANAGER
		    for(i=0; i<MAX_LFSMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, lfs_manager[i].cmdname);
		#endif
	
		#ifndef __DISABLE_MSSMANAGER
		    for(i=0; i<MAX_MSSMANAGER_PROGS; ++i)
		    	hashmap_remove(cmdmap, mss_manager[i].cmdname);
	    #endif
	    
	    return NULL;
	}
	
	if(!cmdmap)
    {
    	dim_typ i;
		cmdmap = hashmap_new();

		for(i=0; i<MAX_PROGRAMMI; ++i)
			hashmap_put(cmdmap, main_menu[i].cmdname, &main_menu[i]);
		
		#ifndef __DISABLE_MULTIUSER		
			for(i=0; i<MAX_MULTIUSER_PROGS; ++i)
				hashmap_put(cmdmap, multiuser_prog[i].cmdname, &multiuser_prog[i]);
		#endif
					
	    for(i=0; i<MAX_ADVCALC_PROGS; ++i)
	    	hashmap_put(cmdmap, adv_calc[i].cmdname, &adv_calc[i]);
	    	
	    #ifndef _DISABLE_CRYPTOGRAPHICHASH
		    for(i=0; i<MAX_CRYPTOGRAPHICHASH_PROGS; ++i)
		    	hashmap_put(cmdmap, cryptographic_hash_prog[i].cmdname, &cryptographic_hash_prog[i]);
		#endif
		
		#ifndef _DISABLE_SPECIALMATRICES
		    for(i=0; i<MAX_SPECMAT_PROGS; ++i)
		    	hashmap_put(cmdmap, specmat_prog[i].cmdname, &specmat_prog[i]);
		#endif
	
	    for(i=0; i<MAX_ALGEBRA_OPERATIONS; ++i)
	    	hashmap_put(cmdmap, alg_operations[i].cmdname, &alg_operations[i]);
	
		#ifndef __DISABLE_SYSTEM
		    for(i=0; i<MAX_SYSTEM_PROGS; ++i)
		    	hashmap_put(cmdmap, system_prog[i].cmdname, &system_prog[i]);
		#endif
	
		#ifndef __DISABLE_VARLISTMANAGER
		    for(i=0; i<MAX_ENVSMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, envs_manager[i].cmdname, &envs_manager[i]);
		#endif
	
		#ifndef __DISABLE_MATMANAGER
		    for(i=0; i<MAX_MATMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, mat_manager[i].cmdname, &mat_manager[i]);
		#endif
	            
	    #ifndef __DISABLE_LOGSMANAGER
		    for(i=0; i<MAX_LOGSMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, logs_manager[i].cmdname, &logs_manager[i]);
		#endif
	
		#ifndef __DISABLE_SYSLOGMANAGER
		    for(i=0; i<MAX_SYSLOGMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, syslog_manager[i].cmdname, &syslog_manager[i]);
		#endif
	
		#ifndef __DISABLE_LAYOUTSMANAGER
		    for(i=0; i<MAX_LAYOUTSMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, layouts_manager[i].cmdname, &layouts_manager[i]);
	    #endif
	
	    for(i=0; i<MAX_SETTINGS; ++i)
	    	hashmap_put(cmdmap, change_settings[i].cmdname, &change_settings[i]);
		
		#ifndef __DISABLE_COLSMANAGER
		    for(i=0; i<MAX_COLSMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, cols_manager[i].cmdname, &cols_manager[i]);
	    #endif
	
		#ifndef __DISABLE_LFSMANAGER
		    for(i=0; i<MAX_LFSMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, lfs_manager[i].cmdname, &lfs_manager[i]);
		#endif
	
		#ifndef __DISABLE_MSSMANAGER
		    for(i=0; i<MAX_MSSMANAGER_PROGS; ++i)
		    	hashmap_put(cmdmap, mss_manager[i].cmdname, &mss_manager[i]);
	    #endif
    }
    
    sprog * found = NULL;
    (void) hashmap_get(cmdmap, cmdname, (void **) &found);
	return found;

}

 inline int  cmpfunc(const void * a, const void * b)
{
   return ( *(ityp*)a - *(ityp*)b );
}

 inline int mpfr_cmpfunc(const void * a, const void * b)
{
   return mpfr_cmp(*(mpfr_t*)a, *(mpfr_t*)b);
}

__MATHSUITE inline void  _showUsage(sprog * const prog)
{
    msprintf(COLOR_ERROR, "\nUSAGE: ");
    msprintf(COLOR_SYSTEM, "%s %s;", prog->cmdname, prog->usage);
    return;
}

__MATHSUITE void  printUsage(sprog * const prog)
{
    _showUsage(prog);
    // prog->program_function(0, NULL);
    return;
}

__MATHSUITE void  prepareToExit(void)
{
	if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
		if(file_exists(PHONY_FILE))
		{
			remove(PHONY_FILE);
			return;
		}
		else
			(void) checkForFHErrors(PHONY_FILE, "w");
			
	session * const currentSession = accessCurrentSession();
	
	if(currentSession->MLSystem.list)
		delAllMatrixFromMatrixList();
		
	if(currentSession->MLSystem.lmpMatrix && currentSession->MLSystem.lmpMatrix->matrix)	
		matrixFree(&(currentSession->MLSystem.lmpMatrix->matrix), currentSession->MLSystem.lmpMatrix->dim);
		
	flushAllMemoizersBuffers();
    _backupColFile();
    
    #ifndef __DISABLE_MULTIUSER
    if(access(server_mode) && access(lastSessionSocket) != INVALID_SOCKET)
		return;
	#endif

    if(isSett(BOOLS_ITEMSAUTOSAVING))
        for(dim_typ i=0; i<MAX_LISTS; ++i)
            updAll(i);
            
    searchProgram(SEARCHPROGRAM_FREESTRING);
    msprintf(COLOR_CREDITS, EXIT_MESSAGE);
    PRINTL();
    return;
}

__MATHSUITE inline void resetLmpMatrix()
{
	accessCurrentLmpMatrix = calloc(1, sizeof(matrixObj));
    errMem(accessCurrentLmpMatrix, VSPACE);
    matrixObj * const currentLmpMatrix = accessCurrentLmpMatrix;
	currentLmpMatrix->matrix = NULL;
    currentLmpMatrix->dim[ROWS] = STARTING_LMP_ROWS;
    currentLmpMatrix->dim[COLUMNS] = STARTING_LMP_COLUMNS;
    return;
}

#ifndef __DISABLE_SERVER
	__MATHSUITE bool getAddressInfo(char *string, struct sockaddr * addr, socklen_t len)
	{
		
		// no reverse lookup in getnameinfo
		int niflags = NI_NUMERICSERV | NI_NUMERICHOST;
		char IP[INET6_ADDRSTRLEN];
		char port[PORT_STRLEN];
		
		// display local address of the socket
		if(getnameinfo(addr, len, IP, INET6_ADDRSTRLEN, port, PORT_STRLEN, NI_NUMERICSERV | NI_NUMERICHOST))
		{
			printErr(ERROR_INPUTOUTPUT, "getnameinfo()");
			return false;
		}
	
		sprintf(string, "%s:%s", IP, port);
		return true;
	}
	
	__MATHSUITE inline void closeUnsetSocket(fd_set * set, int * sock)
	{
		session * sockSession = NULL;
		char socketString[MAX_SOCKETLENGTH];
		register int backupSocket = access(lastSessionSocket);
		itoa(*sock, socketString, 10);
		msyprintf(COLOR_SYSTEM, "Peer closed connection");
		
		if(hashmap_get(access(sessions), socketString, (void**) &sockSession) == MAP_OK)
		{
			access(lastSessionSocket) = *sock;
			userObj ** user = &accessCurrentUser;
			session * const currentSession = accessCurrentSession();
			
			if((*user)->level != LEVEL_GUEST)
			{
				access(lastSessionSocket) = backupSocket;
				msyprintf(COLOR_SYSTEM, ": %s.\n\n", access(curLayout)->database.identifierMode ? (*user)->email : (*user)->username);
				access(lastSessionSocket) = *sock;
				
				if(isSett(BOOLS_ITEMSAUTOSAVING))
				{
					updAll(ENVS);
					updAll(MATRICES);
				}
				
				delAllMatrixFromMatrixList();
				currentSession->MLSystem.list = NULL;
				
				// here you have to free lmpMatrix
				if(currentSession->MLSystem.lmpMatrix)
					matrixFree(&(currentSession->MLSystem.lmpMatrix->matrix), currentSession->MLSystem.lmpMatrix->dim);
					
				currentSession->MLSystem.lmpMatrix = NULL;
				currentSession->MLSystem.cur_item = currentSession->MLSystem.itemsno = 0;
				
				free(*user);
				*user = NULL;
					
			}	
			else
				msyprintf(COLOR_SYSTEM, ".\n\n");	
			
			free(currentSession);
			hashmap_remove(access(sessions), socketString);	
			access(lastSessionSocket) = backupSocket;	
		}
		else
			msyprintf(COLOR_SYSTEM, ".\n\n");
				
		closesocket(*sock);
		FD_CLR(*sock, set);
		*sock = INVALID_SOCKET;
		return;
	}
	
	__MATHSUITE inline sel_typ getFamily(const register sel_typ family)
	{
		return family == 4 ? AF_INET : family == 6 ? AF_INET6 : AF_UNSPEC;
	}
#endif

__MATHSUITE void removeDirEnt(const char dirname[static MAX_PATH_LENGTH], const register dim_typ extNo, char extensions[static extNo][MAX_EXTENSION_LENGTH])
{
	dim_typ i;
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir (dirname)))
	{
		if(extensions)
			while((ent = readdir (dir)))
			{
				char * ext = strrchr(ent->d_name, EXTENSION_DOT_CHAR);
				if(ext)
				{
					++ ext;
					for(i=0; i<sizeof(extensions); ++i)
						if(!strcmp(ext, extensions[i]))
							(void) remove(ent->d_name);
				}
			}
		else
			for( ; (ent = readdir (dir)); (void) remove(ent->d_name));
		closedir (dir);
	}
	return;
}

__MATHSUITE void removeCriticalFiles(const bool mode)
{
	#ifdef WINOS
		if(!(mode & RCF_SKIP)) 
		{
			removeDirEnt(MUTANT_DIR"\\", 2, (char [2][MAX_EXTENSION_LENGTH]) { "c", "h" } );
			(void) remove(MUTANT_DIR"\\mathSuite.ico");
			(void) remove(MUTANT_DIR"\\mathSuite_private.rc");
		}
		if(mode & (RCF_EXPREVAL|RCF_DEEP))
			return;
		if(mode & RCF_DEEP)
		{
			(void) remove(MUTANT_DIR"\\obj\\expreval.o");
			(void) remove(MUTANT_DIR"\\obj\\bessel.o");
			(void) remove(MUTANT_DIR"\\obj\\values.o");
			removeDirEnt(MUTANT_DIR"\\obj\\jburkardt\\", 0, NULL);
		}
		else
		{
			if(mode & RCF_EXPREVAL)
				(void) remove(MUTANT_DIR"\\obj\\expreval.o");
			(void) remove(MUTANT_DIR"\\obj\\mengine.o");
			(void) remove(MUTANT_DIR"\\obj\\algebra.o");
			(void) remove(MUTANT_DIR"\\obj\\main.o");
			(void) remove(MUTANT_DIR"\\obj\\settings.o");
			(void) remove(MUTANT_DIR"\\obj\\programs.o");
			(void) remove(MUTANT_DIR"\\obj\\lists_manager.o");
			(void) remove(MUTANT_DIR"\\obj\\mss_manager.o");
			(void) remove(MUTANT_DIR"\\obj\\envs_manager.o");
			(void) remove(MUTANT_DIR"\\obj\\mat_manager.o");
			(void) remove(MUTANT_DIR"\\obj\\cols_manager.o");
			(void) remove(MUTANT_DIR"\\obj\\layouts_manager.o");
			(void) remove(MUTANT_DIR"\\obj\\logs_manager.o");
		}
	#else
		if(!(skip & RCF_SKIP)) 
		{
			removeDirEnt(MUTANT_DIR"/", 2, (char [2][MAX_EXTENSION_LENGTH]) { "c", "h" } );
			(void) remove(MUTANT_DIR"/mathSuite.ico");
			(void) remove(MUTANT_DIR"/mathSuite_private.rc");
		}
		if(mode & (RCF_EXPREVAL|RCF_DEEP))
			return;
		if(mode & RCF_DEEP)
		{
			(void) remove(MUTANT_DIR"/obj/expreval.o");
			(void) remove(MUTANT_DIR"/obj/bessel.o");
			(void) remove(MUTANT_DIR"/obj/values.o"); 
			removeDirEnt(MUTANT_DIR"/obj/jburkardt/", 0, NULL);
		}
		else
		{
			if(mode & RCF_EXPREVAL)
				(void) remove(MUTANT_DIR"/obj/expreval.o");
			(void) remove(MUTANT_DIR"/obj/mengine.o");
			(void) remove(MUTANT_DIR"/obj/algebra.o");
			(void) remove(MUTANT_DIR"/obj/main.o");
			(void) remove(MUTANT_DIR"/obj/settings.o");
			(void) remove(MUTANT_DIR"/obj/programs.o");
			(void) remove(MUTANT_DIR"/obj/lists_manager.o");
			(void) remove(MUTANT_DIR"/obj/mss_manager.o");
			(void) remove(MUTANT_DIR"/obj/envs_manager.o");
			(void) remove(MUTANT_DIR"/obj/mat_manager.o");
			(void) remove(MUTANT_DIR"/obj/cols_manager.o");
			(void) remove(MUTANT_DIR"/obj/layouts_manager.o");
			(void) remove(MUTANT_DIR"/obj/logs_manager.o");
		}
	#endif
	return;
}

static inline unsigned char * _MDC2(const void *inp, size_t bytes, unsigned char *md)
{
	return MDC2((const unsigned char *) inp, bytes, md);
}

/*
static inline unsigned char * _MD2(const void *inp, size_t bytes, unsigned char *md)
{
	return MD2((const unsigned char *) inp, bytes, md);
}
*/

static inline unsigned char * _MD4(const void *inp, size_t bytes, unsigned char *md)
{
	return MD4((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _MD5(const void *inp, size_t bytes, unsigned char *md)
{
	return MD5((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _RIPEMD160(const void *inp, size_t bytes, unsigned char *md)
{
	return RIPEMD160((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _SHA1(const void *inp, size_t bytes, unsigned char *md)
{
	return SHA1((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _SHA224(const void *inp, size_t bytes, unsigned char *md)
{
	return SHA224((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _SHA256(const void *inp, size_t bytes, unsigned char *md)
{
	return SHA256((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _SHA384(const void *inp, size_t bytes, unsigned char *md)
{
	return SHA384((const unsigned char *) inp, bytes, md);
}

static inline unsigned char * _SHA512(const void *inp, size_t bytes, unsigned char *md)
{
	return SHA512((const unsigned char *) inp, bytes, md);
}


__MATHSUITE void hexhash(const register sel_typ type, const char *str, char *hashRep)
{
	static const dim_typ lens[MAX_HEXHASHES] =
	{
		MDC2_DIGEST_LENGTH,
		// MD2_DIGEST_LENGTH,
		MD4_DIGEST_LENGTH,
		MD5_DIGEST_LENGTH,
		RIPEMD160_DIGEST_LENGTH,
		WHIRLPOOL_DIGEST_LENGTH,
		SHA_DIGEST_LENGTH,
		SHA224_DIGEST_LENGTH,
		SHA256_DIGEST_LENGTH,
		SHA384_DIGEST_LENGTH,
		SHA512_DIGEST_LENGTH
	};
	
	static unsigned char * (* hashFunctions[MAX_HEXHASHES])(const void *, size_t, unsigned char *) =
	{
		_MDC2,
		// _MD2,
		_MD4,
		_MD5,
		_RIPEMD160,
		WHIRLPOOL,
		_SHA1,
		_SHA224,
		_SHA256,
		_SHA384,
		_SHA512
	};
	
	register dim_typ len = lens[type];
	unsigned char hash[len];
	hashFunctions[type]((const void *) str, strlen(str), hash);
	strcpy(hashRep, NULL_CHAR);
	
	for(sel_typ i = 0; i < len; ++i)
		sprintf(hashRep, "%s%02x", hashRep, (unsigned char) hash[i]);
	return;
}

__MATHSUITE void  _handleCmdLine(char * str)
{	
	dim_typ i = 2;
	userObj * user = accessCurrentSession()->user;
	char *cmdtab[MAX_ARGS];
	const register dim_typ len = strlen(str);
	
	sprog * prog;
	
	if(access(mss))
    {
        // char mss_apnt[MAX_FILE_LINES+SIGN_STRING];
        // sprintf(mss_apnt, CMD_BCALC" %s", str);
        // strcpy(str, mss_apnt);
        cmdtab[1] = str;
        prog = &main_menu[MAIN_BCALCULATOR];
    }
	else
	{
		cmdtab[0]=strtok(str,BLANK_STRING);
		prog = searchProgram(cmdtab[0]);
	}
	
	// const sprog * const prog = searchProgram(cmdtab[0]);
	
	if(prog)
    {
    	#ifndef __DISABLE_DATABASE
			if(user->level < prog->level)
				printErr(ERROR_PERMISSIONDENIED, "You haven't permission to execute: %s command.\nRequired Level: %s", prog->cmdname, suite_c.usertypes_names[prog->level]);
			else
		#endif
			{
				if(!access(mss))
				{
					if(prog->argc == ONEARG_ARGC)
						if(len > strlen(cmdtab[0])+1)
							cmdtab[1] = &str[strlen(cmdtab[0])+1];
						else
						{
							printErr(ERROR_INPUTOUTPUT, "You have to insert the single Argument required by the %s command", prog->cmdname);
							return;
						}
					else
						for(i=0; cmdtab[i] != NULL; cmdtab[++i] = strtok(NULL, BLANK_STRING));
						
					if(access(server_mode))
					{
						
						if(prog->argc != ONEARG_ARGC)
						{
							if(prog->argc == INVALID_ARGC)
							{
								printErr(ERROR_INPUTOUTPUT, "This is a blocking I/O command and cannot be executed by a Client");
								return;
							}
							
							if(i-1 < prog->argc)
							{
								printErr(ERROR_INPUTOUTPUT, "This is a semi-blocking I/O command.\nMinimum Number of Arguments: %hu", prog->argc);
								return;
							}
						}
							
					}
				}
				
        		prog->program_function(i-1, &cmdtab[1]);
        	}
	}
	else
    {
	
		if(!strcmp(cmdtab[0], ECHO_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_USER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to do \""ECHO_CMD"\" command", suite_c.usertypes_names[LEVEL_USER]);
				return;
			}
			#endif
			if(len > strlen(ECHO_CMD)+1) 
				msprintf(COLOR_USER, "\n%s\n", &str[strlen(ECHO_CMD)+1]);
			else
				printErr(ERROR_INPUTOUTPUT, "You have to enter at least one param to "ECHO_CMD" Command");
			return;
		}
		else if(!strcmp(cmdtab[0], ERROR_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_USER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to do \""ERROR_CMD"\" command", suite_c.usertypes_names[LEVEL_USER]);
				return;
			}
			#endif
			if(len > strlen(ERROR_CMD)+1) 
				printErr(ERROR_INPUTOUTPUT, &str[strlen(ERROR_CMD)+1]);
			else
				printErr(ERROR_INPUTOUTPUT, "You have to enter at least one param to "ECHO_CMD" Command");
			return;
		}
		else if((!strcmp(cmdtab[0], PRINT_CMD)) || !strcmp(cmdtab[0], SYSPRINT_CMD))
		{
			const bool which = !strcmp(cmdtab[0], PRINT_CMD);
			const char *cmdname = which ?PRINT_CMD:SYSPRINT_CMD;
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_USER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to do \"%s\" command", suite_c.usertypes_names[LEVEL_USER], cmdname);
				return;
			}
			#endif
			if(len > strlen(cmdname)+1) 
			{
				cmdtab[1] = &str[strlen(cmdname)+1];
				char *text = cmdtab[1];
		    	char append[MAX_BUFSIZ];
		    	char backup[MAX_BUFSIZ];
				size_t enter_point, exit_point, i;
				mpfr_t tmp;
				mpfr_init(tmp);
		    	for(enter_point=0 ; text[enter_point] != '\0' ; ++enter_point)
		    		if(text[enter_point] == OPENING_BRACKET && text[enter_point-1]  != ESCAPE_CHARACTER)
		    		{
			    		strcpy(append, NULL_CHAR);
			    		for(exit_point=0 ; text[1+enter_point+exit_point] != CLOSING_BRACKET && text[2+enter_point+exit_point] != ESCAPE_CHARACTER ; ++exit_point)
			    			sprintf(append, "%s%c", append, text[1+enter_point+exit_point]);
			    		strcpy(backup, &text[enter_point+exit_point+2]);
			    		(void) parse(append, &tmp);
			    		mpfr_sprintf(append, "%Rf", tmp);
			    		for (exit_point=0; (*(text + enter_point+exit_point) = *(append + exit_point)) != '\0'; ++ exit_point);
			    		for(i=0 ; (*(text + enter_point+exit_point+i) = *(backup + i)) != '\0'; ++i);
			    	}
		    	mpfr_clear(tmp);
		    	if(which)
		  			msprintf(COLOR_DEFAULT, "\n%s\n", text);
		    	else
		    		msyprintf(COLOR_SYSTEM, "\n%s\n", text);
			}
			else
				printErr(ERROR_INPUTOUTPUT, "You have to enter at least one param to %s Command", cmdname);
			return;
		}
		else if(!strcmp(cmdtab[0], DISABLEFEATURES_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_DEVELOPER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to do \""DISABLEFEATURES_CMD"\" command", suite_c.usertypes_names[LEVEL_DEVELOPER]);
				return;
			}
			#endif
			char *buf = len > strlen(DISABLEFEATURES_CMD) ? &str[strlen(DISABLEFEATURES_CMD)+1] : NULL_CHAR;
			if(buf)
			{
				cmdtab[1] = calloc(1024, sizeof(char));
				sprintf(cmdtab[1], " %s", buf);
			}
			char buffer[MAX_BUFLEN];
			char *token = replace(cmdtab[1], BLANK_STRING, BLANK_STRING"-D__DISABLE_");
			++ token;
			msyprintf(COLOR_SYSTEM, MUTANT_CREATING);
			sprintf(buffer, MAKE_EXECUTABLE" all OPTS=%s && "MUTANT_REMOVEDIRSCMD, token);
			// system("cd "MUTANT_DIR" && rm -f obj/expreval.o && rm -f obj/mengine.o && rm -f obj/algebra.o && rm -f obj/mengine.o && rm -f obj/main.o && rm -f obj/settings.o && rm -f obj/programs.o && rm -f obj/lists_manager.o && rm -f obj/mss_manager.o && rm -f obj/envs_manager.o && rm -f obj/mat_manager.o && rm -f obj/cols_manager.o && rm -f obj/layouts_manager.o && rm -f obj/logs_manager.o");
			removeCriticalFiles(RCF_EXPREVAL);
			// system(MUTANT_REMOVEDIRSFASTCMD);
			system(EXTRACTSRC_CMD);
			system(buffer);
			free(token);
			if(isSett(BOOLS_RUNMUTCODEAFTERCOMP))
			{
				msprintf(COLOR_USER, MUTANT_WELCOME);
				system(MUTANT_FILENAME);
			}	
			if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
				exit(EXIT_MUTANTCODE);
			return;
		}
		else if(!strcmp(cmdtab[0], DISABLEFEATURESFAST_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_DEVELOPER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to do \""DISABLEFEATURESFAST_CMD"\" command", suite_c.usertypes_names[LEVEL_DEVELOPER]);
				return;
			}
			#endif
			char *buf = len > strlen(DISABLEFEATURESFAST_CMD) ? &str[strlen(DISABLEFEATURESFAST_CMD)+1] : NULL_CHAR;
			if(buf)
			{
				cmdtab[1] = calloc(1024, sizeof(char));
				sprintf(cmdtab[1], " %s", buf);
			}
			char buffer[MAX_BUFLEN];
			char *token = replace(cmdtab[1], BLANK_STRING, BLANK_STRING"-D__DISABLE_");
			++ token;
			msyprintf(COLOR_SYSTEM, MUTANT_CREATING);
			sprintf(buffer, MAKE_EXECUTABLE" all OPTS=%s && "MUTANT_REMOVEDIRSCMD, token);
			// system("cd "MUTANT_DIR" && rm -f obj/mengine.o && rm -f obj/algebra.o && rm -f obj/mengine.o && rm -f obj/main.o && rm -f obj/settings.o && rm -f obj/programs.o && rm -f obj/lists_manager.o && rm -f obj/mss_manager.o && rm -f obj/envs_manager.o && rm -f obj/mat_manager.o && rm -f obj/cols_manager.o && rm -f obj/layouts_manager.o && rm -f obj/logs_manager.o");
			removeCriticalFiles(RCF_NONE);
			// system(MUTANT_REMOVEDIRSFASTCMD);
			system(EXTRACTSRC_CMD);
			system(buffer);
			free(token);
			if(isSett(BOOLS_RUNMUTCODEAFTERCOMP))
			{
				msprintf(COLOR_USER, MUTANT_WELCOME);
				system(MUTANT_FILENAME);
			}
			if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
				exit(EXIT_MUTANTCODE);
			return;
		}
		else if(!strcmp(cmdtab[0], DISABLEDEEP_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_OWNER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to execute \""REBUILD_CMD"\" command", suite_c.usertypes_names[LEVEL_OWNER]);
				return;
			}
			#endif
			msyprintf(COLOR_ERROR, "WARNING: \""DISABLEDEEP_CMD"\" command is extremely dangerous!\nIt may cause potentially unrepairable damages to the Mutant Version!");
			(void) getch();
			char *buf = len > strlen(DISABLEDEEP_CMD) ? &str[strlen(DISABLEDEEP_CMD)+1] : NULL_CHAR;
			if(buf)
			{
				cmdtab[1] = calloc(1024, sizeof(char));
				sprintf(cmdtab[1], " %s", buf);
			}
			char buffer[MAX_BUFLEN];
			char *token = replace(cmdtab[1], BLANK_STRING, BLANK_STRING"-D__DISABLEDEEP_");
			++ token;
			msyprintf(COLOR_SYSTEM, MUTANT_CREATING);
			sprintf(buffer, MAKE_EXECUTABLE" all OPTS=%s && "MUTANT_REMOVEDIRSCMD, token);
			// system("cd "MUTANT_DIR" && rm -f obj/expreval.o && rm -f obj/bessel.o && rm -f obj/values.o && rm -f obj/jburkardt/*");
			removeCriticalFiles(RCF_DEEP);
			// system(MUTANT_REMOVEDIRSFASTCMD);
			system(EXTRACTSRC_CMD);
			system(buffer); 
			free(token);
			if(isSett(BOOLS_RUNMUTCODEAFTERCOMP))
			{
				msprintf(COLOR_USER, MUTANT_WELCOME);
				system(MUTANT_FILENAME);
			}
			if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
				exit(EXIT_MUTANTCODE);
			return;
		}
		else if(!strcmp(cmdtab[0], MUTABLE_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_DEVELOPER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to execute \""MUTABLE_CMD"\" command", suite_c.usertypes_names[LEVEL_DEVELOPER]);
				return;
			}
			#endif
			cmdtab[1] = len > strlen(MUTABLE_CMD) ? &str[strlen(MUTABLE_CMD)+1] : NULL_CHAR;
			sprintf(cmdtab[1], " %s", cmdtab[1]);
			char buffer[MAX_BUFLEN];
			char *token = replace(cmdtab[1], BLANK_STRING, BLANK_STRING"-f");
			++ token;
			msyprintf(COLOR_SYSTEM, MUTANT_CREATING);
			sprintf(buffer, MAKE_EXECUTABLE" all OPTS=%s && "MUTANT_REMOVEDIRSCMD, token);
			removeCriticalFiles(RCF_EXPREVAL | RCF_DEEP);
			// system(MUTANT_REMOVEDIRSFASTCMD);
			system(EXTRACTSRC_CMD);
			system(buffer);
			free(token);
			if(isSett(BOOLS_RUNMUTCODEAFTERCOMP))
			{
				msprintf(COLOR_USER, MUTANT_WELCOME);
				system(MUTANT_FILENAME);
			}
			if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
				exit(EXIT_MUTANTCODE);
			return;
		}
		else if(!strcmp(cmdtab[0], INLINE_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_DEVELOPER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to execute \""INLINE_CMD"\" command", suite_c.usertypes_names[LEVEL_DEVELOPER]);
				return;
			}
			#endif
			if(len > strlen(PRINT_CMD))
			{
				cmdtab[1] = &str[strlen(INLINE_CMD)+1];
				FILE *fp = fopen(MSS_TMPFILENAMSRC, "w");
				fprintf(fp, "#include <stdio.h>\n");
				fprintf(fp, "#include <stdlib.h>\n");
				char *command = cmdtab[1];
				if(strrchr(command, OPENING_BRACKET))
				{
					char *headers_string = strtok(command, OPENING_BRACKETSTRING);
					char *headers[MAX_HEADERS];
					// headers_string = strtok(NULL, "{");
					headers_string = strtok(headers_string, CLOSING_BRACKETSTRING);
					command = strtok(NULL, "}");
					headers[0]=strtok(headers_string, ",");
					for(i=0; headers[i] != NULL; headers[++i]=strtok(NULL,","))
						fprintf(fp, "#include <%s>\n", headers[i]);
				}
				fprintf(fp, "int main(int argc, char *argv[]) {\n");
				fprintf(fp, command);
				fprintf(fp, "\nreturn 0;\n}");
				fclose(fp);
				msyprintf(COLOR_SYSTEM, "\nExecuting Inline Code.\n");
				system(MSS_COMPILER" "MSS_TMPFILENAMSRC" -o "MSS_TMPFILENAMEXE" && "MSS_TMPFILENAMEXE);
				msyprintf(COLOR_SYSTEM, "\nInline Code Terminated.\n");
				(void) remove(MSS_TMPFILENAMSRC);
				(void) remove(MSS_TMPFILENAMEXE);
			}
			else
				printErr(ERROR_INPUTOUTPUT, "You have to enter at least one statement to "INLINE_CMD" Command");
			return;
		}
		else if(!strcmp(cmdtab[0], REBUILD_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_OWNER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to execute \""REBUILD_CMD"\" command", suite_c.usertypes_names[LEVEL_OWNER]);
				return;
			}
			#endif
			cmdtab[1] = len > strlen(REBUILD_CMD) ? &str[strlen(REBUILD_CMD)+1] : NULL_CHAR;
			char buffer[MAX_BUFLEN];
			msyprintf(COLOR_SYSTEM, MUTANT_CREATING);
			sprintf(buffer, MAKE_EXECUTABLE" %s all && "MUTANT_REMOVEDIRSCMD, cmdtab[1]);
			system(MAKE_EXECUTABLE" clean");
			removeCriticalFiles(RCF_EXPREVAL | RCF_DEEP);
			// system(MUTANT_REMOVEDIRSFASTCMD);
			system(EXTRACTSRC_CMD);
			system(buffer);
			if(isSett(BOOLS_RUNMUTCODEAFTERCOMP))
			{
				msprintf(COLOR_USER, MUTANT_WELCOME);
				system(MUTANT_FILENAME);
			}
			if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
				exit(EXIT_MUTANTCODE);
			return;
		}
		else if(!strcmp(cmdtab[0], MAKE_CMD))
		{
			#ifndef __DISABLE_DATABASE
			if(user->level < LEVEL_DEVELOPER)
			{
				printErr(ERROR_PERMISSIONDENIED, "You have to be a %s to execute \""MAKE_CMD"\" command", suite_c.usertypes_names[LEVEL_DEVELOPER]);
				return;
			}
			#endif
			cmdtab[1] = len > strlen(MAKE_CMD) ? &str[strlen(MAKE_CMD)+1] : NULL_CHAR;
			char buffer[MAX_BUFLEN];
			msyprintf(COLOR_SYSTEM, MUTANT_CREATING);
			sprintf(buffer, MAKE_EXECUTABLE" %s && "MUTANT_REMOVEDIRSCMD, cmdtab[1]);
			removeCriticalFiles(RCF_EXPREVAL | RCF_DEEP);
			// system(MUTANT_REMOVEDIRSFASTCMD);
			system(EXTRACTSRC_CMD);
			system(buffer);
			if(isSett(BOOLS_RUNMUTCODEAFTERCOMP))
			{
				msprintf(COLOR_USER, MUTANT_WELCOME);
				system(MUTANT_FILENAME);
			}
			if(isSett(BOOLS_EXITAFTERMUTCODEXEC))
				exit(EXIT_MUTANTCODE);
			return;
		}
		
		unsigned char hashRep[HEXHASH_LEN+1] = NULL_CHAR;
    	hexhashSHA512(cmdtab[0], hashRep);
    
	    #ifndef __DISABLE_DATABASE
	    if(user->level >= LEVEL_ADMIN && !strcmp(hashRep, DEVMODE_CMD))
	    {
	    	userObj * user2 = NULL;
	    	const bool result = user->level == LEVEL_ADMIN ? readUserFromUserName(USERNAME_DEVELOPER, &user2) && user2 : false;
	    	if(result)
	    	{
	    		free(user);
	    		user = user2;
	    	}
	    	else
	    		user->level = LEVEL_DEVELOPER;
	    	if(user->level == LEVEL_DEVELOPER)
	    		msyprintf(COLOR_SYSTEM, "\n%s Mode has been successfully enabled.\n", suite_c.usertypes_names[LEVEL_DEVELOPER]);
			return;
	    }
	    
	    if(user->level >= LEVEL_DEVELOPER && !strcmp(hashRep, OWNERMODE_CMD))
	    {
	    	userObj * user2 = NULL;
	    	if(readUserFromUserName(USERNAME_OWNER, &user) && user)
	    	{
	    		free(user);
	    		user = user2;
	    	}
	    	else
	    		user->level = LEVEL_OWNER;
	    	msyprintf(COLOR_PURPLE, "\n%s Mode has been successfully enabled.\n", suite_c.usertypes_names[LEVEL_OWNER]);
			return;
	    }
	    #endif
	    
	    printErr(1, "SubProgram hasn't been found in Program Macro-List");
	}
	
    return;
}

__MATHSUITE bool  _execScriptFiles(const char path[static MAX_PATH_LENGTH])
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
                
            // for(cmdtab[0]=strtok(str,BLANK_STRING),i=0; cmdtab[i] != NULL; cmdtab[++i] = strtok(NULL, BLANK_STRING));
            
			_handleCmdLine(str);
        }

        fclose(fp);
        return true;
    }

    return false;
}

__MATHSUITE bool  _lfLoader(const char path[static MAX_PATH_LENGTH])
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

__MATHSUITE bool  _lfCreate(const char path[static MAX_PATH_LENGTH])
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
		__MATHSUITE  static inline void itoa(const int val, char string[])
		{
			sprintf(string, "%d", val);
			return;
		}
	#endif

	 static inline void ftoa(char *string, const float value, const fsel_typ prec)
	{
		sprintf(string, "%.*f", prec, value);
		return;
	}

	 XMLCALL inline xmlDoc *   xmlInit(const char file_name[static XML_FILENAMES_LENGTH], xmlXPathContext ** xpathCtx)
	{
		xmlInitParser();
		LIBXML_TEST_VERSION;
		xmlDoc * doc = xmlParseFile( file_name );
		(*xpathCtx) = xmlXPathNewContext( doc );
		return doc;
	}

	 XMLCALL inline void   xmlExit(const char file_name[static XML_FILENAMES_LENGTH], xmlDoc ** doc, xmlXPathObject ** xpathObj, xmlXPathContext ** xpathCtx)
	{
		xmlSaveFileEnc(file_name, (*doc), XML_ENCODING);
		xmlFreeDoc((*doc));
    	xmlCleanupParser();
    	(*doc) = (xmlDoc*) NULL;
    	(*xpathObj) = (xmlXPathObject*) NULL;
    	(*xpathCtx) = (xmlXPathContext*) NULL;
    	return;
	}

	 XMLCALL inline bool   xmlWriteInt(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const int value)
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

	 XMLCALL inline bool   xmlWriteBool(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const bool value)
	{
		static const char bools_identifiers[2][SIGN_STRING] =
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

	 XMLCALL inline bool   xmlWriteFloat(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const float value)
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

	 XMLCALL inline bool   xmlWriteString(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, const char string[static MAX_XML_FIELDSTRINGS])
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

	 XMLCALL inline bool   xmlGetInt(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, int * value)
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

	 XMLCALL inline bool   xmlGetBool(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, bool * value)
	{
		static const char bools_identifiers[2][SIGN_STRING] =
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

	 XMLCALL inline bool   xmlGetFloat(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, float * value)
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

	 XMLCALL inline bool   xmlGetString(xmlXPathObject ** xpathObj, xmlXPathContext * xpathCtx, const char * nodeAddress, char string[static MAX_XML_FIELDSTRINGS])
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

	__MATHSUITE XMLCALL void  getProgramSettings(dim_typ which_layout)
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
			
		#ifndef __DISABLE_MULTIUSER
		
		#ifndef __DISABLE_SERVER
			if(!xmlWriteString(&xpathObj, xpathCtx, "/settings/serverSettings/server", cur_layout->server.server))
				printErr(5, "In Writing /settings/serverSettings/server field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/serverSettings/port", cur_layout->server.port))
				printErr(5, "In Writing /settings/serverSettings/port field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/serverSettings/family", cur_layout->server.family))
				printErr(5, "In Writing /settings/serverSettings/family field");	
			if(!xmlWriteBool(&xpathObj, xpathCtx, "/settings/serverSettings/resolution", cur_layout->server.resolution))
		    	printErr(5, "In Writing /settings/serverSettings/resolution field");	
		    if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/serverSettings/nListeningSockets", cur_layout->server.nListeningSockets))
				printErr(5, "In Writing /settings/serverSettings/nListeningSockets field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/serverSettings/nAcceptableConnections", cur_layout->server.nAcceptableConnections))
				printErr(5, "In Writing /settings/serverSettings/nAcceptableConnections field");	
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/serverSettings/bufferLength", cur_layout->server.bufferLength))
				printErr(5, "In Writing /settings/serverSettings/bufferLength field");
		#endif
		#ifndef __DISABLE_DATABASE
			if(!xmlWriteString(&xpathObj, xpathCtx, "/settings/databaseSettings/server", cur_layout->database.server))
				printErr(5, "In Writing /settings/databaseSettings/server field");	
			if(!xmlWriteString(&xpathObj, xpathCtx, "/settings/databaseSettings/username", cur_layout->database.username))
				printErr(5, "In Writing /settings/databaseSettings/username field");
			if(!xmlWriteString(&xpathObj, xpathCtx, "/settings/databaseSettings/password", cur_layout->database.password))
				printErr(5, "In Writing /settings/databaseSettings/password field");
			if(!xmlWriteString(&xpathObj, xpathCtx, "/settings/databaseSettings/database", cur_layout->database.database))
				printErr(5, "In Writing /settings/databaseSettings/database field");
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/databaseSettings/maxUsers", cur_layout->database.maxUsers))
				printErr(5, "In Writing /settings/databaseSettings/maxUsers field");	
			if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/databaseSettings/defaultLevel", cur_layout->database.defaultLevel))
				printErr(5, "In Writing /settings/databaseSettings/defaultLevel field");
			if(!xmlWriteBool(&xpathObj, xpathCtx, "/settings/databaseSettings/identifierMode", cur_layout->database.identifierMode))
		    	printErr(5, "In Writing /settings/databaseSettings/identifierMode field");
		#endif
		#endif

		if(!xmlWriteInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxMatrixListItems", cur_layout->max_matrix_list_items))
			printErr(5, "In Writing /settings/matricesOptions/maxMatrixListItems field");
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
	    	char name[MIN_STRING<<2] = NULL_CHAR;
			char strboolized[MIN_STRING<<1] = NULL_CHAR;
			strboolize(suite_c.bools_names[i], strboolized);
	    	sprintf(name, "/settings/booleanKeys/%s", strboolized);
	        if(!xmlWriteBool(&xpathObj, xpathCtx, name, (cur_layout->bools & suite_c.bools[i].bmask) == suite_c.bools[i].bmask))
	        	printErr(5, "In Writing %s field", name);
	    }

	    xmlExit(tmp->path, &doc, &xpathObj, &xpathCtx);
	    return;
	}

	__MATHSUITE XMLCALL void  resetProgramSettings(layoutObj * const tmp, const char path[static MAX_PATH_LENGTH])
	{
		dim_typ i;
		xmlXPathContext * xpathCtx = (xmlXPathContext*) NULL;
	    xmlXPathObject * xpathObj = (xmlXPathObject*) NULL;
	    xmlDoc * doc = xmlInit(path, &xpathCtx);
	    
	    tmp->database.con = NULL;
			
		char ex_char[1];
		if(xmlGetString(&xpathObj, xpathCtx, "/settings/generalSettings/exitChar", ex_char))
			tmp->exit_char = ex_char[0];
		else
		{
			printErr(5, "In Parsing /settings/generalSettings/exitChar field");
			tmp->exit_char = EXIT_CHAR;
		}
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
		
		#ifndef __DISABLE_MULTIUSER
		char tmp_string[400];
		#ifndef __DISABLE_SERVER
			if(xmlGetString(&xpathObj, xpathCtx, "/settings/serverSettings/server", tmp_string))
				strcpy(tmp->server.server, tmp_string);
			else
			{
				printErr(5, "In Parsing /settings/serverSettings/server field");
				itoa(SERVER_ADDRESS, tmp->server.server, 10);
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/serverSettings/port", &tmp_int))
				tmp->server.port = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/serverSettings/port field");
				tmp->server.port = SERVER_PORT;
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/serverSettings/family", &tmp_int))
				tmp->server.family = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/serverSettings/family field");
				tmp->server.family = SERVER_FAMILY;
			}
			bool tmp_bool = false;
			if(xmlGetBool(&xpathObj, xpathCtx, "/settings/serverSettings/resolution", &tmp_bool))
	        	tmp->server.resolution = tmp_bool;
	        else
			{
				printErr(5, "In Parsing /settings/serverSettings/resolution field");
				tmp->server.resolution = DEFAULT_RESOLUTION;
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/serverSettings/nListeningSockets", &tmp_int))
				tmp->server.nListeningSockets = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/serverSettings/nListeningSockets field");
				tmp->server.nListeningSockets = DEFAULT_LISTENINGSOCKETS;
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/serverSettings/nAcceptableConnections", &tmp_int))
				tmp->server.nAcceptableConnections = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/serverSettings/nAcceptableConnections field");
				tmp->server.nAcceptableConnections = MAX_ACCEPTABLECONNECTIONS;
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/serverSettings/bufferLength", &tmp_int))
				tmp->server.bufferLength = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/serverSettings/bufferLength field");
				tmp->server.bufferLength = CLIENTSERVER_BUFSIZE;
			}
		#endif
			
		#ifndef __DISABLE_DATABASE
			if(xmlGetString(&xpathObj, xpathCtx, "/settings/databaseSettings/server", tmp_string))
				strcpy(tmp->database.server, tmp_string);
			else
			{
				printErr(5, "In Parsing /settings/databaseSettings/server field");
				strcpy(tmp->database.server, DEFAULT_DBSERVER);
			}
			if(xmlGetString(&xpathObj, xpathCtx, "/settings/databaseSettings/username", tmp_string))
				strcpy(tmp->database.username, tmp_string);
			else
			{
				printErr(5, "In Parsing /settings/databaseSettings/username field");
				strcpy(tmp->database.username, DEFAULT_DBUSERNAME);
			}
			if(xmlGetString(&xpathObj, xpathCtx, "/settings/databaseSettings/password", tmp_string))
				strcpy(tmp->database.password, tmp_string);
			else
			{
				printErr(5, "In Parsing /settings/databaseSettings/password field");
				strcpy(tmp->database.password, DEFAULT_DBPASSWORD);
			}
			if(xmlGetString(&xpathObj, xpathCtx, "/settings/databaseSettings/database", tmp_string))
				strcpy(tmp->database.database, tmp_string);
			else
			{
				printErr(5, "In Parsing /settings/databaseSettings/database field");
				strcpy(tmp->database.password, DEFAULT_DATABASE);
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/databaseSettings/defaultLevel", &tmp_int))
				tmp->database.defaultLevel = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/databaseSettings/defaultLevel field");
				tmp->database.defaultLevel = DEFAULT_LEVEL;
			}
			if(xmlGetInt(&xpathObj, xpathCtx, "/settings/databaseSettings/maxUsers", &tmp_int))
				tmp->database.maxUsers = tmp_int;
			else
			{
				printErr(5, "In Parsing /settings/databaseSettings/maxUsers field");
				tmp->database.maxUsers = MAX_USERS;
			}
			if(xmlGetBool(&xpathObj, xpathCtx, "/settings/databaseSettings/identifierMode", &tmp_bool))
	        	tmp->database.identifierMode = tmp_bool;
	        else
			{
				printErr(5, "In Parsing /settings/databaseSettings/identifierMode field");
				tmp->database.identifierMode = DEFAULT_IDENTIFIER;
			}
		#endif
		#endif

		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxMatrixListItems", &tmp_int))
			tmp->max_matrix_list_items = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxMatrixListItems field");
			tmp->max_matrix_list_items = MAX_MATRIXLISTITEMS;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxRows", &tmp_int))
			tmp->matrix_max_rows = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxRows field");
			tmp->matrix_max_rows = MAX_ROWS;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/matricesOptions/maxColumns", &tmp_int))
			tmp->matrix_max_columns = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/matricesOptions/maxColumns field");
			tmp->matrix_max_columns = MAX_COLUMNS;
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
			tmp->basecalc_minbase = BCALC_BASECHANGE_MINBASE;
		}
		if(xmlGetInt(&xpathObj, xpathCtx, "/settings/baseConversions/maxBase", &tmp_int))
			tmp->basecalc_maxbase = tmp_int;
		else
		{
			printErr(5, "In Parsing /settings/baseConversions/maxBase field");
			tmp->basecalc_maxbase = BCALC_BASECHANGE_MAXBASE;
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

	    for(i=tmp->bools=0; i<MAX_BOOL_SETTINGS; ++i)
	    {
	    	tmp_bool=false;
			char name[MIN_STRING<<2] = NULL_CHAR;
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

	__MATHSUITE XMLCALL void  _colFileLoader(const char path[static MAX_PATH_LENGTH])
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
	        	msprintf(_COLOR_ERROR, "ERROR: In Parsing /colors/defaultColor field");
	        	COLOR_DEFAULT = DEFAULT_COLOR;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/errorsColor", &tmp_int))
	        	COLOR_ERROR = tmp_int;
	        else
	        {
	        	msprintf(_COLOR_ERROR, "ERROR: In Parsing /colors/errorsColor field");
	        	COLOR_ERROR = _COLOR_ERROR;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/creditsColor", &tmp_int))
	        	COLOR_CREDITS = tmp_int;
	        else
	        {
	        	msprintf(_COLOR_ERROR, "ERROR: In Parsing /colors/creditsColor field");
	        	COLOR_CREDITS = _COLOR_CREDITS;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/userColor", &tmp_int))
	        	COLOR_USER = tmp_int;
	        else
	        {
	        	msprintf(_COLOR_ERROR, "ERROR: In Parsing /colors/userColor field");
	        	COLOR_USER = _COLOR_USER;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/systemColor", &tmp_int))
	        	COLOR_SYSTEM = tmp_int;
	        else
	        {
	        	msprintf(_COLOR_ERROR, "ERROR: In Parsing /colors/systemColor field");
	        	COLOR_SYSTEM = _COLOR_SYSTEM;
	        }
	        if(xmlGetInt(&xpathObj, xpathCtx, "/colors/authorColor", &tmp_int))
	        	COLOR_AUTHOR = tmp_int;
	        else
	        {
	        	msprintf(_COLOR_ERROR, "ERROR: In Parsing /colors/authorColor field");
	        	COLOR_AUTHOR = _COLOR_AUTHOR;
	        }

			xmlExit(access(colors_path), &doc, &xpathObj, &xpathCtx);
	        return;
	    }

	    __MATHSUITE XMLCALL void  _backupColFile(void)
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
	__WINCALL void clrscr()
	{
	    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
	    COORD coord = {0, 0};
	    DWORD count;
	    CONSOLE_SCREEN_BUFFER_INFO csbi;
	    
	    GetConsoleScreenBufferInfo(hStdOut, &csbi);
	    FillConsoleOutputCharacter(hStdOut, ' ', csbi.dwSize.X * csbi.dwSize.Y, coord, &count);
	    SetConsoleCursorPosition(hStdOut, coord);
	}

     __WINCALL HWND WINAPI   GetConsoleWindowNT()
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

     __WINCALL bool   windowsFileHandler(char *name, const char *ext_pattern, const char default_ext[static MAX_EXTENSION_LENGTH], bool mode)
    {
        // HWND hwnd;
        OPENFILENAME ofn;
        bool result = false; 
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

	int inet_pton(int af, const char *src, void *dst)
	{
		struct sockaddr_storage ss;
		int size = sizeof(ss);
		char src_copy[INET6_ADDRSTRLEN+1];
		
		ZeroMemory(&ss, sizeof(ss));
		//
		strncpy (src_copy, src, INET6_ADDRSTRLEN+1);
		src_copy[INET6_ADDRSTRLEN] = 0;
		
		if (WSAStringToAddress(src_copy, af, NULL, (struct sockaddr *)&ss, &size) == 0) 
			switch(af) 
			{
				case AF_INET:
					*(struct in_addr *)dst = ((struct sockaddr_in *)&ss)->sin_addr;
					return 1;
				case AF_INET6:
					*(struct in6_addr *)dst = ((struct sockaddr_in6 *)&ss)->sin6_addr;
					return 1;
			}
			
		return 0;
	}
	
	const char *inet_ntop(int af, const void *src, char *dst, int size)
	{
		struct sockaddr_storage ss;
		unsigned long s = size;
		
		ZeroMemory(&ss, sizeof(ss));
		ss.ss_family = af;
		
		switch(af) 
		{
			case AF_INET:
				((struct sockaddr_in *)&ss)->sin_addr = *(struct in_addr *)src;
				break;
			case AF_INET6:
				((struct sockaddr_in6 *)&ss)->sin6_addr = *(struct in6_addr *)src;
				break;
			default:
				return NULL;
		}
	
		return (WSAAddressToString((struct sockaddr *)&ss, sizeof(ss), NULL, dst, &s) == 0) ? dst : NULL;
	}
	
#else
	 __MATHSUITE   int getch( )
    {
    	int ch;
        struct termios oldt, newt;
        tcgetattr( STDIN_FILENO, &oldt );
        newt = oldt;
        newt.c_lflag &= ~( ICANON | ECHO );
        tcsetattr( STDIN_FILENO, TCSANOW, &newt );
        ch = getchar();
        tcsetattr( STDIN_FILENO, TCSANOW, &oldt );
        return ch;
    }
#endif

 inline const char * const   getFilename(const char path[static MAX_PATH_LENGTH])
{
    const char * const stkr = strrchr(path, '\\');
    return stkr ? stkr+1 : path;
}

 inline void  updInfo(const char server[static SERVER_LENGTH+5])
{

    char title[MAXX_STRING];
	
	if(access(server_mode) && server)
		sprintf(title, PROGRAM_NAME" - [ USRLOG: %s ] - [ SYSLOG: %s ] - [ %s ] - Listening on: %s",
    		getItemsListNo(LOGS) ? getFilename(listNo(access(lists)[LOGS].cur_item, LOGS)->path) : NONE_STRING, access(sysLog) ? getFilename(access(sysLogPath)) : NONE_STRING, getItemsListNo(LAYOUTS) ? getFilename(listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path) : NULL_CHAR, server);
    else
		sprintf(title, PROGRAM_NAME" - [ %s ] - [ %s ] - [ USRLOG: %s ] - [ SYSLOG: %s ] - [ %s ]", getItemsListNo(ENVS) ? getFilename(listNo(access(lists)[ENVS].cur_item, ENVS)->path) : NULL_CHAR, getItemsListNo(MATRICES) ? getFilename(listNo(access(lists)[MATRICES].cur_item, MATRICES)->path) : NULL_CHAR,
    		getItemsListNo(LOGS) ? getFilename(listNo(access(lists)[LOGS].cur_item, LOGS)->path) : NONE_STRING, access(sysLog) ? getFilename(access(sysLogPath)) : NONE_STRING, getItemsListNo(LAYOUTS) ? getFilename(listNo(access(lists)[LAYOUTS].cur_item, LAYOUTS)->path) : NULL_CHAR);

	if(server && server[0] == UPDINFO_UNSTABLECHAR)
		strcat(title, " - UNSTABLE ENVIRONMENT");
	
	#ifdef WINOS
	    if(!SetConsoleTitle(title))
	        printErr(22, "SetConsoleTitle failed with error: %lu", GetLastError());
	#else
		printf("%c]0;%s%c", '\033', title, '\007');
	#endif

    return;
}

 void   SetColor(const sel_typ ForgC)
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

 static inline void  _min_mpfr_cmpfunc(mpfr_t rop, mpfr_t a, mpfr_t b)
{
    mpfr_set(rop, MPFR_MIN(a,b), MPFR_RNDN);
    return;
}

 static inline void  _max_mpfr_cmpfunc(mpfr_t rop, mpfr_t a, mpfr_t b)
{
    mpfr_set(rop, MPFR_MAX(a,b), MPFR_RNDN);
    return;
}

 static inline ityp  _min_cmpfunc(const register ityp a, const register ityp b)
{
    return MIN(a,b); 
}

 static inline ityp  _max_cmpfunc(const register ityp a, const register ityp b)
{
    return MAX(a,b);
}

__MATHSUITE void MINMAX(ityp *rop, ityp *rop2, unsigned dim, ityp vector[static dim], dim_typ *idx, dim_typ *idx2)
{
	if(dim == 1)
	{
		if(rop)
			*rop = vector[0];
		if(rop2)
			*rop2 = vector[0];
		return;
	}
	if(dim == 2)
	{
		if(rop)	
			*rop = MIN(vector[0], vector[1]);
		if(rop2)
			*rop2 = MAX(vector[0], vector[1]);
		return;
	}
	
	dim_typ i;
	const register dim_typ floor_idx = dim>>1;
	ityp min1, min2, max1, max2;

	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		MINMAX(rop ? &min1 : NULL, rop2 ? &max1 : NULL, floor_idx, vector, NULL, NULL);
		#pragma omp section
		MINMAX(rop ? &min2 : NULL, rop2 ? &max2 : NULL, dim-floor_idx, vector + floor_idx, NULL, NULL);
	}
	
	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		if(rop) 
		{
			*rop = MIN(min1, min2);
			if(idx)
				*idx = ((ityp*)lfind(rop, vector, &dim, sizeof(ityp), cmpfunc)) - vector;
		}
		#pragma omp section
		if(rop2)
		{
			*rop2 = MAX(max1, max2);
			if(idx2)
				*idx2 = ((ityp*)lfind(rop2, vector, &dim, sizeof(ityp), cmpfunc)) - vector;
		}
	}
	
	return;
}

__MATHSUITE void  MPFR_MINMAX(mpfr_t rop, mpfr_t rop2, mpfr_t dim, mpfr_t vector[static mpfr_get_ui(dim, MPFR_RNDN)], mpfr_t idx, mpfr_t idx2)
{
	if((!rop) && !rop2)
		return;
		
	if(mpfr_cmp_ui(dim,1) == 0)
	{
		if(rop) 
			mpfr_set(rop, vector[0], MPFR_RNDN);
		if(rop2)
			mpfr_set(rop2, vector[0], MPFR_RNDN);
		return;
	}
	
	if(mpfr_cmp_ui(dim,2) == 0)
	{
		if(rop)
			mpfr_set(rop, MPFR_MIN(vector[0],vector[1]), MPFR_RNDN);
		if(rop2)
			mpfr_set(rop2, MPFR_MAX(vector[0],vector[1]), MPFR_RNDN);
		return;
	}
	
	unsigned ndim = mpfr_get_ui(dim, MPFR_RNDN);
	mpfr_t floor_idx, ceil_idx, min1, min2, max1, max2;
	mpfr_inits(floor_idx, ceil_idx, NULL);
	mpfr_set_ui(floor_idx, mpfr_get_ui(dim, MPFR_RNDN)>>1, MPFR_RNDN);
	mpfr_sub(ceil_idx, dim, floor_idx, MPFR_RNDN);
	
	if(rop)
		mpfr_inits(min1, min2, NULL);
	if(rop2)
		mpfr_inits(max1, max2, NULL);
		
	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		MPFR_MINMAX(rop ? min1 : NULL, rop2 ? max1 : NULL, floor_idx, vector, NULL, NULL);
		#pragma omp section
		MPFR_MINMAX(rop ? min2 : NULL, rop2 ? max2 : NULL, ceil_idx, vector + mpfr_get_ui(floor_idx, MPFR_RNDN), NULL, NULL);
	}

	#pragma omp parallel sections num_threads(2)
	{
		#pragma omp section
		if(rop)
		{
			mpfr_set(rop, MPFR_MIN(min1, min2), MPFR_RNDN);
			if(idx)
				mpfr_set_si(idx, ((mpfr_t*)lfind(rop, vector, &ndim, sizeof(mpfr_t), mpfr_cmpfunc))-vector, MPFR_RNDN);
			mpfr_clears(min1, min2, NULL);
		}
	
		#pragma omp section
		if(rop2)
		{
			mpfr_set(rop2, MPFR_MAX(max1, max2), MPFR_RNDN);
			if(idx2)
				mpfr_set_si(idx2, ((mpfr_t*)lfind(rop2, vector, &ndim, sizeof(mpfr_t), mpfr_cmpfunc))-vector, MPFR_RNDN);
			mpfr_clears(max1, max2, NULL);
		}
	}
	
	mpfr_clears(floor_idx, ceil_idx, NULL); 
	return;
}

__MATHSUITE  ityp   getDiffTime(struct timeval * t1)
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

__MATHSUITE inline bool   isDomainForbidden(mpfr_t val, bool mode)
{
    if(TYPE_DOMAIN(val))
    {
        msprintf(COLOR_ERROR, "\n%sPUT OVERFLOW ERROR.\n\n", mode ? "OUT":"IN");
        return true;
    }
    return false;
}

 inline void    free_foreach(mpfr_t **matrix, const dim_typ algebra_units, const register dim_typ dim[2], bool mode)
{
    if(mode)
        return;
    #pragma omp parallel for num_threads(algebra_units)
    for(dim_typ i=0; i<algebra_units; ++i)
        matrixFree(&matrix[i], dim);
    return;
}

 inline bool   checkErrMem(const void * pntr)
{
    if(!pntr)
    {
        printErr(12, "An error occurred during Heap Dynamic Memory Allocation");
        msyprintf(COLOR_SYSTEM, "\n(Sub)Program Terminating...\n\n");
    	updInfo(UPDINFO_UNSTABLESTRING);
        #ifdef WINOS
            (void) getch();
        #endif
        return true;
    }
    return false;
}

__MATHSUITE bool   matrixAlloc(mpfr_t **matrix, const register dim_typ dim[2])
{
    if(!(*matrix))
        (*matrix) = NULL;

    (*matrix) = calloc(dim[ROWS]*dim[COLUMNS], sizeof(mpfr_t));
    errMem((*matrix), false);
    
   	for(dim_typ i=0; i<dim[ROWS]*dim[COLUMNS]; ++i)
		mpfr_init_set_ui((*matrix)[i], 0, MPFR_RNDN);

    return true;
}

__MATHSUITE void   _matrixFree(mpfr_t **matrix, const register dim_typ dim[static 2], bool mode)
{
    if(mode || !(*matrix))
        return;
        
    for(dim_typ i=0; i<dim[ROWS]*dim[COLUMNS]; ++i)
		mpfr_clear((*matrix)[i]);  

    free((*matrix));
    (*matrix) = NULL; // to avoid dangling references,
    // even if all this pointer passed to this function
    // are generally locally allocated into static void functions.
    return;
}

__MATHSUITE bool   equalMatrix(mpfr_t **matrix1, mpfr_t *matrix2, const register dim_typ dim[static 2], const bool rlc)
{
	if(rlc && !matrixAlloc(matrix1, dim))
		return false;

    dim_typ i, j;
    // Phisically equalling matrix1 values to matrix2 ones
    for(i=0; i<dim[ROWS]; ++i)
        for(j=0; j<dim[COLUMNS]; ++j)
            mpfr_set(*((*matrix1) + dim[COLUMNS]*i + j), *(matrix2 + dim[COLUMNS]*i + j), MPFR_RNDN);

    return true;
}

__MATHSUITE inline void   _flushMemoizersBuffers(sel_typ mode)
{

    free(access(sysMem)[mode].memoizer);
    access(sysMem)[mode].memoizer = NULL;
    access(sysMem)[mode].current_max_index = 0;
    msyprintf(COLOR_SYSTEM, "%s Function Memoizer has been properly flushed.\n\n", suite_c.memoizers_names[mode]);
    return;
}

__MATHSUITE inline void   flushAllMemoizersBuffers(void)
{
	for(dim_typ i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
		if(access(sysMem)[i].current_max_index)
			_flushMemoizersBuffers(i);

	return;
}

__MATHSUITE dim_typ   selectListItem(dim_typ dim, bool mode, const char *string, const char list[static dim][MIN_STRING])
{
    dim_typ item;
    dim_typ i;

    PRINTN();
    msprintf(COLOR_CREDITS, string);
    msprintf(COLOR_CREDITS, ":\n");
    PRINTL();

    if(mode)
    {
        for(i=0; i<dim; ++i)
            printf("- %hu: %s;\n", i, list[i]); // ext_math.functions[i].name);
        printf("- %hu: Back to previous SubProgram.\n", i);

        PRINTL();

        mpfr_t item2;
		
        while(requires(item2, NULL, NULL_CHAR, NULL_CHAR, PARSER_NOSETTINGS) || isNullVal(item2) || mpfr_cmp_ui(item2, (item = mpfr_get_ui(item2, MPFR_RNDN))) || item < 0 || item > dim)
        {
        	mpfr_clear(item2);
			printErr(1, "Invalid inserted value");
		}

		mpfr_clear(item2);
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

#ifndef __DISABLE_DATABASE

__MATHSUITE bool writeVar(const char *name, const char *tab, mpfr_t val)
{
	
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(tab) > TAB_MAXNAMELENGTH || strlen(name) > VAR_MAXNAMELENGTH)
		return false;
	
	char *tmptab = replace(tab, "\\", "\\\\");
	if(tmptab)
		tab = tmptab;
		
	MYSQL_RES *rx;
	char query[QUERY_LENGTH];
	sprintf(query, "CALL "SP_SELECTNAMEFROMTABBYTABNAME"(\"%s\")", tab);
	
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	userObj * user = accessCurrentUser; 
	
	if(!mysql_fetch_row(rx))
	{
		sprintf(query, "CALL "SP_CREATETAB"(@InsertId, \"%s\")",tab);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_free_result(rx);
			mysql_close(con);
			return false;
		}
		
		mysql_next_result(con);
		strcpy(query, "SELECT @InsertId");
		
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		if ((rx = mysql_store_result(con)) == NULL) 
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		register int mid = -1;
		MYSQL_ROW linsidRow = mysql_fetch_row(rx);
		
		if(linsidRow && linsidRow[0]) 
			mid = atoi(linsidRow[0]);
			
		mysql_free_result(rx);
		sprintf(query, "CALL "SP_CREATEUSERTAB"(%d,%d)",user->iduser,mid);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_free_result(rx);
			mysql_close(con);
			return false;
		}
		
		if(name == NULL)
		{
			mysql_free_result(rx);
			return true;
		}
	}
	
	mysql_free_result(rx);
	sprintf(query, "CALL "SP_SELECTVALFROMVARTABVARTABBYVARNAMEANDTABNAME"(\"%s\", \"%s\")",name,tab);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	if(!mysql_fetch_row(rx))
	{
		mysql_free_result(rx);
		sprintf(query, "CALL "SP_SELECTIDTABFROMUSERTABTABBYTABNAMEANDIDUSER"(\"%s\",%d)", tab, user->iduser);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}

		if ((rx = mysql_store_result(con)) == NULL) 
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
	
		if(mysql_fetch_row(rx)) // if(!mysql_fetch_row(rx))
		{	
			mysql_free_result(rx);
			mpfr_sprintf(query, "CALL "SP_CREATEVAR"(@InsertId, \"%s\",\"%Rf\")",name,val);
			mysql_next_result(con);
			if (mysql_query(con, query))
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
				mysql_close(con);
				return false;
			}
			
			mysql_next_result(con);
			strcpy(query, "SELECT @InsertId");
			
			if (mysql_query(con, query))
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
				mysql_close(con);
				return false;
			}
			
			if ((rx = mysql_store_result(con)) == NULL) 
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
				mysql_close(con);
				return false;
			}
			
			register int mid = -1;
			MYSQL_ROW linsidRow = mysql_fetch_row(rx);
			
			if(linsidRow && linsidRow[0]) 
				mid = atoi(linsidRow[0]);
			
			mysql_free_result(rx);
			sprintf(query, "CALL "SP_CREATETABVAR"(\"%s\",%d)",tab,mid);
			mysql_next_result(con);
			if (mysql_query(con, query))
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
				mysql_close(con);
				return false;
			}
			
			return true;
		}
		
		mysql_free_result(rx);
		return true;
	}

	if(user->level >= LEVEL_ADMIN)
		mpfr_sprintf(query, "CALL "SP_UPDATEVARTABVARTAB"(\"%Rf\",\"%s\",\"%s\")",val,name,tab);
	else
		mpfr_sprintf(query, "CALL "SP_UPDATEVARTABVARTABUSERTAB"(\"%Rf\",\"%s\",\"%s\",%d)",val,name,tab,user->iduser);
		
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE bool readVar(const char *name, const char *tab, mpfr_t rop)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(tab) > TAB_MAXNAMELENGTH || strlen(name) > VAR_MAXNAMELENGTH)
		return false;
	
	char *tmptab = replace(tab, "\\", "\\\\");
	if(tmptab)
		tab = tmptab;
		
	char query[QUERY_LENGTH];
	userObj * user = accessCurrentUser;
	
	if(user->level >= LEVEL_ADMIN)
		sprintf(query, "CALL "SP_SELECTVALFROMVARTABVARTABBYVARNAMEANDTABNAME"(\"%s\",\"%s\")",name,tab);
	else
		sprintf(query, "CALL "SP_SELECTVALFROMVARTABVARTABUSERTABBYVARNAMEANDTABNAMEANDIDUSER"(\"%s\",\"%s\",%d)",name,tab,user->iduser);
		
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_ROW row;
	mpfr_init(rop);
	
	if((row = mysql_fetch_row(rx)) && row[0])
		mpfr_strtofr(rop, row[0], NULL, 0, MPFR_RNDN);
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE int findVarList(const char *tab)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return INVALID_ID;
	}
	
	if(strlen(tab) > TAB_MAXNAMELENGTH)
		return false;
	
	char *tmptab = replace(tab, "\\", "\\\\");
	if(tmptab)
		tab = tmptab;
		
	char query[QUERY_LENGTH];
	userObj * const user = accessCurrentUser;
	
	if(user->level >= LEVEL_ADMIN)
		sprintf(query, "CALL "SP_SELECTIDTABFROMTABBYTABNAME"(\"%s\")", tab);
	else
		sprintf(query, "CALL "SP_SELECTIDTABFROMUSERTABTABBYTABNAMEANDIDUSER"(\"%s\",%d)", tab, user->iduser);
	
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return INVALID_ID;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return INVALID_ID;
	}
	
	MYSQL_ROW row = mysql_fetch_row(rx);	
	const register int idfound = row ? atoi(row[0]) : INVALID_ID; 
	mysql_free_result(rx);	
	return idfound;
}

__MATHSUITE bool readVarList(const char *tab, dim_typ *varlistno, char var[][VAR_MAXNAMELENGTH], char val[][VAR_MAXVALLENGTH])
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(tab) > TAB_MAXNAMELENGTH)
		return false;
	
	char *tmptab = replace(tab, "\\", "\\\\");
	if(tmptab)
		tab = tmptab;
		
	char query[QUERY_LENGTH];
	userObj * const user = accessCurrentUser;
	
	if(user->level >= LEVEL_ADMIN)
		sprintf(query, "CALL "SP_SELECTVARNAMEANDVARVALFROMVARTABVARTABBYTABNAME"(\"%s\")",tab);
	else
		sprintf(query, "CALL "SP_SELECTVARNAMEANDVARVALFROMVARTABVARTABUSERTABBYTABNAMEANDIDUSER"(\"%s\",%d)",tab, user->iduser);

	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_ROW row = mysql_fetch_row(rx);
	
	for(*varlistno=0; row; ++ *varlistno)
	{
		strcpy(var[*varlistno], row[0]);
		strcpy(val[*varlistno], row[1]);
		row = mysql_fetch_row(rx);
	}
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE bool deleteVarList(const char *tab)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(tab) > TAB_MAXNAMELENGTH)
		return false;
	
	char *tmptab = replace(tab, "\\", "\\\\");
	if(tmptab)
		tab = tmptab;
		
	char query[QUERY_LENGTH];
	userObj * const user = accessCurrentUser;
	const bool permissionGranted = user->level >= LEVEL_ADMIN;
	
	if(permissionGranted)
		sprintf(query, "CALL "SP_DELETETABVARFROMTABVARTABBYTABNAME"(\"%s\")",tab);
	else
		sprintf(query, "CALL "SP_DELETETABVARFROMTABVARTABUSERTABBYTABNAMEANDIDUSER"(\"%s\",%d)",tab,user->iduser);
	
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	if(permissionGranted)
		sprintf(query, "CALL "SP_DELETEUSERTABFROMUSERTABTABBYTABNAME"(\"%s\")",tab);
	else
		sprintf(query, "CALL "SP_DELETEUSERTABFROMUSERTABTABBYTABNAMEANDIDUSER"(\"%s\",%d)",tab,user->iduser);
	
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	sprintf(query, "CALL "SP_SELECTUSERTABFROMUSERTABTABBYTABNAMEANDIDUSER"(\"%s\",%d)", tab, user->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
		
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}

	if(!mysql_fetch_row(rx))
	{
		sprintf(query, "CALL "SP_DELETETAB"(\"%s\")",tab);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_free_result(rx);
			mysql_close(con);
			return false;
		}
	}
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE bool writeMat(const char *name, const char *row)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(name) > MAT_MAXNAMELENGTH)
		return false;
	
	char *tmpname = replace(name, "\\", "\\\\");
	if(tmpname)
		name = tmpname;
		
	char query[QUERY_LENGTH];
	sprintf(query, "CALL "SP_SELECTIDMATFROMMATBYMATNAME"(\"%s\")",name);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	userObj * const user = accessCurrentUser; 
	
	if(!mysql_fetch_row(rx))
	{
		mysql_free_result(rx);
		sprintf(query, "CALL "SP_CREATEMAT"(@InsertId, \"%s\")",name);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		mysql_next_result(con);
		strcpy(query, "SELECT @InsertId");
		
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		if ((rx = mysql_store_result(con)) == NULL) 
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		register int mid = -1;
		MYSQL_ROW linsidRow = mysql_fetch_row(rx);
		
		if(linsidRow && linsidRow[0]) 
			mid = atoi(linsidRow[0]);
			
		mysql_free_result(rx);
		sprintf(query, "CALL "SP_CREATEROW"(\"%s\",%d)",row,mid);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		sprintf(query, "CALL "SP_CREATEUSERMAT"(%d,%d)",user->iduser,mid);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		return true;
	}
	
	query[0] = BLANK_CHAR;

	if(user->level >= LEVEL_ADMIN)
		sprintf(query, "CALL "SP_CREATEROWBYMATNAME"(\"%s\",\"%s\")", row, name);
	else
	{
		mysql_free_result(rx);
		sprintf(query, "CALL "SP_SELECTIDMATFROMMATUSERMATBYMATNAMEANDIDUSER"(\"%s\",%d)", name, user->iduser);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		if ((rx = mysql_store_result(con)) == NULL) 
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
			mysql_close(con);
			return false;
		}
		
		MYSQL_ROW dbrow = mysql_fetch_row(rx);
		query[0] = BLANK_CHAR;
		
		if(dbrow)
			sprintf(query, "CALL "SP_CREATEROW"(\"%s\",%d)", row, atoi(dbrow[0]));
	}
	mysql_next_result(con);
	if (query[0] != BLANK_CHAR && mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE int findMat(const char *name)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return INVALID_ID;
	}
	
	if(strlen(name) > MAT_MAXNAMELENGTH)
		return false;
	
	char *tmpname = replace(name, "\\", "\\\\");
	if(tmpname)
		name = tmpname;
		
	char query[QUERY_LENGTH];
	userObj * const user = accessCurrentUser;
	
	if(user->level >= LEVEL_ADMIN)
		sprintf(query, "CALL "SP_SELECTIDMATFROMMATBYMATNAME"(\"%s\")", name);
	else
		sprintf(query, "CALL "SP_SELECTIDMATFROMMATUSERMATBYMATNAMEANDIDUSER"(\"%s\",%d)", name, user->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return INVALID_ID;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return INVALID_ID;
	}
	
	MYSQL_ROW row = mysql_fetch_row(rx);
	const register int idfound = row ? atoi(row[0]) : INVALID_ID;
	mysql_free_result(rx);	
	return idfound;
}

__MATHSUITE bool readMat(const char *name, char mat[])
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(name) > MAT_MAXNAMELENGTH)
		return false;
	
	char *tmpname = replace(name, "\\", "\\\\");
	if(tmpname)
		name = tmpname;
		
	char query[QUERY_LENGTH];
	userObj * const user = accessCurrentUser;
	
	if(user->level >= LEVEL_ADMIN)
		sprintf(query, "CALL "SP_SELECTROWFROMMATROWBYMATNAME"(\"%s\")",name);
	else
		sprintf(query, "CALL "SP_SELECTROWFROMMATROWUSERMATBYMATNAMEANDIDUSER"(\"%s\",%d)",name,user->iduser);	
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	} 
	
	dim_typ i;
	MYSQL_ROW row;

	register bool atLeastOneRow = false;
	char token[MAX_BUFSIZ];
	
	for(row = mysql_fetch_row(rx); row; row = mysql_fetch_row(rx))
	{
		sprintf(token, "%s;", row[0]);
		strcat(mat, token);
		if(!atLeastOneRow)
			atLeastOneRow = true;
	}
	
	mysql_free_result(rx);
	return atLeastOneRow; 
}

__MATHSUITE bool deleteMat(const char *name, const bool skip)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(name) > MAT_MAXNAMELENGTH)
		return false;
	
	char *tmpname = replace(name, "\\", "\\\\");
	if(tmpname)
		name = tmpname;
		
	char query[QUERY_LENGTH];
	userObj * const user = accessCurrentUser; 
	const bool permissionGranted = user->level >= LEVEL_ADMIN;
	
	if(permissionGranted)
		sprintf(query, "CALL "SP_DELETEROWFROMMATROWBYMATNAME"(\"%s\")",name);
	else	
		sprintf(query, "CALL "SP_DELETEROWFROMMATROWUSERMATBYMATNAMEANDIDUSER"(\"%s\",%d)",name,user->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	if(skip)
		return true;
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
		
	if(permissionGranted)
		sprintf(query, "CALL "SP_DELETEUSERMATFROMMATUSERMATBYMATNAME"(\"%s\")",name);
	else	
		sprintf(query, "CALL "SP_DELETEUSERMATFROMMATUSERMATBYMATNAMEANDIDUSER"(\"%s\",%d)",name,user->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	sprintf(query, "CALL "SP_SELECTUSERMATFROMMATUSERMATBYMATNAMEANDIDUSER"(\"%s\",%d)", name, user->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
		
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	if(!mysql_fetch_row(rx))
	{
		sprintf(query, "CALL "SP_DELETEMAT"(\"%s\")",name);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_free_result(rx);
			mysql_close(con);
			return false;
		}
	}
	
	mysql_free_result(rx);
	return true; 
}

/* 
__MATHSUITE bool readUser(const register int uid, char user [static 6][USER_MAXFIELDLENGTH])
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
		
	char query[QUERY_LENGTH];
	sprintf(query, "SELECT * FROM user WHERE iduser=%d", uid); // " AND usermat.idmat=row.idmat AND usermat.iduser=%d;",name,accessCurrentUser->iduser);
	
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}

	MYSQL_ROW row;

	if((row = mysql_fetch_row(rx)))
	{
		strcpy(user[0], row[0]);
		strcpy(user[1], row[1]);
		strcpy(user[2], row[2]);
		strcpy(user[3], row[3]);
		strcpy(user[4], row[4]);
		strcpy(user[5], row[5]);
	}
	
	mysql_free_result(rx);
	return true; 
}
*/

__MATHSUITE bool readUser(const register int uid, userObj ** user)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
		
	char query[QUERY_LENGTH];
	sprintf(query, "CALL "SP_SELECTUSER"(%d)", uid); // " AND usermat.idmat=row.idmat AND usermat.iduser=%d;",name,accessCurrentUser->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}

	MYSQL_ROW row;
	
	if((row = mysql_fetch_row(rx)))
	{
		(*user) = calloc(1, sizeof(userObj));
		errMem((* user), false);
		(*user)->iduser = uid;
		strcpy((*user)->name, row[1]);
		strcpy((*user)->surname, row[2]);
		strcpy((*user)->username, row[3]);
		strcpy((*user)->email, row[4]);
		strcpy((*user)->password, row[5]);
		(*user)->level = atoi(row[6]);
	}
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE bool readUserFromUserName(const char name[static USER_MAXIDENTIFIERLENGTH], userObj ** user)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(name) > USER_MAXUSERNAMELENGTH)
		return false;

	char query[QUERY_LENGTH];
	sprintf(query, "CALL "SP_SELECTUSERBYUSERNAME"(\"%s\")", name); // " AND usermat.idmat=row.idmat AND usermat.iduser=%d;",name,accessCurrentUser->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}

	MYSQL_ROW row;
	
	if((row = mysql_fetch_row(rx)))
	{
		(*user) = calloc(1, sizeof(userObj));
		errMem((* user), false);
		(*user)->iduser = atoi(row[0]);
		strcpy((*user)->name, row[1]);
		strcpy((*user)->surname, row[2]);
		strcpy((*user)->username, name);
		strcpy((*user)->email, row[4]);
		strcpy((*user)->password, row[5]);
		(*user)->level = atoi(row[6]);
	}
	
	mysql_free_result(rx);
	return true; 
}

__MATHSUITE bool readUserFromEmail(const char email[static USER_MAXIDENTIFIERLENGTH], userObj ** user)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(strlen(email) > USER_MAXEMAILLENGTH)
		return false;
		
	char query[QUERY_LENGTH];
	sprintf(query, "CALL "SP_SELECTUSERBYEMAIL"(\"%s\")", email); // " AND usermat.idmat=row.idmat AND usermat.iduser=%d;",name,accessCurrentUser->iduser);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	MYSQL_RES *rx;
			
	if ((rx = mysql_store_result(con)) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}

	MYSQL_ROW row;
	
	if((row = mysql_fetch_row(rx)))
	{
		
		(*user) = calloc(1, sizeof(userObj));
		errMem((* user), false);
		(*user)->iduser = atoi(row[0]);
		strcpy((*user)->name, row[1]);
		strcpy((*user)->surname, row[2]);
		strcpy((*user)->username, row[3]);
		strcpy((*user)->email, email);
		strcpy((*user)->password, row[5]);
		(*user)->level = atoi(row[6]);
	}
	
	mysql_free_result(rx);
	return true; 
}


__MATHSUITE bool writeUser(userObj *restrict user, const bool mode)
{
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return false;
	}
	
	if(!user)
		return false;
	
	const bool asrt = mode == WRITEUSER_INSERT;
	char query[QUERY_LENGTH];
	
	if(strlen(user->name) > USER_MAXNAMELENGTH || strlen(user->surname) > USER_MAXSURNAMELENGTH || strlen(user->username) > USER_MAXUSERNAMELENGTH || strlen(user->email) > USER_MAXEMAILLENGTH || strlen(user->password) > USER_MAXPASSWORDLENGTH)
		return false;
	
	if(asrt == WRITEUSER_INSERT)
		sprintf(query, "CALL "SP_CREATEUSER"(\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d)", user->name, user->surname, user->username, user->email, user->password, user->level);
	else
		sprintf(query, "CALL "SP_UPDATEUSER"(%d,\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%d)", user->iduser, user->name, user->surname, user->username, user->email, user->password, user->level);

	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return false;
	}
	
	return true;
}

__MATHSUITE inline matrixObj * getMatrixFromMatrixList(const register dim_typ which_item)
{
	session * const currentSession = accessCurrentSession();
	
	if(currentSession->MLSystem.itemsno && currentSession->MLSystem.list)
	{
		mNodelist * ref = currentSession->MLSystem.list;
	
		for(register dim_typ count=0; ref && count < currentSession->MLSystem.itemsno; ++ count)
		{
			if(count == which_item)	
				return ref->matrix;
			ref = ref->next;
		}
		
	}
	
	return NULL;	
}

__MATHSUITE bool pushMatrixList(matrixObj * matrix)
{

	mNodelist * newItem = malloc(sizeof(mNodelist));
	errMem(newItem, false);
	newItem->matrix = malloc(sizeof(matrixObj));
	
	if(!newItem->matrix)
	{
		free(newItem);
		return false;
	}
	
	if(!equalMatrix(&(newItem->matrix->matrix), matrix->matrix, matrix->dim, EQUALMATRIX_REALLOC))
	{
		free(newItem->matrix);
		return false;
	}
	
	session * const currentSession = accessCurrentSession();
	const bool maximumReached = currentSession->MLSystem.itemsno == access(curLayout)->max_matrix_list_items;
	
	if(maximumReached)
		delTailMatrixList();
		
	newItem->matrix->dim[ROWS] = matrix->dim[ROWS];
	newItem->matrix->dim[COLUMNS] = matrix->dim[COLUMNS];
	
	if(currentSession->MLSystem.itemsno == 0 || !currentSession->MLSystem.list)
		if(!(currentSession->MLSystem.list = calloc(1, sizeof(mNodelist))))
		{
			matrixFree(&(newItem->matrix->matrix), newItem->matrix->dim);
			free(newItem);
			return false;
		}
		
	if(!maximumReached)	
		++ currentSession->MLSystem.itemsno;

	currentSession->MLSystem.cur_item = 0;
	newItem->next = currentSession->MLSystem.itemsno == 1 ? NULL : currentSession->MLSystem.list;
	currentSession->MLSystem.list = newItem;
	return true;
}

__MATHSUITE void delAllMatrixFromMatrixList(void)
{
	session * const currentSession = accessCurrentSession();
	if(currentSession->MLSystem.itemsno == 0 || !currentSession->MLSystem.list)	
		return;
		
	mNodelist * tmp;
	
	for(mNodelist * ref = currentSession->MLSystem.list; ref; ref = tmp)
	{
		tmp = ref->next;
		free(ref);
		-- currentSession->MLSystem.itemsno;
	}
	
	currentSession->MLSystem.cur_item = 0;
	return; 
}

__MATHSUITE void delTailMatrixList(void)
{
	for(mNodelist * ref = accessCurrentMatrixList; ref; ref = ref->next)
		if(ref->next == NULL) // tail of the list
			free(ref);
		
	return; 
}

__MATHSUITE inline session * accessCurrentSession()
{
	char socketString[MAX_SOCKETLENGTH];
	session * currentSession = NULL;
	
	itoa(access(lastSessionSocket), socketString, 10);
	
	if(hashmap_get(access(sessions), socketString, (void**) &currentSession) == MAP_MISSING)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "FATAL ERROR: Some error occurred during trying to access a Session Item");
		(void) getch();
		exit(HASHMAP_ERROR);
		return NULL;
	}
	
	return currentSession; 
}


__MATHSUITE void dbEstablishConnection(const sel_typ argc, char *argv[static 4])
{
	if(access(curLayout)->database.con)
	{
		mysql_close(access(curLayout)->database.con);
		msyprintf(COLOR_USER, "Previous DB connections has been closed.");
	}
	
	const char *server = argc ? argv[0] : access(curLayout)->database.server;
	const char *username = argc > 1 ? argv[1] : access(curLayout)->database.username;
	const char *password = argc > 2 ? argv[2] : access(curLayout)->database.password;
	const char *database = argc > 3 ? argv[3] : access(curLayout)->database.database;	
	
	if((access(curLayout)->database.con = mysql_init(NULL)) == NULL)
	{
		 printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_init()");
		 return;
	}
	
	MYSQL * const con = access(curLayout)->database.con;
	
	if (mysql_real_connect(con, server, username, password, database, 0, NULL, 0) == NULL) 
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_real_connect(): %s",mysql_error(con));
		mysql_close(con);
		return;
	}
	
	userObj * user2 = NULL;
	session * const currentSession = accessCurrentSession();
	const register int iduser = currentSession->user ? currentSession->user->iduser : MASTER_SESSION;
	
	if(readUser(iduser, &user2) && user2)
	{
		userObj * user = currentSession->user;
		
		if(user)
			free(user);
		
		if(currentSession->MLSystem.list)
		{
			delAllMatrixFromMatrixList();
			free(currentSession->MLSystem.list);
		}
		
		currentSession->user = user2;
		currentSession->lastItem[ENVS] = access(lists)[ENVS].cur_item;
		currentSession->lastItem[MATRICES] = access(lists)[MATRICES].cur_item;
		currentSession->MLSystem.list = malloc(sizeof(mNodelist));
		errMem(currentSession->MLSystem.list, VSPACE);
		// sentinel values. I know this is a bad matrix.
		currentSession->MLSystem.list->matrix = NULL;
		currentSession->MLSystem.list->next = NULL;

		
	}
	return;
}
#endif

__MATHSUITE void  viewProgramSettings(dim_typ which_layout)
{
    msprintf(COLOR_SYSTEM, "\nCurrent Settings Layout:\n\n");

	dim_typ i;
	char password[DBPASSWORD_LENGTH];
    nodelist * const tmp = listNo(which_layout, LAYOUTS);
    layoutObj * const cur_layout = ((layoutObj *)(tmp->data));

    PRINTL();

	msyprintf(COLOR_SYSTEM, "- Exit CHAR: %c;\n", cur_layout->exit_char);
    msyprintf(COLOR_SYSTEM, "- PROGRAM PRECISION: %hu;\n", cur_layout->precision);
    msyprintf(COLOR_SYSTEM, "- STABILIZER FACTOR: %hu;\n", cur_layout->stabilizer_factor);
    msyprintf(COLOR_SYSTEM, "- Min Stirling NUMBER: %hu;\n", cur_layout->min_stirling_number);
    msyprintf(COLOR_SYSTEM, "- ALGEBRA: %s;\n", suite_c.algebra_elements_names[cur_layout->algebra]);
    msyprintf(COLOR_SYSTEM, "- Outlier Constant: %.*f;\n", DEFAULT_PRECISION, cur_layout->outlier_constant);
    PRINTN();
    
    msyprintf(COLOR_SYSTEM, "MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n", mpfr_get_version (), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR, MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
	PRINTN();
    
    #ifndef __DISABLE_MULTIUSER
    
    #ifndef __DISABLE_SERVER
	    msyprintf(COLOR_SYSTEM, "- Server: %s;\n", cur_layout->server.server);
	    msyprintf(COLOR_SYSTEM, "- Port: %d;\n", cur_layout->server.port);
	    msyprintf(COLOR_SYSTEM, "- Family: IPv%s;\n", cur_layout->server.family == AF_INET ? "4":"6");
	    msyprintf(COLOR_SYSTEM, "- Resolution: %s;\n", cur_layout->server.resolution ? "Enabled":"Disabled");
	    msyprintf(COLOR_SYSTEM, "- Max Listening Sockets: %hu;\n", cur_layout->server.nListeningSockets);
	    msyprintf(COLOR_SYSTEM, "- Max Acceptable Connections: %hu;\n", cur_layout->server.nAcceptableConnections);
	    msyprintf(COLOR_SYSTEM, "- Buffer Length: %d;\n", cur_layout->server.bufferLength);
		PRINTN();
	#endif
    
    #ifndef __DISABLE_DATABASE
		msyprintf(COLOR_SYSTEM, "MySQL client version: %s\n", mysql_get_client_info());
	    msyprintf(COLOR_SYSTEM, "- DB Server: %s;\n", cur_layout->database.server);
	    msyprintf(COLOR_SYSTEM, "- DB Username: %s;\n", cur_layout->database.username);
	    strcpy(password, cur_layout->database.password);
	    asteriskize(password);
	    msyprintf(COLOR_SYSTEM, "- DB Password: %s;\n", password);
	    msyprintf(COLOR_SYSTEM, "- Database: %s;\n", cur_layout->database.database);
	    msyprintf(COLOR_SYSTEM, "- Max USERS: %hu;\n", cur_layout->database.maxUsers);
	    msyprintf(COLOR_SYSTEM, "- Default LEVEL: %hu;\n", cur_layout->database.defaultLevel);
		msyprintf(COLOR_SYSTEM, "- IDENTIFIER Mode: %s;\n", suite_c.identifier_names[cur_layout->database.identifierMode]);
		PRINTN();
	#endif
	#endif

	msyprintf(COLOR_SYSTEM, "- Matrices MAX MatrixList ITEMS: %hu;\n", cur_layout->max_matrix_list_items);
    msyprintf(COLOR_SYSTEM, "- Matrices MAX ROWS: %hu;\n", cur_layout->matrix_max_rows);
    msyprintf(COLOR_SYSTEM, "- Matrices MAX COLUMNS: %hu;\n", cur_layout->matrix_max_columns);
    msyprintf(COLOR_SYSTEM, "- Matrices BLOCK SIZE: %hu;\n", cur_layout->block_size);
    msyprintf(COLOR_SYSTEM, "- Matrices Min OSMM Dimension: %hu;\n", cur_layout->min_osmm_dim);
    msyprintf(COLOR_SYSTEM, "- Matrices Min Strassen Dimension: %hu;\n", cur_layout->min_strassen_dim);
    msyprintf(COLOR_SYSTEM, "- Max EIGENVALUES FINDER Iterations: %hu;\n", cur_layout->max_eigvalues_iterations);
    msyprintf(COLOR_SYSTEM, "- Max DSVD Iterations: %hu;\n", cur_layout->max_dsvd_iterations);
	PRINTN();

	for(i=0; i<MAX_MEMOIZABLE_FUNCTIONS; ++i)
	{
		char str[SIGN_STRING] = NULL_CHAR;
		strcpy(str, suite_c.memoizers_names[i]);
		toupper_s(str);
		msyprintf(COLOR_SYSTEM, "- Max %s Memoizable Index: %hu;\n", str, cur_layout->max_memoizable_indices[i]);
	}

	PRINTN();
    msyprintf(COLOR_SYSTEM, "- MIN Processable BASE %hu;\n", cur_layout->basecalc_minbase);
    msyprintf(COLOR_SYSTEM, "- MAX Processable BASE: %hu;\n", cur_layout->basecalc_maxbase);
    msyprintf(COLOR_SYSTEM, "- MAX PROCESSABLE Decimal2Bin NUM: %d;\n", cur_layout->max_changebase_binary_convnum);

    msyprintf(COLOR_SYSTEM, "- NEWTON DIFFTABLES Minimum Dimension: %hu;\n", cur_layout->min_newton_difftables_dim);
    msyprintf(COLOR_SYSTEM, "- NEWTON DIFFTABLES Maximum Dimension: %hu;\n", cur_layout->max_newton_difftables_dim);

    msyprintf(COLOR_SYSTEM, "- MIN Processable ROMAN NUMBER: %hu;\n", cur_layout->min_roman_number);
    msyprintf(COLOR_SYSTEM, "- MAX Processable ROMAN  NUMBER: %hu;\n", cur_layout->max_roman_number);

    msyprintf(COLOR_SYSTEM, "- Pascal's Triangle MIN ROWS: %hu;\n", cur_layout->pascal_triangle_min_rows);
    msyprintf(COLOR_SYSTEM, "- Pascal's Triangle MAX ROWS: %hu.\n", cur_layout->pascal_triangle_max_rows);

	PRINTN();
    msprintf(COLOR_USER, "\n*** Boolean Settings ***\n\n");
    PRINTL();

    for(i=0; i<MAX_BOOL_SETTINGS; ++i)
        msyprintf(COLOR_SYSTEM, "- %s: %s;\n", suite_c.bools_names[i], (cur_layout->bools & suite_c.bools[i].bmask) == suite_c.bools[i].bmask ? "Enabled":"Disabled");

    PRINTL();

    return;
}

__MATHSUITE void  setProgramSettings(dim_typ which_layout)
{
    nodelist * item_data = listNo(which_layout, LAYOUTS);
    layoutObj * const tmp = malloc(sizeof(layoutObj));
    errMem(tmp, VSPACE);
    resetProgramSettings(tmp, item_data->path);
    item_data->data = tmp;
    return;
}

__MATHSUITE inline bool   catchPause()
{
    signal(SIGINT, sigexit);
    if(access(sigresult))
    {
        printf("\nPress any key to continue Program Action.\n");
        printf("or %c to stop Execution.\n", access(curLayout)->exit_char);
        access(sigresult) = false;
        if(getch() == access(curLayout)->exit_char) return true;
    }

    return false;
}

__MATHSUITE inline void  logPrint(logObj * const which_log)
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

__MATHSUITE inline void  logWrite(FILE *fpnt, logObj * const which_log)
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
__MATHSUITE void  getVarListEx(FILE *fpnt, dim_typ which_env)
{
    void *cookie = NULL;
    EXPRTYPE vval;
    exprValList * const tmp = ((exprValList *)((exprType*)(listNo(which_env, ENVS)->data))->var_list);
    char *vname;

    if(fpnt == stdout)
        msprintf(COLOR_USER, "\nVariable list items:\n\n");

	mpfr_init(vval);
    cookie = exprValListGetNext(tmp, &vname, vval, NULL, NULL);
    // cookie = exprValListGetNext(vlist, &vname, &vval, NULL, NULL);

	if(fpnt)
	{
	    do
	    {
	
	        if(fpnt != stdout && (!strcmp(vname, DEFAULT_ENVS_ANSVALNAME)||!strcmp(vname, "global")))
	            break;
	            
	        char str[MIN_STRING];
	        mpfr_sprintf(str,"%s = %Rf;\n", vname, /*access(curLayout)->precision,*/ vval);
	        fprintf2(fpnt, str);
	        // fputs(str, fpnt);
	        cookie = exprValListGetNext(tmp/*suite.exprVars->var_list*/, &vname, vval, NULL, cookie);
	    }
	    while(cookie);
	}
	#ifndef __DISABLE_DATABASE
		else if(access(curLayout)->database.con)
		{
			char tab[MAX_PATH_LENGTH];
			strcpy(tab, listNo(which_env, ENVS)->path);
			(void) strtok(tab, EXTENSION_DOT);
			do
		    {
		
		        if((!strcmp(vname, DEFAULT_ENVS_ANSVALNAME)||!strcmp(vname, "global")))
		            break;
	
				writeVar(vname, tab, vval);
		        cookie = exprValListGetNext(tmp/*suite.exprVars->var_list*/, &vname, vval, NULL, cookie);
		    }
		    while(cookie);
		}
	#endif

	mpfr_clear(vval);
    return;
}

__MATHSUITE const bool   requires(mpfr_t val, const char *cmd_string, const char *string, const char *res_string, const bool options)
{
	mpfr_t fake_val;
    char buf[MAX_BUFSIZ];
    exprObj *e = INIT_OBJLIST;
    EXPRERRTYPE err;
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
        msprintf(COLOR_CREDITS, string);
        PRINT2N();
        (void) gets(buf);
    }

    access(exitHandle) = INVALID_EXITHANDLE;
    // EXTREMELY IMPORTANT
    // to remember to de-comment this if you shall enable override signal handler
   	// signal(SIGINT, equalSpecialMatrix);

    // To check whether User wants to exit from requiringINPUT or Not.
    if(access(sigresult) || (!strcmp(buf, MATRIXSET_COMMAND)) || (!strcmp(buf, MATRIXBACK_COMMAND)) || !buf[0])
    {
    	
        if(!strcmp(buf, MATRIXSET_COMMAND))
            access(exitHandle) = EXITHANDLE_SETCMD;

        else if(!strcmp(buf, MATRIXBACK_COMMAND))
            access(exitHandle) = EXITHANDLE_BACKCMD;

        if(!buf[0])
            access(exitHandle) = EXITHANDLE_EXIT;

        access(sigresult) = false;
        mpfr_init_set_d(val, NULL_VAL, MPFR_RNDN);
        return false;
    }

    if(buf[strlen(buf)-1] != TERMINATING_CHAR)
        strcat(buf, TERMINATING_STRING);

    PRINTN();

    // Creates expression
    err = exprCreate(&e, access(func_list), access(exprVars)->var_list, access(const_list), NULL, 0);
    if(err != EXPR_ERROR_NOERROR)
        {
        msyprintf(COLOR_SYSTEM, "Expr Creation Error.\n");
        longjmp(jumper, err);
        }

    /* Parse expr */
    err = exprParse(e, buf);
    if(err != EXPR_ERROR_NOERROR)
        {
        exprGetErrorPosition(e, &start, &end);
        msyprintf(COLOR_SYSTEM, "Parse Error (%d,%d).\n", start, end);
        longjmp(jumper, err);
        }

    /* Eval expression */
    const bool assert = (isSett(BOOLS_SHOWDIFFTIME) && (options & PARSER_SHOWDIFFTIME) == PARSER_SHOWDIFFTIME);

    if(assert)
    	gettimeofday(&tvBegin, NULL);
    	
    // mpfr_init(fake_val);
    err = exprEval(e, fake_val);

    if(err != EXPR_ERROR_NOERROR)
        msyprintf(COLOR_ERROR, "Eval Error:\n%s.\n", suite_c.exprErrorString[err == EXPR_ERROR_UNKNOWN ? 0 : err+1  ]);


    if((options & PARSER_SHOWRESULT) == PARSER_SHOWRESULT)
    {
		char buf[MAX_BUFSIZ];
        mpfr_sprintf(buf, "\n%s: %Rf.\n\n", res_string, fake_val);
        msprintf(COLOR_USER, buf);

    }

    if(assert)
    {
        PRINTL();
        msprintf(COLOR_SYSTEM, "Average Time: %.*f;\n", SHOWTIME_PRECISION, getDiffTime(&tvBegin));
        PRINTL();
    }

    if(getItemsListNo(ENVS) != STARTING_ENVSNO && access(exprVars)->e_ANS && isSett(BOOLS_SAVERESULTS) && (options & PARSER_SAVERESULT) == PARSER_SAVERESULT)
    	mpfr_set(*(access(exprVars)->e_ANS), fake_val, MPFR_RNDN);
        
	// CRITICAL POINT!!! To be implemented later...
	///refreshExprEvalLists();
    // freeExprEvalLists();

    if(getItemsListNo(ENVS) != STARTING_ENVSNO)
    {
        if(isSett(BOOLS_SHOWVARLIST) && (options & PARSER_SHOWVARLIST) == PARSER_SHOWVARLIST)
            getVarList(stdout);
        if(lazy_exec && isSett(BOOLS_ITEMSAUTOSAVING))
            saveItem(access(lists)[ENVS].cur_item, ENVS);
    }

	mpfr_init_set(val, fake_val, MPFR_RNDN);
	mpfr_clear(fake_val);
    exprFree(e);
    /// mpfr_clear(val);
    // VOLATILE ARRAY SYSTEM (VAS)
    
    for(uint64_t i=0; i<access(PRMSystem).currentIndex; ++i)
		if(accessContainer(i).pnt && accessContainer(i)._volatile)
			del(i);
			
    return false;
}

__MATHSUITE bool  insertDims(dim_typ *righe, dim_typ *colonne)
{
    const dim_typ old_dims[2] =
    {
        (*righe),
        (*colonne)
    };

	mpfr_t tmp, tmp2;

    msprintf(COLOR_CREDITS, "\nEnter ROWS and COLUMNS as expected:\n\[ROWS]\n[COLUMNS].\n");
    PRINTHOWTOBACKMESSAGE();

    while(requires(tmp, NULL, NULL_CHAR, "Inserted ROWS", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, ((*righe) = mpfr_get_ui(tmp, MPFR_RNDN))) ||
			requires(tmp2, NULL, NULL_CHAR, "Inserted COLUMNS", PARSER_SHOWRESULT) || isNullVal(tmp2) || mpfr_cmp_ui(tmp2, ((*colonne) = mpfr_get_ui(tmp2, MPFR_RNDN))) || (*colonne) < 1 || (*colonne) > access(curLayout)->matrix_max_columns || (*righe) < 1 || (*colonne) < 1 || (*righe) > access(curLayout)->matrix_max_rows || (*colonne) > access(curLayout)->matrix_max_columns)
    {
    	mpfr_clears(tmp, tmp2, NULL); 
        CLEARBUFFER();
        if(exitHandleCheck) // if(tmp3[ROWS] == NULL_VAL || tmp3[COLUMNS] == NULL_VAL)
        {
            (*righe) = old_dims[ROWS];
            (*colonne) = old_dims[COLUMNS];
            return false;
        }
        printErr(33, "Invalid [ROWS COLUMNS] format.\nYou have to insert non-negative ROWS and COLUMNS,\n\
and must be respectively less than: %hu and %hu", access(curLayout)->matrix_max_rows, access(curLayout)->matrix_max_columns);
    }
    
    mpfr_clears(tmp, tmp2, NULL); 
    return true;
}

__MATHSUITE bool  insertDim(dim_typ *dim, bool mode)
{
    const dim_typ old_dim = (*dim);
    dim_typ max_dim;

    PRINTN();

    if(mode > COLUMNS)
    {
        max_dim = MAX_RIGHE_PER_COLONNE;
        msprintf(COLOR_CREDITS, "Enter Square Matrix DIMENSION.");
    }
    else
    {
        max_dim = mode ? access(curLayout)->matrix_max_columns : access(curLayout)->matrix_max_rows;
        msprintf(COLOR_CREDITS, "Enter Matrix %s.", mode ? "COLUMNS":"ROWS");
    }
    
    mpfr_t tmp;

	PRINTN();
    PRINTHOWTOBACKMESSAGE();
	
    while(requires(tmp, NULL, NULL_CHAR, "Inserted DIMENSION is", PARSER_SHOWRESULT) || isNullVal(tmp) || mpfr_cmp_ui(tmp, ((*dim) = mpfr_get_ui(tmp, MPFR_RNDN))) || (*dim) < 1 || (*dim) > max_dim)
    {
    	mpfr_clear(tmp);
        CLEARBUFFER();
        if(exitHandleCheck)
        {
            (*dim) = old_dim;
            return false;
        }
        printErr(33, "Invalid inserted Value.\nMust be a non-negative integer less than: %hu", max_dim);
    }
    
	mpfr_clear(tmp);
    return true;
}

// La seguente funzione sarebbe stata la funzione Handler
// del segnale SIGINT. Il problema e' che dovrebbe essere
// sempre e continuamente richiamata la funzione signal
// per far sì che funzioni correttamente. Meglio evitare

/*
static inline matrixObj * const getMatrixType(int type)
{
	switch(type)
	{
		case SPECIALMATRIX_CURRENTMATRIX:
			return access(curMatrix);
		case SPECIALMATRIX_LMPMATRIX:
			return accessCurrentLmpMatrix;
		case SPECIALMATRIX_MATRIXLIST:
			return getMatrixFromMatrixList(accessCurrentSession()->MLSystem.cur_item+type); // popCurrentMatrixList;
		default:
			return NULL;
	}
	return NULL;
}
*/

static inline matrixObj * const getMatrixType(int type)
{
	
	if(type == SPECIALMATRIX_CURRENTMATRIX)
		return access(curMatrix);
	else if(type == SPECIALMATRIX_LMPMATRIX)
		return accessCurrentLmpMatrix;
		
	return getMatrixFromMatrixList(type >= MAX_SPECIALMATRICES ? type-MAX_SPECIALMATRICES : accessCurrentSession()->MLSystem.cur_item);
}

/*
__MATHSUITE void  equalSpecialMatrix(int dest, int source)
{
	static const char * matrixTypes[3] =
    {
    	"Current Matrix",
    	"Last Matrix Printed",
    	"Current MatrixList Head"
    };
    
    const register dim_typ real_dest = dest;
	const register dim_typ real_source = source;
    
    if(dest >= MAX_SPECIALMATRICES)	
    	dest = MAX_SPECIALMATRICES;
    
    if(source >= MAX_SPECIALMATRICES)
    	source = MAX_SPECIALMATRICES;
	
	access(sigresult) = true;
	
	if(dest == SPECIALMATRIX_LMPMATRIX)
	{
		printErr(ERROR_FORBIDDENOPERATION, "You cannot write the Last Matrix Printed");
		return;
	}
	
	if(source == dest)
	{
		printErr(ERROR_FORBIDDENOPERATION, "You cannot equal the same matrix");
		return;
	}
    
    const char * sourceString = matrixTypes[source];
    const char * destString = matrixTypes[dest];
    matrixObj * const sourceMat = getMatrixType(source);
    matrixObj * const destMat = getMatrixType(dest);
    
	if((source == SPECIALMATRIX_CURRENTMATRIX || dest == SPECIALMATRIX_CURRENTMATRIX) && getItemsListNo(MATRICES) == STARTING_MATNO)
		return;
		
	if(dest == SPECIALMATRIX_MATRIXLIST)
	{
		if(sourceMat && !pushMatrixList(sourceMat))
	    {
	    	printErr(ERROR_NOTENOUGHSPACE, "Some error occurred during equalling Current MatrixList Head");
	    	return;
	    }	
	}	
	else
	{
		if(sourceMat && destMat && destMat->matrix && sourceMat->matrix)
	    {
	        if(isSett(BOOLS_ITEMSAUTOSAVING) && !saveItem(access(lists)[MATRICES].cur_item, MATRICES))
	            return;
	            
	        matrixFree(&destMat->matrix, destMat->dim);
	        
	    	if(!equalMatrix(&destMat->matrix, sourceMat->matrix, sourceMat->dim, EQUALMATRIX_REALLOC))
	    		printErr(ERROR_NOTENOUGHSPACE, "Some error occurred during equalling %s to %s", destString, sourceString);
	    		
	        destMat->dim[ROWS] = sourceMat->dim[ROWS];
	        destMat->dim[COLUMNS] = sourceMat->dim[COLUMNS];
	    }
	    else
	    {
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "%s, %s or both are not available", destString, sourceString);
			return;
		}
	}
	
   	msyprintf(COLOR_SYSTEM, "\n%s has been equalled\nto %s.\n\n", destString, sourceString);
    return;
}
*/

__MATHSUITE void  equalSpecialMatrix(int dest, int source)
{
	static const char * matrixTypes[3] =
    {
    	"Current Matrix",
    	"Last Matrix Printed",
    	"Current MatrixList Matrix"
    };
    
    const register dim_typ real_dest = dest;
	const register dim_typ real_source = source;
    
    if(dest > SPECIALMATRIX_MATRIXLIST)	
    	dest = SPECIALMATRIX_MATRIXLIST;
    
    if(source > SPECIALMATRIX_MATRIXLIST)
    	source = SPECIALMATRIX_MATRIXLIST;
	
	access(sigresult) = true;
	
	if(dest == SPECIALMATRIX_LMPMATRIX)
	{
		printErr(ERROR_FORBIDDENOPERATION, "You cannot write the Last Matrix Printed");
		return;
	}
	
	if(real_dest == real_source)
	{
		printErr(ERROR_FORBIDDENOPERATION, "You cannot equal the same matrix");
		return;
	}
    
    const char * sourceString = matrixTypes[source];
    const char * destString = matrixTypes[dest];
    matrixObj * const sourceMat = getMatrixType(real_source);
    
	if((source == SPECIALMATRIX_CURRENTMATRIX || dest == SPECIALMATRIX_CURRENTMATRIX) && getItemsListNo(MATRICES) == STARTING_MATNO)
		return;
		
	if(real_dest == SPECIALMATRIX_MATRIXLIST)
	{
		printf("\npushing\n"); 
		if(sourceMat && !pushMatrixList(sourceMat))
	    {
	    	printErr(ERROR_NOTENOUGHSPACE, "Some error occurred during equalling Current MatrixList Head");
	    	return;
	    }	
	}	
	else
	{
		printf("\nequalling\n");
		matrixObj * const destMat = getMatrixType(real_dest);
		
		if(sourceMat && destMat && destMat->matrix && sourceMat->matrix)
	    {
	        if(isSett(BOOLS_ITEMSAUTOSAVING) && !saveItem(access(lists)[MATRICES].cur_item, MATRICES))
	            return;
	            
	        matrixFree(&destMat->matrix, destMat->dim);
	        
	    	if(!equalMatrix(&destMat->matrix, sourceMat->matrix, sourceMat->dim, EQUALMATRIX_REALLOC))
	    		printErr(ERROR_NOTENOUGHSPACE, "Some error occurred during equalling %s to %s", destString, sourceString);
	    		
	        destMat->dim[ROWS] = sourceMat->dim[ROWS];
	        destMat->dim[COLUMNS] = sourceMat->dim[COLUMNS];
	    }
	    else
	    {
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "%s, %s or both are not available", destString, sourceString);
			return;
		}
	}
	
	if((dest == SPECIALMATRIX_MATRIXLIST || source == SPECIALMATRIX_MATRIXLIST) && (dest != real_dest || source != real_source))
	{
		char str[MAX_STRING] = NULL_CHAR;
		
		if(dest == real_dest)
			sprintf(str, "\n%s ", destString);
		else
			sprintf(str, "\nThe %hu MatrixList Matrix ", real_dest-MAX_SPECIALMATRICES);
		
		if(source == real_source)
			sprintf(str, "%shas been equalled\nto %s.\n\n", str, sourceString);
		else
			sprintf(str, "%shas been equalled\nto the %hu MatrixList Matrix.\n\n", str, real_source-MAX_SPECIALMATRICES);
			
		msyprintf(COLOR_SYSTEM, str);		
	}
	else
   		msyprintf(COLOR_SYSTEM, "\n%s has been equalled\nto %s.\n\n", destString, sourceString);
   		
    return;
}

__MATHSUITE inline void  sigexit(int unused)
{
    access(sigresult) = true;
    return;
}

__MSSHELL_WRAPPER_ __MATHSUITE void showNewtonDifferenceTable(dim_typ n, mpfr_t x[static n], mpfr_t y[access(curLayout)->max_newton_difftables_dim][access(curLayout)->max_newton_difftables_dim], bool mode)
{
	char buf[MAX_BUFSIZ];
    dim_typ i, j;
    msprintf(COLOR_SYSTEM, "\n***********%s Difference Table ***********\n", mode ? "FORWARD" : "BACKWARD");
    
    //display Difference Table
    for(i=0;i<n;++i)
    {
        mpfr_sprintf(buf, "\t"OUTPUT_CONVERSION_FORMAT, x[i]);
        msprintf(COLOR_USER, buf);

        for(j=0; (mode ? (j<(n-i)):(j<=i)) ;++j)
        {
            mpfr_sprintf(buf, "\t"OUTPUT_CONVERSION_FORMAT, y[i][j]);
            msprintf(COLOR_USER, buf);
    	}

        PRINTN();
    }

    return;
}

__MATHSUITE bool  insertNMMatrix(mpfr_t **matrix, const register dim_typ dim[static 2])
{
    if(!matrixAlloc(matrix, dim))
        return false;

    volatile char tmp=0;
    dim_typ i, j;
    dim_typ start_col_index;


    // we must seek for backtracking...
    // BY Inserting somewhere a 'return false' statement.
    // (FINAL BACKTRACKING RESPONSE)
    for(i=start_col_index=0; tmp != -1 && i<dim[ROWS]; ++i)
        for(j=start_col_index; tmp != -1 && j<dim[COLUMNS]; ++j)
        {
            while((tmp = insertElement((*matrix), (dim_typ2){i, j}, dim[COLUMNS], false, MPFR_SET_ELEMENTS)) != 1 && tmp != -2)
                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                    if(access(curMatrix)->dim[ROWS] != dim[ROWS] || access(curMatrix)->dim[COLUMNS] != dim[COLUMNS])
                        printErr(1, "You cannot use Current Matrix because\nit doesn't have %hu Rows and %hu Columns", dim[ROWS], dim[COLUMNS]);
                    else
                    {
                        (void) equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim, EQUALMATRIX_NOREALLOC);
                        msyprintf(COLOR_SYSTEM, "\nYou're correctly using Current Matrix.\n\n");
                        break;
                    }
                else
                {
                    if(!i)
                    {
                        matrixFree(matrix, dim);
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
                    matrixFree(matrix, dim);
                    return false;
                }
            }

        }

    return true;
}

__MATHSUITE volatile char  insertElement(mpfr_t *restrict matrix, const register dim_typ dim[static 2], const register dim_typ columns, const bool square, const bool mpfr_init_elements)
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
    msprintf(COLOR_CREDITS, "\nEnter [%hu,%hu] Matrix Element.\n", dim[ROWS]+1, dim[COLUMNS]+1);
    bool once_executed = false;
    sel_typ tmp = 1;
    mpfr_t mpftmp; 
    // tmp = true;
    
    char str[MIN_STRING];

    do
    {
        CLEARBUFFER();
        sprintf(str, "[%hu,%hu] Matrix Element correctly inserted", dim[ROWS]+1, dim[COLUMNS]+1);
    	requires(mpftmp, NULL, NULL_CHAR, str, PARSER_SHOWRESULT);
    	
		if(tmp)
	    	if(mpfr_init_elements && !once_executed)
	    	{
	    		once_executed = true;
	    		mpfr_init_set(*(matrix + columns*dim[ROWS] + dim[COLUMNS]), mpftmp, MPFR_RNDN);
	    	}
			else
				mpfr_set(*(matrix + columns*dim[ROWS] + dim[COLUMNS]), mpftmp, MPFR_RNDN);
				
		if((tmp = !isNullVal(mpftmp)))
			mpfr_clear(mpftmp);
        CLEARBUFFER();

        if(access(exitHandle) == EXITHANDLE_SETCMD)
            return -1;

        if(access(exitHandle) == EXITHANDLE_BACKCMD)
            return -2;

        if(access(exitHandle) == EXITHANDLE_EXIT && (!dim[ROWS]) && (!dim[COLUMNS]))
            return -3;

        if(mpfr_zero_p(*(matrix + columns*dim[ROWS] + dim[COLUMNS])) && dcheck && INVERSE_OPS && ((__pmode__ >= ALGOPS_MATRIXMULTIPLICATION && __pmode__ <= ALGOPS_DOTPRODUCT )|| __pmode__ == ALGOPS_SCALARDIVISIONMATRIX))
           printErr(33, "You cannot enter a 0 because program is performing a Division somewhere");

    }
    while((tmp != 1 && !(dim[ROWS]) && !(dim[COLUMNS])) || (!(*(matrix + columns*dim[ROWS] + dim[COLUMNS])) && dcheck && INVERSE_OPS && ((__pmode__ >= ALGOPS_DOTPRODUCT && __pmode__ <= ALGOPS_SCALARDIVISIONMATRIX)
    ||__pmode__ == ALGOPS_SCALARDIVISIONMATRIX))||(tmp && isDomainForbidden(*(matrix + columns*dim[ROWS] + dim[COLUMNS]), INPUT)));
	
    CLEARBUFFER();
    return tmp;
}

__MATHSUITE inline volatile bool  checkBackTracking(volatile char tmp, dim_typ *colonna)
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

__MATHSUITE inline volatile  sel_typ checkBackTracking2(volatile char tmp, dim_typ *riga, dim_typ *colonna, dim_typ *start_col_index, dim_typ colonne)
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

__MATHSUITE bool   enterMatrix(mpfr_t **matrix, dim_typ *righe, dim_typ *colonne, bool square, bool view)
{
    volatile char tmp = 1;
    dim_typ start_col_index = 0;

	msprintf(COLOR_CREDITS, "\n\nEnter%s Matrix", square ? " Quad" : "");
    msprintf(COLOR_CREDITS, " by following procedure you'll see.\n");
    msprintf(COLOR_CREDITS, "And when you reach desired rows %s columns dimensions, press ENTER.\n\n",square ? "=":"and");

    (*matrix) = malloc(sizeof(mpfr_t)<<1);
    errMem((*matrix), false);

    (*righe) = 1;

    dim_typ analog_rows = 2;
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
                matrixFree(matrix, ((dim_typ2){*righe, *colonne}));
                return false;
            }
        }

        analog_columns = (*colonne)+1;

        (*matrix) = realloc((*matrix), sizeof(mpfr_t)*analog_rows*analog_columns);
        errMem((*matrix), false);

        if((tmp = insertElement((*matrix), (dim_typ2){0, (*colonne)}, *colonne, square, MPFR_INIT_ELEMENTS)) == -1 && getItemsListNo(MATRICES) != STARTING_MATNO)
        {
            if(square && access(curMatrix)->dim[ROWS] != access(curMatrix)->dim[COLUMNS])
            {
                printErr(1, "You cannot use Current Matrix because\nit isn't a Square one");
                (*colonne) --;
            }
            else
            {
            	matrixFree(matrix, ((dim_typ2){*righe, *colonne}));
                if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim, EQUALMATRIX_REALLOC))
                    return false;
                msyprintf(COLOR_SYSTEM, "\nYou are correctly using Current Matrix.\n\n");
                (*righe) = access(curMatrix)->dim[ROWS];
                (*colonne) = access(curMatrix)->dim[COLUMNS];
                return true;
            }
        }
    }

    tmp = 1;
    (*colonne) --;

    (*matrix) = realloc((*matrix), sizeof(mpfr_t)*analog_rows*(*colonne));
    errMem((*matrix), false);


    for(*righe = 1; (square ? *righe < *colonne : (tmp && __pmode__ != ALGOPS_DOTPRODUCT)); ++(*righe))
    {
        dim_typ i;
        analog_rows = (*righe)+1;

        (*matrix) = realloc((*matrix), sizeof(mpfr_t)*analog_rows*(*colonne));
        errMem((*matrix), false);

        for(i=start_col_index; tmp && i<*colonne; ++i)
        {
            /// if(tmp == -2) return false; // simple backtracking on second or major rows...
            // BACKTRACKING ON SECOND or MAJOR ROWS EXPERIMENTAL FORMULA

            while((tmp = insertElement((*matrix), (dim_typ2){*righe, i}, *colonne, square, MPFR_INIT_ELEMENTS)) != 1 && tmp != -2 && i)
                if(getItemsListNo(MATRICES) != STARTING_MATNO && tmp == -1)
                    if(square && access(curMatrix)->dim[ROWS] != access(curMatrix)->dim[COLUMNS])
                        printErr(1, "You cannot use Current Matrix because\nit isn't a Square one");
                    else
                    {
                    	matrixFree(matrix, ((dim_typ2){*righe, *colonne}));
                        if(!equalMatrix(matrix, access(curMatrix)->matrix, access(curMatrix)->dim, EQUALMATRIX_REALLOC))
                            return false;
                        msyprintf(COLOR_SYSTEM, "\nYou're correctly using Current Matrix.\n\n");
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
               		matrixFree(matrix, ((dim_typ2){*righe, *colonne}));
                    return false;
                }
            }
        }

    }

    if(!(square) && __pmode__ != ALGOPS_DOTPRODUCT)
       (*righe) --;


    if(view)
    {
        msprintf(COLOR_SYSTEM, "\nInserted [%hu X %hu] Matrix is:\n\n", *righe, *colonne);
        printMatrix(stdout, matrix, ((dim_typ2){*righe, *colonne}));
    }

    return true;
}

