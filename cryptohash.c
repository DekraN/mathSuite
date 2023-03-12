#ifndef __DISABLE_CRYPTOGRAPHICHASH


#include "dutils.h"


__MSSHELL_WRAPPER_ static void  mdc2(const sel_typ argc, char ** argv);
/*
__MSSHELL_WRAPPER_ static void  md2(const sel_typ argc, char ** argv);
*/
__MSSHELL_WRAPPER_ static void  md4(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  md5(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  ripemd160(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  whirlpool(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  sha1(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  sha224(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  sha256(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  sha384(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void  sha512(const sel_typ argc, char ** argv);

sprog cryptographic_hash_prog[MAX_CRYPTOGRAPHICHASH_PROGS] =
{
	[CRYPTOGRAPHICHASH_MDC2] =
	{
		CMD_MDC2,
		NAME_MDC2,
		USAGE_MDC2,
		#ifndef __DISABLE_DATABASE
		LEVEL_MDC2,
		#endif
		mdc2,
		ARGC_MDC2,
		AUTOMATIC,
		CHILD
	},
	/*
	[CRYPTOGRAPHICHASH_MD2] =
	{
		CMD_MD2,
		NAME_MD2,
		USAGE_MD2,
		#ifndef __DISABLE_DATABASE
		LEVEL_MD2,
		#endif
		md2,
		AUTOMATIC,
		CHILD
	},
	*/
	[CRYPTOGRAPHICHASH_MD4] =
	{
		CMD_MD4,
		NAME_MD4,
		USAGE_MD4,
		#ifndef __DISABLE_DATABASE
		LEVEL_MD4,
		#endif
		md4,
		ARGC_MD4,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_MD5] =
	{
		CMD_MD5,
		NAME_MD5,
		USAGE_MD5,
		#ifndef __DISABLE_DATABASE
		LEVEL_MD5,
		#endif
		md5,
		ARGC_MD5,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_RIPEMD160] =
	{
		CMD_RIPEMD160,
		NAME_RIPEMD160,
		USAGE_RIPEMD160,
		#ifndef __DISABLE_DATABASE
		LEVEL_RIPEMD160,
		#endif
		ripemd160,
		ARGC_RIPEMD160,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_WHIRLPOOL] =
	{
		CMD_WHIRLPOOL,
		NAME_WHIRLPOOL,
		USAGE_WHIRLPOOL,
		#ifndef __DISABLE_DATABASE
		LEVEL_WHIRLPOOL,
		#endif
		whirlpool,
		ARGC_WHIRLPOOL,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_SHA1] =
	{
		CMD_SHA1,
		NAME_SHA1,
		USAGE_SHA1,
		#ifndef __DISABLE_DATABASE
		LEVEL_SHA1,
		#endif
		sha1,
		ARGC_SHA1,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_SHA224] =
	{
		CMD_SHA224,
		NAME_SHA224,
		USAGE_SHA224,
		#ifndef __DISABLE_DATABASE
		LEVEL_SHA224,
		#endif
		sha224,
		ARGC_SHA224,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_SHA256] =
	{
		CMD_SHA256,
		NAME_SHA256,
		USAGE_SHA256,
		#ifndef __DISABLE_DATABASE
		LEVEL_SHA256,
		#endif
		sha256,
		ARGC_SHA256,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_SHA384] =
	{
		CMD_SHA384,
		NAME_SHA384,
		USAGE_SHA384,
		#ifndef __DISABLE_DATABASE
		LEVEL_SHA384,
		#endif
		sha384,
		ARGC_SHA384,
		AUTOMATIC,
		CHILD
	},
	[CRYPTOGRAPHICHASH_SHA512] =
	{
		CMD_SHA512,
		NAME_SHA512,
		USAGE_SHA512,
		#ifndef __DISABLE_DATABASE
		LEVEL_SHA512,
		#endif
		sha512,
		ARGC_SHA512,
		AUTOMATIC,
		CHILD
	}
};

__MSSHELL_WRAPPER_ static void  mdc2(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want MDC-2 to hash.\n\n");
		(void) gets(buf);
	}
	
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[MDC2_HASHREPLENGTH+1];
		hexhashMDC2(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nMDC-2 hash:\n%s\n\n", hashRep);
	}
	
}

/*
__MSSHELL_WRAPPER_ static void  md2(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ];
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want MD-2 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[MD2_HASHREPLENGTH+1];
		hexhashMD2(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nMD-2 hash:\n%s\n\n", hashRep);
	}
}
*/

__MSSHELL_WRAPPER_ static void  md4(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want MD-4 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[MD4_HASHREPLENGTH+1];
		hexhashMD4(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nMD-4 hash:\n%s\n\n", hashRep);
	}
	
}

__MSSHELL_WRAPPER_ static void  md5(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want MD-5 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[MD5_HASHREPLENGTH+1];
		hexhashMD5(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nMD-5 hash:\n%s\n\n", hashRep);
	}
	
}

__MSSHELL_WRAPPER_ static void  ripemd160(const sel_typ argc, char ** argv)
{
	
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want RIPEMD-160 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[RIPEMD160_HASHREPLENGTH+1];
		hexhashRIPEMD160(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nRIPEMD-160 hash:\n%s\n\n", hashRep);
	}
}

__MSSHELL_WRAPPER_ static void  whirlpool(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want WHIRLPOOL to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[WHIRLPOOL_HASHREPLENGTH+1];
		hexhashWHIRLPOOL(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nWHIRLPOOL hash:\n%s\n\n", hashRep);
	}
	
}

__MSSHELL_WRAPPER_ static void  sha1(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want SHA-1 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[SHA1_HASHREPLENGTH+1];
		hexhashSHA1(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nSHA-1 hash:\n%s\n\n", hashRep);
	}
	
	return;
}

__MSSHELL_WRAPPER_ static void  sha224(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want SHA-224 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[SHA224_HASHREPLENGTH+1];
		hexhashSHA224(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nSHA-224 hash:\n%s\n\n", hashRep);
	}
	
	return;
}

__MSSHELL_WRAPPER_ static void  sha256(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want SHA-256 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[SHA256_HASHREPLENGTH+1];
		hexhashSHA256(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nSHA-256 hash:\n%s\n\n", hashRep);
	}
	
	return;
}

__MSSHELL_WRAPPER_ static void  sha384(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want SHA-384 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[SHA384_HASHREPLENGTH+1];
		hexhashSHA384(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nSHA-384 hash:\n%s\n\n", hashRep);
	}
	
	return;
}

__MSSHELL_WRAPPER_ static void  sha512(const sel_typ argc, char ** argv)
{
	char buf[MAX_BUFSIZ] = NULL_CHAR;
	
	if(argv && argv[0])
		strcpy(buf, argv[0]);
	else
	{
		// strcpy(buf, NULL_CHAR);
		msyprintf(COLOR_CREDITS, "Enter a string which you want SHA-512 to hash.\n\n");
		(void) gets(buf);
	}
		
	const register dim_typ len = strlen(buf);
	
	if(len < 1 || len>=MAX_BUFSIZ)
		printErr(ERROR_INPUTOUTPUT, "Invalid string. You have to insert a String whose length is between 0 and " STR(MAX_BUFSIZ) " characters");
	else
	{
		unsigned char hashRep[SHA512_HASHREPLENGTH+1];
		hexhashSHA512(buf, hashRep);	
		msyprintf(COLOR_USER, "\n\nSHA-512 hash:\n%s\n\n", hashRep);
	}	
	
	return;
}

#endif
