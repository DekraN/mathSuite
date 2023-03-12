#include "dutils.h"

#ifndef __DISABLE_MULTIUSER

#ifndef __DISABLE_SERVER
	__MSSHELL_WRAPPER_ static void serverMode(const sel_typ argc, char * argv [static 7]);
#endif
#ifndef __DISABLE_DATABASE
__MSSHELL_WRAPPER_ static void connectDatabase(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void login(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void _register(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void editUser(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void permit(const sel_typ argc, char ** argv);
__MSSHELL_WRAPPER_ static void remPermit(const sel_typ argc, char ** argv);
#endif


sprog multiuser_prog[MAX_MULTIUSER_PROGS] =
{
	#ifndef __DISABLE_SERVER
    [MULTIUSER_SERVERMODE] =
    {
    	CMD_SERVERMODE,
    	NAME_SERVERMODE, 
    	USAGE_SERVERMODE,
    	#ifndef __DISABLE_DATABASE
    	LEVEL_SERVERMODE,
    	#endif
    	serverMode,
    	ARGC_SERVERMODE,
    	AUTOMATIC,
    	CHILD
    }
    #endif
    #ifndef __DISABLE_DATABASE
    ,
    [MULTIUSER_DATABASE] = 
    {
    	CMD_DATABASE,
    	NAME_DATABASE,
    	USAGE_DATABASE,
    	LEVEL_DATABASE,
    	connectDatabase,
    	ARGC_DATABASE,
    	AUTOMATIC,
    	CHILD
    },
    [MULTIUSER_LOGIN] =
	{
		CMD_LOGIN,
		NAME_LOGIN,
		USAGE_LOGIN,
		LEVEL_LOGIN,
		login,
		ARGC_LOGIN,
		AUTOMATIC,
		CHILD
	},
	[MULTIUSER_REGISTER] =
	{
		CMD_REGISTER,
		NAME_REGISTER,
		USAGE_REGISTER,
		LEVEL_REGISTER,
		_register,
		ARGC_REGISTER,
		AUTOMATIC,
		CHILD
	},
	[MULTIUSER_EDITUSER] =
	{
		CMD_EDITUSER,
		NAME_EDITUSER,
		USAGE_EDITUSER,
		LEVEL_EDITUSER,
		editUser,
		ARGC_EDITUSER,
		AUTOMATIC,
		CHILD
	},
	[MULTIUSER_REMPERMIT] =
	{
		CMD_PERMIT,
		NAME_PERMIT,
		USAGE_PERMIT,
		LEVEL_PERMIT,
		permit,
		ARGC_PERMIT,
		AUTOMATIC,
		CHILD
	},
	[PROG_REMPERMIT] =
	{
		CMD_REMPERMIT,
		NAME_REMPERMIT,
		USAGE_REMPERMIT,
		LEVEL_REMPERMIT,
		remPermit,
		ARGC_REMPERMIT,
		AUTOMATIC,
		CHILD
	}
	#endif
};

#ifndef __DISABLE_SERVER

	/*
	struct args
	{
		dim_typ j;
		ssize_t n;
		fd_set rset;
		fd_set allset;
		char * buf;
		int *conn_set;
		mpfr_t * response;
	};

	static void * serveClient(void * argvoid)
	{
		struct args * arg = (struct args *) argvoid;
		arg->buf[arg->n] ='\0';
		msprintf(COLOR_USER, "Received message: '%s'\n", arg->buf);
		//reply to client
		(void) parse(arg->buf, arg->response);
		mpfr_sprintf(arg->buf, "%Rf", arg->response);
		if (send(arg->conn_set[arg->j], arg->buf, strlen(arg->buf), 0) == -1)
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "send()");
			// client closed connection, perform cleaning tasks
			closeUnsetSocket(&(arg->allset), &(arg->conn_set[arg->j]));
		}
		return NULL;
	}

	
	__MSSHELL_WRAPPER_ void serverMode(const sel_typ argc, char ** argv)
	{
		dim_typ i;
		#ifdef WINOS
			WSADATA wsaData;
			(void) WSAStartup(MAKEWORD(2, 2), &wsaData);
		#endif
		
		const dim_typ nListeningSockets = argc > 4 ? atoi(argv[4]) : DEFAULT_LISTENINGSOCKETS;
		const bool resolve = argc > 3 ? atoi(argv[3]) : DEFAULT_RESOLUTION;
		struct addrinfo *res, *p;
		struct addrinfo hints = { 0 };
		int lsock[nListeningSockets];
		
		// prepare select() sets
		fd_set rset, allset;
	
		FD_ZERO(&rset);
		FD_ZERO(&allset);
		// FD_SET(-1, &rset);
		
		hints.ai_family = argc > 2 ? getFamily(atoi(argv[2])) : AF_UNSPEC; // hints.ai_family = argc > 2 ? getFamily((atoi(argv[2]))) : AF_INET; //IPv4 family;
		hints.ai_flags = AI_PASSIVE;
		if(!resolve)
			hints.ai_flags |= AI_NUMERICHOST;
		hints.ai_socktype = SOCK_STREAM;
		
		char address[INET6_ADDRSTRLEN];
		char baseport[PORT_STRLEN];
		char port[PORT_STRLEN];
		
		if(argc)
			strcpy(address, argv[0]); 
		else
			itoa(INADDR_LOOPBACK, address, 10);
			
		if(argc > 1)
			strcpy(port, argv[1]); 
		else
			itoa(SERVER_PORT, port, 10);
		
		allset = rset;
		strcpy(baseport, port);
		
		register bool atLeastOneSocket = false;
		
		for(i=0; i<nListeningSockets; ++i)
		{
			itoa(atoi(baseport)+i, port, 10);
			if(getaddrinfo(address, port, &hints, &res))
			{
				printErr(ERROR_INPUTOUTPUT, "getaddrinfo()");
				continue;
			}
			
			for (p = res; p != NULL ; p = p->ai_next)
			{
				lsock[i] = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
				
				if (lsock[i] == -1)
				{
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "socket()");
					continue;	
				}
		
				if (bind(lsock[i], p->ai_addr, p->ai_addrlen) == -1)
				{
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "bind()");
					closesocket(lsock[i]);
					continue;
				}
				break;
			}
			
			if (listen(lsock[i], SERVER_BACKLOG))
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "listen()");
				closesocket(lsock[i]);
				continue;
			}
			
			char servername[INET_ADDRSTRLEN];
			
			if(!getAddressInfo(servername, p->ai_addr, p->ai_addrlen))
				continue;
				
			msprintf(COLOR_USER, "\nServer address: %s", servername);
			FD_SET(lsock[i], &rset);
			if(!atLeastOneSocket)
				atLeastOneSocket = true;
			freeaddrinfo(res);
		}
		
		if(!atLeastOneSocket)
			return;
		
		const dim_typ nAcceptableConnections = argc > 5 ? atoi(argv[5]) : MAX_ACCEPTABLECONNECTIONS;
		int conn_set[nAcceptableConnections];
		
		#pragma omp parallel for num_threads(nAcceptableConnections)
		for(i=0; i<nAcceptableConnections; ++i)
			 conn_set[i] = -1;
			 
		allset = rset;
		dim_typ j;
		int peerfd;
		mpfr_t response;
		SOCKADDR_STORAGE peer_addr;
		int len = sizeof(peer_addr);
		bool quit = false; //regola il loop infinito nel server
		bool connected = false; //regola la gestione della connessione col client
		int result;
		// const struct timeval * const inftime = calloc(1, sizeof(struct timeval));
	
		while (!quit) 
		{
			// perform select()
			rset = allset;
			const register dim_typ totalSocketsSize = nListeningSockets + nAcceptableConnections;
			ityp totalSockets[totalSocketsSize];
			#pragma omp parallel for num_threads(nListeningSockets)
			for(i=0; i<nListeningSockets; ++i)
				totalSockets[i] = lsock[i];
			#pragma omp parallel for num_threads(nAcceptableConnections)
			for(i=0; i<nAcceptableConnections; ++i)
				totalSockets[nListeningSockets+i] = conn_set[i];
				
			const register int maxFd = MAX(_MAX(nListeningSockets, (ityp*)lsock), _MAX(nAcceptableConnections, (ityp*)conn_set)); //    _MAX(totalSocketsSize, totalSockets);
			
			if ((result = select(maxFd+1, &rset, NULL, NULL, NULL)) == SOCKET_ERROR)
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "select()");
				continue;
			}
			else if (!result) // timeout on select(), should never get here
				continue;
				
			// check listening sockets
			for (i = 0; i < nListeningSockets; ++i)
				if (FD_ISSET(lsock[i], &rset))
				{
					if((peerfd = accept(lsock[i], (struct sockaddr *)&peer_addr, &len)) == -1)
						printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "accept()");
					else
					{
						char clientaddr[INET_ADDRSTRLEN] = "";
						if(!getAddressInfo(clientaddr, (struct sockaddr *)&peer_addr, sizeof(peer_addr)))
							return;
						printf("\nAccepted a new TCP connection from %s\n", clientaddr);
						
						bool is_registered = false;
						
						for (j = 0; j < nAcceptableConnections; ++j)
							if (conn_set[j] == -1)
							{
								conn_set[j] = peerfd;
								is_registered = true;
								break;
							}
							
						if(is_registered)
							FD_SET(peerfd, &allset);
						else
						{
							printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "Reached Maximum Client queue size");
							for(j=0; j<nListeningSockets; ++j)
								closesocket(lsock[i]);
						}
					}
					-- result;
				}
				
			char buf[SERVER_BUFSIZE];
			ssize_t n = 0;
		
			// check active connections
			for(j=0; result > 0; ++j)
			{
				if ((n = recv(conn_set[j], buf, SERVER_BUFSIZE-1, 0)) == -1)
				{
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "recv()");
					closeUnsetSocket(&allset, &conn_set[j]);

				}
				else if (n==0)
				{
					msyprintf(COLOR_SYSTEM, "Peer closed connection.\n");
					// client closed connection, perform cleaning tasks
					closeUnsetSocket(&allset, &conn_set[j]);
				}
				else
				{
					pthread_attr_t thread_attr;
					pthread_attr_init(&thread_attr);
					pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED);
					pthread_t thread;
					struct args arg;
					arg.j = j;
					arg.n = n;
					arg.rset = rset;
					arg.allset = allset;
					arg.conn_set = conn_set;
					arg.buf = buf;
					arg.response = &response;
					pthread_create(&thread, &thread_attr, serveClient, &arg);
					pthread_attr_destroy(&thread_attr);
				}
				-- result;
			}
		}
		
		mpfr_clear(response);
		#pragma omp parallel for num_threads(nListeningSockets)
		for(i=0; i<nListeningSockets; ++i)
			closesocket(lsock[i]);
		WSACleanup();
		return;
	}
	*/
	
	__MSSHELL_WRAPPER_ void serverMode(const sel_typ argc, char * argv[static 7])
	{	
		dim_typ i;
		register int oldSocket = access(lastSessionSocket);
		const int bufferLength = argc > 6 ? atoi(argv[6]) : access(curLayout)->server.bufferLength;
		
		access(serverBuffer) = calloc(bufferLength, sizeof(char));
		errMem(access(serverBuffer), VSPACE);	
		access(server_mode) = true;
		
		#ifndef __DISABLE_DATABASE
			userObj * currentUser = accessCurrentUser;
			const sel_typ old_level = currentUser->level;
			currentUser->level = LEVEL_GUEST;
		#endif
		
		#ifdef WINOS
			WSADATA wsaData;
			(void) WSAStartup(MAKEWORD(2, 2), &wsaData);
		#endif
		
		const dim_typ nListeningSockets = argc > 4 ? atoi(argv[4]) : access(curLayout)->server.nListeningSockets; // DEFAULT_LISTENINGSOCKETS;
		const bool resolve = argc > 3 ? atoi(argv[3]) : access(curLayout)->server.resolution; // DEFAULT_RESOLUTION;
		struct addrinfo *res, *p;
		struct addrinfo hints = { 0 };
		int lsock[nListeningSockets];
		
		// prepare select() sets
		fd_set rset, allset;
	
		FD_ZERO(&rset);
		FD_ZERO(&allset);
		// FD_SET(-1, &rset);
		
		hints.ai_family = getFamily(argc > 2 ? atoi(argv[2]) : access(curLayout)->server.family); // access(curLayout)->server.family; // AF_UNSPEC; // hints.ai_family = argc > 2 ? getFamily((atoi(argv[2]))) : AF_INET; //IPv4 family;
		hints.ai_flags = AI_PASSIVE;
		
		if(!resolve)
			hints.ai_flags |= AI_NUMERICHOST;
			
		hints.ai_socktype = SOCK_STREAM;
		
		char address[SERVER_LENGTH];
		char baseport[PORT_STRLEN];
		char port[PORT_STRLEN];
		
		strcpy(address, argc ? argv[0] : access(curLayout)->server.server); 
			
		if(argc > 1)
			strcpy(port, argv[1]); 
		else
			itoa(access(curLayout)->server.port, port, 10);
		
		allset = rset;
		strcpy(baseport, port);
		
		register bool atLeastOneSocket = false;
		
		for(i=0; i<nListeningSockets; ++i)
		{
			itoa(atoi(baseport)+i, port, 10);
			if(getaddrinfo(address, port, &hints, &res))
			{
				printErr(ERROR_INPUTOUTPUT, "getaddrinfo()");
				continue;
			}
			
			for (p = res; p != NULL ; p = p->ai_next)
			{
				lsock[i] = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
				
				if (lsock[i] == -1)
				{
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "socket()");
					continue;	
				}
		
				if (bind(lsock[i], p->ai_addr, p->ai_addrlen) == -1)
				{
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "bind()");
					closesocket(lsock[i]);
					continue;
				}
				break;
			}
			
			if (listen(lsock[i], SERVER_BACKLOG))
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "listen()");
				closesocket(lsock[i]);
				continue;
			}
			
			char servername[INET_ADDRSTRLEN];
			
			if(!getAddressInfo(servername, p->ai_addr, p->ai_addrlen))
				continue;
				
			msprintf(COLOR_USER, "\nServer address: %s", servername);
			
			if(!i)
			{	
				char buf[SERVER_LENGTH+5];
				sprintf(buf, "%s+%hu", servername, nListeningSockets-1);
				updInfo(buf);
			}
				
			FD_SET(lsock[i], &rset);
			if(!atLeastOneSocket)
				atLeastOneSocket = true;
			freeaddrinfo(res);
		}
		
		if(!atLeastOneSocket)
			return;
		
		const dim_typ nAcceptableConnections = argc > 5 ? atoi(argv[5]) : access(curLayout)->server.nAcceptableConnections; // MAX_ACCEPTABLECONNECTIONS;		
		int conn_set[nAcceptableConnections];
		
		#pragma omp parallel for num_threads(nAcceptableConnections)
		for(i=0; i<nAcceptableConnections; ++i)
			 conn_set[i] = -1;
			 
		allset = rset;
		dim_typ j;
		int peerfd;
		mpfr_t response;
		SOCKADDR_STORAGE peer_addr;
		int len = sizeof(peer_addr);
		bool quit = false; //regola il loop infinito nel server
		bool connected = false; //regola la gestione della connessione col client
		int result;
		// const struct timeval * const inftime = calloc(1, sizeof(struct timeval));
	
		while (!quit) 
		{
			// perform select()
			rset = allset;
			ityp maxlsock, maxaccptcon;
			_MAX(&maxlsock, nListeningSockets, (ityp*)lsock);
			_MAX(&maxaccptcon, nAcceptableConnections, (ityp*)conn_set);
			const bool identifierMode = access(curLayout)->database.identifierMode;
			const char * identifierString = suite_c.identifier_names[identifierMode];
			const register int maxFd = MAX(maxlsock, maxaccptcon); // _MAX(totalSocketsSize, totalSockets);
			
			if ((result = select(maxFd+1, &rset, NULL, NULL, NULL)) == SOCKET_ERROR)
			{
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "select()");
				continue;
			}
			else if (!result) // timeout on select(), should never get here
				continue;
				
			char socketString[MAX_SOCKETLENGTH];
				
			// check listening sockets
			for (i = 0; i < nListeningSockets; ++i)
				if (FD_ISSET(lsock[i], &rset))
				{
					if((peerfd = accept(lsock[i], (struct sockaddr *)&peer_addr, &len)) == -1)
						printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "accept()");
					else
					{
						bool is_registered = false;
						char clientaddr[INET_ADDRSTRLEN] = "";
						if(!getAddressInfo(clientaddr, (struct sockaddr *)&peer_addr, sizeof(peer_addr)))
							return;
						
						printf("\nAccepted a new TCP connection from %s\n\n", clientaddr);
			        	
						for (j = 0; j < nAcceptableConnections; ++j)
							if (conn_set[j] == -1)
							{
								conn_set[j] = peerfd;
								is_registered = true;
								break;
							}
							
						// supposed to be EXTREMELY FAST..
						session * const newSession = calloc(1, sizeof(session));
						newSession->user = currentUser;
						itoa(peerfd, socketString, 10);
						(void) hashmap_put(access(sessions), socketString, newSession);
							
						if(is_registered)
							FD_SET(peerfd, &allset);
						else
						{
							printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "Reached Maximum Client queue size");
							for(j=0; j<nListeningSockets; ++j)
								closesocket(lsock[i]);
						}
					}
					-- result;
				}
				
			char buf[SERVER_BUFSIZE];
			ssize_t n = 0;
		
			// check active connections
			for(j=0; result > 0; ++j)
				if ( (conn_set[j] != -1)  && FD_ISSET(conn_set[j], &rset))
				{
					if ((n = recv(conn_set[j], buf, SERVER_BUFSIZE-1, 0)) == -1)
					{
						printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "recv()");
						closeUnsetSocket(&allset, &conn_set[j]);
					}
					else if (n==0)
						// client closed connection, perform cleaning tasks
						closeUnsetSocket(&allset, &conn_set[j]);
					else 
					{
						
						/*
						int flag = 0;
						(void) setsockopt(conn_set[j], IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
						*/
						buf[n] ='\0';
						
						#ifndef __DISABLE_DATABASE
						
						char buf_copy[SERVER_BUFSIZE];
						strcpy(buf_copy, buf);
						char *cmdname[3] = { NULL_CHAR };
						cmdname[0] = strtok(buf_copy, BLANK_STRING);
						
						if(strcmp(cmdname[0], CMD_LOGIN) && strcmp(cmdname[0], CMD_REGISTER))
						{
							itoa(conn_set[j], socketString, 10);
							session * unknownSession = NULL;
							(void) hashmap_get(access(sessions), socketString, (void **) &unknownSession);
							msprintf(COLOR_USER, "Received message: '%s' from User with %s: %s\n\n", buf, identifierString, identifierMode ? unknownSession->user->email : unknownSession->user->username);
						}
						else
						{
							cmdname[1] = strtok(NULL, BLANK_STRING);
							if(cmdname[1])
							{
								if((cmdname[2] = strtok(NULL, BLANK_STRING)))
								{
									asteriskize(cmdname[2]);
									msprintf(COLOR_USER, "Received message: '%s %s %s'\n\n", cmdname[0], cmdname[1], cmdname[2]);
								}
							}
							
						}
						#else
							msprintf(COLOR_USER, "Received message: '%s'\n\n", buf);
						#endif
						//reply to client
						// (void) parse(buf, &response);
						// access(lastSessionSocket) = accessCurrentSession()->sessionSocket = conn_set[j];
						
						// printf("user is: %s\n\n", accessCurrentUser->name);
						if((access(lastSessionSocket) = conn_set[j]) != oldSocket) // if user is logged
						{
							session * const currentSession = accessCurrentSession();
							if(access(lists)[ENVS].cur_item != currentSession->lastItem[ENVS])
								passToItem(currentSession->lastItem[ENVS], ENVS, PTI_NONE);
							if(access(lists)[MATRICES].cur_item != currentSession->lastItem[MATRICES])
								passToItem(currentSession->lastItem[MATRICES], MATRICES, PTI_NONE);
															
						}

						_handleCmdLine(buf);
						// PRINTSIGNAL();
						// mpfr_sprintf(buf, "%Rf", response);
						/*
						flag = 1;
						(void) setsockopt(conn_set[j], IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
						(void) send(conn_set[j], NULL, 0, 0);
						*/
						// printf(access(serverBuffer)); // testing only
						if (send(conn_set[j], access(serverBuffer), strlen(access(serverBuffer)), 0) == -1)
						{
							printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "send()");
							// client closed connection, perform cleaning tasks
							closeUnsetSocket(&allset, &conn_set[j]);
						}
						
						strcpy(access(serverBuffer), NULL_CHAR); 
						
					}
					-- result;
				}
			
			
			access(lastSessionSocket) = INVALID_SOCKET;
		}
		
		mpfr_clear(response);
		#pragma omp parallel for num_threads(nListeningSockets)
		for(i=0; i<nListeningSockets; ++i)
			closesocket(lsock[i]);
		WSACleanup();
		currentUser->level = old_level;
		return;
	}
#endif

#ifndef __DISABLE_DATABASE
__MSSHELL_WRAPPER_ void connectDatabase(const sel_typ argc, char ** argv)
{
	dbEstablishConnection(argc, argv);
	return;
}

static inline int hashmap_cmpfunc(void* a, void *b)
{
	return ((session*)b)->user->iduser == ((userObj*)a)->iduser;
}

__MSSHELL_WRAPPER_ static void login(const sel_typ argc, char ** argv)
{
	if(!access(curLayout)->database.con)
    	dbEstablishConnection(0, NULL);

    userObj * user = NULL;
    const bool identifierMode = access(curLayout)->database.identifierMode;
	const char * identifierString = suite_c.identifier_names[identifierMode];
	bool (* const identifierFunc)(const char [static USER_MAXIDENTIFIERLENGTH],  userObj **) = identifierMode ? readUserFromEmail : readUserFromUserName;
	
	if(argc < 2 || (!identifierFunc(argv[0], &user)) || !user)
    {
    	printErr(ERROR_INPUTOUTPUT, "You have to insert an existing %s and Password", identifierString);    
		_showUsage(&multiuser_prog[MULTIUSER_LOGIN]);	
		PRINTN(); 
    	return;
    }
    
    
    unsigned char hashRep[HEXHASH_LEN+1] = NULL_CHAR;
    hexhashSHA512(argv[1], hashRep);
	
	if(!strcmp(hashRep, user->password))
	{
		register dim_typ i;
		register dim_typ candidate=access(curLayout)->database.maxUsers;
		char lastSessionSocket[5];
		itoa(access(lastSessionSocket), lastSessionSocket, 10);
		
		session * const currentSession = accessCurrentSession();
		session * tmpSession = NULL;
		
		if(access(server_mode))
		{
			if(hashmap_iterate(access(sessions), hashmap_cmpfunc, user) == MAP_OK)
			{		
				currentSession->MLSystem.list = calloc(1, sizeof(mNodelist));
				errMem(currentSession->MLSystem.list, free(user));
				
				currentSession->MLSystem.cur_item = STARTING_ITEMSNO;
				currentSession->MLSystem.list->matrix = NULL;
				currentSession->MLSystem.list->next = NULL;
				currentSession->MLSystem.lmpMatrix = calloc(1, sizeof(matrixObj));
				errMem(currentSession->MLSystem.lmpMatrix, free(user));
				currentSession->MLSystem.lmpMatrix->matrix = NULL;
				currentSession->MLSystem.lmpMatrix->dim[ROWS] = STARTING_LMP_ROWS;
				currentSession->MLSystem.lmpMatrix->dim[COLUMNS] = STARTING_LMP_COLUMNS;
			}
			else
			{
				printErr(ERROR_PERMISSIONDENIED, "User %d with %s: %s has already connected to the server", user->iduser, identifierString, argv[0]);
				free(user);
				return;
			}
			
		}
		else
		{
			userObj ** user = &accessCurrentUser;
			free(*user);
			*user = NULL;
		}
		
		currentSession->user = user;
		currentSession->lastItem[ENVS] = access(lists)[ENVS].cur_item;
		currentSession->lastItem[MATRICES] = access(lists)[MATRICES].cur_item;
		
		MYSQL * const con = access(curLayout)->database.con;
	
		if(!con)
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		else
		{
			MYSQL_ROW row;
			MYSQL_RES *rx;
			char query[QUERY_LENGTH];
			sprintf(query, "CALL "SP_SELECTDEFAULTTABS"(%d)", user->iduser);
			mysql_next_result(con);
			if (mysql_query(con, query))	
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			else
			{	
				if ((rx = mysql_store_result(con)) == NULL) 
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
				else
				{
					for(row = mysql_fetch_row(rx); row; row = mysql_fetch_row(rx))
					{
						strcat(row[0], EXTENSION_DOT""DEFAULT_VARLIST_FILE_EXTENSION); 
						listInsertProc(row[0], ENVS);
					}
					mysql_free_result(rx);	
				}
			}
				
			sprintf(query, "CALL "SP_SELECTDEFAULTMATS"(%d)", user->iduser);
			mysql_next_result(con);
			if (mysql_query(con, query))	
				printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			else
			{
				if ((rx = mysql_store_result(con)) == NULL) 
					printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
				else
				{
					for(row = mysql_fetch_row(rx); row; row = mysql_fetch_row(rx))
					{
						strcat(row[0], EXTENSION_DOT""DEFAULT_MATRIX_FILE_EXTENSION); 
						listInsertProc(row[0], MATRICES);
					}
					mysql_free_result(rx);	
				}	
			}
		}
		msprintf(COLOR_USER, "\nYou have been logged in as %s: %s.\nWelcome %s %s!\n", identifierString, argv[0], user->name, user->surname);
	}
	else
	{
		free(user);
		printErr(ERROR_PERMISSIONDENIED, "Incorrect password. Try again");
		_showUsage(&multiuser_prog[MULTIUSER_LOGIN]);	
		PRINTN(); 
	}
	return;
}

__MSSHELL_WRAPPER_ static void _register(const sel_typ argc, char ** argv)
{
	if(!access(curLayout)->database.con)
    	dbEstablishConnection(0, NULL);

    userObj * user = NULL;
    
    if(argc < 4)
    {
    	printErr(ERROR_INPUTOUTPUT, "You have to insert at least an Username, an e-mail, and two matching Passwords");    
		_showUsage(&multiuser_prog[MULTIUSER_REGISTER]);
		PRINTN(); 
		return;
    }
    
    if(strcmp(argv[2], argv[3]))
    {
    	printErr(ERROR_INPUTOUTPUT, "The two inserted passwords don't match. Try again");
    	return;
    }
    
    if(readUserFromUserName(argv[0], &user) && user) 
    {
    	free(user);
    	printErr(ERROR_PERMISSIONDENIED, "Already existing Username: %s. Try again", argv[0]);    
		_showUsage(&multiuser_prog[MULTIUSER_REGISTER]);
		PRINTN(); 
    	return;
    }
    
    if(readUserFromEmail(argv[1], &user) && user) 
    {
    	free(user);
    	printErr(ERROR_PERMISSIONDENIED, "Already existing Email: %s. Try again", argv[1]);    
		_showUsage(&multiuser_prog[MULTIUSER_REGISTER]);
		PRINTN(); 
    	return;
    }
    
    user = calloc(1, sizeof(userObj));
    errMem(user, VSPACE);
    const char * identifierString = suite_c.identifier_names[access(curLayout)->database.identifierMode];
    unsigned char hashRep[HEXHASH_LEN+1] = NULL_CHAR;
    
    hexhashSHA512(argv[2], hashRep);
    strcpy(user->name, argc > 4 ? argv[4] : NULL_CHAR);
	strcpy(user->surname, argc > 5 ? argv[5] : NULL_CHAR);
	strcpy(user->username, argv[0]);
	strcpy(user->email, argv[1]);
	strcpy(user->password, hashRep);
	user->level = access(curLayout)->database.defaultLevel;
	
	if(writeUser(user, WRITEUSER_INSERT))
	{
		session * thisSession = NULL;
		session * const currentSession = accessCurrentSession();
		char lastSessionSocket[MAX_SOCKETLENGTH];
		itoa(access(lastSessionSocket), lastSessionSocket, 10);

		if(access(server_mode))
		{
			currentSession->MLSystem.list = calloc(1, sizeof(mNodelist));
			errMem(currentSession->MLSystem.list, free(user));
			
			currentSession->MLSystem.cur_item = STARTING_ITEMSNO;
			currentSession->MLSystem.list->matrix = NULL;
			currentSession->MLSystem.list->next = NULL;
			currentSession->MLSystem.lmpMatrix = calloc(1, sizeof(matrixObj));
			errMem(currentSession->MLSystem.lmpMatrix, free(user));
			currentSession->MLSystem.lmpMatrix->matrix = NULL;
			currentSession->MLSystem.lmpMatrix->dim[ROWS] = STARTING_LMP_ROWS;
			currentSession->MLSystem.lmpMatrix->dim[COLUMNS] = STARTING_LMP_COLUMNS;
		}
		else
		{
			userObj ** user = &accessCurrentUser;
			free(*user);
			*user = NULL;
		}
		
		accessCurrentUser = user;
		msprintf(COLOR_USER, "You have been successfully registered (and logged) as %s: %s.\n", identifierString, argv[0]);
	}
	else
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "It could not be possible to add a new %s: %s", identifierString, argv[0]);
		_showUsage(&multiuser_prog[MULTIUSER_REGISTER]);	
		PRINTN(); 
	}
	return;
}

__MSSHELL_WRAPPER_ static void editUser(const sel_typ argc, char ** argv)
{
	if(!access(curLayout)->database.con)
    	dbEstablishConnection(0, NULL);

    userObj * user = NULL;
    
    static bool (* const identifierFuncs[2])(const char [static USER_MAXIDENTIFIERLENGTH],  userObj **) =
	{
		readUserFromUserName,
		readUserFromEmail
	};
	
	
	const bool identifierMode = access(curLayout)->database.identifierMode;
	const char * identifierString = suite_c.identifier_names[identifierMode];
	const char * nonIdentifierString = suite_c.identifier_names[!identifierMode];
	
    if(argc < 5 || (!identifierFuncs[identifierMode](argv[0], &user)) || !user)
    {
    	printErr(ERROR_INPUTOUTPUT, "You have to insert at least an existing %s, %s, a Password and two matching Passwords", identifierString, nonIdentifierString);    
		_showUsage(&multiuser_prog[MULTIUSER_EDITUSER]);	
		PRINTN(); 
    	return;
    }
    
    if(strcmp(identifierMode ? user->username : user->email, argv[1]))
    {
    	userObj * user2 = NULL;
		if(identifierFuncs[!identifierMode](argv[1], &user2) && user2) 
	    {
	    	free(user2);
	    	printErr(ERROR_PERMISSIONDENIED, "Already existing %s: %s. Try again", nonIdentifierString, argv[1]);    
			_showUsage(&multiuser_prog[MULTIUSER_REGISTER]);
			PRINTN(); 
	    	return;
	    }
    }
    
    if(strcmp(argv[3], argv[4]))
    {
    	free(user);
    	printErr(ERROR_INPUTOUTPUT, "The two inserted passwords don't match. Try again");
    	return;
    }
    
    unsigned char hashRep[HEXHASH_LEN+1] = NULL_CHAR;
    hexhashSHA512(argv[2], hashRep);
    
    if(strcmp(hashRep, user->password))
	{
		printErr(ERROR_PERMISSIONDENIED, "Incorrect password. Try again");
		_showUsage(&multiuser_prog[MULTIUSER_EDITUSER]);	
		free(user);
		PRINTN();
		return;
	}
    
    hexhashSHA512(argv[3], hashRep);
    strcpy(identifierMode ? user->username : user->email, argv[1]);
    strcpy(user->name, argc > 5 ? argv[5] : NULL_CHAR);
	strcpy(user->surname, argc > 6 ? argv[6] : NULL_CHAR);
	// strcpy(user->username, argv[0]);
	strcpy(user->password, hashRep);
	(void) writeUser(user, WRITEUSER_UPDATE);
	msprintf(COLOR_USER, "You have been successfully modified the %s: %s.\n", identifierString, argv[0]);
	
	/*
	else
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "It could not be possible to modify the new Username: %s", user->username);
		_showUsage(&multiuser_prog[MULTIUSER_EDITUSER]);	
		PRINTN(); 
	}
	*/
	
	free(user);
	return;
}

__MSSHELL_WRAPPER_ static void permit(const sel_typ argc, char ** argv)
{
	if(!access(curLayout)->database.con)
    	dbEstablishConnection(0, NULL);
    	
    const bool identifierMode = access(curLayout)->database.identifierMode;
    const char * identifierString = suite_c.identifier_names[identifierMode];
    bool (* const identifierFunc)(const char [static USER_MAXIDENTIFIERLENGTH],  userObj **) = identifierMode ? readUserFromEmail : readUserFromUserName;
    	
	if(argc < 2)
	{
		printErr(ERROR_INPUTOUTPUT, "You have to insert a valid Item's Name, a valid %s\nand possibly a valid Item's Type", identifierString);
		_showUsage(&multiuser_prog[MULTIUSER_REMPERMIT]);
		PRINTN(); 
		return;
	}
	
	userObj * user = NULL;
	const bool type = ENVS+(argc > 2 && atoi(argv[2]));
	const register int item_id = type ? findMat(argv[0]) : findVarList(argv[0]);
	const char *typestring = type ? "mat":"tab";
	
	if(item_id == INVALID_ID)
	{
		printErr(ERROR_NOSUCHFILEORDIR, "%s: %s has not been found on Database", suite_c.listsnames[type], argv[0]);
		_showUsage(&multiuser_prog[MULTIUSER_REMPERMIT]);
		PRINTN(); 
		return;
	}	
	
	if((!identifierFunc(argv[1], &user)) || !user)
	{
		printErr(ERROR_NOSUCHFILEORDIR, "User with %s: %s has not been found on Database", identifierString, argv[1]);
		_showUsage(&multiuser_prog[MULTIUSER_REMPERMIT]);
		PRINTN(); 
		return;
	}
	
	const register int iduser = user->iduser;
	free(user);
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return;
	}
	
	char query[QUERY_LENGTH];
	userObj * currentUser = accessCurrentUser;
	
	if(user->level < LEVEL_ADMIN)
	{
		if(type)
			sprintf(query, "CALL "SP_SELECTIDUSERFROMUSERMATBYIDUSERANDIDMAT"(%d,%d)", currentUser->iduser, item_id);
		else
			sprintf(query, "CALL "SP_SELECTIDUSERFROMUSERTABBYIDUSERANDIDTAB"(%d,%d)", currentUser->iduser, item_id);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return;
		}
		
		MYSQL_RES *rx;
				
		if ((rx = mysql_store_result(con)) == NULL) 
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
			mysql_close(con);
			return;
		}
		
		if(!mysql_fetch_row(rx))
		{
			printErr(ERROR_PERMISSIONDENIED, "You haven't permission to do this");
			mysql_free_result(rx);
			return;
		}
		mysql_free_result(rx);
	}

	if(type)
		sprintf(query, "CALL "SP_CREATEUSERTAB"(%d,%d)", iduser, item_id);
	else
		sprintf(query, "CALL "SP_CREATEUSERMAT"(%d,%d)", iduser, item_id);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return;
	}
	
	msprintf(COLOR_USER, "\nIt has been given the permission to edit the %s:\n%s\nto the User with %s: %s.\n", suite_c.listsnames[type], argv[0], identifierString, argv[1]);
	return;
}

__MSSHELL_WRAPPER_ static void remPermit(const sel_typ argc, char ** argv)
{
	if(!access(curLayout)->database.con)
    	dbEstablishConnection(0, NULL);
	
	const bool identifierMode = access(curLayout)->database.identifierMode;
    const char * identifierString = suite_c.identifier_names[identifierMode];
    bool (* const identifierFunc)(const char [static USER_MAXIDENTIFIERLENGTH],  userObj **) = identifierMode ? readUserFromEmail : readUserFromUserName;
    
	if(argc < 2)
	{
		printErr(ERROR_INPUTOUTPUT, "You have to insert a valid Item's Name, a valid %s\nand possibly a valid Item's Type", identifierString);
		_showUsage(&multiuser_prog[MULTIUSER_REMPERMIT]);
		PRINTN(); 
		return;
	}
	
	userObj * user = NULL;
	const bool type = ENVS+(argc > 2 && atoi(argv[2]));
	const register int item_id = type ? findMat(argv[0]) : findVarList(argv[0]);
	const char *typestring = type ? "mat":"tab";
	
	if(item_id == INVALID_ID)
	{
		printErr(ERROR_NOSUCHFILEORDIR, "%s: %s has not been found on Database", suite_c.listsnames[type], argv[0]);
		_showUsage(&multiuser_prog[MULTIUSER_REMPERMIT]);
		PRINTN(); 
		return;
	}	
	
	if((!identifierFunc(argv[1], &user)) || !user)
	{
		printErr(ERROR_NOSUCHFILEORDIR, "User with %s: %s has not been found on Database", identifierString, argv[1]);
		_showUsage(&multiuser_prog[MULTIUSER_REMPERMIT]);
		PRINTN(); 
		return;
	}
	
	const register int iduser = user->iduser;
	free(user);
	MYSQL * const con = access(curLayout)->database.con;
	
	if(!con)
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "DB is not connected!");
		return;
	}
	
	char query[QUERY_LENGTH];
	userObj * const currentUser = accessCurrentUser; 
	
	if(currentUser->level < LEVEL_ADMIN)
	{
		if(type)
			sprintf(query, "CALL "SP_SELECTIDUSERFROMUSERMATBYIDUSERANDIDMAT"(%d,%d)", currentUser->iduser, item_id);
		else
			sprintf(query, "CALL "SP_SELECTIDUSERFROMUSERTABBYIDUSERANDIDTAB"(%d,%d)", currentUser->iduser, item_id);
		mysql_next_result(con);
		if (mysql_query(con, query))
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
			mysql_close(con);
			return;
		}
		
		MYSQL_RES *rx;
				
		if ((rx = mysql_store_result(con)) == NULL) 
		{
			printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_store_result(): %s",mysql_error(con));
			mysql_close(con);
			return;
		}
		
		if(!mysql_fetch_row(rx))
		{
			printErr(ERROR_PERMISSIONDENIED, "You haven't permission to do this");
			mysql_free_result(rx);
			return;
		}
		mysql_free_result(rx);
	}

	if(type)
		sprintf(query, "CALL "SP_DELETEUSERTAB"(%d,%d)", iduser, item_id);
	else
		sprintf(query, "CALL "SP_DELETEUSERMAT"(%d,%d)", iduser, item_id);
	mysql_next_result(con);
	if (mysql_query(con, query))
	{
		printErr(ERROR_RESOURCETEMPORARYUNAVAILABLE, "mysql_query(): %s",mysql_error(con));
		mysql_close(con);
		return;
	}
	
	msprintf(COLOR_USER, "\nIt has been removed the permission to edit the %s:\n%s\nto the User with %s: %s.\n", suite_c.listsnames[type], argv[0], identifierString, argv[1]);
	return;
}

#endif
#endif
