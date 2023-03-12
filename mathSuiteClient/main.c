// mathSuite v10 Client Version
// 09/05/2016 Test


#include <stdio.h>
#include <stdlib.h>

#include <winsock2.h>
#include <ws2tcpip.h>

// #pragma comment(lib, "ws2_32.lib")

#define EXIT_FAILURE 1
#define MAX_BUFLEN 10240
#define PORT_STRLEN 6
#define EXIT_CMD "exit"
#define WELLKNOWN_PORT 49152
#define WELLKNOWN_SERVER INADDR_LOOPBACK
#define INADDR_LOOPBACKSTRING "127.0.0.1"
#define DEFAULT_RESOLUTION 0
#define EXIT_MESSAGE "\n\n\nThank you for using this program.\n"
#define PRINTN() printf("\n");
#define PRINTL() printf("\n________________________________________________________________________________\n\n") // enhanced formatting


static void printErr(const char *err_string)
{
	printf("\n\nERROR: %s.\n\n", err_string);
	system("PAUSE");
	exit(EXIT_FAILURE);
	return;
}

static char getAddressInfo(char *string, struct sockaddr * addr, socklen_t len)
{
	
	// no reverse lookup in getnameinfo
	int niflags = NI_NUMERICSERV | NI_NUMERICHOST;
	char IP[INET6_ADDRSTRLEN];
	char port[PORT_STRLEN];
 
	// display local address of the socket
	if(getnameinfo(addr, len, IP, INET6_ADDRSTRLEN, port, PORT_STRLEN, NI_NUMERICSERV | NI_NUMERICHOST))
	{
		printErr("getnameinfo()");
		return 0;
	}

	sprintf(string, "%s:%s", IP, port);
	return 1;
}

static inline char getFamily(const register char family)
{
	return family == 4 ? AF_INET : family == 6 ? AF_INET6 : AF_UNSPEC;
}

int main(int argc, char *argv[])
{
	printf("Welcome to mathSuite v10 Client Version.\n");
	printf("Type an existing command or write \"info\" to see Commands List.\n\n");
	
	WSADATA wsaData;
	char buf[MAX_BUFLEN];
	int sockfd;
	int n;
	
	(void) WSAStartup(MAKEWORD(2, 2), &wsaData);
	
	// SOCKADDR_STORAGE addr = { 0 };
	
	const char resolve = argc > 4 ? atoi(argv[4]) : DEFAULT_RESOLUTION;
	struct addrinfo *res, *p;
	struct addrinfo hints = { 0 };
	
	hints.ai_family = argc > 3 ? getFamily((atoi(argv[3]))) : AF_UNSPEC; //IPv4 family;
	if(!resolve)
		hints.ai_flags = AI_NUMERICHOST;
		
	hints.ai_socktype = SOCK_STREAM;
	
	char address[INET6_ADDRSTRLEN];
	char port[PORT_STRLEN];
	
	if(argc > 1)
		strcpy(address, argv[1]); 
	else
		itoa(WELLKNOWN_SERVER, address, 10);
		
	if(argc > 2)
		strcpy(port, argv[2]); 
	else
		itoa(WELLKNOWN_PORT, port, 10);
	
	if(getaddrinfo(address, port, &hints, &res))
	{
		printErr("getaddrinfo()");
		return;
	}
	
	for (p = res; p != NULL ; p = p->ai_next)
	{
		sockfd = socket(p->ai_family, p->ai_socktype, p->ai_protocol);
		
		if (sockfd == -1)
		{
			printErr("socket()");
			continue;	
		}

		// whether the socket is TCP or UDP, it will be connected to destination address
		if (connect(sockfd, p->ai_addr, p->ai_addrlen) == -1)
		{
			printErr("connect()");
			closesocket(sockfd);
			continue;
		}
		break;
	}

	/*
	int flag;
	(void) setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
	*/
	
	// const char *line = TERMINATING_STRING;

	for((void) gets(buf) ; strcmp(buf, EXIT_CMD) ; (void) gets(buf))
	{
		/*
		flag = 0;
		(void) setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
		*/
		if((n = send(sockfd, buf, strlen(buf), 0)) == -1)
			printErr("send() error");
		/*
		flag = 1;
		(void) setsockopt(sockfd, IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
		(void) send(sockfd, NULL, 0, 0);
		*/
			
		// printf("\n\nINFO: Client queued %d bytes to server.\n", (int) n);
	
		/*
		for((n = recv(sockfd, buf, MAX_BUFLEN-1, 0)), buf[n] ='\0'; strcmp(&buf[strlen(buf)-strlen(line)], line); (n = recv(sockfd, buf, MAX_BUFLEN-1, 0)), buf[n] ='\0')
		{
			PRINTL();
			printf(buf);
			if(n == -1)
			{
				printErr("recv()");
				break;
			}
			else if(!n)
			{
				printErr("Server closed connection");	
				break;
			}
			PRINTL();
		}
		*/
		
		
		if((n = recv(sockfd, buf, MAX_BUFLEN-1, 0)) == -1)
			printErr("recv() error");
		else if(!n)
			printErr("Server closed connection");
		
		
		buf[n] = '\0';	
		printf(buf);
		PRINTL();
		
		// printf("\nSERVER RESPONSE: %s\n\n", buf);		
		
	}
	
	printf(EXIT_MESSAGE);
	freeaddrinfo(res);
	closesocket(sockfd);
	WSACleanup();
	system("PAUSE");
	return EXIT_SUCCESS;
}
