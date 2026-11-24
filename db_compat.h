#ifndef DB_COMPAT_H
#define DB_COMPAT_H

#ifdef USE_MYSQL
    #include <mysql/mysql.h>
#else
    #include <sqlite3.h>
    #include <stdbool.h>
    #include <stddef.h>

    typedef char **MYSQL_ROW;

    typedef struct {
        char ***rows;       /* rows[row_idx][col_idx] */
        int num_rows;
        int num_cols;
        int current_row;
    } MYSQL_RES;

    typedef struct {
        sqlite3 *db;
        char last_error[512];
        sqlite3_int64 last_insert_id;
        MYSQL_RES *pending_res;
    } MYSQL;

    MYSQL *compat_mysql_init(MYSQL *mysql);
    MYSQL *compat_mysql_real_connect(MYSQL *mysql, const char *host, const char *user,
                                     const char *passwd, const char *db, unsigned int port,
                                     const char *unix_socket, unsigned long clientflag);
    void compat_mysql_close(MYSQL *mysql);
    const char *compat_mysql_error(MYSQL *mysql);
    int compat_mysql_query(MYSQL *mysql, const char *q);
    MYSQL_RES *compat_mysql_store_result(MYSQL *mysql);
    MYSQL_ROW compat_mysql_fetch_row(MYSQL_RES *result);
    void compat_mysql_free_result(MYSQL_RES *result);
    int compat_mysql_next_result(MYSQL *mysql);

    #define mysql_init          compat_mysql_init
    #define mysql_real_connect  compat_mysql_real_connect
    #define mysql_close         compat_mysql_close
    #define mysql_error         compat_mysql_error
    #define mysql_query         compat_mysql_query
    #define mysql_store_result  compat_mysql_store_result
    #define mysql_fetch_row     compat_mysql_fetch_row
    #define mysql_free_result   compat_mysql_free_result
    #define mysql_next_result   compat_mysql_next_result
#endif

/* Fallback per sistemi RHEL 8 / OpenSSL 1.1.1+ / 3.x dove MDC-2 è disabilitato */
/*
#include <openssl/opensslconf.h>
#include <openssl/evp.h>

#if defined(OPENSSL_NO_MDC2) || !defined(MDC2_DIGEST_LENGTH)
    #ifndef MDC2_DIGEST_LENGTH
        #define MDC2_DIGEST_LENGTH 16
    #endif

    static inline unsigned char *MDC2(const unsigned char *d, size_t n, unsigned char *md) {
        static unsigned char static_md[MDC2_DIGEST_LENGTH];
        unsigned char full_digest[EVP_MAX_MD_SIZE];
        unsigned int out_len = 0;

        if (!md) md = static_md;

        EVP_MD_CTX *ctx = EVP_MD_CTX_new();
        if (ctx) {
            if (EVP_DigestInit_ex(ctx, EVP_sha256(), NULL) == 1 &&
                EVP_DigestUpdate(ctx, d, n) == 1 &&
                EVP_DigestFinal_ex(ctx, full_digest, &out_len) == 1) {
                memcpy(md, full_digest, MDC2_DIGEST_LENGTH);
            }
            EVP_MD_CTX_free(ctx);
        }
        return md;
    }
#else
    #include <openssl/mdc2.h>
#endif
*/

#endif /* DB_COMPAT_H */
