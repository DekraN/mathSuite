#ifndef USE_MYSQL

#include "db_compat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_SP_ARGS 10

/* Schema SQLite incorporato per auto-inizializzazione del file .db */
static const char *INIT_SCHEMA_SQL =
    "PRAGMA foreign_keys = ON;"
    "CREATE TABLE IF NOT EXISTS mat (idmat INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT NOT NULL UNIQUE);"
    "CREATE TABLE IF NOT EXISTS row (idrow INTEGER PRIMARY KEY AUTOINCREMENT, row TEXT, idmat INTEGER NOT NULL, FOREIGN KEY (idmat) REFERENCES mat(idmat) ON DELETE CASCADE);"
    "CREATE TABLE IF NOT EXISTS tab (idtab INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT NOT NULL UNIQUE);"
    "CREATE TABLE IF NOT EXISTS var (idvar INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT NOT NULL, val TEXT NOT NULL);"
    "CREATE TABLE IF NOT EXISTS tabvar (idtab INTEGER NOT NULL, idvar INTEGER NOT NULL, PRIMARY KEY (idtab, idvar), FOREIGN KEY (idtab) REFERENCES tab(idtab) ON DELETE CASCADE, FOREIGN KEY (idvar) REFERENCES var(idvar) ON DELETE CASCADE);"
    "CREATE TABLE IF NOT EXISTS user (iduser INTEGER PRIMARY KEY AUTOINCREMENT, name TEXT NOT NULL, surname TEXT NOT NULL, username TEXT NOT NULL UNIQUE, email TEXT NOT NULL UNIQUE, password TEXT NOT NULL, level INTEGER NOT NULL DEFAULT 1);"
    "CREATE TABLE IF NOT EXISTS usermat (iduser INTEGER NOT NULL, idmat INTEGER NOT NULL, \"default\" INTEGER NOT NULL DEFAULT 0, PRIMARY KEY (iduser, idmat), FOREIGN KEY (iduser) REFERENCES user(iduser) ON DELETE CASCADE, FOREIGN KEY (idmat) REFERENCES mat(idmat) ON DELETE CASCADE);"
    "CREATE TABLE IF NOT EXISTS usertab (iduser INTEGER NOT NULL, idtab INTEGER NOT NULL, \"default\" INTEGER NOT NULL DEFAULT 0, PRIMARY KEY (iduser, idtab), FOREIGN KEY (iduser) REFERENCES user(iduser) ON DELETE CASCADE, FOREIGN KEY (idtab) REFERENCES tab(idtab) ON DELETE CASCADE);"
    "CREATE VIEW IF NOT EXISTS mat_row AS SELECT row.row AS row_row, mat.name AS mat_name, row.idrow AS row_idrow FROM row JOIN mat ON row.idmat = mat.idmat;"
    "CREATE VIEW IF NOT EXISTS mat_usermat AS SELECT mat.idmat AS mat_idmat, mat.name AS mat_name, usermat.iduser AS usermat_iduser, usermat.idmat AS usermat_idmat, usermat.\"default\" AS usermat_default FROM mat JOIN usermat ON mat.idmat = usermat.idmat;"
    "CREATE VIEW IF NOT EXISTS mat_row_usermat AS SELECT mat.idmat AS mat_idmat, mat.name AS mat_name, row.idrow AS row_idrow, row.idmat AS row_idmat, row.row AS row_row, usermat.iduser AS usermat_iduser, usermat.idmat AS usermat_idmat, usermat.\"default\" AS usermat_default FROM mat JOIN usermat ON usermat.idmat = mat.idmat JOIN row ON mat.idmat = row.idmat;"
    "CREATE VIEW IF NOT EXISTS tabvar_tab AS SELECT tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab, tab.idtab AS tab_idtab, tab.name AS tab_name FROM tabvar JOIN tab ON tab.idtab = tabvar.idtab;"
    "CREATE VIEW IF NOT EXISTS tabvar_tab_usertab AS SELECT tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab, tab.idtab AS tab_idtab, tab.name AS tab_name, usertab.idtab AS usertab_idtab, usertab.iduser AS usertab_iduser, usertab.\"default\" AS usertab_default FROM tabvar JOIN tab ON tab.idtab = tabvar.idtab JOIN usertab ON usertab.idtab = tab.idtab;"
    "CREATE VIEW IF NOT EXISTS usertab_tab AS SELECT usertab.iduser AS usertab_iduser, usertab.idtab AS usertab_idtab, usertab.\"default\" AS usertab_default, tab.name AS tab_name FROM usertab JOIN tab ON tab.idtab = usertab.idtab;"
    "CREATE VIEW IF NOT EXISTS var_tabvar AS SELECT var.idvar AS var_idvar, tabvar.idvar AS tabvar_idvar, var.val AS var_val, var.name AS var_name FROM var JOIN tabvar ON var.idvar = tabvar.idvar JOIN tab ON tab.idtab = tabvar.idtab;"
    "CREATE VIEW IF NOT EXISTS var_tabvar_tab AS SELECT var.idvar AS var_idvar, tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab, var.val AS var_val, var.name AS var_name, tab.idtab AS tab_idtab, tab.name AS tab_name FROM var JOIN tabvar ON var.idvar = tabvar.idvar JOIN tab ON tab.idtab = tabvar.idtab;"
    "CREATE VIEW IF NOT EXISTS var_tabvar_tab_usertab AS SELECT var.idvar AS var_idvar, tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab, var.val AS var_val, var.name AS var_name, tab.idtab AS tab_idtab, tab.name AS tab_name, usertab.iduser AS usertab_iduser, usertab.\"default\" AS usertab_default FROM var JOIN tabvar ON var.idvar = tabvar.idvar JOIN tab ON tab.idtab = tabvar.idtab JOIN usertab ON usertab.idtab = tab.idtab;"
    "INSERT OR IGNORE INTO user (iduser, name, surname, username, email, password, level) VALUES "
    "(0, 'MASTER', 'SESSION', 'root', 'root@mathsuite.org', 'c916114f019f8b185d533b3a4fffdb8bebf9c81c8ceb6d832388adac0cadfb162f1422c8b599514479cf8695416b59b747caa061f9424cc10b10397b9787f338', 2),"
    "(1, 'MASTER', 'SESSION', 'developer', 'developer@mathsuite.org', '71bc535a84a94e2c95a7a0d7de5b9c2b2feb6f688ebad343543b8efc79d76c6116fe0b0dce79ef69c15d1d374c030ee3f5a798518567baf5a3615782d82fd198', 3),"
    "(2, 'MASTER', 'SESSION', 'owner', 'owner@mathsuite.org', '592bf0c5506db50b5241500fb57e538e21c3165c991691069e33b07882ad9f69abeb30e57e3746c07bff58daf6178aa247f49feba651d980c987f4a92440f0a4', 4);";

static char *compat_strdup(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char *p = (char *)malloc(len);
    if (p) memcpy(p, s, len);
    return p;
}

static void set_error(MYSQL *mysql, const char *msg) {
    if (!mysql) return;
    snprintf(mysql->last_error, sizeof(mysql->last_error), "%s", msg ? msg : "Unknown error");
}

static MYSQL_RES *alloc_res(int num_cols) {
    MYSQL_RES *res = (MYSQL_RES *)calloc(1, sizeof(MYSQL_RES));
    if (res) res->num_cols = num_cols;
    return res;
}

static void append_row(MYSQL_RES *res, char **cols) {
    if (!res) return;
    char ***tmp = (char ***)realloc(res->rows, sizeof(char **) * (res->num_rows + 1));
    if (!tmp) return;
    res->rows = tmp;
    res->rows[res->num_rows++] = cols;
}

/* Esegue una query SQLite parametrizzata e salva eventuale resultset in mysql->pending_res */
static int exec_sql(MYSQL *mysql, const char *sql, int argc, const char **argv, int expect_rows, int capture_insert_id) {
    sqlite3_stmt *stmt = NULL;
    if (mysql->pending_res) {
        compat_mysql_free_result(mysql->pending_res);
        mysql->pending_res = NULL;
    }

    if (sqlite3_prepare_v2(mysql->db, sql, -1, &stmt, NULL) != SQLITE_OK) {
        set_error(mysql, sqlite3_errmsg(mysql->db));
        return 1;
    }

    for (int i = 0; i < argc; i++) {
        if (argv[i] == NULL) {
            sqlite3_bind_null(stmt, i + 1);
        } else {
            sqlite3_bind_text(stmt, i + 1, argv[i], -1, SQLITE_TRANSIENT);
        }
    }

    if (expect_rows) {
        int cols = sqlite3_column_count(stmt);
        MYSQL_RES *res = alloc_res(cols);
        int rc;
        while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
            char **row = (char **)calloc(cols, sizeof(char *));
            for (int c = 0; c < cols; c++) {
                const unsigned char *val = sqlite3_column_text(stmt, c);
                row[c] = val ? compat_strdup((const char *)val) : NULL;
            }
            append_row(res, row);
        }
        if (rc != SQLITE_DONE) {
            set_error(mysql, sqlite3_errmsg(mysql->db));
            compat_mysql_free_result(res);
            sqlite3_finalize(stmt);
            return 1;
        }
        mysql->pending_res = res;
    } else {
        int rc = sqlite3_step(stmt);
        if (rc != SQLITE_DONE && rc != SQLITE_ROW) {
            set_error(mysql, sqlite3_errmsg(mysql->db));
            sqlite3_finalize(stmt);
            return 1;
        }
        if (capture_insert_id) {
            mysql->last_insert_id = sqlite3_last_insert_rowid(mysql->db);
        }
    }

    sqlite3_finalize(stmt);
    mysql->last_error[0] = '\0';
    return 0;
}

/* Parser degli argomenti di CALL spName(arg1, arg2, ...) con gestione apici singoli e escape */
static int parse_call(const char *query, char *sp_name, size_t sp_size, char args[MAX_SP_ARGS][1024], int *arg_count) {
    const char *p = query;
    while (isspace((unsigned char)*p)) p++;
    if (strncasecmp(p, "CALL", 4) != 0) return 0;
    p += 4;
    while (isspace((unsigned char)*p)) p++;

    /* Estrae nome procedure (eventualmente con backtick) */
    if (*p == '`') p++;
    size_t i = 0;
    while (*p && *p != '(' && *p != '`' && !isspace((unsigned char)*p) && i < sp_size - 1) {
        sp_name[i++] = *p++;
    }
    sp_name[i] = '\0';
    while (*p && *p != '(') p++;
    if (*p != '(') return 0;
    p++; /* salta '(' */

    *arg_count = 0;
    while (*p) {
        while (isspace((unsigned char)*p)) p++;
        if (*p == ')') break;
        if (*arg_count >= MAX_SP_ARGS) break;

        char *dst = args[*arg_count];
        size_t d = 0;

        if (*p == '\'' || *p == '"') {
            char quote = *p++;
            while (*p && d < 1023) {
                if (*p == '\\' && *(p + 1)) {
                    p++;
                    dst[d++] = *p++;
                } else if (*p == quote && *(p + 1) == quote) {
                    dst[d++] = quote;
                    p += 2;
                } else if (*p == quote) {
                    p++;
                    break;
                } else {
                    dst[d++] = *p++;
                }
            }
            dst[d] = '\0';
        } else {
            while (*p && *p != ',' && *p != ')' && d < 1023) {
                dst[d++] = *p++;
            }
            while (d > 0 && isspace((unsigned char)dst[d - 1])) d--;
            dst[d] = '\0';
        }
        (*arg_count)++;

        while (isspace((unsigned char)*p)) p++;
        if (*p == ',') p++;
    }
    return 1;
}

/* Dispatcher delle Stored Procedure di mathsuite.sql */
static int dispatch_sp(MYSQL *mysql, const char *sp, char args[MAX_SP_ARGS][1024], int argc) {
    const char *p[MAX_SP_ARGS];
    for (int i = 0; i < argc; i++) p[i] = args[i];

    /* --- MATRICES --- */
    if (strcmp(sp, "createMat") == 0 && argc >= 2) {
        /* args[0] = @InsertId, args[1] = _matName */
        return exec_sql(mysql, "INSERT INTO mat(name) VALUES(?);", 1, &p[1], 0, 1);
    }
    if (strcmp(sp, "createRow") == 0 && argc >= 2) {
        return exec_sql(mysql, "INSERT INTO row(row, idmat) VALUES(?, ?);", 2, p, 0, 0);
    }
    if (strcmp(sp, "createRowByMatName") == 0 && argc >= 2) {
        return exec_sql(mysql, "INSERT INTO row(row, idmat) VALUES(?, (SELECT idmat FROM mat WHERE name = ?));", 2, p, 0, 0);
    }
    if (strcmp(sp, "createUserMat") == 0 && argc >= 2) {
        return exec_sql(mysql, "INSERT INTO usermat(iduser, idmat) VALUES(?, ?);", 2, p, 0, 0);
    }
    if (strcmp(sp, "deleteMat") == 0 && argc >= 1) {
        return exec_sql(mysql, "DELETE FROM mat WHERE name = ?;", 1, p, 0, 0);
    }
    if (strcmp(sp, "deleteRowFromMatRowByMatName") == 0 && argc >= 1) {
        return exec_sql(mysql, "DELETE FROM row WHERE idmat IN (SELECT idmat FROM mat WHERE name = ?);", 1, p, 0, 0);
    }
    if (strcmp(sp, "deleteRowFromMatRowUserMatByMatNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql,
            "DELETE FROM row WHERE idmat IN ("
            "SELECT mat.idmat FROM mat JOIN usermat ON mat.idmat = usermat.idmat WHERE mat.name = ? AND usermat.iduser = ?);",
            2, p, 0, 0);
    }
    if (strcmp(sp, "deleteUserMat") == 0 && argc >= 2) {
        return exec_sql(mysql, "DELETE FROM usermat WHERE iduser = ? AND idmat = ?;", 2, p, 0, 0);
    }
    if (strcmp(sp, "deleteUserMatFromMatUserMatByMatName") == 0 && argc >= 1) {
        return exec_sql(mysql, "DELETE FROM usermat WHERE idmat IN (SELECT idmat FROM mat WHERE name = ?);", 1, p, 0, 0);
    }
    if (strcmp(sp, "deleteUserMatFromMatUserMatByMatNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "DELETE FROM usermat WHERE idmat IN (SELECT idmat FROM mat WHERE name = ?) AND iduser = ?;", 2, p, 0, 0);
    }
    if (strcmp(sp, "selectDefaultMats") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT mat_name FROM mat_usermat WHERE usermat_iduser = ? AND usermat_default = 1;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectIdMatFromMatByMatName") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT idmat FROM mat WHERE name = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectIdMatFromMatUserMatByMatNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT mat_idmat FROM mat_usermat WHERE mat_name = ? AND usermat_iduser = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectIdUserFromUserMatByIdUserAndIdMat") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT iduser FROM usermat WHERE iduser = ? AND idmat = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectRowFromMatRowByMatName") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT row_row FROM mat_row WHERE mat_name = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectRowFromMatRowUserMatByMatNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT row_row FROM mat_row_usermat WHERE mat_name = ? AND usermat_iduser = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectUserMatFromMatUserMatByMatNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT usermat.* FROM usermat JOIN mat ON usermat.idmat = mat.idmat WHERE mat.name = ? AND usermat.iduser = ?;", 2, p, 1, 0);
    }

    /* --- VARLISTS (TAB & VAR) --- */
    if (strcmp(sp, "createTab") == 0 && argc >= 2) {
        /* args[0] = @InsertId, args[1] = _name */
        return exec_sql(mysql, "INSERT INTO tab(name) VALUES(?);", 1, &p[1], 0, 1);
    }
    if (strcmp(sp, "createVar") == 0 && argc >= 3) {
        /* args[0] = @InsertId, args[1] = _name, args[2] = _val */
        return exec_sql(mysql, "INSERT INTO var(name, val) VALUES(?, ?);", 2, &p[1], 0, 1);
    }
    if (strcmp(sp, "createTabVar") == 0 && argc >= 2) {
        return exec_sql(mysql, "INSERT INTO tabvar(idtab, idvar) VALUES((SELECT idtab FROM tab WHERE name = ?), ?);", 2, p, 0, 0);
    }
    if (strcmp(sp, "createUserTab") == 0 && argc >= 2) {
        return exec_sql(mysql, "INSERT INTO usertab(iduser, idtab) VALUES(?, ?);", 2, p, 0, 0);
    }
    if (strcmp(sp, "deleteTab") == 0 && argc >= 1) {
        return exec_sql(mysql, "DELETE FROM tab WHERE name = ?;", 1, p, 0, 0);
    }
    if (strcmp(sp, "deleteTabVarFromTabVarTabByTabName") == 0 && argc >= 1) {
        return exec_sql(mysql, "DELETE FROM tabvar WHERE idtab IN (SELECT idtab FROM tab WHERE name = ?);", 1, p, 0, 0);
    }
    if (strcmp(sp, "deleteTabVarFromTabVarTabUserTabByTabNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql,
            "DELETE FROM tabvar WHERE idtab IN ("
            "SELECT tab.idtab FROM tab JOIN usertab ON usertab.idtab = tab.idtab WHERE tab.name = ? AND usertab.iduser = ?);",
            2, p, 0, 0);
    }
    if (strcmp(sp, "deleteUserTab") == 0 && argc >= 2) {
        return exec_sql(mysql, "DELETE FROM usertab WHERE iduser = ? AND idtab = ?;", 2, p, 0, 0);
    }
    if (strcmp(sp, "deleteUserTabFromUserTabTabByTabName") == 0 && argc >= 1) {
        return exec_sql(mysql, "DELETE FROM usertab WHERE idtab IN (SELECT idtab FROM tab WHERE name = ?);", 1, p, 0, 0);
    }
    if (strcmp(sp, "deleteUserTabFromUserTabTabByTabNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "DELETE FROM usertab WHERE idtab IN (SELECT idtab FROM tab WHERE name = ?) AND iduser = ?;", 2, p, 0, 0);
    }
    if (strcmp(sp, "selectDefaultTabs") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT tab_name FROM usertab_tab WHERE usertab_iduser = ? AND usertab_default = 1;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectIdTabFromTabByTabName") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT idtab FROM tab WHERE name = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectIdTabFromUserTabTabByTabNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT usertab_idtab FROM usertab_tab WHERE tab_name = ? AND usertab_iduser = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectIdUserFromUserTabByIdUserAndIdTab") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT iduser FROM usertab WHERE iduser = ? AND idtab = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectNameFromTabByTabName") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT name FROM tab WHERE name = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectUserTabFromUserTabTabByTabNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT usertab.* FROM usertab JOIN tab ON tab.idtab = usertab.idtab WHERE tab.name = ? AND usertab.iduser = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectValFromVarTabVarTabByVarNameAndTabName") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT var_val FROM var_tabvar_tab WHERE var_name = ? AND tab_name = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "selectValFromVarTabVarTabUserTabByVarNameAndTabNameAndIdUser") == 0 && argc >= 3) {
        return exec_sql(mysql, "SELECT var_val FROM var_tabvar_tab_usertab WHERE var_name = ? AND tab_name = ? AND usertab_iduser = ?;", 3, p, 1, 0);
    }
    if (strcmp(sp, "selectVarNameAndVarValFromVarTabVarTabByTabName") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT var_name, var_val FROM var_tabvar_tab WHERE tab_name = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectVarNameAndVarValFromVarTabVarTabUserTabByTabNameAndIdUser") == 0 && argc >= 2) {
        return exec_sql(mysql, "SELECT var_name, var_val FROM var_tabvar_tab_usertab WHERE tab_name = ? AND usertab_iduser = ?;", 2, p, 1, 0);
    }
    if (strcmp(sp, "updateVarTabVarTab") == 0 && argc >= 3) {
        /* args: _val, _varName, _tabName */
        return exec_sql(mysql,
            "UPDATE var SET val = ? WHERE name = ? AND idvar IN ("
            "SELECT tabvar.idvar FROM tabvar JOIN tab ON tab.idtab = tabvar.idtab WHERE tab.name = ?);",
            3, p, 0, 0);
    }
    if (strcmp(sp, "updateVarTabVarTabUserTab") == 0 && argc >= 4) {
        /* args: _val, _varName, _tabName, _idUser */
        return exec_sql(mysql,
            "UPDATE var SET val = ? WHERE name = ? AND idvar IN ("
            "SELECT tabvar.idvar FROM tabvar JOIN tab ON tab.idtab = tabvar.idtab "
            "JOIN usertab ON usertab.idtab = tab.idtab WHERE tab.name = ? AND usertab.iduser = ?);",
            4, p, 0, 0);
    }

    /* --- USERS --- */
    if (strcmp(sp, "createUser") == 0 && argc >= 6) {
        return exec_sql(mysql, "INSERT INTO user(name, surname, username, email, password, level) VALUES(?, ?, ?, ?, ?, ?);", 6, p, 0, 1);
    }
    if (strcmp(sp, "selectUser") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT * FROM user WHERE iduser = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectUserByEmail") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT * FROM user WHERE email = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "selectUserByUserName") == 0 && argc >= 1) {
        return exec_sql(mysql, "SELECT * FROM user WHERE username = ?;", 1, p, 1, 0);
    }
    if (strcmp(sp, "updateUser") == 0 && argc >= 7) {
        /* args: _idUser, _name, _surname, _username, _email, _password, _level */
        const char *reordered[7] = { p[1], p[2], p[3], p[4], p[5], p[6], p[0] };
        return exec_sql(mysql, "UPDATE user SET name = ?, surname = ?, username = ?, email = ?, password = ?, level = ? WHERE iduser = ?;", 7, reordered, 0, 0);
    }

    char err[256];
    snprintf(err, sizeof(err), "Compat layer: Stored Procedure sconosciuta '%s'", sp);
    set_error(mysql, err);
    return 1;
}

/* --- Implementazione API MySQL --- */

MYSQL *compat_mysql_init(MYSQL *mysql) {
    if (!mysql) {
        mysql = (MYSQL *)calloc(1, sizeof(MYSQL));
    } else {
        memset(mysql, 0, sizeof(MYSQL));
    }
    return mysql;
}

MYSQL *compat_mysql_real_connect(MYSQL *mysql, const char *host, const char *user,
                                 const char *passwd, const char *db, unsigned int port,
                                 const char *unix_socket, unsigned long clientflag) {
    (void)host; (void)user; (void)passwd; (void)port; (void)unix_socket; (void)clientflag;
    if (!mysql) return NULL;

    char db_path[512];
    if (db && db[0] != '\0') {
        if (strstr(db, ".db") || strstr(db, ".sqlite") || strchr(db, '/') || strchr(db, '\\')) {
            snprintf(db_path, sizeof(db_path), "%s", db);
        } else {
            snprintf(db_path, sizeof(db_path), "%s.db", db);
        }
    } else {
        snprintf(db_path, sizeof(db_path), "mathsuite.db");
    }

    if (sqlite3_open(db_path, &mysql->db) != SQLITE_OK) {
        set_error(mysql, sqlite3_errmsg(mysql->db));
        if (mysql->db) {
            sqlite3_close(mysql->db);
            mysql->db = NULL;
        }
        return NULL;
    }

    char *errmsg = NULL;
    if (sqlite3_exec(mysql->db, INIT_SCHEMA_SQL, NULL, NULL, &errmsg) != SQLITE_OK) {
        set_error(mysql, errmsg);
        sqlite3_free(errmsg);
        sqlite3_close(mysql->db);
        mysql->db = NULL;
        return NULL;
    }

    return mysql;
}

void compat_mysql_close(MYSQL *mysql) {
    if (!mysql) return;
    if (mysql->pending_res) {
        compat_mysql_free_result(mysql->pending_res);
        mysql->pending_res = NULL;
    }
    if (mysql->db) {
        sqlite3_close(mysql->db);
        mysql->db = NULL;
    }
}

const char *compat_mysql_error(MYSQL *mysql) {
    return mysql ? mysql->last_error : "Invalid MYSQL handle";
}

int compat_mysql_query(MYSQL *mysql, const char *q) {
    if (!mysql || !mysql->db || !q) return 1;

    /* 1. Intercetta "SELECT @InsertId" usato subito dopo createVar / createMat / createTab */
    const char *trimmed = q;
    while (isspace((unsigned char)*trimmed)) trimmed++;

    if (strncasecmp(trimmed, "SELECT @InsertId", 16) == 0) {
        if (mysql->pending_res) {
            compat_mysql_free_result(mysql->pending_res);
            mysql->pending_res = NULL;
        }
        MYSQL_RES *res = alloc_res(1);
        char **row = (char **)calloc(1, sizeof(char *));
        char id_buf[32];
        snprintf(id_buf, sizeof(id_buf), "%lld", (long long)mysql->last_insert_id);
        row[0] = compat_strdup(id_buf);
        append_row(res, row);
        mysql->pending_res = res;
        mysql->last_error[0] = '\0';
        return 0;
    }

    /* 2. Intercetta chiamate a Stored Procedure "CALL ..." */
    char sp_name[128];
    char args[MAX_SP_ARGS][1024];
    int argc = 0;
    if (parse_call(trimmed, sp_name, sizeof(sp_name), args, &argc)) {
        return dispatch_sp(mysql, sp_name, args, argc);
    }

    /* 3. Fallback per eventuali query SQL standard (SELECT vs INSERT/UPDATE/DELETE) */
    int is_select = (strncasecmp(trimmed, "SELECT", 6) == 0 || strncasecmp(trimmed, "PRAGMA", 6) == 0);
    return exec_sql(mysql, trimmed, 0, NULL, is_select, !is_select);
}

MYSQL_RES *compat_mysql_store_result(MYSQL *mysql) {
    if (!mysql) return NULL;
    MYSQL_RES *res = mysql->pending_res;
    mysql->pending_res = NULL;
    return res;
}

MYSQL_ROW compat_mysql_fetch_row(MYSQL_RES *result) {
    if (!result || result->current_row >= result->num_rows) return NULL;
    return result->rows[result->current_row++];
}

void compat_mysql_free_result(MYSQL_RES *result) {
    if (!result) return;
    if (result->rows) {
        for (int r = 0; r < result->num_rows; r++) {
            if (result->rows[r]) {
                for (int c = 0; c < result->num_cols; c++) {
                    free(result->rows[r][c]);
                }
                free(result->rows[r]);
            }
        }
        free(result->rows);
    }
    free(result);
}

int compat_mysql_next_result(MYSQL *mysql) {
    (void)mysql;
    /* In MySQL restituisce >0 se ci sono altri resultset, -1 quando non ce ne sono più */
    return -1;
}

#include <openssl/opensslv.h>
#if OPENSSL_VERSION_NUMBER >= 0x30000000L
#include <openssl/provider.h>

__attribute__((constructor))
static void compat_load_openssl_providers(void) {
    OSSL_PROVIDER_load(NULL, "default");
    OSSL_PROVIDER_load(NULL, "legacy");
}
#endif

#endif /* !USE_MYSQL */
