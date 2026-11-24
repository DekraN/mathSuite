PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS mat (
    idmat INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS row (
    idrow INTEGER PRIMARY KEY AUTOINCREMENT,
    row TEXT,
    idmat INTEGER NOT NULL,
    FOREIGN KEY (idmat) REFERENCES mat(idmat) ON DELETE CASCADE ON UPDATE NO ACTION
);
CREATE INDEX IF NOT EXISTS idmat_fk_idx ON row(idmat);

CREATE TABLE IF NOT EXISTS tab (
    idtab INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS var (
    idvar INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    val TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS tabvar (
    idtab INTEGER NOT NULL,
    idvar INTEGER NOT NULL,
    PRIMARY KEY (idtab, idvar),
    FOREIGN KEY (idtab) REFERENCES tab(idtab) ON DELETE CASCADE ON UPDATE NO ACTION,
    FOREIGN KEY (idvar) REFERENCES var(idvar) ON DELETE CASCADE ON UPDATE NO ACTION
);
CREATE INDEX IF NOT EXISTS idvar_fk2_idx ON tabvar(idvar);

CREATE TABLE IF NOT EXISTS user (
    iduser INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    surname TEXT NOT NULL,
    username TEXT NOT NULL UNIQUE,
    email TEXT NOT NULL UNIQUE,
    password TEXT NOT NULL,
    level INTEGER NOT NULL DEFAULT 1
);

CREATE TABLE IF NOT EXISTS usermat (
    iduser INTEGER NOT NULL,
    idmat INTEGER NOT NULL,
    "default" INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (iduser, idmat),
    FOREIGN KEY (iduser) REFERENCES user(iduser) ON DELETE CASCADE ON UPDATE NO ACTION,
    FOREIGN KEY (idmat) REFERENCES mat(idmat) ON DELETE CASCADE ON UPDATE NO ACTION
);
CREATE INDEX IF NOT EXISTS usermat_fk2_idx ON usermat(idmat);

CREATE TABLE IF NOT EXISTS usertab (
    iduser INTEGER NOT NULL,
    idtab INTEGER NOT NULL,
    "default" INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY (iduser, idtab),
    FOREIGN KEY (iduser) REFERENCES user(iduser) ON DELETE CASCADE ON UPDATE NO ACTION,
    FOREIGN KEY (idtab) REFERENCES tab(idtab) ON DELETE CASCADE ON UPDATE NO ACTION
);
CREATE INDEX IF NOT EXISTS usertab_fk2_idx ON usertab(idtab);

-- Viste (identiche a MySQL per compatibilità)
CREATE VIEW IF NOT EXISTS mat_row AS
SELECT row.row AS row_row, mat.name AS mat_name, row.idrow AS row_idrow
FROM row JOIN mat ON row.idmat = mat.idmat;

CREATE VIEW IF NOT EXISTS mat_usermat AS
SELECT mat.idmat AS mat_idmat, mat.name AS mat_name,
       usermat.iduser AS usermat_iduser, usermat.idmat AS usermat_idmat, usermat."default" AS usermat_default
FROM mat JOIN usermat ON mat.idmat = usermat.idmat;

CREATE VIEW IF NOT EXISTS mat_row_usermat AS
SELECT mat.idmat AS mat_idmat, mat.name AS mat_name,
       row.idrow AS row_idrow, row.idmat AS row_idmat, row.row AS row_row,
       usermat.iduser AS usermat_iduser, usermat.idmat AS usermat_idmat, usermat."default" AS usermat_default
FROM mat
JOIN usermat ON usermat.idmat = mat.idmat
JOIN row ON mat.idmat = row.idmat;

CREATE VIEW IF NOT EXISTS tabvar_tab AS
SELECT tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab,
       tab.idtab AS tab_idtab, tab.name AS tab_name
FROM tabvar JOIN tab ON tab.idtab = tabvar.idtab;

CREATE VIEW IF NOT EXISTS tabvar_tab_usertab AS
SELECT tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab,
       tab.idtab AS tab_idtab, tab.name AS tab_name,
       usertab.idtab AS usertab_idtab, usertab.iduser AS usertab_iduser, usertab."default" AS usertab_default
FROM tabvar
JOIN tab ON tab.idtab = tabvar.idtab
JOIN usertab ON usertab.idtab = tab.idtab;

CREATE VIEW IF NOT EXISTS usertab_tab AS
SELECT usertab.iduser AS usertab_iduser, usertab.idtab AS usertab_idtab,
       usertab."default" AS usertab_default, tab.name AS tab_name
FROM usertab JOIN tab ON tab.idtab = usertab.idtab;

CREATE VIEW IF NOT EXISTS var_tabvar AS
SELECT var.idvar AS var_idvar, tabvar.idvar AS tabvar_idvar,
       var.val AS var_val, var.name AS var_name
FROM var
JOIN tabvar ON var.idvar = tabvar.idvar
JOIN tab ON tab.idtab = tabvar.idtab;

CREATE VIEW IF NOT EXISTS var_tabvar_tab AS
SELECT var.idvar AS var_idvar, tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab,
       var.val AS var_val, var.name AS var_name,
       tab.idtab AS tab_idtab, tab.name AS tab_name
FROM var
JOIN tabvar ON var.idvar = tabvar.idvar
JOIN tab ON tab.idtab = tabvar.idtab;

CREATE VIEW IF NOT EXISTS var_tabvar_tab_usertab AS
SELECT var.idvar AS var_idvar, tabvar.idvar AS tabvar_idvar, tabvar.idtab AS tabvar_idtab,
       var.val AS var_val, var.name AS var_name,
       tab.idtab AS tab_idtab, tab.name AS tab_name,
       usertab.iduser AS usertab_iduser, usertab."default" AS usertab_default
FROM var
JOIN tabvar ON var.idvar = tabvar.idvar
JOIN tab ON tab.idtab = tabvar.idtab
JOIN usertab ON usertab.idtab = tab.idtab;

-- Seed utenti di sistema (iduser 0, 1, 2 come nel dump)
INSERT OR IGNORE INTO user (iduser, name, surname, username, email, password, level) VALUES
(0, 'MASTER', 'SESSION', 'root', 'root@mathsuite.org', 'c916114f019f8b185d533b3a4fffdb8bebf9c81c8ceb6d832388adac0cadfb162f1422c8b599514479cf8695416b59b747caa061f9424cc10b10397b9787f338', 2),
(1, 'MASTER', 'SESSION', 'developer', 'developer@mathsuite.org', '71bc535a84a94e2c95a7a0d7de5b9c2b2feb6f688ebad343543b8efc79d76c6116fe0b0dce79ef69c15d1d374c030ee3f5a798518567baf5a3615782d82fd198', 3),
(2, 'MASTER', 'SESSION', 'owner', 'owner@mathsuite.org', '592bf0c5506db50b5241500fb57e538e21c3165c991691069e33b07882ad9f69abeb30e57e3746c07bff58daf6178aa247f49feba651d980c987f4a92440f0a4', 4);
