#!/usr/bin/python
__author__ = 'p0054421'

# -*- coding: utf-8 -*-

import sqlite3 as lite
import sys

import pydicom

con = None
dbb = '/home/p0054421/MEGA/calcul/LCS_tractor/dicom/sofa_db'
dicc = '/home/p0054421/Downloads/testDicom/test.dcm'


imm = pydicom.read_file(dicc)
try:
    con = lite.connect(dbb)

    cur = con.cursor()
    cur.execute('SELECT SQLITE_VERSION()')

    data = cur.fetchone()

    print "SQLite version: %s" % data

except lite.Error, e:

    print "Error %s:" % e.args[0]
    sys.exit(1)

finally:

    if con:
        con.close()
