import sqlite3
import threading


class SQLiteBase:

    def __init__(self, infile):
        self.dbfilepath = infile
        self.conn_by_threadid = dict()

    def getConnCursor(self):

        curThreadID = threading.get_ident()

        if curThreadID in self.conn_by_threadid:
            conn = self.conn_by_threadid[curThreadID]
            return conn, conn.cursor()

        print("Creating new db connection for thread", curThreadID)

        conn = sqlite3.connect(self.dbfilepath)
        c = conn.cursor()

        self.conn_by_threadid[curThreadID] = conn

        return conn, c