#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Background and base for making database queries.
"""

import os
import time
import urllib

import logging
logger = logging.getLogger(__name__)

try:
    import psycopg2
except ImportError:
    psycopg2 = None


class Database(object):
    """
    Database mixin provides a standard process for connecting to a database
    """
    CREDENTIALS_FILENAME = '~/.config/tsig/tsig-dbinfo'
    CREDENTIALS_URL = 'http://tessellate.mit.edu/tsig/tsig-dbinfo'

    def __init__(self, dbinfo_loc=None,
                 dbhost=None, dbname=None, dbport=5432,
                 dbuser=None, dbpass=None,
                 dbtable=None):

        # get database credentials and other connection parameters
        self.dbinfo = dict()
        if dbhost is not None:
            self.dbinfo['dbhost'] = dbhost
        if dbname is not None:
            self.dbinfo['dbname'] = dbname
        if dbport is not None:
            self.dbinfo['dbport'] = dbport
        if dbuser is not None:
            self.dbinfo['dbuser'] = dbuser
        if dbpass is not None:
            self.dbinfo['dbpass'] = dbpass
        if dbtable is not None:
            self.dbinfo['dbtable'] = dbtable

        # override database parameters with info from file or url
        self.dbinfo.update(self.read_dbinfo(dbinfo_loc))

        # FIXME: eliminate this hack - it belongs in a query results object
        self.column_names = []

    @classmethod
    def read_dbinfo(cls, path=None):
        """
        Try to read from the specified location, which could be a file or a
        ULR.  If the specified location does not work, fail fast.

        If nothing is specified, then fallback to first trying the default
        file location, then try the web server and cache result.
        """
        default_filename = cls.CREDENTIALS_FILENAME
        default_url = cls.CREDENTIALS_URL
        if not path or '://' in path:
            url = path or default_url
            filename = os.path.abspath(os.path.expanduser(default_filename))
            try:
                if not os.path.isfile(filename):
                    try:
                        os.makedirs(os.path.dirname(filename))
                    except os.error:
                        pass
                    testfile = urllib.URLopener()
                    testfile.retrieve(url, filename)
            except Exception as err:
                logging.error("Download of credentials failed: "
                              "%s (%s)" % (str(err), url))
        else:
            filename = os.path.abspath(os.path.expanduser(path))

        # now read the parameters
        dbinfo = dict()
        try:
            with open(filename, 'r') as f:
                for line in f:
                    name, value = line.split('=')
                    dbinfo[name.strip()] = value.strip()
        except IOError:
            raise Exception("No %s database credentials at %s" % (
                cls.__name__, filename))
        return dbinfo

    @staticmethod
    def check_dbinfo(dbinfo):
        missing = []
        for k in ['dbhost', 'dbname', 'dbuser', 'dbpass', 'dbtable']:
            if k not in dbinfo:
                missing.append(k)
        if missing:
            raise Exception("Missing database parameters: %s" %
                            ','.join(missing))

    def _do_query(self, sqlcmd, *args, **kw):
        if psycopg2 is None:
            logger.error("PostgreSQL python module psycopg2 is required.")
            return None

        logger.debug("query: %s" % sqlcmd)
        t0 = time.time()

        conn = psycopg2.connect(
            host=self.dbinfo['dbhost'], database=self.dbinfo['dbname'],
            port=self.dbinfo['dbport'],
            user=self.dbinfo['dbuser'], password=self.dbinfo['dbpass'],
            connect_timeout=5)

        cur = conn.cursor()
        t1 = time.time()
        logger.debug("query start")
        cur.execute(sqlcmd, *args, **kw)
        t2 = time.time()
        logger.debug("query complete")
        if 'select' in sqlcmd.lower():
            result = cur.fetchall()
            # FIXME: redesign this
            self.column_names = [desc[0] for desc in cur.description]
        else:
            result = conn.commit()

        t3 = time.time()
        logger.debug("data retrieved")
        cur.close()
        conn.close()
        t4 = time.time()
        logger.info("query: setup=%.3f execute=%.3f fetch=%.3f close=%.3f"
                    % (t1 - t0, t2 - t1, t3 - t2, t4 - t3))
        return result

