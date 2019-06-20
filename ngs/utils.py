#!/usr/bin/env python
#--coding:utf-8--
"""
utils.py
2019-06-19: function timer added 
"""
import time, logging, sys, os
from datetime import datetime
#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4, 2.75)
mpl.rcParams["font.size"] = 10.0
import seaborn as sns
sns.set_style("white")
import pylab
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

__author__ = "CAO Yaqiang"
__date__ = "2019-05-28"
__email__ = "caoyaqiang0410@gmail.com"

def timer(func,REPEATS=5):
    """
    Timer as decorator
    REPEATS: int, repeat time to run target function. 
    """
    def wrapper(*args, **kwargs):
        start = time.time()
        for _ in range(REPEATS):
            v = func(*args, **kwargs)
        t = (time.time() - start ) / REPEATS
        print("{} time elapsed: {}".format(func.__name__ , t))
        return v
    return wrapper

def getLogger(fn=os.getcwd() + "/" + os.path.basename(__file__) + ".log"):
    """
    Return the logger with output to stdout and log file. 
    """
    #get the current time
    date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
    #set up logging, both write log info to console and log file
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(name)-6s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=fn,
        filemode='a')
    logger = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.NOTSET)
    return logger


def callSys(cmds, logger):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        logger.info(c)
        try:
            os.system(c)
        except:
            logger.error(c)


def cFlush(r):
    """
    One line flush.
    """
    sys.stdout.write("\r%s" % r)
    sys.stdout.flush()


class PET(object):
    #cA is the center of left read
    __slots__ = [
        "chromA", "chromB", "startA", "startB", "endA", "endB", "strandA",
        "strandB", "cA", "cB", "distance", "cis"
    ]

    def __init__(self, d):
        """
        d is line = line.split( "\n" )[ 0 ].split( "\t" ) from BEDPE file 
        """
        self.chromA = d[0]
        self.startA = int(d[1])
        self.endA = int(d[2])
        self.strandA = d[8]
        self.chromB = d[3]
        self.startB = int(d[4])
        self.endB = int(d[5])
        self.strandB = d[9]
        self.cis = False
        if self.chromA == self.chromB:
            self.cis = True
            self.cA = (self.startA + self.endA) / 2
            self.cB = (self.startB + self.endB) / 2
            self.distance = abs(self.cB - self.cA)
