#!/usr/bin/env python
#--coding:utf-8--
"""
utils.py
2019-06-19: function timer added 
2019-06-24: PET modified to more powerful
"""
import time, logging, sys, os
from datetime import datetime
#global settings
import matplotlib as mpl
mpl.use("pdf")
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["figure.figsize"] = (4 * 0.8, 2.75 * 0.8)
mpl.rcParams["font.size"] = 8.0
import seaborn as sns
sns.set_style("white")
import pylab
import brewer2mpl
colors = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"


def isTool(name):
    """
    Check if a tool is on PATh and marked as executable.

    Parameters
    ---
    name: str

    Returns
    ---
    True or False
    """
    from distutils.spawn import find_executable
    if find_executable(name) is not None:
        return True
    else:
        return False


def timer(func, REPEATS=5):
    """
    TIMER for estimate the running time of a funciton, can be used as decorator.

    Parameters
    ---
    func: funciton
    REPEATS: int, repeat time to run target function. 

    Usage
    ---
    @timer
    def run():
        for i in range(1,10000):
            pass
    """

    def wrapper(*args, **kwargs):
        start = time.time()
        for _ in range(REPEATS):
            v = func(*args, **kwargs)
        t = (time.time() - start) / REPEATS
        print("{} time elapsed: {}".format(func.__name__, t))
        return v

    return wrapper


def getLogger(fn=os.getcwd() + "/" + os.path.basename(__file__) + ".log"):
    """
    Setting up the logger systems.

    Parameters
    ----
    fn: str, file name to store the logging infor, default is genreated with time and script name

    Returns
    ----
    logging.loger 
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

    Parameters
    ---
    cmds: list, commands to run in the bash
    logger: logging.loger
    """
    for c in cmds:
        logger.info(c)
        try:
            os.system(c)
        except:
            logger.error(c)


def cFlush(r):
    """
    One line flush to show the programming process.

    Parameters
    ---
    r: str
    """
    sys.stdout.write("\r%s" % r)
    sys.stdout.flush()


class PET(object):
    """
    Paired-end tags / PETs object.
    """
    __slots__ = [
        "chromA",
        "chromB",
        "startA",
        "startB",
        "endA",
        "endB",
        "strandA",
        "strandB",
        "cA",
        "cB",
        "distance",
        "cis",
        "length",
        "mid",
        "mapq",
        "start",
        "end",
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
        self.mapq = int(d[7])
        if self.chromA == self.chromB:
            self.cis = True
            #adjust the left end and right end to make sure left is alwasy small than right
            if self.startA + self.endA > self.startB + self.endB:
                self.startA, self.startB = self.startB, self.startA
                self.endA, self.endB = self.endB, self.endA
                self.strandA, self.strandB = self.strandB, self.strandA
            self.cA = (self.startA + self.endA) / 2
            self.cB = (self.startB + self.endB) / 2
            self.distance = abs(self.cB - self.cA)  #used for long-range PETs
            self.length = self.endB - self.startA  #fragment length, used for short-range PETs
            self.mid = (self.startA + self.endB) / 2  #middle of the fragment
            self.start = self.startA
            self.end = self.endB
        else:
            self.distance = None
            self.length = None
            self.mid = None

