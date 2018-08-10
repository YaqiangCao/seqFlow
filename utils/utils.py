#!/usr/bin/env python2.7
#--coding:utf-8--
"""
utils.py
"""
import time, logging, sys, os
from datetime import datetime

__author__ = "CAO Yaqiang"
__date__ = "2018-07-12"
__email__ = "caoyaqiang0410@gmail.com"


def getlogger(fn=os.getcwd() + "/" + os.path.basename(__file__) + ".log"):
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


def call_sys(cmds):
    """
    Call systematic commands without return.
    """
    for c in cmds:
        logger.info(c)
        try:
            os.system(c)
        except:
            logger.error(c)
