#!/usr/bin/python

# Ludo - CesbBIO (2019) #
import sys
import os
import numpy as np
import datetime

#from osgeo import ogr

class infolog():
    """ Really simple and usfull log manager """
    def __init__(self):
      self.level = 'INFO'
      self.listLevels = ["DEBUG","INFO","WARNING","ERROR","FATAL"] 
      self.Lev = {"DEBUG":"DEBUG",
		  "INFO":"\033[32mINFO\033[0m",
		  "WARNING":"\033[33mWARNING\033[0m",
		  "ERROR":"\033[31mERROR\033[0m",
		  "FATAL":"\033[31mFATAL\033[0m"}
    
    def set_level(self,level):
      if(level in self.listLevels):
        self.level = level

    def get_level(self):
      return self.level

    def msg(self,message,msglevel="INFO"):
      if(msglevel in self.listLevels):
        t = datetime.datetime.now()
        text = "[%s][%s] %s"%(t.strftime('%Y-%m-%d %H:%M:%S'),self.Lev[msglevel],message)
        msgidx = self.listLevels.index(msglevel)
        idx = self.listLevels.index(self.level)
        if( msgidx >= idx):
          print(text)
      else:
        print("[ERROR] Level not defined to print: %s"%(message))


if __name__ == '__main__':

    print("test")
    log = infolog()
    log.set_level("WARNING")
    print log.get_level()
    log.msg("debug","DEBUG")
    log.msg("info")
    log.msg("warning","WARNING")
    log.msg("error","ERROR")
    log.msg("fatal","FATAL")
    
