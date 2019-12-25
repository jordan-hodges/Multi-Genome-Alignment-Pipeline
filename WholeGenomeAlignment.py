import os, sys, glob, json, gzip, tarfile, time, random, Logging

import BlastPrep, AlignmentPrep, ProcessAlignments
from tools import *

config = json.load(open('config.JSON'))

