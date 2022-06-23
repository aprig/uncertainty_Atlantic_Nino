###################################################################################
# Load all useful modules (see /home/gula/Desktop/python_v2/Modules/)
###################################################################################


import conda
import sys
import os


import string
import cmocean ## Color map 
import numpy.ma as ma ## mask array
import matplotlib.pyplot as plt ## to make figures
import matplotlib.dates as pltd ## Dates on the plots
import scipy.signal as sc   ## for detrending the data
import numpy as np
import time as tm
import numpy.ma as ma
import xarray as xr
import scipy.stats as scstats
import scipy.interpolate as sc_interp
import matplotlib.dates as mdates
import Atools as Atools
import matplotlib as mpl
import glob
import pandas as pd
import cftime as cft
from datetime import date
from dateutil.relativedelta import relativedelta
from datetime import datetime
from datetime import timedelta
from pathlib import Path
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import gridspec
from scipy import stats
from mpl_toolkits.basemap import Basemap ## To make maps