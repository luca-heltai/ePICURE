from distutils.core import setup, Extension
from setuptools import setup, Extension
import numpy as np


module1 = Extension('_Bas',
                    sources = ['utilities/_Bas.c'])

setup (name = 'NURBS',
       version = '1.0',
       description = 'NURBS',
       include_dirs = [np.get_include()],
       ext_modules = [module1])
