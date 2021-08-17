from setuptools import Extension, setup

module = Extension('spkmeans',  sources = ['spkmeansmodule.c'])

setup(name='spkmeans',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])