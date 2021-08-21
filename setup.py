from setuptools import Extension, setup

module = Extension('spkmeansModule',  sources = ['spkmeansmodule.c', 'spkmeans.c'])

setup(name='spkmeansModule',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])