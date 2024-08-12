import sys
from setuptools import setup, find_packages, Extension


setup(
    name='jobcraft',
    version='0.1',
    #packages=find_packages(['jobcraft','jobcraft.*']),
    packages=find_packages(),
    author='Martin Vondrak',
)
