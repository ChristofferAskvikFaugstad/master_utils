from setuptools import setup

setup(name='utils',
   version='1.1',
   author='Christoffer Askvik Faugstad',
   packages=['.'] # Directory where package is located
   )

# Installed with 
# pip install -e .
# This makes utils avaliablee everywhere
# To uninstall
# pip uninstall .
# Everything have to happen in main folder