from setuptools import setup

setup(
   name='lssd',
   version='0.1.1',
   description='A Python package for data analysis in luminescence sediment provenance.',
   author='João Bueno',
   author_email='jbueno@usp.br',
   packages=['lssd'],  #same as name
   install_requires=['numpy', 'matplotlib', 'pandas'], #external packages as dependencies
)
