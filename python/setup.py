import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "mwaspbiepy"

list_files=[]
list_files.append('../src/em_muller_trans_wrap.f90')

FLIBS = os.getenv('MWASP_LIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None,FLIBS))
FLIBS.append('../lib-static/libmwasp.a')
if platform == "darwin":
    FLIBS.append('-L/usr/local/lib')
if platform == "linux" or platform == "linux2":
    FLIBS.append('-lopenblas')

mwasp = []
mwasp.append('em_solver_wrap_mem')
mwasp.append('em_solver_open_geom')
mwasp.append('em_solver_wrap')
mwasp.append('em_solver_wrap_postproc')
mwasp.append('em_sol_exact')
mwasp.append('em_plot_surf_current_vtk')


ext_helm = Extension(
    name='mwaspbie',
    sources=list_files,
    f2py_options=['only:']+mwasp+[':'],
    extra_f90_compile_args=["-std=legacy"],
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Leslie Greengard, Manas Rachh, Felipe Vico",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines for solving electromagnetic transmission problems",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "pytest"
    ],
    ext_modules=[ext_helm],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)
