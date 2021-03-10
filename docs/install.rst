Installation
============

Obtaining fmm3dbie
******************

The source code can be downloaded from https://github.com/mrachh/miniwasp 


Dependencies
************

This library is supported for unix/linux, and Mac OSX.

For the basic libraries

* Fortran compiler, such as ``gfortran`` packaged with GCC
* GNU make
* `FMM3D <https://github.com/flatironinstitute/FMM3D>`_
* `fmm3dbie <https://github.com/fastalgorithms/fmm3dbie>`_
* Blocked linear algebra routines and LAPACK (default used: Netlib BLAS
  on Linux machines, and framework accelerate on MacOS)

Optional:

* for building Python wrappers you will need ``python`` and ``pip`` 

Quick install instructions
*********************************************

Make sure you have dependencies downloaded, and `cd` into your fmm3dbie
directory. 

-  For linux, run ``make python``.
-  For Mac OSX, run ``cp make.inc.macos.gnu make.inc`` followed by ``make python``.

This should compile the static library
in ``lib-static/``, and build the python package ``mwaspbie`` in
``python``.

In order to link against the dynamic library, you will have to update
the ``LD_LIBRARY_PATH`` environment
variable on Linux and ``DYLD_LIBRARY_PATH`` environment variable on Mac OSX
to the installation directory.
You may then link to the miniwasp library using the ``-L$(FMMBIE_INSTALL_DIR) -lfmm3dbie -L$(FMM_INSTALL_DIR) -lfmm3d -L$(MWASPBIE_INSTALL_DIR) -lmwasp`` 
option.

.. note :: 
   On MacOSX, /usr/local/lib is included by default in the
   DYLD_LIBRARY_PATH.


To verify successful compilation of the program, run ``make test``
which compiles some fortran test drivers in ``test/`` linked against
the static library, after which it
runs the test programs. 

If ``make test`` fails, see more detailed instructions below. 

If there is an error in testing on a standard set-up,
please file a bug report as a New Issue at https://github.com/mrachh/miniwasp/issues

Building Python wrappers
****************************

First make sure you have python3 and pip3 installed. 

You may then execute ``make python`` (after copying over the
operating system specific make.inc.* file to make.inc) which calls
pip for the install. 

See ``python/example1.py`` to see
usage examples for the Python wrappers.


A few words about Python environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There can be confusion and conflicts between various versions of Python and installed packages. It is therefore a very good idea to use virtual environments. Here's a simple way to do it (after installing python-virtualenv)::

  Open a terminal
  virtualenv -p /usr/bin/python3 env1
  . env1/bin/activate

Now you are in a virtual environment that starts from scratch. All pip installed packages will go inside the env1 directory. (You can get out of the environment by typing ``deactivate``)


Tips for installing dependencies
**********************************

On Ubuntu linux
~~~~~~~~~~~~~~~~

On Ubuntu linux (assuming python3 as opposed to python)::

  sudo apt-get install make build-essential gfortran libopenblas-dev 


On Fedora/CentOS linux
~~~~~~~~~~~~~~~~~~~~~~~~

On a Fedora/CentOS linux system, these dependencies can be installed as 
follows::

  sudo yum install make gcc gcc-c++ gcc-gfortran libgomp openblas-devel 

.. _mac-inst:

On Mac OSX
~~~~~~~~~~~~~~~~~~~~~~~~

First setup Homebrew as follows. If you don't have Xcode, install
Command Line Tools by opening a terminal (from /Applications/Utilities/)
and typing::

  xcode-select --install

Then install Homebrew by pasting the installation command from
https://brew.sh

Then do::
  
  brew install gcc openblas 
  

Tips for installing optional dependencies
******************************************

Installing python3 and pip3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Ubuntu linux
##################

::

  sudo apt-get install python3 python3-pip


On Mac OSX
############

Make sure you have homebrew installed. See `Tips for installing dependencies -> On Mac OSX <install.html#mac-inst>`__ 

::
  
  brew install python3

Then use `make python3` instead of `make python`. You will only need to
do this in case the default version of `python` and `pip` is not >=3.0 


