##################################
# CUT STRUCTURE WITH A PLANE     #
# Silva 28-10-19                 #
##################################

This script needs a couple of custom functions to work.
These are defined in the mudules found in util. You need to have these in your python path.
To achieve this, this add this line in the .bashrc:
  export PYTHONPATH="/path/to/here/util:$PYTHONPATH" 

Also, you need to have a couple of external python packages:
 - Atomistic Simulation Environment (ASE). To install it use (sudo) apt-get install python-ase (debian-based linux) or pip install --upgrade --user ase (MaxOs) (https://wiki.fysik.dtu.dk/ase/install.html)
 - NumPy. To install python -m pip install --user numpy scipy (https://scipy.org/install.html)
 - argparser and logging should already be part of Python3
