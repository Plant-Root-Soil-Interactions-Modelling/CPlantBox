# sudo apt install python3-rpy2

# In R: update.packages(), and install.packages('ggplot2'), remove.packages()
# Notes:
# first, install.packages('devtools'), read messages, and install missing ubuntu libraries (sudo apt install xxx)
# if version conflicts, just reinstall packages with install.packages()
# package "sf" needs: sudo apt install libgdal-dev
#

#
# Pyhton environments:
#
# sudo pip3 install virtualenv
# virtualenv myRenv # will create a directory myRenv with python and libraries
# source myRenv/bin/activate
# pip3 install stuff

# in are use the environtment with:
# use_virtualenv("~/workspace/granar/R/myRenv")
# but did not work, (specific Python 3 version required)

#
# for conda, see https://granar.github.io/granar_examples/
#

from rpy2 import robjects

pi = robjects.r['pi']
print(pi)
print(type(pi))
