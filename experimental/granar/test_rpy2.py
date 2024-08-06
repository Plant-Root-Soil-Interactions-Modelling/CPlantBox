# sudo apt install python3-rpy2

# In R: update.packages(), and install.packages('ggplot2'), remove.packages()

from rpy2 import robjects

pi = robjects.r['pi']
print(pi)
print(type(pi))
