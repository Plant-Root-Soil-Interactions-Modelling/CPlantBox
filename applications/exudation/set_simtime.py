import fileinput
from os import fdopen
import sys

for ii in range(0,2): 
    for jj in range(0,2): 
        path = "simulations/plant"+str(ii+1)+"/exudate"+str(jj+1)
        
        for kk in range(10,40): 
            file_in = path+"/example_exudate.py"
            file_out = path+"/day"+str(kk+1)+"/example_exudate.py"

            fin = open(file_in,"rt")
            fout = open(file_out,"wt")
            for line in fin:
                fout.write(line.replace("for ii in range(0,1):", "for ii in range(0,"+str(kk+1)+"):"))
            fin.close()
            fout.close()
            

