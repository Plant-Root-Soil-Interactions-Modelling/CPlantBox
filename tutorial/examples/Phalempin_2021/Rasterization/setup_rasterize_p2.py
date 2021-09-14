import fileinput
from os import fdopen
import sys
from shutil import copyfile

        
for kk in range(0,20): 
    copyfile("p2_rasterize.py", "phase2/rasterize"+str(kk+1)+".py")
    file_in = "p2_rasterize.py"
    file_out = "phase2/rasterize"+str(kk+1)+".py"
    
    fin = open(file_in,"rt")
    fout = open(file_out,"wt")
    for line in fin:
        fout.write(line.replace("mm =1", "mm ="+str(kk+1)).replace("rs.setSeed(1)", "rs.setSeed("+str(kk+1)+")"))
    fin.close()
    fout.close()
            

