#PSP_readDataFile.py  
from __future__ import print_function   
import csv
import numpy as np
                
def scanDataFile(file, delimiter):
    reader = csv.reader(open(file, "rt"), delimiter=delimiter)
    nrRows = 0
    for row in reader:
        if (nrRows == 0): 
            nrCols = len(row) 
        else:
            if (len(row) != nrCols): 
                #wrong file: nr of fields not coherent
                return (nrRows, nrCols, False)
        nrRows += 1
    return (nrRows, nrCols, True)
        
# multi-columns or single row data (only numeric)                
def readDataFile(file, nrHeaderFields, delimiter, isPrintScreen):
    nrRows, nrCols, isFileOk = scanDataFile(file, delimiter)
    if (isFileOk == False): return (nrRows, False)
    
    if (nrRows == 1):
        A = np.zeros((nrRows, nrCols-nrHeaderFields))
    else:
        A = np.zeros((nrRows-nrHeaderFields, nrCols))
        
    reader = csv.reader(open(file, "r"), delimiter=delimiter)
    i = 0    
    for row in reader:
        if (isPrintScreen): print(row)
        if (nrRows == 1):
            #single row
            for j in range(nrHeaderFields, len(row)):
                A[i, j-nrHeaderFields] = float(row[j])
        else:
            #multi-columns
            if (i >= nrHeaderFields):
                for j in range(0, len(row)):
                    A[(i-nrHeaderFields), j] = float(row[j])
        i += 1 
          
    return(A, True)

# multi-columns generic data (text, date, etc.)
def readGenericDataFile(file, nrHeaderRows, delimiter, isPrintScreen):
    nrRows, nrCols, isFileOK = scanDataFile(file, delimiter)
    if (isPrintScreen): print ('nrRows =', nrRows, ' nrCols =', nrCols)
    
    if (isFileOK == False): return (nrRows, False)
    myReader = csv.reader(open(file, "rt"), delimiter=delimiter)
    
    A = []
    i = 0
    for myRow in myReader:
        if (isPrintScreen): print(myRow)
        if (i >= nrHeaderRows):
            A.append(myRow)
        i += 1     
    return(A, True)