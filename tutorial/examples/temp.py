import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

volNew = np.array([0.00073057, 0.00126148, 0.00217822, 0.00376115, 0.00649442,
       0.01121398, 0.0193633 , 0.03343482, 0.05773228])
NC = 10
matrix_size = NC -1
sub_diag_values = -1
main_diag_values = 1  # Zeros on the main diagonal
matrix = np.diag(np.full(matrix_size-1,sub_diag_values), k=1) + np.diag(np.full(matrix_size,main_diag_values), k=0)
matrix[-1,] = volNew