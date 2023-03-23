import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# Define the number of rows and columns in the matrix
n = 100

# Define the number of bands in the matrix
num_bands = 5

# Create a random sparse band matrix with num_bands bands
diagonal = np.random.rand(n)
off_diagonal = np.random.rand(num_bands-1, n-1)
matrix = sp.diags([off_diagonal[num_bands-2-i,:] for i in range(num_bands-1)] + [diagonal] + [off_diagonal[i,:] for i in range(num_bands-1)], list(range(1-num_bands, num_bands)), shape=(n,n)).tocsr()

# Create a random right-hand side vector as a numpy array
b_np = np.random.rand(n)

# Convert the numpy array to a sparse matrix
#b = sp.csr_matrix(b_np)
b = sp.csr_matrix(b_np.reshape(-1, 1))

# Solve the linear system using the spsolve function
x = spla.spsolve(matrix, b)

# Print the solution vector
print(x)
