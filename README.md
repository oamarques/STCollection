# STCollection

STCollection contains matrices that have been used for testing LAPACK's symmetric tridiagonal eigensolvers and bidiagonal SVD algorithms. 
The collection includes symmetric tridiagonal matrices and (upper) bidiagonal matrices. 

The original set of matrices in the collection has been described in 
*Marques, Demmel, Voemel, and Parlett, A Testing Infrastructure for Symmetric Tridiagonal Eigensolvers, ACM TOMS, 35:2008.*
Matrices in the collection have been used in the experiments described in 
*Demmel, Marques, Parlett, and Voemel, Performance and Accuracy of LAPACK's Symmetric Tridiagonal Eigensolvers, SIAM J. Sci. Comput., 30:2008,* and
*Marques, Demmel, and Vasconcelos, Bidiagonal SVD Computation via an Associated Tridiagonal Eigenproblem, ACM TOMS, 46:2020.*
Difficult cases and cases that expose bugs in LAPACK are continuously added to the collection.

**Matrices**

The matrices are stored in the directory DATA, in files with extension dat. In each of these files, the matrix dimension is given in the first row 
and then (row index, diagonal entry, off-diagonal entry) tuples in the subsequent rows. Each dat file (matrix) has a corresponding file with extension 
eig (eigenvalue distribution of the matrix) or sv (singular values of the matrix). The eigenvalue distributions are showed in files with extension 
jpeg, and also in log scale (negative eigenvalues are represented in magenta.)

**Testing Infrastructure**

The folder stetester contains a set of Fortran subroutines that have been used to test LAPACK's symmetric tridiagonal eigensolvers. See the README and
the README.in files in that folder for compilation instructions and a description of the functionalities provided by the testing infrastructure. The infrastructure 
can be also used to generate a variety of matrices (in addition to the ones in the DATA folder).
