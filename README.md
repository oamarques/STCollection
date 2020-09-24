# STCollection

STCollection contains matrices that have been used for testing LAPACK's symmetric tridiagonal eigensolvers and bidiagonal SVD algorithms. 
The collection includes symmetric tridiagonal matrices and (upper) bidiagonal matrices. 

The original set of matrices in the collection was described in Marques, Demmel, Voemel, and Parlett, A Testing Infrastructure for Symmetric 
Tridiagonal Eigensolvers, ACM TOMS, 35:2008, and the collection was used in the experiments described in Demmel, Marques, Parlett, and Voemel, 
Performance and Accuracy of LAPACK's Symmetric Tridiagonal Eigensolvers, SIAM J. Sci. Comput., 30:2008. Difficult cases, and cases that expose 
occasional bugs in LAPACK, are continuously added to the collection.

The matrices are stored in the directory DATA, in files with extension dat. In each of these files, the matrix dimension is given in the first row 
and then tuples (row index, diagonal entry, off-diagonal entry) in the subsequent rows. Each dat file (matrix) has a corresponding file with extension 
eig (eigenvalue distribution of the matrix) or sv (singular values of the matrices). The eigenvalue distributions are showed in files with extension 
jpeg, and also in log scale (negative eigenvalues are represented with a different color.)
