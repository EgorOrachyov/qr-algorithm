# QR-algorithm

QR-decomposition based QR-algorithm for eigenvalues evaluation
of symmetric matrix with real values with OpenMP directives for
parallelization of computations for multi-core systems. 
QR-decomposition, used in the algorithm based on modified Gram-Schmidt process.

## Algorithm

Listing 1. QR-algorithm pseudocode:
```
1. function QR-algorithm(in A, in k, out v)
2.  X = A
3.  for i in range 1:k
4.   QR-decomposition(X, Q, R)
5.   X = R * Q
6.  endfor
7.  v = {X[i,i] for i in range 1:n}
8. endfunction
```

Listing 2. QR-decomposition based on MGS algorithm pseudocode:
```
1.  function QR-decomposition(in A, out Q, out R)
2.   V = A
3.   R = E(nxn)
4.   for i in range 1:n
5.    Q[i] = V[i] / |V[i]|
6.    R[i,i] *= |V[i]|
7.    for j in range i+1:n
8.     R[j,i] += <Q[i], V[j]>
9.     V[j] = V[j] - R[j,i] * Q[i]
10.   endfor 
11.  endfor
12. endfunction
```

## Matrices structure

- A - dense matrix of real values in column-major order
- Q - dense matrix of real values in column-major order
- R - dense matrix of real values in row-major order
- X - dense matrix of real values in column-major order

## OpenMP parallelization

- For-each loop lines 7 - 10 in listing 2.
- Matrix-matrix multiplication line 5 in listing 1.

## License

This project licensed under MIT license. Full [license text](https://github.com/EgorOrachyov/qr-algorithm/blob/main/LICENSE.md).

## Also

This task done as part of `High-performace computing I` university course.
