# Gaussian-Elimination-Matlab

## Overview
This set of MATLAB functions provides implementations for Gaussian Elimination to solve linear equation systems (LES) and compute the inverse of a square matrix. The three functions are named gaussel1, gaussel2, and gaussel3, each building upon the previous one to offer additional features.
- Tests are available within the source codes 

## Functions

1). gaussel1.m - Gaussian Elimination for Linear Equation Systems

- Input Parameters: Coefficient matrix (A) and the right-side vector (b) of LES.

- Output Argument: Solution vector (x).

- Functionality:
  
  . Utilizes MATLAB row-operations for algorithm organization.

  . Handles cases where GE can't be solved without row or column swap.

  . Warns the user in the case of an underdetermined LES, providing a base solution.

  . Optionally displays matrices A(i) during computation.

  . Supports LES with multiple right sides.

2). gaussel2.m - Gaussian Elimination with Partial and Whole Pivoting

- Additional Input Parameter: Boolean parameter to choose between partial and whole pivoting.

- Output: Solution vector (x).

- Functionality:
  
  . Extends gaussel1 to support partial and whole pivoting.

  . Allows the user to choose the pivoting method.

  . Switches to whole pivoting if partial pivoting gets stuck.

  . Informs the user about the switch in pivoting method.

  . Optionally displays matrices A(i) during computation.

  . Takes into account the impact of pivoting on the solution.

3). gaussel3.m - Gaussian Elimination for Matrix Inversion

- Input Parameter: Square matrix for which the inverse is to be computed.

- Output: Inverse matrix.

- Functionality:
  
  . Checks input arguments before computation.

  . Computes the determinant of the matrix.

  . Calls gaussel1 to handle multiple right-sides during computation.
