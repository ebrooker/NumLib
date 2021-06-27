# NumLib
Numerical Library: Demonstration of various numerical methods/recipes.

# Notes
Currently, this package is being written in modern Fortran (90 and newer). There may be extensions to Python, C/C++, and/or MATLAB as time permits.

Additionally, "package" is a bit of a misnomer right now, considering there is no compilation into a library at the moment.

# Install
The package can be compiled by navigating to the NumLib/bin/ directory and either using the CMake file or the compile.sh file (which also executes a test program).

It is recommended having a relatively newer compiler to capture more modern Fortran language elements (such as F2003 and F2008).
  - Compilation has been tested on GNU 9.1.0 and 9.2.0 version gfortran compilers.


# Methods
Numerical methods supported:

  - Linear Algebra
    - Inner and Outer Products (workarounds for intrinsics issues)
    - Gaussian Elimination (Basic, Row Echelon From, and Reduced REF)
    - Gauss-Jordan Inverse using RREF([A|I]) = [I|B] --> B = inv(A)
    - Thomas solver for tridiagonal matrix systems
    - Back substitution for solving Ax=b with Gaussia Elimination
    
  - Numerical differentiation
    - Forward differencing (1st/2nd order stencils up to 4th derivatives)
    - Backward differencing (1st/2nd order stencils up to 4th derivatives)
    - Central differencing (2nd/4th order stencils up to 4th derivatives)
    - Richardson Extrapolation

  - Numerical Integration
    - Ordinary Differential Equations
      - Explicit Runge Kutta Methods
        - Forward Euler, Midpoint, Heun's
        - Ralston's, SSPRK3, Generic RK4

  - Unconstrained Optimization
    - Line Search w/ Backtracking for minimization
      - Steepest Descent, Inexact Newton, Broyden-Fletcher-Goldfarb-Shanno
    - Univariate Rootfinding
      - Newton's, Halley's, Steffenson's
      - Muller's, Ridder's, Secant, Bisection
      - False Position, ITP, Brent-Dekker
      - Inverse Quadratic Interpolation

Numerical methods in development:
  - Unconstrained Optimization
    - Controlled Random Search for minimization (stochastic method)
  - Linear Algebra
    - Forward Substitution
    - LU, LUP, LUPQ, LDU, and Cholesky-variant decompositions
