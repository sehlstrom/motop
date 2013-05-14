# motop
mtop is a Matlab package that enables topology optimization using the [CALFEM Matlab library](http://sourceforge.net/projects/calfem/ "CALFEM Matlab library").

## Introduction
The general aim of topology optimization is to determine the optimal placement of a given material in space. In other words, the goal is to determine which points should be filled with material and which points should be voids. When applying finite element discretization to the problem, we can think of the geometrical representation as pixels of black and white representing solid and void areas. Given a design domain  Ω of finite elements, we thus seek to find a optimal subset Ωmat of elements that should be filled with material.

In general, a topology optimization problem involves solving *x* in the following equation:

    |  min     f(x,y{x})
    | x,y{x}
    |
    |          | behavioral constrains on y{x}
    |  s.t.    | equilibrium constraint
    |          | design constrains on x

where:
  * *f* is the objective function that is to be minimized; it usually measures some feature of the structure such as it's stiffness,
  * *x* is a vector of design parameters; usually *x* represents some material property or a geometrical property, e.g. density or thickness,
  * *y{x}* is a state variable describing some response of the structure, e.g. the stress.

