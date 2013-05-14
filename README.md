# motop

mtop is a Matlab package that enables topology optimization using the CALFEM matlab library (see http://sourceforge.net/projects/calfem/).

## Introduction
In general, a topoplogy optimization problem involves solving *x* in the following equation:

    |  min     f(x,y{x})
    | x,y{x}
    |
    |          | behavioral constrains on y{x}
    |  s.t.    | equilibrium constraint
    |          | design constrains on x

where:
  * *f* is the objective function that is to be minimized; it usually meassures some feature of the structre such as it's stiffness,
  * *x* is a vector of design parameters; usually *x* represents some material property or a geometrical property, e.g. density or thickness,
  * *y{x}* describes some response of the stucture, e.g. the stress.

