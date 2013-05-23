# Contents

## Filter functions
Filter functions, 'F*.m', have the following syntax:

    x = F(H, Hs, x)
    [x, dF] = F(H, Hs, x, dF)

where 'H' and 'Hs' can be computed using the __FUsetup.m__ function.

Functions:
* __FDensity.m__ density filter; operates on both x and dF
* __FSensitivity.m__ sensitivity filter; operates on dF

## Young's Modulus Interpolation Functions
All Young's modulus interpolation functions, 'E*.m', have the following syntax:

    E = E(varargin)
    [E, dE] = E(varargin)
    [E, dE, ddE] = E(varargin)

Functions:
* __ELin.m__ linear interpolation of Young's modulus
* __EModSIMP.m__ ModSIMP interpolation of Young's modulus
* __ERAMP.m__ RAMP interpolation of Young's modulus