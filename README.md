# Solver Library Project

This repo contains a small solver library for first-order ordinary differential equations (ODEs) of the form

```math
y' = f(y,t)
```

where `y` is a **vector** of unknowns, `y'` its time-derivative and `t` the time.

## Vector Class

This vector class mimics the concept of a mathematical vector with a given number of components and works for `double` and `std::complex<double>` numbers. It supports the following operations:

- Addition: add another vector to the vector itself
- Scaling: scale every component of the vector by a scalar factor
- Accessing a specific element by index
- Printing: print the vector content in a human-readable format

Unit tests were written for the implementation and can be found in the 'test' folder.

## Functions

- A function that scales the input vector by a fixed constant `a`. The time argument is unused. Mathematically speaking: `f(y,t) = a*y`
- A function that takes an arbitrary number of other functions alongside a scaling value for each of them. The function returns the sum of all function     values scaled with the factors. Mathematically speaking:
  `f(y,t) = a_1 * f_1(y,t) + a_2 * f_2(y,t) + ...`

Unit tests were written for the implementation and can be found in the 'test' folder.

## Solvers

The ODEs are solved for a start vector `y_0` at start time `t_0` until a final time `t_f`. The time interval `[t_0,t_f]` is subdivided into a number of time steps of size `h`.
The following numerical schemes were implemented to solve ODEs:

- explicit Euler (EE):
    ```math
    y_{n+1} = y_n + h * f(y,t)
    ```
- 4-stage Runge-Kutta (RK4):

  ![](https://wikimedia.org/api/rest_v1/media/math/render/svg/0b7865da10afe692831b0cdb375f9e41021c5da2)

  ![](https://wikimedia.org/api/rest_v1/media/math/render/svg/4b038c70313036aabe58cdc5d6ec6ecb098dbb70)

- Two-step Adams Bashforth (AB2), with RK4 and Euler as first iteration:

  ![](https://wikimedia.org/api/rest_v1/media/math/render/svg/e7c19f4ef1a113146a6c50f5e32e6f4b07765c4d)

A2B can't be used for the first iteration from `t_0` to `t_0+h`, since it requires two previous solution vectors. For the first step, one can use either EE or RK4.

These schemes are repeatedly applied starting from `t_0` and `y_0` until the final time `t_f` is reached and the final solution vector `y_f` is obtained.

Testing:

Unit tests were written for all schemes and sample scenarios for the solvers implementation and can be found in the 'test' folder:
    - `y' = -y` with `y(0) = 1`
    - `y' = sin(t)` with `y(0) = 1`
    - multidimensional vector
