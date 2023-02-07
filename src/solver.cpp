#include "solver.h"

int
main()
{
  //// 1. Define time and step size parameters for the solver:

  double t0 = 0;   // initial time
  double tf = 1;   // final time
  double h  = 0.1; // step size (h <= tf-t0)

  //// 2. Define method of integration for the solver (method_adams irrelevant for RK4 and Euler)):

  // std::string method = "RK4"; // Runge-Kutta 4th order
  // std::string method = "Euler"; // Explicit Euler
  std::string method = "Adams"; // Adams-Bashforth (method_adams must be defined)

  // std::string method_adams = "RK4"; // method for first iteration of Adams-Bashforth
  std::string method_adams = "Euler"; // method for first iteration of Adams-Bashforth

  //// 3. Define initial function vector y0 for the solver (you can change the vector dimension):

  ODE::Vector<double> y0(3);                   // vector dimension
  double              y0_array[3] = {1, 2, 3}; // initial function vector
  y0.set_elements(y0_array);                   // set initial function vector

  class Input_Function : public ODE::Function<double>
  {
  public:
    void
    value(ODE::Vector<double>& output, const ODE::Vector<double>& input, double time)
    {
      // Define the input function here:
      // E.g.: y' = y*t*0.1
      double res_array[3] = {input[0] * time * 0.1,
                             input[1] * time * 0.1,
                             input[2] * time * 0.1}; // input function

      output.set_elements(res_array); // set input function
    };                                // end of value function
  };                                  // end of class Input_Function

  Input_Function input_function; // create input function object

  //// 4. Call the solver:

  ODE::interface<double>(&input_function, t0, tf, h, y0, method, method_adams);

  return 0;
}