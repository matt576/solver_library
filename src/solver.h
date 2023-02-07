#ifndef SRC_SOLVER_H
#define SRC_SOLVER_H

#include "function.h"

namespace ODE
{
  template <typename Number>
  class Solver // The solver class
  {
  public:
    Function<Number>* f; // The pointer of the function to be integrated. Function is an abstract
                         // class, which does not allow for it to passed as an argument directly.
    double         t0;   // Initial time
    double         tf;   // Final time
    double         h;    // Step size
    Vector<Number> y;    // The solution vector
    Vector<Number> y0;   // The initial condition vector

    virtual void
    solve(Function<Number>* f, double t0, double tf, double h, Vector<Number> y0)
    {
      this->f  = f;              // The function is set
      this->t0 = t0;             // The initial time is set
      this->tf = tf;             // The final time is set
      this->h  = h;              // The step size is set
      this->y0 = y0;             // The initial condition vector is set
    };                           // The solve function
    virtual ~Solver() = default; // The virtual destructor
  };
  namespace Euler // Explicit Euler namespace
  {
    template <typename Number>
    class EulerSolver : public Solver<Number> // The explicit Euler solver class
    {
    public:
      Function<Number>* f;   // The pointer of the function to be integrated
      double            t0;  // The initial time
      double            tf;  // The final time
      double            h;   // Step size
      Vector<Number>    y;   // The solution vector
      Vector<Number>    y0;  // The initial condition vector
      Vector<Number>    y00; // The vector of the previous step (used in the Adams-Bashforth scheme)

      EulerSolver(Function<Number>* f, double t0, double tf, double h, Vector<Number> y0)
      {
        this->f  = f;                  // The function is set
        this->t0 = t0;                 // The initial time is set
        this->tf = tf;                 // The final time is set
        this->h  = h;                  // The step size is set
        this->y0 = y0;                 // The initial condition vector is set
        this->solve(f, t0, tf, h, y0); // The solve function is called
      };

      void
      solve(Function<Number>* f, double t0, double tf, double h, Vector<Number>& y0)
      {
        int            n = (tf - t0) / h; // Number of iterations
        Vector<Number> y = y0; // The solution vector is set to the initial condition vector
        double         t = t0; // The initial time is set

        for (int i = 0; i < n; i++) // The loop over the iterations
          {
            y00 = y;                          // The previous step is set to the current step
            ODE::Vector<Number> res(y0.size); // The result vector is created with the same size as
                                              // the initial condition vector

            (*(this->f)).value(res, y00, t0 + i * h); // The function is evaluated at the current
                                                      // time and the current solution vector
            ODE::Scaler<Number> scaler(h);            // The scaler object is created
            scaler.value(res, res, t0 + i * h); // The result vector is scaled by the step size

            y = y + res; // The solution vector is updated

            std::cout << "iteration: " << i << std::endl;
            std::cout << "y: ";
            y.print_array();
          };
        this->y   = y;   // The solution vector is set
        this->y00 = y00; // The previous step is set
        std::cout << std::endl << "Solution Euler: ";
        (this->y).print_array(); // The solution vector is printed
      };
    };
  } // namespace Euler

  namespace RungeKutta4
  {
    template <typename Number>
    class k1 : public Function<Number>
    {
    public:
      Function<Number>* f_ptr;

      k1(Function<Number>* f_ptr)
      {
        this->f_ptr = f_ptr;
      };
      void
      value(Vector<Number>& output, const Vector<Number>& input, double time)
      {
        (*(this->f_ptr)).value(output, input, time);
      }
    };

    template <typename Number>
    class k2 : public Function<Number>
    {
    public:
      Function<Number>* f_ptr;
      Vector<Number>    k1;
      double            h;

      k2(Function<Number>* f_ptr, double h, Vector<Number> k1)
      {
        this->f_ptr = f_ptr;
        this->h     = h;
        this->k1    = k1;
      };
      void
      value(Vector<Number>& output, const Vector<Number>& input, double time)
      {
        Scaler<Number> scaler(this->h / 2);
        Vector<Number> d_vector(input.size);
        scaler.value(d_vector, this->k1, time);

        Number input_argument_array[input.size];

        for (int i = 0; i < input.size; i++)
          {
            input_argument_array[i] = input[i] + d_vector[i];
          }
        Vector<Number> input_argument(input.size);
        input_argument.set_elements(input_argument_array);

        (*(this->f_ptr)).value(output, input_argument, time + (this->h) / 2);
      }
    };

    template <typename Number>
    class k3 : public Function<Number>
    {
    public:
      Function<Number>* f_ptr;
      Vector<Number>    k2;
      double            h;

      k3(Function<Number>* f_ptr, double h, Vector<Number>& k2)
      {
        this->f_ptr = f_ptr;
        this->h     = h;
        this->k2    = k2;
      };
      void
      value(Vector<Number>& output, const Vector<Number>& input, double time)
      {
        Scaler<Number> scaler(this->h / 2);
        Vector<Number> d_vector(input.size);
        scaler.value(d_vector, this->k2, time);

        Number input_argument_array[input.size];

        for (int i = 0; i < input.size; i++)
          {
            input_argument_array[i] = input[i] + d_vector[i];
          }
        Vector<Number> input_argument(input.size);
        input_argument.set_elements(input_argument_array);

        (*(this->f_ptr)).value(output, input_argument, time + (this->h) / 2);
      }
    };

    template <typename Number>
    class k4 : public Function<Number>
    {
    public:
      Function<Number>* f_ptr;
      Vector<Number>    k3;
      double            h;

      k4(Function<Number>* f_ptr, double h, Vector<Number>& k3)
      {
        this->f_ptr = f_ptr;
        this->h     = h;
        this->k3    = k3;
      };
      void
      value(Vector<Number>& output, const Vector<Number>& input, double time)
      {
        Scaler<Number> scaler(this->h);
        Vector<Number> d_vector(input.size);
        scaler.value(d_vector, this->k3, time);

        Number input_argument_array[input.size];

        for (int i = 0; i < input.size; i++)
          {
            input_argument_array[i] = input[i] + d_vector[i];
          }
        Vector<Number> input_argument(input.size);
        input_argument.set_elements(input_argument_array);

        (*(this->f_ptr)).value(output, input_argument, time + (this->h));
      }
    };

    template <typename Number> // The RungeKuttaSolver class is derived from the Solver class
    class RungeKuttaSolver : public Solver<Number>
    {
    public:
      Function<Number>* f;   // The function pointer is set
      double            t0;  // The initial time is set
      double            tf;  // The final time is set
      double            h;   // The step size is set
      Vector<Number>    y;   // The solution vector is set
      Vector<Number>    y0;  // The initial condition is set
      Vector<Number>    y00; // The previous step is set
      RungeKuttaSolver(Function<Number>* f,
                       double            t0,
                       double            tf,
                       double            h,
                       Vector<Number>&   y0) // The constructor is defined
      {
        this->f  = f;                  // The function pointer is set
        this->t0 = t0;                 // The initial time is set
        this->tf = tf;                 // The final time is set
        this->h  = h;                  // The step size is set
        this->y0 = y0;                 // The initial condition is set
        this->solve(f, t0, tf, h, y0); // The solve function is called
      };
      void
      solve(Function<Number>* f,
            double            t0,
            double            tf,
            double            h,
            Vector<Number>&   y0) // The solve function is defined
      {

        int n = (tf - t0) / h; // The number of steps is calculated

        Vector<Number> y = y0; // The solution vector is set to the initial condition
        for (int i = 0; i < n; ++i)
          {
            y00 = y; // The previous step is set to the current solution

            k1<Number>     k1(f);              // The k1 function is defined
            Vector<Number> k1_value(y0.size);  // The k1 value vector is defined
            k1.value(k1_value, y, t0 + i * h); // The k1 value is calculated

            // std::cout << "k1_value: ";
            // k1_value.print_array();

            k2<Number>     k2(f, h, k1_value); // The k2 function is defined
            Vector<Number> k2_value(y0.size);  // The k2 value vector is defined
            k2.value(k2_value, y, t0 + i * h); // The k2 value is calculated

            // std::cout << "k2_value: ";
            // k2_value.print_array();

            k3<Number>     k3(f, h, k2_value); // The k3 function is defined
            Vector<Number> k3_value(y0.size);  // The k3 value vector is defined
            k3.value(k3_value, y, t0 + i * h); // The k3 value is calculated

            // std::cout << "k3_value: ";
            // k3_value.print_array();

            k4<Number>     k4(f, h, k3_value); // The k4 function is defined
            Vector<Number> k4_value(y0.size);  // The k4 value vector is defined
            k4.value(k4_value, y, t0 + i * h); // The k4 value is calculated

            // std::cout << "k4_value: ";
            // k4_value.print_array();

            double         coefficients[4] = {1, 2, 2, 1};  // The coefficients are defined
            Vector<double> coefficients_vector(4);          // The coefficients vector is defined
            coefficients_vector.set_elements(coefficients); // The coefficients vector is set

            Function<Number>*   newp[4] = {&k1, &k2, &k3, &k4}; // The function pointers are defined
            Linear_Comb<Number> linear_combination(
              newp, coefficients_vector); // The linear combination function is defined
            Vector<Number> linear_combination_value(
              y0.size); // The linear combination value vector is defined

            linear_combination.value(linear_combination_value,
                                     y,
                                     t0 + i * h); // The linear combination value is calculated
            // std::cout << "linear_combination_value: ";
            // linear_combination_value.print_array();

            Scaler<Number> scaler(h / 6);         // The scaler function is defined
            Vector<Number> scaler_value(y0.size); // The scaler value vector is defined
            scaler.value(scaler_value,
                         linear_combination_value,
                         t0 + i * h); // The scaler value is calculated

            // std::cout << "scaler_value: ";
            // scaler_value.print_array();

            y = y + scaler_value; // The solution vector is updated

            std::cout << "iteration: " << i << std::endl;
            std::cout << "y: ";
            y.print_array();
          };
        this->y = y;
        this->print_solution(); // The solution is printed
      };

      void
      print_solution()
      {
        std::cout << std::endl << "Solution RK4: ";
        (this->y).print_array();
      }
      ~RungeKuttaSolver() = default; // The destructor is defined
    };
  } // namespace RungeKutta4

  namespace AdamsBashforth
  {
    template <typename Number>
    class AdamsBashforthSolver : public Solver<Number>
    {
    public:
      Function<Number>* f;      // The function pointer
      double            t0;     // The initial time
      double            tf;     // The final time
      double            h;      // The step size
      Vector<Number>    y;      // The solution vector at timestep 0
      Vector<Number>    y1;     // The solution vector at timestep 1
      Vector<Number>    y2;     // The solution vector at timestep 2
      Vector<Number>    y0;     // The initial condition
      std::string       method; // The method

      AdamsBashforthSolver(Function<Number>* f,
                           double            t0,
                           double            tf,
                           double            h,
                           Vector<Number>&   y0,
                           std::string       method)
      {
        this->f      = f;                      // The function pointer is set
        this->t0     = t0;                     // The initial time is set
        this->tf     = tf;                     // The final time is set
        this->h      = h;                      // The step size is set
        this->y0     = y0;                     // The initial condition is set
        this->method = method;                 // The method is set
        this->solve(f, t0, tf, h, y0, method); // The solve function is called
      };
      void
      solve(Function<Number>* f,
            double            t0,
            double            tf,
            double            h,
            Vector<Number>    y0,
            std::string       method)
      {

        int            n = (tf - t0) / h; // The number of timesteps is calculated
        Vector<Number> y = y0;            // The solution vector is set

        // The method for the first timestep is defined
        if (method == "RK4")
          {
            ODE::RungeKutta4::RungeKuttaSolver<Number> rk4(f, t0, 1 * h, h, y0);
            y2 = rk4.y;
          }
        else if (method == "Euler")
          {
            ODE::Euler::EulerSolver<Number> euler(f, t0, 1 * h, h, y0);
            y2 = euler.y;
          }
        else
          {
            std::cout << "Invalid method" << std::endl;
          }

        // The solution steps are updated
        y1 = y2; // y_1 is updated from the other method
        y  = y0; // y is chosen as the inital value

        for (int i = 1; i < n; i++)
          {
            Vector<Number> fx(y0.size);   // The function value vector is defined
            Vector<Number> fx_1(y0.size); // The function value vector is defined at timestep +1

            (*f).value(fx, y, t0 + i * h); // The function value is calculated
            (*f).value(fx_1,
                       y1,
                       t0 + (i + 1) * h); // The function value is calculated at timestep +1

            Vector<Number> part_1(y0.size);   // The part 1 of the solution vector is defined
            Scaler<Number> scaler(3 * h / 2); // The scaler function is defined
            scaler.value(part_1,
                         fx_1,
                         t0 + (i + 1) * h); // The scaler value is calculated to the part 1

            Vector<Number> part_2(y0.size);        // The part 2 of the solution vector is defined
            Scaler<Number> scaler2(-h / 2);        // The scaler function is defined
            scaler2.value(part_2, fx, t0 + i * h); // The scaler value is calculated to the part 2

            y2 = y1 + part_1 + part_2; // The solution vector is updated at timestep +2
            y  = y1;                   // The solution vector is updated at timestep 0
            y1 = y2;                   // The solution vector is updated at timestep +1

            std::cout << "iteration: " << i << std::endl;
            std::cout << "y: ";
            y2.print_array(); // The solution vector is printed
          }
        this->y2 = y2; // The solution vector is set
        std::cout << "Solution AB: ";
        (this->y2).print_array();
      };
      ~AdamsBashforthSolver() = default; // The destructor is defined
    };

  } // namespace AdamsBashforth

  // The interface function is defined
  template <typename Number>
  void
  interface(Function<Number>* f,      // The function pointer
            double            t0,     // The initial time
            double            tf,     // The final time
            double            h,      // The step size
            Vector<Number>&   y0,     // The initial condition
            std::string       method, // The method
            std::string       method_adams) // The method for Adams-Bashforth
  {
    // The method is defined and solved
    if (method == "RK4")
      {
        RungeKutta4::RungeKuttaSolver<Number> rk4_solver(f, t0, tf, h, y0); 
      }
    else if (method == "Euler")
      {
        Euler::EulerSolver<Number> euler_solver(f, t0, tf, h, y0);
      }
    else if (method == "Adams")
      {
        // The method for Adams-Bashforth is defined and solved
        if (method_adams == "Euler")
          {
            AdamsBashforth::AdamsBashforthSolver<double> ab_solver(f, t0, tf, h, y0, "Euler");
          }
        else if (method_adams == "RK4")
          {
            AdamsBashforth::AdamsBashforthSolver<double> ab_solver(f, t0, tf, h, y0, "RK4");
          }
      }
  }
} // namespace ODE

#endif // SRC_SOLVER_H