#ifndef SRC_FUNCTION_H
#define SRC_FUNCTION_H

#include "vector.h"

namespace ODE
{
  /**
   * Interface for vector-valued functions.
   * A concrete function needs to derive from this interface and override the value() method.
   *
   * @tparam Number The type of number. `double` and `std::complex` are allowed.
   */
  template <typename Number>
  class Function
  {
  public:
    /**
     * Virtual destructor.
     */
    virtual ~Function() = default;

    /**
     * Return the value of the function computed from Vector @p input and @p time in Vector @p
     * output.
     */
    virtual void
    value(Vector<Number>& output, const Vector<Number>& input, double time) = 0;
  };

  template <typename Number>
  class Scaler : public Function<Number> // The scaler class
  {
  public:
    Number scalar;        // The scalar
    Scaler(Number scalar) // Constructor
    {
      this->scalar = scalar; // The scalar is set
    }
    void
    value(Vector<Number>& output, const Vector<Number>& input, double time) // The value function
    {
      Number res[(input).size]; // The result vector is created
      for (int i = 0; i < (input).size; i++)
        {
          Number input_number = input[i]; // The i-th element of the input vector is taken
          Number output_number =
            input_number *
            this->scalar; // The i-th element of the input vector is multiplied by the scalar
          res[i] = output_number; // The i-th element of the output vector is set
        }
      output.set_elements(res); // The output vector is set
    }
  };
  template <typename Number>
  class Linear_Comb : public Function<Number> // The linear combination class
  {
  public:
    Vector<Number>     scalars;        // The scalars vector
    Function<Number>** first_function; // The pointer of the pointer of the first function
    int                size;           // The size of the functions array

    Linear_Comb(Function<Number>* functions[], Vector<Number> scalars) // Constructor
    {
      this->scalars = scalars;              // The scalars vector is set
      this->size    = sizeof(functions[0]); // The size of the functions array is set

      Function<Number>**
        new_array[scalars.size]; // The array of the pointers to the pointers of the functions
      for (int j = 0; j < scalars.size; j++)
        {
          new_array[j] = &(functions[j]); // The pointers of the pointers of the functions are set
        }
      this->first_function =
        new_array[0]; // The pointer of the pointer of the first function is set
    }
    void
    value(Vector<Number>& output, const Vector<Number>& input, double time) // The value function
    {
      Vector<Number> res(input.size); // The result vector is created
      for (int i = 0; i < scalars.size; ++i)
        {
          Vector<Number> temp(input.size); // The temporary vector is created
          (**(this->first_function + i))
            .value(temp, input, time);       // The value of the i-th function is computed
          Scaler<Number> scaler(scalars[i]); // The scaler is created
          scaler.value(temp, temp, time);    // The scaler is applied to the temporary vector
          res = res + temp;                  // The temporary vector is added to the result vector
        };
      output = res; // The result vector is set
    };
  };

  template <typename Number>
  class Function_Test : public ODE::Function<Number> // The function test class
  {
  public:
    double return_value;

    Function_Test(Number return_value) // Constructor
    {
      this->return_value = return_value; // The return value is set
    }
    void
    value(ODE::Vector<double>& output, const ODE::Vector<double>& input, double time)
    {
      double res_array[3] = {1 * this->return_value,
                             2 * this->return_value,
                             3 * this->return_value}; // Example of a function
      output.set_elements(res_array);                 // The output vector is set
    };
  };
} // namespace ODE

#endif // SRC_FUNCTION_H