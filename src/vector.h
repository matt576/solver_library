#ifndef SRC_VECTOR_H
#define SRC_VECTOR_H

#include <complex>
#include <iostream>
namespace ODE
{
  /**
   * The numerical Vector class.
   *
   * @tparam Number The type of number. `double` and `std::complex` are allowed.
   */
  template <typename Number>
  class Vector
  {
  public:
    int     size;    // Size of the vector
    Number* ptr_arr; // Pointer of the first element

    /**
     * Constructor of the vector.
     *
     * @param size The size of the vector.
     */
    Vector(int size = 1)
    {
      Number array[size];     // Array of the vector is created
      this->size = size;      // Size of the vector is set
      ptr_arr    = &array[0]; // Pointer of the first element is set
    };

    Vector
    operator+(Vector& another_vector) // Defining the operator + for summing vectors
    {
      Vector res(this->size); // Result vector is created

      Number* new_array = new Number[this->size]; // Array of the result vector is created
      for (int i = 0; i < this->size; i++)
        {
          new_array[i] = *(this->ptr_arr + i) + *(another_vector.ptr_arr + i);
        }

      res.ptr_arr = &new_array[0]; // Pointer of the first element of the result vector is set
      res.size    = this->size;    // Size of the result vector is set

      return res;
    };

    Number
    operator[](int i) const
    {
      return ptr_arr[i]; // Defining the operator [] for accessing the elements of the vector
    }

    void
    set_elements(Number arr[]) // Defining the method for setting the elements of the vector
    {
      Number* new_array = new Number[size];
      for (int j = 0; j < this->size; j++)
        {
          new_array[j] = arr[j]; // Elements of the vector are set
        }
      this->ptr_arr = &new_array[0]; // Pointer of the first element of the vector is set
    };

    void
    print_array() // Defining the method for printing the elements of the vector
    {
      std::cout << "The array is: [";
      for (int i = 0; i < this->size - 1; i++)
        {
          std::cout << *(this->ptr_arr + i) << " "; // Elements of the vector are printed
        }
      std::cout << *(this->ptr_arr + this->size - 1) << "]" << std::endl;
    }

    Vector<Number>
    scale(Number a) // Defining the method for scaling the vector
    {
      Vector<Number> res(this->size);                    // Result vector is created
      Number*        new_array = new Number[this->size]; // Array of the result vector is created
      for (int i = 0; i < this->size; i++)
        {
          new_array[i] = *(this->ptr_arr + i) * a; // Elements of the result vector are set
        }
      res.ptr_arr = &new_array[0]; // Pointer of the first element of the result vector is set
      res.size    = this->size;    // Size of the result vector is set
      return res;
    };
  };
}; // namespace ODE

#endif // SRC_VECTOR_H
