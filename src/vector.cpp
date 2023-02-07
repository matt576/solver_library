#include "function.h"

int
main()
{
  //// 1. Declaring a Vector object:
  ODE::Vector<double> v(3);
  double              my_vector[3] = {4, 2, 3};
  v.set_elements(my_vector);

  //// 2. Printing the vector:
  v.print_array();

  //// 3. Scaling a vector:
  ODE::Vector<double> v22(3);
  v22 = v.scale(10); // v22 = 10 * v
  v22.print_array();

  //// 4. Accessing the elements of a vector:
  double a = v[1];
  std::cout << a << std::endl;

  //// 5. Adding two vectors:
  v = v + v;
  v.print_array();

  return 0;
}