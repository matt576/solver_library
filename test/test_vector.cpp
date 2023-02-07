#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN // This tells doctest to provide a main() - only do this
                                           // in one cpp file per test executable!

#include "../src/solver.h"
#include "doctest.h"

#include <sstream>
#include <string>

TEST_CASE("Vector Class Addition")
{
  ODE::Vector<double> v(3);
  ODE::Vector<double> w(3);
  ODE::Vector<double> z(3);
  ODE::Vector<double> y(3);

  double my_vector[3]   = {1, 2, 3};
  double my_vector_2[3] = {0, 1, 2};
  double my_vector_3[3] = {1, 3, 5};

  v.set_elements(my_vector);
  w.set_elements(my_vector_2);
  y.set_elements(my_vector_3);

  z = v + w;

  for (int i = 0; i < 3; ++i)
    CHECK(z[i] == y[i]);
};

TEST_CASE("Vector Addition with Complex Numbers")
{
  ODE::Vector<std::complex<double>> vj(3);
  ODE::Vector<std::complex<double>> wj(3);
  ODE::Vector<std::complex<double>> zj(3);
  ODE::Vector<std::complex<double>> yj(3);

  std::complex<double> my_vector_1[3] = {0, 0, 0};
  my_vector_1[0]                      = std::complex<double>(1.0, 1.0);
  my_vector_1[1]                      = std::complex<double>(1.0, 2.0);
  my_vector_1[2]                      = std::complex<double>(3.0, 4.0);

  std::complex<double> my_vector_2[3] = {0, 0, 0};
  my_vector_2[0]                      = std::complex<double>(1.0, 1.0);
  my_vector_2[1]                      = std::complex<double>(1.0, 1.0);
  my_vector_2[2]                      = std::complex<double>(1.0, 1.0);

  std::complex<double> my_vector_3[3] = {0, 0, 0};
  my_vector_3[0]                      = std::complex<double>(2.0, 2.0);
  my_vector_3[1]                      = std::complex<double>(2.0, 3.0);
  my_vector_3[2]                      = std::complex<double>(4.0, 5.0);

  vj.set_elements(my_vector_1);
  wj.set_elements(my_vector_2);
  yj.set_elements(my_vector_3);

  zj = vj + wj;

  for (int i = 0; i < 3; ++i)
    CHECK(zj[i] == yj[i]);
};

TEST_CASE("Vector Class Scaling")
{
  ODE::Vector<double> v(3);
  ODE::Vector<double> y(3);
  double              my_vector[3]   = {1, 2, 3};
  double              my_vector_2[3] = {2, 4, 6};
  v.set_elements(my_vector);
  y.set_elements(my_vector_2);

  v = v.scale(2);

  for (int i = 0; i < 3; ++i)
    CHECK(v[i] == y[i]);
};

TEST_CASE("Vector Scaling with Complex")
{
  ODE::Vector<std::complex<double>> vj(3);
  ODE::Vector<std::complex<double>> yj(3);

  std::complex<double> my_vector_1[3] = {0, 0, 0};
  my_vector_1[0]                      = std::complex<double>(1.0, 1.0);
  my_vector_1[1]                      = std::complex<double>(1.0, 2.0);
  my_vector_1[2]                      = std::complex<double>(3.0, 4.0);

  std::complex<double> my_vector_2[3] = {0, 0, 0};
  my_vector_2[0]                      = std::complex<double>(2.0, 2.0);
  my_vector_2[1]                      = std::complex<double>(2.0, 4.0);
  my_vector_2[2]                      = std::complex<double>(6.0, 8.0);

  vj.set_elements(my_vector_1);
  yj.set_elements(my_vector_2);

  vj = vj.scale(2);

  for (int i = 0; i < 3; ++i)
    CHECK(vj[i] == yj[i]);
};

TEST_CASE("Vector Class Index")
{
  ODE::Vector<double> v(3);
  double              my_vector[3] = {1, 2, 3};
  v.set_elements(my_vector);

  double a = v[1];

  CHECK(a == 2);
};

TEST_CASE("Vector Index with Complex")
{
  ODE::Vector<std::complex<double>> vj(3);

  std::complex<double> my_vector_1[3] = {0, 0, 0};
  my_vector_1[0]                      = std::complex<double>(1.0, 1.0);
  my_vector_1[1]                      = std::complex<double>(1.0, 2.0);
  my_vector_1[2]                      = std::complex<double>(3.0, 4.0);

  vj.set_elements(my_vector_1);

  std::complex<double> a = vj[1];

  CHECK(a == std::complex<double>(1.0, 2.0));
};

TEST_CASE("Vector Printing")
{
  ODE::Vector<double> v(3);
  double              my_vector[3] = {1, 2, 3};
  v.set_elements(my_vector);

  std::stringstream buffer;
  std::streambuf*   old = std::cout.rdbuf(buffer.rdbuf());

  v.print_array();

  std::cout.rdbuf(old);

  std::string expected_output = "The array is: [1 2 3]\n";
  std::string actual_output   = buffer.str();
  REQUIRE(actual_output == expected_output);
};

TEST_CASE("Vector Printing with Complex")
{
  ODE::Vector<std::complex<double>> vj(3);
  std::complex<double>              my_vector_1[3] = {0, 0, 0};
  my_vector_1[0]                                   = std::complex<double>(1.0, 1.0);
  my_vector_1[1]                                   = std::complex<double>(1.0, 2.0);
  my_vector_1[2]                                   = std::complex<double>(3.0, 4.0);
  vj.set_elements(my_vector_1);

  std::stringstream buffer;
  std::streambuf*   old = std::cout.rdbuf(buffer.rdbuf());

  vj.print_array();

  std::cout.rdbuf(old);

  std::string expected_output = "The array is: [(1,1) (1,2) (3,4)]\n";
  std::string actual_output   = buffer.str();
  REQUIRE(actual_output == expected_output);
};

TEST_CASE("Scaler Function")
{
  ODE::Vector<double> v(3);
  ODE::Vector<double> v_scaled(3);
  double              my_vector[3]   = {8, 4, 6};
  double              my_vector_2[3] = {80, 40, 60};
  v.set_elements(my_vector);
  v_scaled.set_elements(my_vector_2);

  ODE::Vector<double> v2(3);
  ODE::Scaler<double> s(10);
  s.value(v2, v, 2);

  for (int i = 0; i < 3; ++i)
    CHECK(v2[i] == v_scaled[i]);
};

TEST_CASE("Testing Linear Combination Function")
{
  ODE::Vector<double> scalars_test(3);
  double              scalars_test_array[3] = {1, 2, 5};
  scalars_test.set_elements(scalars_test_array);

  ODE::Vector<double> sol(3);
  double              sol_array[3] = {200, 400, 600};
  sol.set_elements(sol_array);

  std::cout << "scalars_test" << std::endl;
  scalars_test.print_array();

  ODE::Function_Test<double> s1(10); // f_1 = [10, 20, 30] , 10*1 + 20*2 + 30*5 = 200
  ODE::Function_Test<double> s2(20); // f_2 = [20, 40, 60] , 20*1 + 40*2 + 60*5 = 400
  ODE::Function_Test<double> s3(30); // f_3 = [30, 60, 90] , 30*1 + 60*2 + 90*5 = 600

  ODE::Function<double>* functions[3] = {&s1, &s2, &s3}; // 3 functions of 3 elements each

  ODE::Linear_Comb<double> l(functions, scalars_test);

  ODE::Vector<double> l_value(3);
  l.value(l_value, scalars_test, 1);
  l_value.print_array();

  for (int i = 0; i < 3; ++i)
    CHECK(l_value[i] == sol[i]);
};

TEST_CASE("Solver tests for `y' = sin(t)` with `y(0) = 1`")
{
  double t0 = 0;
  double tf = 1;
  double h  = 0.1;

  ODE::Vector<double> y0_test_2(1);
  double              y0_test_2_array[1] = {1};
  y0_test_2.set_elements(y0_test_2_array);

  class Test_2_Function : public ODE::Function<double>
  {
  public:
    void
    value(ODE::Vector<double>& output, const ODE::Vector<double>& input, double time)
    {
      double res_array[1] = {std::sin(time)};
      output.set_elements(res_array);
    };
  };
  SUBCASE("Runge Kutta 4")
  {
    Test_2_Function test_2_function;

    // std::stringstream buffer; // for testing output with prints
    // std::streambuf*   old = std::cout.rdbuf(buffer.rdbuf());

    ODE::RungeKutta4::RungeKuttaSolver<double> rk4_solver_test_2(
      &test_2_function, t0, tf, h, y0_test_2);

    CHECK(rk4_solver_test_2.y[0] < 1.4599698 + 0.1);
    CHECK(rk4_solver_test_2.y[0] > 1.4599698 - 0.1);
    // 1.4599698 is the exact solution, 0.1 tolerance is used to account for floating point errors

    // std::cout.rdbuf(old);

    // std::string expected_output = "The array is: [1.4597]\n";
    // std::string actual_output   = buffer.str();
    // REQUIRE(actual_output == expected_output);
  };

  SUBCASE("Euler")
  {
    Test_2_Function test_2_function;

    ODE::Euler::EulerSolver<double> euler_solver_test_2(&test_2_function, t0, tf, h, y0_test_2);

    CHECK(euler_solver_test_2.y[0] < 1.4599698 + 0.1);
    CHECK(euler_solver_test_2.y[0] > 1.4599698 - 0.1);
  };

  SUBCASE("Adams Bashforth with RK4")
  {
    Test_2_Function test_2_function;

    ODE::AdamsBashforth::AdamsBashforthSolver<double> ab_solver(
      &test_2_function, t0, tf, h, y0_test_2, "RK4");

    CHECK(ab_solver.y2[0] < 1.4599698 + 0.1);
    CHECK(ab_solver.y2[0] > 1.4599698 - 0.1);
  };

  SUBCASE("Adams Bashforth with Euler")
  {
    Test_2_Function test_2_function;

    ODE::AdamsBashforth::AdamsBashforthSolver<double> ab_solver(
      &test_2_function, t0, tf, h, y0_test_2, "Euler");

    CHECK(ab_solver.y2[0] < 1.4599698 + 0.1);
    CHECK(ab_solver.y2[0] > 1.4599698 - 0.1);
  };
};

TEST_CASE("Solver tests for `y' = -y` with `y(0) = 1`")
{
  double t0 = 0;
  double tf = 1;
  double h  = 0.1;

  ODE::Vector<double> y0_test_1(1);
  double              y0_test_1_array[1] = {1};
  y0_test_1.set_elements(y0_test_1_array);

  class Test_1_Function : public ODE::Function<double>
  {
  public:
    void
    value(ODE::Vector<double>& output, const ODE::Vector<double>& input, double time)
    {
      double res_array[1] = {-input[0]};
      output.set_elements(res_array);
    };
  };

  SUBCASE("Runge Kutta 4")
  {
    Test_1_Function test_1_function;

    ODE::RungeKutta4::RungeKuttaSolver<double> rk4_solver_test_1(
      &test_1_function, t0, tf, h, y0_test_1);

    CHECK(rk4_solver_test_1.y[0] < 0.36787944117144233 + 0.1);
    CHECK(rk4_solver_test_1.y[0] > 0.36787944117144233 - 0.1);
  };

  SUBCASE("Euler")
  {
    Test_1_Function test_1_function;

    ODE::Euler::EulerSolver<double> euler_solver_test_1(&test_1_function, t0, tf, h, y0_test_1);

    CHECK(euler_solver_test_1.y[0] < 0.36787944117144233 + 0.1);
    CHECK(euler_solver_test_1.y[0] > 0.36787944117144233 - 0.1);
  };

  SUBCASE("Adams Bashforth with RK4")
  {
    Test_1_Function test_1_function;

    ODE::AdamsBashforth::AdamsBashforthSolver<double> ab_solver(
      &test_1_function, t0, tf, h, y0_test_1, "RK4");

    CHECK(ab_solver.y2[0] < 0.36787944117144233 + 0.1);
    CHECK(ab_solver.y2[0] > 0.36787944117144233 - 0.1);
  };

  SUBCASE("Adams Bashforth with Euler")
  {
    Test_1_Function test_1_function;

    ODE::AdamsBashforth::AdamsBashforthSolver<double> ab_solver(
      &test_1_function, t0, tf, h, y0_test_1, "Euler");

    CHECK(ab_solver.y2[0] < 0.36787944117144233 + 0.1);
    CHECK(ab_solver.y2[0] > 0.36787944117144233 - 0.1);
  };
};

TEST_CASE("Solver tests for multidimensional vector of functions")
{
  double t0 = 0;
  double tf = 1;
  double h  = 0.1;

  ODE::Vector<double> y0(3);
  double              y0_array[3] = {1, 2, 3};
  y0.set_elements(y0_array);

  class Linear_Function : public ODE::Function<double>
  {
  public:
    void
    value(ODE::Vector<double>& output, const ODE::Vector<double>& input, double time)
    {
      double res_array[3] = {input[0] * time * 0.1, input[1] * time * 0.1, input[2] * time * 0.1};
      output.set_elements(res_array);
    };
  };

  SUBCASE("Runge Kutta 4")
  {
    Linear_Function linear_function;

    ODE::RungeKutta4::RungeKuttaSolver<double> rk4_solver(&linear_function, t0, tf, h, y0);

    CHECK(rk4_solver.y[0] < 1.0512711 + 0.1);
    CHECK(rk4_solver.y[0] > 1.0512711 - 0.1);
    CHECK(rk4_solver.y[1] < 2.10254219 + 0.1);
    CHECK(rk4_solver.y[1] > 2.10254219 - 0.1);
    CHECK(rk4_solver.y[2] < 3.15381329 + 0.1);
    CHECK(rk4_solver.y[2] > 3.15381329 - 0.1);
  };

  SUBCASE("Euler")
  {
    Linear_Function linear_function;

    ODE::Euler::EulerSolver<double> euler_solver(&linear_function, t0, tf, h, y0);

    CHECK(euler_solver.y[0] < 1.0512711 + 0.1);
    CHECK(euler_solver.y[0] > 1.0512711 - 0.1);
    CHECK(euler_solver.y[1] < 2.10254219 + 0.1);
    CHECK(euler_solver.y[1] > 2.10254219 - 0.1);
    CHECK(euler_solver.y[2] < 3.15381329 + 0.1);
    CHECK(euler_solver.y[2] > 3.15381329 - 0.1);
  };

  SUBCASE("Adams Bashforth with RK4")
  {
    Linear_Function linear_function;

    ODE::AdamsBashforth::AdamsBashforthSolver<double> ab_solver(
      &linear_function, t0, tf, h, y0, "RK4");

    CHECK(ab_solver.y2[0] < 1.0512711 + 0.1);
    CHECK(ab_solver.y2[0] > 1.0512711 - 0.1);
    CHECK(ab_solver.y2[1] < 2.10254219 + 0.1);
    CHECK(ab_solver.y2[1] > 2.10254219 - 0.1);
    CHECK(ab_solver.y2[2] < 3.15381329 + 0.1);
    CHECK(ab_solver.y2[2] > 3.15381329 - 0.1);
  };

  SUBCASE("Adams Bashforth with Euler")
  {
    Linear_Function linear_function;

    ODE::AdamsBashforth::AdamsBashforthSolver<double> ab_solver(
      &linear_function, t0, tf, h, y0, "Euler");

    CHECK(ab_solver.y2[0] < 1.0512711 + 0.1);
    CHECK(ab_solver.y2[0] > 1.0512711 - 0.1);
    CHECK(ab_solver.y2[1] < 2.10254219 + 0.1);
    CHECK(ab_solver.y2[1] > 2.10254219 - 0.1);
    CHECK(ab_solver.y2[2] < 3.15381329 + 0.1);
    CHECK(ab_solver.y2[2] > 3.15381329 - 0.1);
  };
};