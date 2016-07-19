/*
   non-uniform finite_difference implementation 1, based on: https://github.com/tuxication/nufd
*/


#include <iostream>

#define my_assert2(x, msg) {if(!(x)) std::cout << "Assersion error " << std::endl << msg << std::endl; exit(1); }

#include <vector>
/*
#include <algorithm>
#include <iomanip>
#include <random>
#include <chrono>
*/
#include <cassert>


using namespace std;

std::vector<double> fdcoef(unsigned int mord, unsigned int nord, double x0, const std::vector<double>::const_iterator grid);
std::vector<double> fd(unsigned int m, unsigned int n, const std::vector<double>& grid, const std::vector<double>& u);

// header end



std::vector<double> fdcoef(unsigned int mord, const unsigned int nord, double x0, const std::vector<double>::const_iterator grid) {
  // this routine implements simple recursions for calculating the weights
  // of finite difference formulas for any order of derivative and any order
  // of accuracy on one-dimensional grids with arbitrary spacing.

  // from Bengt Fornberg's article
  // generation of finite difference formulas on arbitrary spaced grids.
  // math. comp., 51(184):699-706, 1988.

  // input:
  // mord       = the order of the derivative
  // nord       = order of accuracy n
  // x0         = point at which to evaluate the coefficients
  // grid[nord] = array containing the grid starting at the lowest bound
  //              use during finite difference scheme

  // output:
  // coef[nord] = coefficients of the finite difference formula
  std::vector<double> coef(nord, 0.0);

  // local variables
  //constexpr
  int nmmin(min(nord,mord));
  double c1, c2, c3, c4, alpha;

  // more precision for weight calculations results
  // in a smaller error on output coefficients
  long double weight[nmmin+1][nord][nord];
  for (int i(0); i<nmmin+1; i++)
    for (int j(0); j<nord; j++)
      for (int k(0); k<nord; k++)
        weight[i][j][k] = 0.0;

  // recursive algorithm implementation
  weight[0][0][0] = 1.0;
  c1 = 1.0;
  for (int nn(1); nn<nord; nn++) {
    c2 = 1.0;
    for (int nu(0); nu<nn; nu++) {
      c3 = grid[nn] - grid[nu];
      c2 = c2*c3;
      c4 = 1.0/c3;
      alpha = grid[nn] - x0;
      weight[0][nn][nu] = c4 * (alpha * weight[0][nn-1][nu]);

      for (int mm(1); mm<nmmin+1; mm++)
        weight[mm][nn][nu] = c4 * (alpha * weight[mm][nn-1][nu] - mm * weight[mm-1][nn-1][nu]);
    }
    alpha = grid[nn-1] - x0;
    weight[0][nn][nn] = c1/c2 * (-alpha*weight[0][nn-1][nn-1]);
    c4 = c1/c2;

    for (int mm(1); mm<nmmin+1; mm++)
      weight[mm][nn][nn] = c4 * (mm * weight[mm-1][nn-1][nn-1] - alpha * weight[mm][nn-1][nn-1]);

    c1 = c2;
  }

  // load the coefficients
  for (int nu(0); nu<nord; nu++)
    coef[nu] = double(weight[mord-1][nord-1][nu]);

  return coef;
}

std::vector<double> fd(unsigned int m, unsigned int n, const std::vector<double>& grid, const std::vector<double>& u) {
  // this routine computes the order m derivatives
  // using n points on an arbitrary grid

  // input:
  // m           1=value, 2=1st diff, 3=2nd diff, 4=3rd diff, ...
  // n           = number of points use in fd schemes
  // grid[ngrid] = array of independent values
  // u[ngrid]    = function values at the grid points

  // output:
  // du[ngrid]   = first derivative values at the grid points
  size_t ngrid(grid.size());
  std::vector<double> du(ngrid, 0.0);
  std::vector<double> coef(n, 0.0);

  // validate the size of the grid and number of points
  // used in the finite difference scheme
  my_assert2(n<=ngrid, "");

  // use to point at the first element used
  // in the finite diff scheme to pass to fdcoef
  auto begin = grid.begin();
  auto end = grid.end();

  // number of forward and backward points
  int fb((n-1)/2);

  // beginning of the grid (forward differences)
  for (int i(0); i<fb; i++) {
    coef = fdcoef(m, n, grid[i], begin);
    for (int j(0); j<n; j++)
      du[i] = du[i] + coef[j]*u[j];
  }

  // middle of the grid (central differences)
  for (int i(fb); i<ngrid-fb; i++) {
    coef = fdcoef(m, n, grid[i], begin+i-fb);
    for (int j(0); j<n; j++)
      du[i] = du[i] + coef[j]*u[i-fb+j];
  }

  // end of grid (backward differences)
  for (size_t i(ngrid-fb); i<ngrid; i++) {
    coef = fdcoef(m, n, grid[i], end-n);
    for (int j(0); j<n; j++)
      du[i] = du[i] + coef[j]*u[ngrid-n+j];
  }

  return du;
}


// the functions fss002 and fss04 were on the original test_fdcoef.f
// found http://cococubed.asu.edu/code_pages/fdcoef.shtml
// Only kept them if someone want to see the general method for fixed order
// There both combine into one algorithm inside fd() with variable diff order and accuracy
std::vector<double> fss002(const std::vector<double>& grid, const std::vector<double>& u) {
  // this routine computes second order accurate first derivatives
  // on an arbitrary grid

  // input:
  // m:          1=value, 2=1st diff, 3=2nd diff
  // n           = number of points use in fd schemes
  // grid[ngrid] = array of independent values
  // u[ngrid]    = function values at the grid points

  // output:
  // du[ngrid]   = first derivative values at the grid points
  size_t ngrid(grid.size());
  std::vector<double> du(ngrid, 0.0);

  // m: 1=value, 2=1st diff, 3=2nd diff
  // n: number of points use in fd schemes
  unsigned int m(2), n(5);
  std::vector<double> coef(n, 0.0);
  auto grid_it = grid.cbegin();

  // forward differences
  coef = fdcoef(m, n, grid[0], grid_it);
  du[0] = coef[0]*u[0] + coef[1]*u[1] + coef[2]*u[2];

  // middle of the grid; central differences
  for (size_t i(1); i<ngrid-1; i++) {
    coef = fdcoef(m, n, grid[i], grid_it+i-1);
    du[i] = coef[0] * u[i - 1] + coef[1] * u[i] + coef[2] * u[i + 1];
  }

  // backward differences
  coef = fdcoef(m, n, grid[ngrid-1], grid_it+ngrid-3);
  du[ngrid-1] = coef[0] * u[ngrid - 3] + coef[1] * u[ngrid - 2] + coef[2] * u[ngrid-1];

  return du;
}

std::vector<double> fss004(const std::vector<double>& grid, const std::vector<double>& u) {
  // this routine computes fourth order accurate second derivatives
  // on an arbitrary grid

  // input:
  // n           = number of points use in fd schemes
  // ngrid       = number of points in the grid,
  // grid[ngrid] = array of independent values
  // u[ngrid]    = function values at the grid points

  // output:
  // du[ngrid]   = second derivative values at the grid points
  size_t ngrid(grid.size());
  std::vector<double> du(ngrid, 0.0);

  // local variables
  // m: 1=value, 2=1st diff, 3=2nd diff
  // n: number of points use in fd schemes
  unsigned int m(3), n(5);
  std::vector<double> coef(n, 0.0);
  auto grid_it = grid.cbegin();

  // first point; one sided
  // evaluation point + 4 to the right
  coef = fdcoef(m, n, grid[0], grid_it);
  for (size_t j(0); j<n; j++)
    du[0] = du[0] + coef[j]*u[j];

  // second point; one sided
  // 1 point to the left + evaluation point + 3 to the right
  coef = fdcoef(m, n, grid[1], grid_it);
  for (size_t j(0); j<n; j++)
    du[1] = du[1] + coef[j]*u[j];

  // middle of the grid; central differences
  for (size_t i(2); i<ngrid-2; i++) {
    coef = fdcoef(m, n, grid[i], grid_it+i-2);

    for (size_t j(0); j<n; j++)
      du[i] = du[i] + coef[j]*u[i+j-2];
  }

  // second to last point; one sided
  // 4 to the left + evaluation point
  coef = fdcoef(m, n, grid[ngrid-2], grid_it+ngrid-n);
  for (size_t j(0); j<n; j++)
    du[ngrid-2] = du[ngrid-2] + coef[j]*u[ngrid+j-n];

  // last point; one sided
  // 3 to the right + evaluation point + 1 point to the right
  coef = fdcoef(m, n, grid[ngrid-1], grid_it+ngrid-n);
  for (size_t j(0); j<n; j++)
    du[ngrid-1] = du[ngrid-1] + coef[j]*u[ngrid+j-n];

  return du;
}
