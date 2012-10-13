/*
    Copyright (C) 2012  Alexander Lapin

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
	Contacts:
	name: Alexander Lapin
	e-mail: lapinra@gmail.com
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using namespace std;

const double MAX_NORM = 1000;
const double coef[] = {-1.1, 0.2, 4.7};

double alpha = coef[0];

double P (double x, double y)
{
  return x + y;
}

double Q (double x, double y)
{
  return (alpha - 1) * y + (alpha - 2) * x - x * x * x - x * x * y;
}

double P_n025 (double x, double y, double tau)
{
  return P(x + 0.5 * tau * P(x, y), y + 0.5 * tau * Q(x, y));
}

double Q_n025 (double x, double y, double tau)
{
  return Q(x + 0.5 * tau * P(x, y), y + 0.5 * tau * Q(x, y));
}


double P_n05 (double x, double y, double tau)
{
  return P(x + 0.5 * tau * P_n025(x, y, tau), y + 0.5 * tau * Q_n025(x, y, tau));
}

double Q_n05 (double x, double y, double tau)
{
  return Q(x + 0.5 * tau * P_n025(x, y, tau), y + 0.5 * tau * Q_n025(x, y, tau));
}


double P_n075 (double x, double y, double tau)
{
  return P(x + tau * P_n05(x, y, tau), y + tau * Q_n05(x, y, tau));
}

double Q_n075 (double x, double y, double tau)
{
  return Q(x + tau * P_n05(x, y, tau), y + tau * Q_n05(x, y, tau));
}


double next_x (double x, double y, double tau)
{
  return x + (1.0 / 6.0) * tau * (P(x, y) + 2 * P_n025(x, y, tau) + 2 * P_n05(x, y, tau) + P_n075(x, y, tau));
}

double next_y (double x, double y, double tau)
{
  return y + (1.0 / 6.0) * tau * (Q(x, y) + 2 * Q_n025(x, y, tau) + 2 * Q_n05(x, y, tau) + Q_n075(x, y, tau));
}

// Intersect segment ({x1, y1}, {x2, y2}) with line {y = 0}
double get_intersection (double x1, double y1, double x2, double y2)
{
  return x1 - y1 * (x2 - x1) / (y2 - y1);
}

// and check that it is intersection with ray {x > 0, y = 0}
bool check_intersection (double x1, double x2, double x)
{
  return (x > 0.0 && x >= min (x1, x2) && x <= max (x1, x2));
}

double norm (double x, double y)
{
  return x * x + y * y;
}

FILE *out;

void find_period (double x0, double y0, double tau, FILE *out)
{
  double x, y, x_pr, y_pr, new_x, new_y;
  double xt;
  double half_tau;
  double x_half_tau, y_half_tau;
  double x_intersection;
  int k;
  double error;

  k = 0;
  error = 0.0;

  x = x_half_tau = x0;
  y = y_half_tau = y0;

  xt = -1.0; // Indicates that no intersections were found before

  if (tau > 1e-2 - 1e-6 && tau < 1e-2 + 1e-6)
    fprintf (out, "%.5lf %.5lf\n", x, y);

  while (true)
    {
      x_pr = x;
      y_pr = y;

      x = next_x (x_pr, y_pr, tau);
      y = next_y (x_pr, y_pr, tau);

      fprintf (out, "%.5lf %.5lf\n", x, y);

      k++;

      new_x = next_x (x_half_tau, y_half_tau, tau * 0.5);
      new_y = next_y (x_half_tau, y_half_tau, tau * 0.5);
      x_half_tau = new_x;
      y_half_tau = new_y;

      new_x = next_x (x_half_tau, y_half_tau, tau * 0.5);
      new_y = next_y (x_half_tau, y_half_tau, tau * 0.5);
      x_half_tau = new_x;
      y_half_tau = new_y;

      error = max (error, norm (x_half_tau - x, y_half_tau - y));

      if (norm (x, y) > MAX_NORM)
        {
          cout << "Not periodic: too far" << endl;
          break;
        }

      if (norm (x, y) < max (4 * error, 1e-16))
        {
          cout << "Not periodic: too close" << endl;
          break;
        }

      x_intersection = get_intersection (x_pr, y_pr, x, y);

      if (check_intersection (x_pr, x, x_intersection))
        {
          
          if (norm (xt - x_intersection, 0.0) < max (4 * error, 1e-16))
            {
              printf ("Period: %.5lf\nError: %.3le\n", k * tau, sqrt (error));
              break;
            }

          xt = x_intersection;
          k = 0;
        }
    }
}

int main ()
{
  cout << "Input alpha:" << endl;
  cin >> alpha;
  
  cout << "Input x0, y0:" << endl;
  double x0, y0;
  cin >> x0 >> y0;

  cout << "Input tau:" << endl;
  double tau;
  cin >> tau;
  
  out = fopen ("output.txt", "w");
  
  find_period (x0, y0, tau, out);
  
  fclose(out);
 
  return 0;
}

