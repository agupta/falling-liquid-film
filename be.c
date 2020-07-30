#include "grid/multigrid1D.h"
#include "poisson.h"
//#include "view.h"

// Physical
#define rhoWater (998.2071)
#define muWater (1.0016e-3)
#define pAir (101325)
#define gamma (72.86e-3)
#define g (9.80665)

// Problem specific
#define h_s (125e-6)
#define theta (pi/4.)

// scales
#define L (h_s)
#define U ((rhoWater * g * sq(h_s) * sin(theta)) / (2 * muWater))
#define T (L/U)
#define P ((muWater * U) / L)

// dimensionless constants
#define Re ((rhoWater * h_s * U)/muWater)
#define Ca ((muWater*U)/gamma)

void output_precomputed(FILE * file) {
  // 17 is a magic number
  fprintf(file, "Scales:\n")
  fprintf(file, "L: %.17g\n", L);
  fprintf(file, "U: %.17g\n", U);
  fprintf(file, "T: %.17g\n", T);
  fprintf(file, "P: %.17g\n\n", P);

  fprintf(file, "Dimensionless quantities:")
  fprintf(file, "Re: %.17g\n", Re);
  fprintf(file, "Ca: %.17g\n", Ca);
}

int main () {
  output_precomputed(stdout);
}
