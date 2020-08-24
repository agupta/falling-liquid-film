#include "grid/multigrid1D.h"
#include "poisson.h"
//#include "view.h"

// Physical
#define rhoWater (998.2071)
#define muWater (1.0016e-3)
#define pAir (101325.)
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

// explicit BE solver loop
void solve_explicit (scalar h, scalar F, double dt) {
  scalar dh[], hx[], hxx[], hxxx[], hxxxx[], Fx[]; // e.g. hx corresponds to \partial_x h
  foreach () {
    // Central difference approximations.
    hx[] = (h[1] - h[-1])/(2.*Delta);
    hxx[] = (h[1] - 2.*h[] + h[-1])/sq(Delta);
    hxxx[] = (h[2] - 2.*h[1] + 2.*h[-1] - h[-2])/(2.*cube(Delta));
    hxxxx[] = (h[2] - 4.*h[1] + 6.*h[] - 4.*h[-1] + h[-2])/(sq(sq(Delta)));
    Fx[] = (F[1] - F[-1])/(2.*Delta);
  }
  foreach()
    dh[] = F[]
           - sq(h[])*(hx[]*(2. - 2./tan(theta)) + hxxx[]/Ca)
           + (cube(h[])/3.)*(2.*hxx[]/tan(theta) - hxxxx[]/Ca)
           + Re*(8.*cube(h[])*hx[]*F[]/3. + 2.*sq(sq(h[]))*Fx[]/3. - 16.*cube(h[])*sq(h[])*sq(hx[])/5. - 8.*sq(cube(h[]))*hxx[]/15.);
  // foreach()
  //   printf("%g %g %g %g %g %g %g %g\n", h[], hx[], hxx[], hxxx[], hxxxx[], Fx[], dh[], dh[]*dt);
  foreach() 
    h[] += dt*dh[];
  boundary ({h});
}

void output_precomputed(FILE * file) {
  // 17 is a magic number
  fprintf(file, "Scales:\n");
  fprintf(file, "L: %.17g\n", L);
  fprintf(file, "U: %.17g\n", U);
  fprintf(file, "T: %.17g\n", T);
  fprintf(file, "P: %.17g\n\n", P);

  fprintf(file, "Dimensionless quantities:\n");
  fprintf(file, "Re: %.17g\n", Re);
  fprintf(file, "Ca: %.17g\n\n\n", Ca);
}

int main () {
  output_precomputed(stdout);
  init_grid(512);
  periodic(right);
  scalar h[], F[];
  foreach() {
    h[] = 1. + 1e-3*sin(2.*pi*x);
    F[] = 1e-1*cos(16.*pi*x);
  }
  boundary ({h, F});

  double dt = 5e-15;
  int i = 0;
  FILE * fp = fopen("explicit.out", "w");
  for (double t = 0; t <= 1e-8; t += dt, i++) {
    if (i % 10000 == 0) {
     foreach()
	     fprintf (fp, "%g %g %g\n", t, x, h[]);
     fputs ("\n", fp);
    }
    if (i % 100000 == 0) {
      fprintf(stdout, "t = %g\n", t);
    }
    //solve_explicit (h, F, dt);
  }
  fclose(fp);
}
