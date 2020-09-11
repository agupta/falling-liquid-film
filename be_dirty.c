#include "grid/multigrid1D.h"
#include "poisson.h"
#include "utils.h"
#include "output.h"
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
#define width (64.)

// scales
#define L (h_s)
#define U ((rhoWater * g * sq(h_s) * sin(theta)) / (2 * muWater))
#define T (L/U)
#define P ((muWater * U) / L)

// dimensionless constants
#define Re ((rhoWater * h_s * U)/muWater)
#define Ca ((muWater*U)/gamma)

// simulation times
#define MAX_TIME (1e0)

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
           - sq(h[])*hx[]*((2. - 2.*hx[]/tan(theta)) + hxxx[]/Ca)
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
  output_precomputed(stderr);
  init_grid(128);
  size(width);
  periodic(right);
  scalar h[], F[];
  foreach() {
    h[] = 1. + 1e-2*sin(2.*pi*x/width);
    F[] = 0.;
  }
  boundary ({h, F});
  
  // here is the point at which we lose mass conservation at 512
  //double dt = 6e-8;
  // here is the point at which we lose mass conservation at 128
  //double dt = 1.6e-5;

  double dt = 1e-5;
  
  // int i = 0;
  double t = 0;
  FILE * fp = fopen("out/explicit.out", "w");

  fprintf(stderr, "start!\n");

  //for (double t = 0; t <= MAX_TIME; t += dt, i++) {
  for (int i = 0; i <= 1000000; t += dt, i++) {
    if (i % 1000 == 0) { // 100 + 1 slices
      foreach()
	      fprintf (fp, "%g %g %g %g\n", t, x, h[], F[]);
      fputs ("\n", fp);
    }
    if (i % 100000 == 0) { // 10 + 1 times
      // fprintf(stdout, "t = %g\n", t); // convenience progress
      fputs ("h >", stdout);
      foreach()
	      fprintf (stdout, "%g ", h[]);
      fputs ("\n", stdout);
      if (i != 100000) { // not the last time
        fputs("F >", stdout);
        foreach()
          scanf("%lf", &F[]);
      }
    }
    solve_explicit (h, F, dt);
  }
  fclose(fp);

  fprintf(stderr, "finish!\n");

  // TODO: print time. see utils.h
}