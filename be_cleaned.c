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

//  scalar dh[], hx[], hxx[], hxxx[], q[], qx[]; // e.g. hx corresponds to \partial_x h

// explicit BE solver loop
void solve_explicit (scalar h, scalar F, double dt) {
  scalar dh[], hx[], hxx[], hxxx[], q[], qx[]; // e.g. hx corresponds to \partial_x h
  foreach () {
    // Central difference approximations.
    hx[] = (h[1] - h[-1])/(2.*Delta);
    hxx[] = (h[1] - 2.*h[] + h[-1])/sq(Delta);
    hxxx[] = (h[2] - 2.*h[1] + 2.*h[-1] - h[-2])/(2.*cube(Delta));
  }
  foreach()
    q[] = 2.*cube(h[])/3. - (cube(h[])/3.)*(2.*hx[]/tan(theta) - hxxx[]/Ca) + Re*(8.*sq(cube(h[]))*hx[]/15. - 2.*sq(sq(h[]))*F[]/3.);
  foreach()
    qx[] = (q[1] - q[-1])/(2.*Delta);
  foreach()
    dh[] = F[] - qx[];
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
  init_grid(256);
  size(width);
  periodic(right);
  scalar h[], F[];
  foreach() {
    h[] = 1. + 1e-2*sin(2.*pi*x/width);
    F[] = 0.;
  }
  boundary ({h, F});
  
  double dt = 1e-7;
  int i = 0;
  FILE * fp = fopen("out/explicit.out", "w");
  for (double t = 0; t <= MAX_TIME; t += dt, i++) {
    if (i % (int)(MAX_TIME/(dt*100)) == 0) { // ~1000 slices
      foreach()
	      fprintf (fp, "%g %g %g\n", t, x, h[]);
      fputs ("\n", fp);
    }
    if (i % (int)(MAX_TIME/(dt*5)) == 0) { // indicate progress 5 + 1 times
      fprintf(stdout, "t = %g\n", t);
    }
    solve_explicit (h, F, dt);
  }
  fclose(fp);
  // TODO: print time. see utils.h
}
