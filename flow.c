#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

#define MAX_LEVEL 5
#define TIME_STEP 1
#define MAX_TIME 300

#define theta pi/4. // angle the ramp is to flat

#define beta 0.5

#define rhoRatio 1./784. // todo: compute from rhoWater, rhoAir
#define muRatio 1./50.

#define Fr 1./sqrt(9.81)
//#define We 1./0.07286
#define We 1./0.0003

int main () {
  init_grid(1 << MAX_LEVEL);

  rho1 = 1.;
  mu1 = 1.;

  rho2 = rhoRatio;
  mu2 = muRatio;

  // Surface tension coefficient
  f.sigma = 1./We;

  size(32);

  run();
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] += sin(theta)/sq(Fr); // along the ramp
  foreach_face(y)
    av.y[] -= cos(theta)/sq(Fr); // perpendicular to ramp
}

event init(t = 0) {
  mask (y > 1 ? top : none);

  // Flow out of the right -> into the left.
  periodic(right);

  // positive -> fluid 1
  // negative -> fluid 2 
  // bottom beta% is oil.
  fraction (f, beta-y + 0.04*sin(4*pi*x));

  // Initially velocity is 0 everywhere.
  foreach () {
    u.x[] = 0.;
    u.y[] = 0.;
  }
  
  // No-slip boundary conditions.
  u.t[top] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
}

event animationFluids(t += TIME_STEP; t <= MAX_TIME) {
  view(tx=-0.5, ty=-0.5);
  clear();
  squares("f", spread = -1, map=cool_warm, linear = true, min=0, max=MAX_LEVEL);
  draw_vof("f", lc = {0.0,0.0,0.0}, lw=1);
  // box();
  save ("Fluids.mp4");
}

event logfile (t += TIME_STEP; t <= MAX_TIME) {
  // log to file
  static FILE * fp = fopen("flow.log", "a");
  fprintf (fp, "%g %d %g\n", t, i, dt);
  fclose(fp);
}