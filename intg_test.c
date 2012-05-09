
#include <stdio.h>

#include "bkgd.h"
#include "bkgd_intg.h"
#include "bkgd_param.h"


#define STEP 0.001

int main(int argc, char **argv) {
  BkgdIntg b_intg;
  BkgdParam parm;
  double x, y, f;

  parm.t = 0.1;
  parm.u = 1.2e-8;
  parm.t_dist = &bkgd_t_dist_gamma;
  parm.gamma_shape = 0.75;
  parm.gamma_scale = 0.133742;
  parm.gamma_c = 1.00116;

  parm.a = 1e-5;
  parm.b = 1.0;
  

  b_intg.w = gsl_integration_workspace_alloc(BKGD_INTG_LIMIT);
  b_intg.f.params = &parm;
  b_intg.p = &parm;

  /* b_intg.f.function = &bkgd_drv1_blk_integrand;*/


  b_intg.f.function = &bkgd_drv2_blk_integrand; 

  parm.r_near = 0.0;
  parm.r_far  = 1.98364e-05;
  
  for(x = parm.a; x <= parm.b; x += STEP) {

    f = bkgd_t_dist_gamma(x, &parm);
    /* f = 0.0; */

    y = bkgd_drv1_blk_integrand(x, &parm);

    fprintf(stderr, "%.10g %g %g\n", f, x, y);
  }


  parm.intg_tab_drv2_blk = 
    interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  /*
  parm.intg_tab_drv1_blk = interp_tab_create(&bkgd_calc_blk_integral, &b_intg);
  */

  return 0;
}
