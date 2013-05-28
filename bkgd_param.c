#include <glib.h>
#include <stdio.h>
#include <string.h>

#include <util/config.h>
#include "bkgd_param.h"
#include "bkgd.h"
#include "interp_tab.h"
#include "bkgd_intg.h"


#define BKGD_PARAM_MIN_T_DIST 1e-6
#define BKGD_PARAM_B_STEP 0.001
#define BKGD_PARAM_A_STEP 0.001

void adjust_t_upbound(BkgdParam *p) {
  double f;

  f = p->t_dist(p->b, p);
  while(f < BKGD_PARAM_MIN_T_DIST && p->b > p->a) {
    p->b -= BKGD_PARAM_B_STEP;
    f = p->t_dist(p->b, p);
  }
  if(p->b < p->a) {
    g_error("adjust_t_upbound: could not find reasonable upper bound on t");
  }
}


void adjust_t_lowbound(BkgdParam *p) {
  double f;

  f = p->t_dist(p->a, p);
  while(f < BKGD_PARAM_MIN_T_DIST && p->a < p->b) {
    p->a += BKGD_PARAM_B_STEP;
    f = p->t_dist(p->a, p);
  }
  if(p->a > p->b) {
    g_error("adjust_t_upbound: could not find reasonable lower bound on t");
  }
}



/**
 * The exponential distribution we use is truncated.
 * It is defined as:
 *
 *   f(t) = C e^(-lambda*t)  for a <= t <= b
 *   f(t) = 0                otherwise
 *
 * We need to find a lambda that gives mean t for this truncated
 * distribution, and C that normalizes the distribution so the
 * integral over its range is 1.
 *
 * Use this definition of mean to find value:
 * 
 */
static void calc_exp_dist_param(BkgdParam *parm) {
  double upper, lower, lambda, mean, c, a, b, err;
  double upper_bound;

  err = -1.0;
  a = parm->a;
  b = parm->b;

  fprintf(stderr, "estimating C and lambda for truncated exponential "
	  "distribution:\n  a=%g <= t <= b=%g, mean t=%g\n", a, b, parm->t);

  /* start with 1/t as upper bound, and move it upwards
   * by doubling as needed
   */
  lower = 0.0;
  upper = 1.5/parm->t;
  mean = (1.0 + a * exp(-upper*a) - b*exp(-upper*b))/upper;
  lambda = upper;
  while(mean > parm->t) {
    lower = upper;
    upper *= 2.0;
    mean = (1.0 + a * exp(-upper*a) - b*exp(-upper*b))/upper;
  }

  /* do binary search for value that gives desired mean, taking into
   * account that function is monotonically decreasing with increasing
   * lambda
   */
  while(err < 0.0 || err > mean*BKGD_PARAM_EPS_REL) {
    /* next point is midpoint between bounds */
    lambda = (lower + upper) / 2;
    
    mean = (1.0 + a * exp(-lambda*a) - b*exp(-lambda*b))/lambda;

    err = fabs(mean - parm->t);

    if(mean > parm->t) {
      /* cur value is too low, so need higher lambda, 
       * set midpoint as new lower bound 
       */
      lower = lambda;
    } else {
      upper = lambda;
    }
  }

  c = lambda / (exp(-lambda*a) - exp(-lambda*b));
  parm->exp_lambda = lambda;
  parm->exp_c = c;
  
  fprintf(stderr, "  final estimates:\n  lambda=%g, c=%g, err=%g\n", 
	  lambda, c, err);


  /* integration gets into trouble once exp(-x) drops below minimum
   * representable value on machine, define upper bound on t to
   * account for this
   */
  upper_bound = BKGD_PARAM_EXP_LOW_BOUND / -lambda;
  if(upper_bound < parm->a) {
    g_error("calc_exp_dist_param: cannot integrate exponential\n"
	    "with mean %g and specified range [%g..%g] because\n"
	    "too close to 0. upper bound=%g\n",
	    parm->t, parm->a, parm->b, upper_bound);
  }

  if(upper_bound < parm->b) {
    fprintf(stderr, "  using b=%g as upper bound on range of integration "
	    "instead of %g\n", upper_bound, parm->b);
    parm->b = upper_bound;
  }
}





/*
 * integrand used to calculate expectation of truncated gamma distr
 */
double gamma_expect_integrand(double x, void *v) {
  double f;

  f = bkgd_t_dist_gamma(x, v);
  return x * f;
}


/**
 * Used to find scale parameter for truncated gamma distribution.
 * This function will return 0, when scale parameter gives desired
 * mean.
 */
double gamma_root_func(void *v) {
  BkgdIntg *bi = v;
  double expect, err;

  /* bi->p->gamma_scale = scale;*/
  
  /* calculate mean */

  bkgd_gsl_integration_wrapper(&bi->f, bi->p->a, bi->p->b, bi->w,
			       bi->p, &expect, &err);

  /* fprintf(stderr, "scale: %g, expect_t: %g, target_t: %g\n", 
   * bi->p->gamma_scale, expect, bi->p->t);
   */

  /* we set this up so that return value is 0 when 
   * the expectation (mean) of distribution is equal to
   * to the desired value 't'
   */
  return expect - bi->p->t;
}


/**
 * Finds normalization constant for truncated Gamma distribution.
 * (since distribution is truncated, it no longer sums to 1 without
 * this constant)
 */
static void set_gamma_dist_norm_const(BkgdParam *parm) {
  BkgdIntg b_intg;
  double intgl, err;

  /* now find normalizing constant to ensure PDF integral = 1 despite
   * truncation
   */
  b_intg.w = gsl_integration_workspace_alloc(BKGD_INTG_LIMIT);
  b_intg.f.params = parm;
  b_intg.p = parm;
  b_intg.f.function = &bkgd_t_dist_gamma;  


  parm->gamma_c = 1.0;
  bkgd_gsl_integration_wrapper(&b_intg.f, parm->a, parm->b, b_intg.w,
			       parm, &intgl, &err);

  parm->gamma_c = 1.0 / intgl;

  /* fprintf(stderr, "a:%g b:%g " */
  /* 	  "gamma_shape:%g gamma_scale:%g gamma_c:%g intgl:%g\n",  */
  /* 	  parm->a, parm->b, parm->gamma_shape, parm->gamma_scale, */
  /* 	  parm->gamma_c, intgl); */

  if(intgl > 1.001) {
    g_error("set_gamma_dist_norm_const: Numerical estimation of "
	    "normalization constant for truncated gamma distribution "
	    "failed. Integral of distribution should be < 1.0 but "
	    "got %g\n", intgl);
  }
  
  
  gsl_integration_workspace_free(b_intg.w);
}




static void calc_gamma_dist_param(BkgdParam *parm) {
  BkgdIntg b_intg;
  double low_brk, up_brk, mid, y, err;

  parm->gamma_c = 1.0;
  b_intg.w = gsl_integration_workspace_alloc(BKGD_INTG_LIMIT);
  b_intg.f.params = parm;
  b_intg.p = parm;

  /* we want to integrate this function to find expectation of
   * truncated gamma distribution, given a shape and scale parameters
   */
  b_intg.f.function = &gamma_expect_integrand;

  fprintf(stderr, "estimating scale for truncated gamma "
	  "distribution:\n  shape=%g, a=%g <= t <= b=%g, mean t=%g\n",
	  parm->gamma_shape, parm->a, parm->b, parm->t);

  /* fit scale parameter by searching for root of function */
  up_brk = 1.0;
  low_brk = 1e-4;
  
  /* find lower bracket that is on left side of root */
  /* fprintf(stderr, "finding low bracket\nlow_brk=%g\n", low_brk); */

  parm->gamma_scale = low_brk;
  set_gamma_dist_norm_const(parm);
  y = gamma_root_func(&b_intg);
  while(y > 0.0) {
    low_brk *= 0.5;
    parm->gamma_scale = low_brk;
    set_gamma_dist_norm_const(parm);
    y = gamma_root_func(&b_intg);
    /* fprintf(stderr, "low_brk=%g, y=%g\n", low_brk, y); */
  }

  /* find upper bracket that is on right side of root */
  /* fprintf(stderr, "finding upper bracket\n"); */
  parm->gamma_scale = up_brk;
  set_gamma_dist_norm_const(parm);
  y = gamma_root_func(&b_intg);
  while(y < 0.0) {
    up_brk *= 2.0;
    parm->gamma_scale = up_brk;
    set_gamma_dist_norm_const(parm);
    y = gamma_root_func(&b_intg);

    if(up_brk > 1e100) {
      g_error("calc_gamma_dist_param: unable to find a reasonable "
	      "scale parameter for truncated gamma distribution "
	      "with shape %g that will give mean(t) of %g. "
	      "Use a different shape parameter or a different mean t.",
	      parm->gamma_shape, parm->t);
    }

    /* fprintf(stderr, "up_brk=%g, y=%g\n", up_brk, y); */
  }

  /* fprintf(stderr, "low_brk=%g, up_brk=%g\n", low_brk, up_brk); */

  /* perform binary search for root, halving size of bracket each time */
  err = up_brk - low_brk;
  mid = (low_brk + up_brk) * 0.5;
  while(err > mid*BKGD_PARAM_EPS_REL) {
    parm->gamma_scale = mid;
    set_gamma_dist_norm_const(parm);
    y = gamma_root_func(&b_intg);
        
    if(y < 0.0) {
      low_brk = mid;
    } else {
      up_brk = mid;
    }

    mid = (low_brk + up_brk) * 0.5;
    err = up_brk - low_brk;
    /* fprintf(stderr, "low_brk=%g, up_brk=%g, err=%g\n",
     * low_brk, up_brk, err); */
  }

  fprintf(stderr, "  gamma_scale=%g, err=%g\n", 
	  parm->gamma_scale, err);

}







static void create_intg_tabs(BkgdParam *parm) {
  BkgdIntg b_intg;

  b_intg.w = gsl_integration_workspace_alloc(BKGD_INTG_LIMIT);
  b_intg.f.params = parm;
  b_intg.p = parm;

  fprintf(stderr, "creating interp tab for blk 2nd deriv integrals\n");
  b_intg.f.function = &bkgd_drv2_blk_integrand;
  parm->intg_tab_drv2_blk = 
    interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  fprintf(stderr, "creating interp tab for site 2nd deriv integrals\n");
  b_intg.f.function = &bkgd_drv2_site_integrand;
  parm->intg_tab_drv2_site = 
    interp_tab_create_1d(&bkgd_calc_site_integral, &b_intg);


  fprintf(stderr, "creating interp tab for blk integrals\n");
  b_intg.f.function = &bkgd_blk_integrand;
  parm->intg_tab_blk = interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  fprintf(stderr, "creating interp tab for site integrals\n");
  b_intg.f.function = &bkgd_site_integrand;
  parm->intg_tab_site = interp_tab_create_1d(&bkgd_calc_site_integral, 
					     &b_intg);

  fprintf(stderr, "creating interp tab for blk 1st deriv integrals\n");
  b_intg.f.function = &bkgd_drv1_blk_integrand;
  parm->intg_tab_drv1_blk = 
    interp_tab_create(&bkgd_calc_blk_integral, &b_intg);

  fprintf(stderr, "creating interp tab for site 1st deriv integrals\n");
  b_intg.f.function = &bkgd_drv1_site_integrand;
  parm->intg_tab_drv1_site = 
    interp_tab_create_1d(&bkgd_calc_site_integral, &b_intg);

  gsl_integration_workspace_free(b_intg.w);
}



/**
 * Retrieves parameters that are used for calculating B values from
 * the config and creates a new BkgdParam structure.
 */
BkgdParam *bkgd_param_new(Config *config) {
  BkgdParam *parm;
  char *dist_type;

  parm = g_new(BkgdParam,1);

  /* t=sh is selection coefficient times heterozygosity coef */
  parm->t = config_get_double(config, "PARAM_T");
  fprintf(stderr, "t=%g\n", parm->t);

  /* u is deleterious mutation rate */
  parm->u = config_get_double(config, "PARAM_U");

  fprintf(stderr, "u=%g\n", parm->u);

  dist_type = config_get_str(config, "PARAM_T_DIST_TYPE");

  if(config_get_boolean(config, "USE_SUM_APPROXIMATION")) {
    parm->apprx_sum = TRUE;
    /* convert B contrib thresh to sum thresh */
    parm->max_sum_thresh = log(1.0-BKGD_PARAM_MAX_SUM_THRESH) / -parm->u;
    fprintf(stderr, "using approx sum upper bound threshold=%g\n",
	    parm->max_sum_thresh);
  } else {
    parm->apprx_sum = FALSE;
    parm->max_sum_thresh = 0;
  }

  if(strcmp(dist_type, "POINT")==0) {
    parm->a = 0;
    parm->b = 0;
    parm->t_dist = NULL;
  }
  else if(strcmp(dist_type, "EXPONENTIAL")==0) {
    parm->t_dist = &bkgd_t_dist_exp;
    parm->a = config_get_double(config, "PARAM_T_DIST_TRUNC_LOW");
    parm->b = config_get_double(config, "PARAM_T_DIST_TRUNC_HIGH");
    calc_exp_dist_param(parm);
  }
  else if(strcmp(dist_type, "GAMMA")==0) {
    parm->t_dist = &bkgd_t_dist_gamma;
    parm->a = config_get_double(config, "PARAM_T_DIST_TRUNC_LOW");
    parm->b = config_get_double(config, "PARAM_T_DIST_TRUNC_HIGH");

    parm->gamma_shape = config_get_double(config, "PARAM_GAMMA_SHAPE");

    calc_gamma_dist_param(parm);
  }
  else {
    g_error("get_bkgd_param: PARAM_T_DIST_TYPE must be one of "
	    "POINT, EXPONENTIAL, GAMMA");
  }

  create_intg_tabs(parm);

  return parm;
}


/**
 * Frees memory associated with provided structure
 */
void bkgd_param_free(BkgdParam *param) {
  /* free interpolation tables if they are defined */
  if(param->intg_tab_blk) {
    interp_tab_free(param->intg_tab_blk);
  }
  if(param->intg_tab_site) {
    interp_tab_free(param->intg_tab_site);
  }
  if(param->intg_tab_drv1_blk) {
    interp_tab_free(param->intg_tab_drv1_blk);
  }
  if(param->intg_tab_drv1_site) {
    interp_tab_free(param->intg_tab_drv1_site);
  }
  if(param->intg_tab_drv2_blk) {
    interp_tab_free(param->intg_tab_drv2_blk);
  }
  if(param->intg_tab_drv2_site) {
    interp_tab_free(param->intg_tab_drv2_site);
  }

  g_free(param);
}
