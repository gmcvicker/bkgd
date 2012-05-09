#include <glib.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

#include <util/config.h>
#include <util/util.h>

#include <gsl/gsl_multimin.h>

#include "bkgd_evo_mdl.h"

#define MAX_ITER 50000
#define MAX_SUBSET_ITER 50

#define MAX_SIMPLEX_SZ 0.005

#define INIT_STEP_SIZE 0.01
#define LINE_MIN_TOL 0.005
#define GRAD_TOL 0.005

#define LL_SMALL_DIFF 0.01
#define LL_LARGE_DIFF 0.1

#define MIN_TYPE_CONJUGATE_PR 1
#define MIN_TYPE_VECTOR_BFGS 2
#define MIN_TYPE_VECTOR_BFGS2 3
#define MIN_TYPE_SIMPLEX 4



/* squares a value */
static double sqr(double x) {
  return x*x;
}



static void write_gsl_vector(FILE *fh, const gsl_vector *x) {
  size_t i;
  char *sep;
  
  sep = "";
  fprintf(fh, "<");
  for(i = 0; i < x->size; i++) {
    fprintf(fh, "%s%g", sep, gsl_vector_get(x, i));
    sep=", ";
  }
  fprintf(fh, ">\n");

}



/**
 * Extracts set of unlocked parameters from the provided model and
 * creates a gsl_vector of parameters that can be passed to a GSL
 * minimizer. 
 */
gsl_vector *model_to_gsl_param(Model *mdl) {
  gsl_vector *x;
  int i, x_idx;

  x = gsl_vector_alloc(mdl->n_unlocked);

  x_idx = 0;
  for(i = 0; i < mdl->n_param; i++) {
    if(!mdl->param[i].locked) {

      if(x_idx >= mdl->n_unlocked) {
	g_error("model_to_gsl_param: expected %d unlocked paramters, "
		"but there are at least %d", mdl->n_unlocked, x_idx+1);
      }

      if(mdl->param[i].up_bound) {
	/* untransform parameters */
	double val = -log(mdl->param[i].max_val-mdl->param[i].val) * 
	  mdl->param[i].scale;
	gsl_vector_set(x, x_idx, val);
      } else {
	gsl_vector_set(x, x_idx, mdl->param[i].val * mdl->param[i].scale);
      }

      x_idx++;
    }
  }

  if(x_idx != mdl->n_unlocked) {
    model_write_param(mdl, stderr);
    g_error("mdl_to_gsl_param: expected %d unlocked parameters, "
	    "but only found %d", mdl->n_unlocked, x_idx);
  }

  return x;
}



/**
 * Extracts gradient consisting of derivatives of LL taken wrt each
 * unlocked model paramter and sets them in the provided gsl vector.
 */
void model_to_gsl_grad(Model *mdl, const gsl_vector *param, gsl_vector *grad) {
  int i, grad_idx;

  if(grad->size != mdl->n_unlocked) {
    g_error("model_to_gsl_grad: expected gradient size (%zd) to "
	    "equal number of unlocked parameters (%d)", grad->size,
	    mdl->n_unlocked);
  }

  grad_idx = 0;
  for(i = 0; i < mdl->n_param; i++) {
    if(!mdl->param[i].locked) {

      if(grad_idx >= mdl->n_unlocked) {
	g_error("model_to_gsl_grad: expected %d unlocked paramters, "
		"but there are at least %d", mdl->n_unlocked, grad_idx+1);
      }

      double drv;

      if(mdl->param[i].up_bound) {
	/* apply transformation for graceful upper bound */
	drv = -exp(mdl->param[i].val)*mdl->param[i].deriv/mdl->param[i].scale;
      } else {
	drv = -mdl->param[i].deriv/mdl->param[i].scale;
      }

      /* flip sign of derivative if we flipped sign of parameter */
      if(!mdl->param[i].allow_neg && gsl_vector_get(param, grad_idx) < 0.0) {
	if(mdl->param[i].val < 0.0) {
	  /* sanity check */
	  g_error("model_to_gsl_grad: expected param value (%g) to be "
		  "positive", mdl->param[i].val);
	}

	drv = -drv;
      }

      gsl_vector_set(grad, grad_idx, drv);
      grad_idx++;
    }
  }

  if(grad_idx != mdl->n_unlocked) {
    model_write_param_ln(mdl, stderr);
    model_write_grad_ln(mdl, stderr);
    g_error("mdl_to_gsl_grad: expected %d unlocked parameters, "
	    "but only found %d", mdl->n_unlocked, grad_idx);
  }

  return;
}





/**
 * Updates a model's unlocked paramters from the values in the
 * provided vector
 */
void model_update_param(Model *mdl, const gsl_vector *x) {
  int i, x_idx;
  double val;

  if(x->size != mdl->n_unlocked) {
    g_error("model_update_param: expected number of unlocked "
	    "model parameters (%d) to equal gsl_param size (%zd)",
	    mdl->n_unlocked, x->size);
  }

  x_idx = 0;
  for(i = 0; i < mdl->n_param; i++) {
    if(!mdl->param[i].locked) {
      /* copy value from gsl vector to unlocked parameter */
      if(mdl->param[i].up_bound) {
	/* apply transformation to limit parameter to upper bound */
	val = mdl->param[i].max_val - 
	  exp(-gsl_vector_get(x, x_idx) / mdl->param[i].scale);
      } else {
	val = gsl_vector_get(x, x_idx) / mdl->param[i].scale;
      }

      if((val < 0.0) && (!mdl->param[i].allow_neg)) {
	/* flip sign if negative */
	val = -val;
      }

      mdl->param[i].val = val;

      x_idx++;
    }
  }

  if(x_idx != x->size) {
    g_error("model_update_param: expected %zd unlocked parameters "
	    "but only found %d", x->size, x_idx);
  }
}



/**
 * Wrapper around LLH function called by GSL minimizer.
 */
double min_f(const gsl_vector *v, void *v_mdl) {
  Model *mdl;
  double ll;

  mdl = v_mdl;

  /* update model parameters */
  model_update_param(mdl, v);

/*   fprintf(stderr, "\nmin_f:\n  "); */
/*   fprintf(stderr, "  gsl_param: "); */
/*   write_gsl_vector(stderr, v); */
/*   fprintf(stderr, "  model_param: "); */
/*   model_write_param_ln(mdl, stderr); */

  /* calculate log-likelihood but not gradient*/
  ll = mdl->llhood(mdl, TRUE, FALSE);

  /*  fprintf(stderr, "  LL=%.10f\n", ll);*/

  /* return negative log-likelihood to minimizer */
  return -ll;
}



/**
 * Wrapper around gradient function called by GSL minimizer.
 */
void min_df(const gsl_vector *v, void *v_mdl, gsl_vector *grad) {
  Model *mdl;

  mdl = v_mdl;

  /* update model parameters */
  model_update_param(mdl, v);

/*   fprintf(stderr, "\nmin_df\n"); */
/*   fprintf(stderr, "  model_param: "); */
/*   model_write_param_ln(mdl, stderr); */
/*   fprintf(stderr, "  gsl_param: "); */
/*   write_gsl_vector(stderr, v); */

  /* calculate gradient */
  mdl->llhood(mdl, FALSE, TRUE);


  /* update gradient that is given back to GSL */
  model_to_gsl_grad(mdl, v, grad);

/*   fprintf(stderr, "  model_grad: "); */
/*   model_write_grad_ln(mdl, stderr); */
/*   fprintf(stderr, "  gsl_grad: "); */
/*   write_gsl_vector(stderr, grad); */

}


/**
 * Wrapper around combined LL / gradient function called by
 * GSL minimizer
 */

void min_fdf(const gsl_vector *v, void *v_mdl, double *f, gsl_vector *grad) {
  Model *mdl;

  mdl = v_mdl;

  fprintf(stderr, "\nmin_fdf\n");

  /* update model parameters */
  model_update_param(mdl, v);

/*   fprintf(stderr, "  model_param: "); */
/*   model_write_param_ln(mdl, stderr); */
/*   fprintf(stderr, "  gsl_param: "); */
/*   write_gsl_vector(stderr, v); */


  /* calculate log-likelihood and gradient */
  *f = -mdl->llhood(mdl, TRUE, TRUE);

  /* update gradient that is given back to GSL */
  model_to_gsl_grad(mdl, v, grad);

/*   fprintf(stderr, "  model_grad: "); */
/*   model_write_grad_ln(mdl, stderr); */
/*   fprintf(stderr, "  gsl_grad: "); */
/*   write_gsl_vector(stderr, grad); */
/*   fprintf(stderr, "  LL=%g", -*f); */

}




/**
 * Outputs model parameter values, iteration number and
 * currelt log-likelihood to provided filehandle
 */
void write_param(FILE *fh, Model *mdl, int iter, double llhood) {
  fprintf(fh, "iter=%d ", iter);
  model_write_param(mdl, fh);
  fprintf(fh, " LL=%.10g\n", llhood);
}





void set_step_size(gsl_vector *x, gsl_vector *ss) {
  int i;

  if(x->size != ss->size) {
    g_error("set_step_size: expected step size vector length (%zd) "
	    "to equal param vector length (%zd)",
	    x->size, ss->size);
  }

  for(i = 0; i < x->size; i++) {
    gsl_vector_set(ss, i, fabs(gsl_vector_get(x,i))*INIT_STEP_SIZE);
  }
}




/*
 * Calls GSL nelder mead simplex minimizer, returns max log likelihood found.
 */
double minimize_simplex(Model *mdl, int max_iter, int min_type) {
  gsl_multimin_fminimizer *minimizer;
  gsl_multimin_function min_func;
  gsl_vector *x, *ss;
  int status;
  double ll, size;
  int iter = 0;

  /* Set initial parameter values */
  x = model_to_gsl_param(mdl);

  /* set initial step sizes */
  ss = gsl_vector_alloc(mdl->n_unlocked);
  set_step_size(x, ss);

   
  /* Initialize method and iterate */
  min_func.n = mdl->n_unlocked;
  min_func.f = min_f;
  min_func.params = mdl;

  minimizer = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, 
					    mdl->n_unlocked);

  gsl_multimin_fminimizer_set(minimizer, &min_func, x, ss);

  fprintf(stderr, "initial parameters: ");
  model_write_param_ln(mdl, stderr);

  ll = -(minimizer->fval);

  do {
    iter++;

    status = gsl_multimin_fminimizer_iterate(minimizer);

    if(status) {
      fprintf(stderr, gsl_strerror(status)); 
      break;
    }
        
    size = gsl_multimin_fminimizer_size(minimizer);
    status = gsl_multimin_test_size(size, MAX_SIMPLEX_SZ);
     

    /* copy current param vals from minimizer back into model and write
     * current parameters to stderr
     */
    ll = -minimizer->fval;
    model_update_param(mdl, minimizer->x);
    write_param(stderr, mdl, iter, ll);
  
  } while (status == GSL_CONTINUE && iter < max_iter);


  write_param(stderr, mdl, iter, -minimizer->fval);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(minimizer);
  
  return ll;
}


/*
 * Calls GSL minimizer function, returns max log likelihood found.
 */
double minimize(Model *mdl, int max_iter, int min_type) {
  gsl_multimin_fdfminimizer *minimizer;
  gsl_multimin_function_fdf min_func;
  gsl_vector *x;
  int status;
  double ll;
  int iter = 0;

  /* Set initial parameter values */
  x = model_to_gsl_param(mdl);
   
  /* Initialize method and iterate */
  min_func.n = mdl->n_unlocked;
  min_func.f = min_f;
  min_func.df = min_df;
  min_func.fdf = min_fdf;
  min_func.params = mdl;

  if(min_type == MIN_TYPE_VECTOR_BFGS) {
    minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, mdl->n_unlocked);
  }
  else if(min_type == MIN_TYPE_VECTOR_BFGS2) {
    minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, mdl->n_unlocked);
  }
  else if(min_type == MIN_TYPE_CONJUGATE_PR) {
    minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, mdl->n_unlocked);
  } 
  else {
    g_error("minimize: unknown minimizer type");
    minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, mdl->n_unlocked);
  }


  gsl_multimin_fdfminimizer_set(minimizer, &min_func, x, 
				INIT_STEP_SIZE, LINE_MIN_TOL);

  fprintf(stderr, "initial parameters: ");
  model_write_param_ln(mdl, stderr);

  ll = -minimizer->f;

  do {
    iter++;
    
    status = gsl_multimin_fdfminimizer_iterate(minimizer);

    if(status) {
      fprintf(stderr, gsl_strerror(status)); 
      fprintf(stderr, "\ngradient: ");
      model_write_grad_ln(mdl, stderr);
      break;
    }
    
    /*
     * size = gsl_multimin_fminimizer_size(minimizer);
     * status = gsl_multimin_test_size(size, MAX_SIMPLEX_SZ);
     */
    status = gsl_multimin_test_gradient (minimizer->gradient, GRAD_TOL);


    /* copy current param vals from minimizer back into model and write
     * current parameters to stderr
     */
    ll = -minimizer->f;
    model_update_param(mdl, minimizer->x);
    write_param(stderr, mdl, iter, ll);
  
  } while (status == GSL_CONTINUE && iter < max_iter);


  write_param(stderr, mdl, iter, -minimizer->f);
  
  gsl_vector_free(x);
  gsl_multimin_fdfminimizer_free(minimizer);
  
  return ll;
}




double try_min_param_subset(Model *mdl, int min_type) {
  GList *rand_unlocked, *cur;
  ModelParam *param;
  int n_new_locked;
  double ll;

  /* Lock a random subset of the parameters in the model */
  n_new_locked = mdl->n_unlocked/2;
  if(n_new_locked >= mdl->n_unlocked) {
    n_new_locked = mdl->n_unlocked - 1;
  }

  if(n_new_locked >= 0) {
     rand_unlocked = model_get_rand_unlocked_param(mdl, n_new_locked);
  } else {
    rand_unlocked = NULL;
  }

  cur = rand_unlocked;
  while(cur != NULL) {
    param = cur->data;
    param->locked = TRUE;
    mdl->n_unlocked--;
    cur = g_list_next(cur);
  }

  fprintf(stderr, "locked %d randomly selected parameters\n", n_new_locked);
  model_write_param_ln(mdl, stderr);

  if(min_type == MIN_TYPE_SIMPLEX) {
    ll = minimize_simplex(mdl, MAX_SUBSET_ITER, min_type);
  } else {
    ll = minimize(mdl, MAX_SUBSET_ITER, min_type);
    
    /*
     * fprintf(stderr, "trying simplex method\n");
     * ll = minimize_simplex(mdl, 100, MIN_TYPE_SIMPLEX);
     */
  }

  /* unlock the parameters that we locked */
  cur = rand_unlocked;
  while(cur != NULL) {
    param = cur->data;
    param->locked = FALSE;
    mdl->n_unlocked++;
    cur = g_list_next(cur);
  }
  g_list_free(rand_unlocked);

  return ll;
}




/*
 *
 * Retrieves output filehandle
 */
FILE *get_output_fh(Config *config) {
  char *dir, *file, *path;
  FILE *fh;

  if(config_has_key(config, "OUTPUT_DIR")) {
    dir = config_get_str(config, "OUTPUT_DIR");
    file = config_get_str(config, "OUTPUT_FILE");
    path = g_strconcat(dir, file, NULL);
    fh = fopen(path, "w");

    if(fh == NULL) {
      g_error("get_output_file: could not open file '%s'", path);
    }

    g_free(path);
  } else {
    fh = stdout;
  }
  return fh;
}


/**
 * Closes output filehandle
 */
void close_output_fh(FILE *fh) {
  if(fh != stdout) {
    fclose(fh);
  }
}





/**
 * Gets the type of minimizer that should be used
 */
int get_min_type(Config *config) {
  char *min_str;

  min_str = config_get_str(config, "MINIMIZER_TYPE");

  if(strcmp(min_str, "CONJUGATE_PR")==0) {
    fprintf(stderr, "Using Polak-Ribiere conjugate gradient minimizer\n");
    return MIN_TYPE_CONJUGATE_PR;
  }
  if(strcmp(min_str, "VECTOR_BFGS")==0) {
    fprintf(stderr, "Using BFGS minimizer\n");
    return MIN_TYPE_VECTOR_BFGS;
  }
  if(strcmp(min_str, "VECTOR_BFGS2")==0) {
    fprintf(stderr, "Using BFGS2 minimizer\n");
    return MIN_TYPE_VECTOR_BFGS2;
  }
  if(strcmp(min_str, "SIMPLEX")==0) {
    fprintf(stderr, "Using Nelder-Mead SIMPLEX minimizer\n");
    return MIN_TYPE_SIMPLEX;
  }

  g_error("unknown minimizer type '%s'", min_str);
  return MIN_TYPE_VECTOR_BFGS2;

}



double run_max_ll(Model *mdl, int min_type, int max_tries, double diff_sz) {
  int improved, n_tries;
  double new_ll, ll;


  fprintf(stderr, "maximizing LL\n");
  improved = TRUE;
  ll = 0;
  while(improved) {
    if(min_type == MIN_TYPE_SIMPLEX) {
      ll = minimize_simplex(mdl, MAX_ITER, min_type);
    } else {
      ll = minimize(mdl, MAX_ITER, min_type);
    }

    improved = FALSE;
    n_tries = 0;
    while(!improved && n_tries < max_tries) {
      /* try locking subset of param and performing extra iterations */
      new_ll = try_min_param_subset(mdl, min_type);

      n_tries++;
      if(new_ll - ll > diff_sz) {
	improved = TRUE;
	fprintf(stderr, "found improvement\n");
      } else {
	fprintf(stderr, "no improvement\n");
      }
    }
    
    if(!improved && min_type != MIN_TYPE_SIMPLEX) {
      /* try one last simplex refinement */
      fprintf(stderr, "trying simplex refinment\n");
      new_ll = minimize_simplex(mdl, MAX_ITER, min_type);

      if(new_ll - ll > diff_sz) {
	improved = TRUE;
	fprintf(stderr, "found improvement\n");
      } else {
	fprintf(stderr, "no improvement\n");
      }
    }
  }

  return ll;
}


/**
 * Main function
 */
int main(int argc, char **argv) {
  Config *config;
  Model *mdl;
  int min_type, max_tries, combine_data_bin, unbinned_refinement;
  double ll;
  FILE *out_fh;
    
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");
  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  min_type = get_min_type(config);
  mdl = bkgd_evo_mdl_new(config);
  bkgd_evo_mdl_read_data(mdl, config);
  out_fh = get_output_fh(config);

  max_tries = config_get_long(config, "MAX_RESTART_TRIES");
 
  combine_data_bin = config_get_boolean(config, "COMBINE_DATA_BINS");
  unbinned_refinement = config_get_boolean(config, "DO_UNBINNED_REFINEMENT");


  if(unbinned_refinement && config_get_boolean(config, "USE_MDIV_CATS")) {
    g_error("cannot perform unbinned refinement if using MDIV categories");
  }

  

  if(combine_data_bin) {
    fprintf(stderr, "first pass: maximizing LL of binned data\n");
    /* do a first pass with binned data */
    bkgd_evo_mdl_bin_data(mdl, config);

    ll = run_max_ll(mdl, min_type, max_tries, LL_LARGE_DIFF);

    /* re-read unbinned data */
    if(unbinned_refinement) {
      bkgd_evo_mdl_free_data(mdl);
      bkgd_evo_mdl_read_data(mdl, config);

      fprintf(stderr, "maximizing LL of unbinned data\n");
      ll = run_max_ll(mdl, min_type, max_tries, LL_SMALL_DIFF);
      fprintf(stderr, "final LL=%.10f\n", ll);
    }
  } else { 
    fprintf(stderr, "maximizing LL of unbinned data\n");
    ll = run_max_ll(mdl, min_type, max_tries, LL_SMALL_DIFF);
    fprintf(stderr, "final LL=%.10f\n", ll);
  }

  write_param(out_fh, mdl, 0, ll);

  close_output_fh(out_fh);
  bkgd_evo_mdl_free(mdl);

  return 0;
}
