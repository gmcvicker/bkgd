#include <glib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <util/util.h>

#include "model.h"
#include "bkgd_evo_mdl.h"




static FILE *get_out_fh(Config *config, char *parm_name) {
  char *path;
  char *dir;
  FILE *fh;

  dir = config_get_str(config, "OUTPUT_DIR");

  path = g_strconcat(dir, parm_name, ".txt", NULL);

  fh = fopen(path, "w");
  if(fh == NULL) {
    g_error("get_out_fh: could not open file '%s'", path);
  }

  g_free(path);

  return fh;
}



void scan_parm_area(FILE *fh, Model *mdl, char *param_name, double start, 
		    double end, int n_step) {  
  int i;
  double orig_val, cur_val, step_sz;
  double ll, deriv;


  fprintf(stderr, "%s [%g-%g]\n", param_name, start, end);

  orig_val = model_get_param_val(mdl, param_name);

  step_sz = (end - start) / (double)n_step;

  /* report LL at original param value as well */
  ll = bkgd_evo_mdl_calc_ll(mdl, TRUE, TRUE);
  fprintf(stderr, "first call: LL=%.10f\n", ll);

  deriv = model_get_param_deriv(mdl, param_name);
  fprintf(fh, "%g %.10g %.10g\n", orig_val, ll, deriv);

  cur_val = start;
  for(i = 0; i < n_step; i++) {
    model_set_param_val(mdl, param_name, cur_val);
    ll = bkgd_evo_mdl_calc_ll(mdl, TRUE, TRUE);
    deriv = model_get_param_deriv(mdl, param_name);
    fprintf(fh, "%g %.10g %.10g\n", cur_val, ll, deriv);
    
    cur_val += step_sz;
  }
  
  /* restore paramter to original value */
  model_set_param_val(mdl, param_name, orig_val);
}




#define PARAM_OFFSET 0.02

void explore_parm(Config *config, char *parm_name, Model *mdl) {
  char *key;
  double val, start, end;
  int n_step;
  FILE *fh;
  

  fprintf(stderr, "parm_name=%s\n", parm_name);

  n_step = config_get_long(config, "N_SCAN_STEP");

  val = model_get_param_val(mdl, parm_name);

  key = g_strconcat(parm_name, "_START", NULL);
  util_str_uc(key);  
  if(config_has_key(config, key)) {
    start = config_get_double(config, key);
  } else {
    start = val - (val * PARAM_OFFSET);
  }
  g_free(key);

  key = g_strconcat(parm_name, "_END", NULL);
  util_str_uc(key);  
  if(config_has_key(config, key)) {
    end = config_get_double(config, key);
  } else {
    end = val + (val * PARAM_OFFSET);
  }
  g_free(key);
  
  n_step = config_get_long(config, "N_SCAN_STEP");

  fh = get_out_fh(config, parm_name);
  scan_parm_area(fh, mdl, parm_name, start, end, n_step);
  fclose(fh);
}




void explore_ll_surface(Config *config, Model *mdl) {
  char **parm_names;
  int n_parm, i;

  parm_names = config_get_str_array(config, "PARAM_NAMES", &n_parm);

  for(i = 0; i < n_parm; i++) {
    explore_parm(config, parm_names[i], mdl);
  }
}



/**
 * Main function
 */
int main(int argc, char **argv) {
  Config *config;
  Model *mdl;
    
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");
  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  fprintf(stderr, "creating model\n");
  mdl = bkgd_evo_mdl_new(config);
  bkgd_evo_mdl_read_data(mdl, config);

  fprintf(stderr, "scanning LL surface\n");
  explore_ll_surface(config, mdl);

  bkgd_evo_mdl_free(mdl);

  return 0;
}
