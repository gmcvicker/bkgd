#include <glib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <util/util.h>

#include "model.h"
#include "bkgd_evo_mdl.h"
#include "branch.h"



static FILE *get_out_fh(Config *config, Branch *br, char *parm_name,
			char *prefix) {
  char *path;
  char *dir;
  FILE *fh;

  dir = config_get_str(config, "OUTPUT_DIR");

  path = g_strconcat(dir, prefix, br->name, "_", parm_name, ".txt", NULL);

  fh = fopen(path, "w");
  if(fh == NULL) {
    g_error("get_out_fh: could not open file '%s'", path);
  }

  g_free(path);

  return fh;
}



void scan_parm_area(Config *config, Model *mdl, char *param_name, 
		    Branch *br, double start,  double end, int n_step) {  
  int i;
  double orig_val, cur_val, step_sz;
  double deriv;
  BkgdEvoMdlParam param;
  BkgdEvoMdlConfig *bem_conf;
  BkgdBin bin;
  FILE *br_fh, *br_prob_fh, *br_prob_i_fh, *br_prob_v_fh;

  br_fh = get_out_fh(config, br, param_name, "br_len_");
  br_prob_fh = get_out_fh(config, br, param_name, "br_pr_");
  br_prob_i_fh = get_out_fh(config, br, param_name, "br_pr_i_");
  br_prob_v_fh = get_out_fh(config, br, param_name, "br_pr_v_");


  bin.B_ex = 0.5;
  bin.B_nex = 0.5;
  bin.lB_ex = log(bin.B_ex);
  bin.lB_nex = log(bin.B_nex);
  bin.cat = 0.03;

  bem_conf = mdl->config;
  
  orig_val = model_get_param_val(mdl, param_name);

  step_sz = (end - start) / (double)n_step;
  

  param_set_zero(&param);

  /* report branch len at original param value as well */
  set_bkgd_evo_param_from_model(mdl, bem_conf, &param);


  param.B = exp(bin.lB_ex  * param.u_ex_scale + 
		bin.lB_nex * param.u_nex_scale);

  bkgd_evo_mdl_calc_mu(&param, mdl, &bin);

  fprintf(stderr, "mu: %g\n", param.mu);
  fprintf(stderr, "mu_a: %g\n", param.mu_a);
  fprintf(stderr, "mu_b: %g\n", param.mu_b);
  fprintf(stderr, "mu_i: %g\n", param.mu_i);
  fprintf(stderr, "mu_v: %g\n", param.mu_v);
  
  param.k_hcg = exp(-0.5 * param.T_hcg / (param.B * param.N_hc));

  br->set_len(br, &param);
  branch_set_prob(br, &param, bem_conf);
  br->set_dlen(br, &bin, &param);
  branch_set_dprob(br, &bin, &param, bem_conf);

  /* report branch len and deriv */
  bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dlen);
  deriv = model_get_param_deriv(mdl, param_name);
  fprintf(br_fh, "%g %.10g %.10g\n", orig_val, br->len, deriv);

  /* report branch probability and derivatives */
  bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dprob);
  deriv = model_get_param_deriv(mdl, param_name);
  fprintf(br_prob_fh, "%g %.10g %.10g\n", orig_val, br->prob, deriv);

  /* report branch transition probability and derivatives */
  bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dprob_i);
  deriv = model_get_param_deriv(mdl, param_name);
  fprintf(br_prob_i_fh, "%g %.10g %.10g\n", orig_val, br->prob_i, deriv);

  /* report branch transition probability and derivatives */
  bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dprob_v);
  deriv = model_get_param_deriv(mdl, param_name);
  fprintf(br_prob_v_fh, "%g %.10g %.10g\n", orig_val, br->prob_v, deriv);


  cur_val = start;
  for(i = 0; i < n_step; i++) {
    /* update param val */
    model_set_param_val(mdl, param_name, cur_val);

    /* calculate branch len and partial derivs */
    set_bkgd_evo_param_from_model(mdl, bem_conf, &param);


    param.B = exp(bin.lB_ex  * param.u_ex_scale + 
		  bin.lB_nex * param.u_nex_scale);

    bkgd_evo_mdl_calc_mu(&param, mdl, &bin);
    
    param.k_hcg = exp(-0.5 * param.T_hcg / (param.B * param.N_hc));

    br->set_len(br, &param);
    branch_set_prob(br, &param, bem_conf);
    br->set_dlen(br, &bin, &param);
    branch_set_dprob(br, &bin, &param, bem_conf);

    /* report branch len and deriv */
    bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dlen);
    deriv = model_get_param_deriv(mdl, param_name);
    fprintf(br_fh, "%g %.10g %.10g\n", cur_val, br->len, deriv);

    /* report branch probability and derivatives */
    bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dprob);
    deriv = model_get_param_deriv(mdl, param_name);
    fprintf(br_prob_fh, "%g %.10g %.10g\n", cur_val, br->prob, deriv);

    /* report branch transition probability and derivatives */
    bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dprob_i);
    deriv = model_get_param_deriv(mdl, param_name);
    fprintf(br_prob_i_fh, "%g %.10g %.10g\n", cur_val, br->prob_i, deriv);

    /* report branch transition probability and derivatives */
    bkgd_evo_mdl_set_mdl_deriv(mdl, &br->dprob_v);
    deriv = model_get_param_deriv(mdl, param_name);
    fprintf(br_prob_v_fh, "%g %.10g %.10g\n", cur_val, br->prob_v, deriv);
    
    cur_val += step_sz;
  }
  
  fclose(br_fh);
  fclose(br_prob_fh);
  fclose(br_prob_i_fh);
  fclose(br_prob_v_fh);
  

  /* restore parameter to original value */
  model_set_param_val(mdl, param_name, orig_val);
}




#define PARAM_OFFSET 0.02

void explore_parm(Config *config, char *parm_name, Branch *br, Model *mdl) {
  char *key;
  double val, start, end;
  int n_step;
  
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
  
  scan_parm_area(config, mdl, parm_name, br, start, end, n_step);

}




void explore_branches(Config *config, Model *mdl) {
  Branch *br;
  int i, j, n_parm;
  BkgdEvoMdlConfig *bem_conf;
  char **parm_names;

  parm_names = config_get_str_array(config, "PARAM_NAMES", &n_parm);

  bem_conf = mdl->config;

  for(i = 0; i < bem_conf->n_branch; i++) {
    br = &bem_conf->branches[i];

    for(j = 0; j < n_parm; j++) {
      fprintf(stderr, "%s\n", parm_names[j]);
      explore_parm(config, parm_names[j], br, mdl);
    }
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

  mdl = bkgd_evo_mdl_new(config);
  bkgd_evo_mdl_read_data(mdl, config);

  explore_branches(config, mdl);

  bkgd_evo_mdl_free(mdl);

  return 0;
}
