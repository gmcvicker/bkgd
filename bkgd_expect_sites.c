
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "branch.h"
#include "cltype.h"
#include "bkgd_evo_mdl.h"
#include "bkgd_data.h"

#define B_INCREMENT 0.001
#define BRANCH_LEN_SMALL 1e-12


void write_obs_counts(FILE *fh, Model *mdl) {
  BkgdEvoMdlParam param;
  BkgdEvoMdlConfig *conf;
  BkgdEvoMdlData *data;
  ColType *cltype;
  int j;
  long i;
   
  conf = mdl->config;
  data = mdl->data;

  param_set_zero(&param);
  set_bkgd_evo_param_from_model(mdl, conf, &param);

  /* write header */
  fprintf(fh, "B");
  for(j = 0; j < conf->n_cltype; j++) {
    fprintf(fh, "\t%s", conf->cltypes[j]->name);
  }
  fprintf(fh, "\n");

  for(i = 0; i < data->n_bin; i++) {
    param.B = exp(data->bin[i].lB_ex  * param.u_ex_scale + 
		  data->bin[i].lB_nex * param.u_nex_scale);

    if(param.B < 0.0 || param.B > 1.0 || isnan(param.B) || isinf(param.B)) {
      g_error("bkgd_evo_mdl_calc_ll: invalid B (%g)\n"
	      "B_ex=%g, B_nex=%g, u_ex_scale=%g, u_nex_scale=%g\n",
	      param.B, data->bin[i].B_ex, data->bin[i].B_nex,
	      param.u_ex_scale, param.u_nex_scale);
    }

    /* write number of observed column type for each bin */
    fprintf(fh, "%g", param.B);
    for(j = 0; j < conf->n_cltype; j++) {
      cltype = conf->cltypes[j];
      cltype->set_n(cltype, &data->bin[i]);
      fprintf(fh, "\t%ld", cltype->n);
    }
    fprintf(fh, "\n");
  }
}




void write_expect_sites(FILE *fh, Model *mdl) {
  BkgdEvoMdlParam param;
  BkgdEvoMdlConfig *conf;
  Branch *br;
  ColType *cltype;
  double b;
  int j;
   
  conf = mdl->config;

  param_set_zero(&param);
  set_bkgd_evo_param_from_model(mdl, conf, &param);

  b = 0.0;

  if(conf->subst_mdl->id == SUBST_MDL_JUKES_CANTOR) {
    param.mu = model_get_param_val(mdl, "mu");
    param.lambda = model_get_param_val(mdl, "lambda");
  } 
  else if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
    param.mu_i = model_get_param_val(mdl, "mu_i");
    param.mu_v = model_get_param_val(mdl, "mu_v");
    param.lambda_i = model_get_param_val(mdl, "lambda_i");
    param.lambda_v = model_get_param_val(mdl, "lambda_v");
  } else {
    g_error("unknown substitution model type");
  }


  /* write header */
  fprintf(fh, "B");
  for(j = 0; j < conf->n_cltype; j++) {
    fprintf(fh, "\t%s", conf->cltypes[j]->name);
  }
  fprintf(fh, "\n");

  while(b <= 1.0) {
    b += B_INCREMENT;

    param.B = b;

    /* calculate probability HC coalescent predates HC/G speciation */
    if(param.T_hcg == 0.0) {
      param.k_hcg = 1.0;
    } else {
      param.k_hcg = exp(-0.5 * param.T_hcg / (param.B * param.N_hc));

      if(param.k_hcg > 1.0 || param.k_hcg < 0.0) {
	g_error("bkgd_evo_mdl_calc_ll: invalid probability of HC "
		"coalescent predating gorilla speciation: k_hcg=%g\n",
		param.k_hcg);
      }
    }

    /* calculate branch lengths and substitution probs for current bin */
    for(j = 0; j < conf->n_branch; j++) {
      br = &conf->branches[j];

      br->set_len(br, &param);

      if(br->len < 0.0 || isnan(br->len) || isinf(br->len)) {
	g_warning("write_expect_sites: invalid length (%g) for branch %s",
		br->len, br->name);
	br->len = BRANCH_LEN_SMALL;
      }

      branch_set_prob(br, &param, conf);
    }

    /* calc prob of observing column type for bin
     * doing non-conserved columns first
     */
    for(j = 0; j < conf->n_cltype; j++) {
      cltype = conf->cltypes[j];
      if(cltype->subst_type != SUBST_TYPE_CONSERVED) {
	cltype_set_prob(cltype, &param, conf);      
      }
    }
    /* now calc prob of conserved col type (depends on other col probs) */
    if(conf->cons_cltype) {
      cltype = conf->cons_cltype;
      cltype_cons_set_prob(cltype, &param, conf);
    }

    /* now write out probabilities for each col type */
    fprintf(fh, "%g", param.B);
    for(j = 0; j < conf->n_cltype; j++) {
      cltype = conf->cltypes[j];
      fprintf(fh, "\t%g", cltype->prob);
    }
    fprintf(fh, "\n");

  }
}




FILE *get_expect_file(Config *config) {
  char *dir, *file, *path;
  FILE *fh;

  dir = config_get_str(config, "OUTPUT_DIR");
  file = config_get_str(config, "EXPECT_OUTPUT_FILE");
  
  path = g_strconcat(dir, file, NULL);

  fh = fopen(path, "w");

  if(fh == NULL) {
    g_error("could not open file '%s'", path);
  }

  g_free(path);

  return fh;
}



FILE *get_obs_file(Config *config) {
  char *dir, *file, *path;
  FILE *fh;

  dir = config_get_str(config, "OUTPUT_DIR");
  file = config_get_str(config, "OBS_OUTPUT_FILE");
  
  path = g_strconcat(dir, file, NULL);

  fh = fopen(path, "w");

  if(fh == NULL) {
    g_error("could not open file '%s'", path);
  }

  g_free(path);

  return fh;
}


/**
 * Main function
 */
int main(int argc, char **argv) {
  Config *config;
  Model *mdl;
  double mean_cat;
  FILE *obs_fh, *exp_fh;
    
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");
  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  obs_fh = get_obs_file(config);
  exp_fh = get_expect_file(config);

  mdl = bkgd_evo_mdl_new(config);
  bkgd_evo_mdl_read_data(mdl, config);
  mean_cat = bkgd_data_mean_cat_val(mdl->data);

  fprintf(stderr, "mean category value: %g\n", mean_cat);
  
  bkgd_evo_mdl_bin_data(mdl, config);

  write_expect_sites(exp_fh, mdl);
  write_obs_counts(obs_fh, mdl);

  bkgd_evo_mdl_free(mdl);


  fclose(obs_fh);
  fclose(exp_fh);

  return 0;
}
