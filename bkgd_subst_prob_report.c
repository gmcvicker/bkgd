
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "branch.h"
#include "cltype.h"
#include "bkgd_evo_mdl.h"
#include "bkgd_data.h"
#include "model.h"



static void update_branch_probs(BkgdEvoMdlConfig *conf, 
                                BkgdEvoMdlParam *param) {
  int i;
  Branch *br;

  for(i = 0; i < conf->n_branch; i++) {
    br = &conf->branches[i]; 
    br->set_len(br, param);

    if(br->len < 0.0 || isnan(br->len) || isinf(br->len)) {
      g_warning("%s:%d: invalid length (%g) for branch %s", __FILE__, __LINE__,
                br->len, br->name);
      br->len = BRANCH_LEN_SMALL;
    }

    branch_set_prob(br, param, conf);
  }
}



static void update_k_hcg(BkgdEvoMdlParam *param) {

  if(param->T_hcg == 0.0) {
    param->k_hcg = 1.0;
  } else {
    param->k_hcg = exp(-0.5 * param->T_hcg / (param->B * param->N_hc));
    
    if(param->k_hcg > 1.0 || param->k_hcg < 0.0) {
      g_error("%s:%d: invalid probability of HC "
	      "coalescent predating gorilla speciation: k_hcg=%g\n",
	      __FILE__, __LINE__, param->k_hcg);
    }
  }
}



void write_subst_probs(FILE *fh, Model *mdl) {
  BkgdEvoMdlParam param;
  BkgdEvoMdlConfig *conf;
  BkgdEvoMdlData *data;
  ColType *cltype;
  int j;
  long i, *cltype_ttls;
  double *single_subst_ttls, *dbl_subst_ttls, single_prp, dbl_prp;
   
  conf = mdl->config;
  data = mdl->data;

  param_set_zero(&param);
  set_bkgd_evo_param_from_model(mdl, conf, &param);
  model_write_param_ln(mdl, stderr);

  cltype_ttls = g_new(long, conf->n_cltype);
  single_subst_ttls = g_new(double, conf->n_cltype);
  dbl_subst_ttls = g_new(double, conf->n_cltype);
  for(i = 0; i < conf->n_cltype; i++) {
    cltype_ttls[i] = 0;
    dbl_subst_ttls[i] = 0.0;
    single_subst_ttls[i] = 0.0;
  }

  for(i = 0; i < data->n_bin; i++) {
    param.B = exp(data->bin[i].lB_ex  * param.u_ex_scale + 
		  data->bin[i].lB_nex * param.u_nex_scale);

    if(param.B < 0.0 || param.B > 1.0 || isnan(param.B) || isinf(param.B)) {
      g_error("%s:%d: invalid B (%g)\n"
	      "B_ex=%g, B_nex=%g, u_ex_scale=%g, u_nex_scale=%g\n",
              __FILE__, __LINE__,
	      param.B, data->bin[i].B_ex, data->bin[i].B_nex,
	      param.u_ex_scale, param.u_nex_scale);
    }
    
    /* update mutation rate */
    bkgd_evo_mdl_calc_mu(&param, mdl, &data->bin[i]);

    /* update probability HC coalescent predates HC/G speciation */
    update_k_hcg(&param);

    /* update branch subsitution probabilities */
    update_branch_probs(conf, &param);

    /* set column type probabilities */
    for(j = 0; j < conf->n_cltype; j++) {
      cltype = conf->cltypes[j];

      if(cltype->subst_type == SUBST_TYPE_CONSERVED) {
	/* could handle these separately */
	continue;
      }

      cltype->set_n(cltype, &data->bin[i]);
      cltype_set_prob(cltype, &param, conf);

      single_prp = cltype->prob_single / 
	(cltype->prob_ttl_double + cltype->prob_single);

      dbl_prp = 1.0 - single_prp;

      dbl_subst_ttls[j] += (double)cltype->n * dbl_prp;
      single_subst_ttls[j] += (double)cltype->n * single_prp;
      cltype_ttls[j] += cltype->n;
    }
  }



  fprintf(fh, "COL.TYPE\tN\tN.SINGLE\tN.DBL\tPRP.SINGLE\tPRP.DBL\n");
	  
  for(i = 0; i < conf->n_cltype; i++) {
    cltype = conf->cltypes[i];

    if(cltype->subst_type != SUBST_TYPE_CONSERVED) {

      fprintf(fh, "%s\t%ld\t%g\t%g\t%g\t%g\n",
	      cltype->name, cltype_ttls[i], single_subst_ttls[i], 
	      dbl_subst_ttls[i], single_subst_ttls[i] / (double)cltype_ttls[i],
	      dbl_subst_ttls[i] / (double)cltype_ttls[i]);
    }
  }

  g_free(cltype_ttls);
  g_free(single_subst_ttls);
  g_free(dbl_subst_ttls);

}




/**
 * Main function
 */
int main(int argc, char **argv) {
  Config *config;
  Model *mdl;
  double mean_cat;
    
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");
  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  mdl = bkgd_evo_mdl_new(config);
  bkgd_evo_mdl_read_data(mdl, config);
  mean_cat = bkgd_data_mean_cat_val(mdl->data);

  fprintf(stderr, "mean category value: %g\n", mean_cat);

  if(config_get_boolean(config, "BIN_DATA")) {
    bkgd_evo_mdl_bin_data(mdl, config);
  }

  write_subst_probs(stdout, mdl);

  bkgd_evo_mdl_free(mdl);


  return 0;
}
