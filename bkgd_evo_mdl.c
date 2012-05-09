#include <glib.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <util/config.h>
#include <util/util.h>
#include <numer/elog.h>

#include <gsl/gsl_multimin.h>

#include "model.h"
#include "bkgd_evo_mdl.h"
#include "branch.h"
#include "cltype.h"
#include "bkgd_data.h"



/**
 * Adds a new branch definition to the provided model
 * configuration structure.
 */
static Branch *add_branch(BkgdEvoMdlConfig *conf, char *name, 
			  void (*set_len)(Branch *, const BkgdEvoMdlParam *),
			  void (*set_dlen)(Branch *, const BkgdBin *, 
					   const BkgdEvoMdlParam *)) {

  Branch *br;

  if(conf->n_branch >= MAX_BRANCH) {
    g_error("add_branch: number of branches exceeds max (%d)", MAX_BRANCH);
  }

  conf->n_branch++;

  br = &conf->branches[conf->n_branch-1];

  br->id = conf->n_branch-1;
  br->name = g_strdup(name);
  br->len = 0.0;
  br->jc = 0.0;
  br->prob = 0.0;

  br->set_len = set_len;
  br->set_dlen = set_dlen;

  param_set_zero(&br->dlen);
  param_set_zero(&br->dprob);

  return br;
}




static void add_cltype(BkgdEvoMdlConfig *conf, ColType *cltype) {
  if(conf->n_cltype >= MAX_COL_TYPE) {
    g_error("add_col_type: number of col types exceeds max (%d)", 
	    MAX_COL_TYPE);
  }
  
  conf->n_cltype++;
  conf->cltypes[conf->n_cltype-1] = cltype;
  
  if(cltype->subst_type == SUBST_TYPE_CONSERVED) {
    if(conf->cons_cltype != NULL) {
      g_error("add_col_type: there is already a column of type CONSERVED but"
	      " should only be one");
    }
    conf->cons_cltype = conf->cltypes[conf->n_cltype-1];
  }

  return;
}



/* returns ptr to first branch with specified name in config 
 *
 */
static Branch *get_branch(BkgdEvoMdlConfig *conf, const char *name) {
  int i;

  for(i = 0; i < conf->n_branch; i++) {
    if(strcmp(conf->branches[i].name, name) == 0) {
      return &conf->branches[i];
    }
  }

  g_error("get_branch: unknown branch %s", name);
  return NULL;
}


static ColType *get_cltype(BkgdEvoMdlConfig *conf, const char *name) {
  int i;

  for(i = 0; i < conf->n_cltype; i++) {
    if(strcmp(conf->cltypes[i]->name, name) == 0) {
      return conf->cltypes[i];
    }
  }

  g_error("get_cltype: unknown cltype %s", name);
  return NULL;
}





static void write_col_probs(FILE *fh, BkgdEvoMdlConfig *conf) {
  int i;
  double ttl_p;

  ttl_p = 0.0;
  for(i = 0; i < conf->n_cltype; i++) {
    fprintf(fh, "pi_%s=%g\n", conf->cltypes[i]->name, conf->cltypes[i]->prob);
    ttl_p += conf->cltypes[i]->prob;
  }

  fprintf(fh, "ttl=%g\n", ttl_p);
}



/**
 * Sets the parameters of Background evolutionary model
 * by retrieving named attributes from the provided model
 * structure.
 */
void set_bkgd_evo_param_from_model(const Model *mdl, 
				   const BkgdEvoMdlConfig *conf,
				   BkgdEvoMdlParam *param) {
  
  param->N_hc = model_get_param_val(mdl, "N_hc");
  param->T_hc = model_get_param_val(mdl, "T_hc");
  param->u_ex_scale = model_get_param_val(mdl, "u_ex_scale");
  param->u_nex_scale = model_get_param_val(mdl, "u_nex_scale");

  if(conf->mdl_type == MDL_TYPE_HCM) {
    param->N_hcm = model_get_param_val(mdl, "N_hcm");
    param->T_hcm = model_get_param_val(mdl, "T_hcm");
  } 
  else if(conf->mdl_type == MDL_TYPE_HCGOM) {
    param->N_hcg = model_get_param_val(mdl, "N_hcg");
    param->N_hcgo = model_get_param_val(mdl, "N_hcgo");
    param->N_hcgom = model_get_param_val(mdl, "N_hcgom");

    param->T_hcg = model_get_param_val(mdl, "T_hcg");
    param->T_hcgo = model_get_param_val(mdl, "T_hcgo");
    param->T_hcgom = model_get_param_val(mdl, "T_hcgom");
  }
  else if(conf->mdl_type == MDL_TYPE_HCOM) {
    param->N_hco = model_get_param_val(mdl, "N_hco");
    param->N_hcom = model_get_param_val(mdl, "N_hcom");

    param->T_hco = model_get_param_val(mdl, "T_hco");
    param->T_hcom = model_get_param_val(mdl, "T_hcom");
  }
  else if(conf->mdl_type == MDL_TYPE_HCGO) {
    param->N_hcg = model_get_param_val(mdl, "N_hcg");
    param->N_hcgo = model_get_param_val(mdl, "N_hcgo");

    param->T_hcg = model_get_param_val(mdl, "T_hcg");
    param->T_hcgo = model_get_param_val(mdl, "T_hcgo");
  }
  else if(conf->mdl_type != MDL_TYPE_HC) {
    g_error("bkgd_evo_param_from_model: unknown model type %d",
	    conf->mdl_type);
  }
    
  /* mu is set later because can be vary for each data bin */
  param->mu = 0.0;
  param->mu_a = 0.0;
  param->mu_b = 0.0;
  param->mu_c = 0.0;

  return;
}





/**
 * sets mutation rate parameters for provided data bin
 */
void bkgd_evo_mdl_calc_mu(BkgdEvoMdlParam *param, Model *mdl, BkgdBin *bin) {
  BkgdEvoMdlConfig *conf;
  double cat;

  conf = mdl->config;

  if(conf->mu_type == MU_TYPE_CAT_QUAD) {
    if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      g_error("calc_mu: CAT_QUAD mutation model not implemented for KIMURA "
	      "subst model");
    }

    /* mutation rate is a quadratic function of category value */
    cat = bin->cat;
    param->mu_a = model_get_param_val(mdl, "mu_a"); 
    param->mu_b = model_get_param_val(mdl, "mu_b"); 
    param->mu_c = model_get_param_val(mdl, "mu_c"); 
    param->mu = (param->mu_a * cat * cat) + (param->mu_b * cat) + param->mu_c;
  }
  else if(conf->mu_type == MU_TYPE_CAT_LIN) {
    cat = bin->cat;
    param->mu_a = model_get_param_val(mdl, "mu_a");
    param->mu_b = model_get_param_val(mdl, "mu_b");
    param->mu_c = 0.0;

    if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      param->mu_i = (param->mu_a * cat);
      param->mu_v = (param->mu_b * cat);
      param->mu = param->mu_i + param->mu_v;

      if(conf->use_double_substs) {
        param->lambda_i = model_get_param_val(mdl, "lambda_i");
        param->lambda_v = model_get_param_val(mdl, "lambda_v");
      } else {
        param->lambda_i = 0.0;
        param->lambda_v = 0.0;
      }
      param->lambda = 0.0;
    } else {
      param->mu = (param->mu_a * cat) + param->mu_b;
      if(conf->use_double_substs) {
        param->lambda = model_get_param_val(mdl, "lambda");
      } else {
        param->lambda = 0.0;
      }
      param->lambda_i = 0.0;
      param->lambda_v = 0.0;
    }
  }
  else {
    if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      param->mu_i = model_get_param_val(mdl, "mu_i");
      param->mu_v = model_get_param_val(mdl, "mu_v");
      param->mu = param->mu_i + param->mu_v;

      if(conf->use_double_substs) {
        param->lambda_i = model_get_param_val(mdl, "lambda_i");
        param->lambda_v = model_get_param_val(mdl, "lambda_v");
      } else {
        param->lambda_i = 0.0;
        param->lambda_v = 0.0;
      }
      param->lambda = 0.0;
    } else {
      param->mu = model_get_param_val(mdl, "mu");
      if(conf->use_double_substs) {
        param->lambda = model_get_param_val(mdl, "lambda");
      } else {
        param->lambda = 0.0;
      }
      param->lambda_i = 0.0;
      param->lambda_v = 0.0;
    }
  }

  if(param->mu < 0.0) {
    param->mu = 0.0;
  }
}




/**
 * Calculates contributions to partial derivatives of LL function wrt
 * each of the model parameters for current data bin. Contributions are
 * added to totals in provided deriv argument.
 */
static void calc_deriv(BkgdEvoMdlParam *ll_deriv,
		       BkgdEvoMdlConfig *conf, 
		       const BkgdEvoMdlParam *param, 
		       const BkgdBin *bin) {
  Branch *br;
  ColType *cltype;
  int i;
  double coef;

  /* calculate branch length and branch substitution probability 
   * partial derivatives
   */
  for(i = 0; i < conf->n_branch; i++) {
    br = &conf->branches[i];
    br->set_dlen(br, bin, param);
    branch_set_dprob(br, bin, param, conf);
  }

  /* calculate column probability partial derivatives */
  for(i = 0; i < conf->n_cltype; i++) {
    cltype = conf->cltypes[i];

    if(cltype->subst_type == SUBST_TYPE_CONSERVED) {
      /* only do columns with substitutions in first step */
      continue;
    }

    cltype_set_dprob(cltype, param, conf);
    
    if(cltype->n > 0) {
      /* add contribution to LL partial derivs */
      coef = cltype->n / cltype->prob;

      if(isnan(coef)) {
	fprintf(stderr, "nan coef, cltype->name=%s, cltype->n=%ld, "
		"cltype->prob=%g\n",cltype->name, cltype->n, cltype->prob);
      }

      param_add(&cltype->dprob, coef, ll_deriv);
    }
  }

  /* now conserved column */
  if(conf->cons_cltype != NULL) {
    cltype = conf->cons_cltype;

    cltype_cons_set_dprob(cltype, param, conf);

    if(cltype->n) {
      /* add contribution to LL partial derivs */
      coef = cltype->n / cltype->prob;
      param_add(&cltype->dprob, coef, ll_deriv);
    }
  }

  return;
}



/* sets derivatives in model */
void bkgd_evo_mdl_set_mdl_deriv(Model *mdl, BkgdEvoMdlParam *deriv) {
  BkgdEvoMdlConfig *conf;

  conf = mdl->config;

  model_set_param_deriv(mdl, "u_ex_scale", deriv->u_ex_scale);
  model_set_param_deriv(mdl, "u_nex_scale", deriv->u_nex_scale);
  model_set_param_deriv(mdl, "N_hc", deriv->N_hc);
  model_set_param_deriv(mdl, "T_hc", deriv->T_hc);

  if(conf->mu_type == MU_TYPE_SINGLE) {
    if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      model_set_param_deriv(mdl, "mu_i", deriv->mu_i);
      model_set_param_deriv(mdl, "mu_v", deriv->mu_v);
      if(conf->use_double_substs) {
        model_set_param_deriv(mdl, "lambda_i", deriv->lambda_i);
        model_set_param_deriv(mdl, "lambda_v", deriv->lambda_v);
      }
    } else {
      model_set_param_deriv(mdl, "mu", deriv->mu);
      if(conf->use_double_substs) {
        model_set_param_deriv(mdl, "lambda", deriv->lambda);
      }
    }
  } 
  else if(conf->mu_type == MU_TYPE_CAT_QUAD) {
    if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      g_error("CAT_QUAD mutation model not yet implemented "
	      "for kimura subst model");
    }

    model_set_param_deriv(mdl, "mu_a", deriv->mu_a);
    model_set_param_deriv(mdl, "mu_b", deriv->mu_b);
    model_set_param_deriv(mdl, "mu_c", deriv->mu_c);
    if(conf->use_double_substs) {
      model_set_param_deriv(mdl, "lambda", deriv->lambda);
    }
  } 
  else if(conf->mu_type == MU_TYPE_CAT_LIN) {
    model_set_param_deriv(mdl, "mu_a", deriv->mu_a);
    model_set_param_deriv(mdl, "mu_b", deriv->mu_b);

    if(conf->use_double_substs) {
      if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
        model_set_param_deriv(mdl, "lambda_i", deriv->lambda_i);
        model_set_param_deriv(mdl, "lambda_v", deriv->lambda_v);
      } else {
        model_set_param_deriv(mdl, "lambda", deriv->lambda);
      }
    }

  }
  else {
    g_error("bkgd_evo_mdl_set_mdl_deriv: unknown mutation model type");
  }

  if(conf->mdl_type == MDL_TYPE_HCM) {
    model_set_param_deriv(mdl, "N_hcm", deriv->N_hcm);
    model_set_param_deriv(mdl, "T_hcm", deriv->T_hcm);
  } 
  else if(conf->mdl_type == MDL_TYPE_HCGOM) {
    model_set_param_deriv(mdl, "N_hcg", deriv->N_hcg);
    model_set_param_deriv(mdl, "N_hcgo", deriv->N_hcgo);
    model_set_param_deriv(mdl, "N_hcgom", deriv->N_hcgom);
    model_set_param_deriv(mdl, "T_hcg", deriv->T_hcg);
    model_set_param_deriv(mdl, "T_hcgo", deriv->T_hcgo);
    model_set_param_deriv(mdl, "T_hcgom", deriv->T_hcgom);

  } else if(conf->mdl_type == MDL_TYPE_HCOM) {
    model_set_param_deriv(mdl, "N_hco", deriv->N_hco);
    model_set_param_deriv(mdl, "N_hcom", deriv->N_hcom);
    model_set_param_deriv(mdl, "T_hco", deriv->T_hco);
    model_set_param_deriv(mdl, "T_hcom", deriv->T_hcom);

  } else if(conf->mdl_type == MDL_TYPE_HCGO) {
    model_set_param_deriv(mdl, "N_hcg", deriv->N_hcg);
    model_set_param_deriv(mdl, "N_hcgo", deriv->N_hcgo);
    model_set_param_deriv(mdl, "T_hcg", deriv->T_hcg);
    model_set_param_deriv(mdl, "T_hcgo", deriv->T_hcgo);
  } else if(conf->mdl_type != MDL_TYPE_HC) {
    g_error("set deriv: unknown model type");
  }
}



typedef struct {
  long n;
  double ttl_prob;
  double ttl_prob_single;
  double ttl_prob_double;
  double ttl_ll;
} ColTypeStat;

/**
 * Calculates the log likelihood (log probability of model given the
 * data) and gradient. If the calc_ll flag is FALSE, the
 * log-likelihood is not calculated and 0.0 is returned. If the
 * calc_grad flag is FALSE the gradient is not calculated and the
 * model partial derivatives are not updated.
 */
double bkgd_evo_mdl_calc_ll(Model *mdl, int calc_ll, int calc_grad) {
  BkgdEvoMdlParam param, ll_deriv;
  Branch *br;
  ColType *cltype;
  BkgdEvoMdlConfig *conf;
  BkgdEvoMdlData *data;
  double ll;
  ColTypeStat *clstat, *cons_clstat;
  long i;
  int j;
   
  conf = mdl->config;
  data = mdl->data;

  if(data == NULL) {
    g_error("bkgd_evo_mdl_calc_ll: data must be set before "
	    "likelihood can be calculated\n");
  }

  param_set_zero(&param);
  param_set_zero(&ll_deriv);
  set_bkgd_evo_param_from_model(mdl, conf, &param);

  ll = 0.0;

  clstat  = g_new(ColTypeStat, conf->n_cltype);
  cons_clstat = NULL;
  for(i = 0; i < conf->n_cltype; i++) {
    clstat[i].n = 0;
    clstat[i].ttl_prob = 0.0;
    clstat[i].ttl_prob_double = 0.0;
    clstat[i].ttl_prob_single = 0.0;
    clstat[i].ttl_ll = 0.0;
  }

  for(i = 0; i < data->n_bin; i++) {
    /* Combine non-exonic and exonic B vals and rescale by new
     * deleterious rate estimates. Use precomputed log B values so
     * that rescaling is more efficient (do not need to take logs
     * first and then re-exponentiate)
     */ 
    param.B = exp(data->bin[i].lB_ex  * param.u_ex_scale + 
		  data->bin[i].lB_nex * param.u_nex_scale);

    if(param.B < 0.0 || param.B > 1.0 || isnan(param.B) || isinf(param.B)) {
      g_error("bkgd_evo_mdl_calc_ll: invalid B (%g)\n"
	      "B_ex=%g, B_nex=%g, u_ex_scale=%g, u_nex_scale=%g\n",
	      param.B, data->bin[i].B_ex, data->bin[i].B_nex,
	      param.u_ex_scale, param.u_nex_scale);
    }

    /* calculate mutation rate for current bin */
    bkgd_evo_mdl_calc_mu(&param, mdl, &data->bin[i]);

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
	g_warning("bkgd_evo_mdl_calc_ll: invalid length (%g) for branch %s",
		br->len, br->name);
	br->len = BRANCH_LEN_SMALL;
      }

      branch_set_prob(br, &param, conf);
    }

    /* calc prob of observing each observed column type for bin */
    for(j = 0; j < conf->n_cltype; j++) {
      cltype = conf->cltypes[j];

      /* do only non-conserved columns first */
      if(cltype->subst_type == SUBST_TYPE_CONSERVED) {
	cons_clstat = &clstat[j];
	continue;
      }

      cltype->set_n(cltype, &data->bin[i]);

      /* compute prob of a single observed column */      
      cltype_set_prob(cltype, &param, conf);
      
      if(cltype->n > 0) {
	/* add contribution to log-likelihood */
	if(cltype->prob > 1.0 || cltype->prob < 0.0) {
	  g_error("bkgd_evo_mdl_calc_ll: invalid probability (%g) "
		  "for column type '%s'", cltype->prob, cltype->name);
	}

	if(cltype->prob < 1e-20) {
	  fprintf(stderr, "LOW prob column %s: p=%g\n", cltype->name, 
		  cltype->prob);
	}

	ll += (double)cltype->n * log(cltype->prob);

	/* for debugging and analysis purposes record some statistics
	 * about each column type
	 */
	clstat[j].n += cltype->n;
	clstat[j].ttl_prob += cltype->prob;
	clstat[j].ttl_prob_double += cltype->prob_ttl_double;
	clstat[j].ttl_prob_single += cltype->prob_single;
	clstat[j].ttl_ll += (double)cltype->n * log(cltype->prob);
      }
    }

    /* now add the contribution from the conserved column */
    if(conf->cons_cltype != NULL) {
      cltype = conf->cons_cltype;
      cltype->set_n(cltype, &data->bin[i]);

      /* compute prob of a single observed column */
      cltype_cons_set_prob(cltype, &param, conf);

      if(cltype->n > 0) {      
	/* add contribution to log-likelihood */
	if(cltype->prob > 1.0 || cltype->prob < 0.0) {
	  g_error("bkgd_evo_mdl_calc_ll: invalid probability (%g) "
		  "for column type '%s'", cltype->prob, cltype->name);
	}

	if(cltype->prob < 1e-20) {
	  fprintf(stderr, "LOW prob column %s: p=%g\n", cltype->name, 
		  cltype->prob);
	}

	ll += (double)cltype->n * log(cltype->prob);

	/* fprintf(stderr, "bin=%ld, CONS prob=%g\n", i, cltype->prob);*/
	cons_clstat->n += cltype->n;
	cons_clstat->ttl_prob += cltype->prob;
	cons_clstat->ttl_prob_double += cltype->prob_ttl_double;
	cons_clstat->ttl_prob_single += cltype->prob_single;
	cons_clstat->ttl_ll += (double)cltype->n * log(cltype->prob);
      }
    }

    if(calc_grad) {
      /* perform derivative calculations **/
      calc_deriv( &ll_deriv, conf, &param,  &data->bin[i]);
    }
  }

  if(calc_grad) {
    /* update model's derivatives */
    bkgd_evo_mdl_set_mdl_deriv(mdl, &ll_deriv);    
    /* model_write_grad_ln(mdl, stderr);*/
  }


/*   for(i = 0; i < conf->n_cltype; i++) { */
/*     double dbl_pp, single_pp, avg_p; */
/*     char *name; */
    
/*     name = conf->cltypes[i]->name; */
/*     avg_p = clstat[i].ttl_prob / (double)clstat[i].n; */
/*     single_pp = clstat[i].ttl_prob_single / clstat[i].ttl_prob; */
/*     dbl_pp = clstat[i].ttl_prob_double / clstat[i].ttl_prob; */

/*     fprintf(stderr, "  %s: n=%ld, LL=%g, p(%s)=%g, " */
/* 	    "p(single|%s)=%g, p(double|%s)=%g\n", */
/*  	    name, clstat[i].n, clstat[i].ttl_ll, name, avg_p, */
/* 	    name, single_pp, name, dbl_pp); */
/*   } */

  g_free(clstat);
  
  return ll;
}







/**
 * Reads model type from configuration
 */
static int get_mdl_type(Config *config) {
  char *mdl_type;

  mdl_type = config_get_str(config, "MODEL_TYPE");
  
  if(strcmp(mdl_type, "HC")==0) {
    return MDL_TYPE_HC;
  }
  if(strcmp(mdl_type, "HCM")==0) {
    return MDL_TYPE_HCM;
  }
  if(strcmp(mdl_type, "HCOM")==0) {
    return MDL_TYPE_HCOM;
  }
  if(strcmp(mdl_type, "HCGOM")==0) {
    return MDL_TYPE_HCGOM;
  }
  if(strcmp(mdl_type, "HCGO")==0) {
    return MDL_TYPE_HCGO;
  }

  g_error("get_mdl_type: model type should be HC or HCM not '%s'",
	  mdl_type);

  return MDL_TYPE_HC;
}


/**
 * Reads the substitution model type from configuration
 */
static SubstModel *subst_mdl_new(Config *config) {
  char *mdl_type;
  SubstModel *smdl;

  smdl = g_new(SubstModel, 1);

  mdl_type = config_get_str(config, "SUBST_MODEL_TYPE");

  smdl->name = g_strdup(mdl_type);
  
  if(strcmp(mdl_type, "JUKES_CANTOR")==0) {
    smdl->id = SUBST_MDL_JUKES_CANTOR;
    smdl->n_subst_type = 1;
    smdl->subst_type_ids = g_new(int, 1);
    smdl->subst_type_ids[0] = SUBST_TYPE_SUBST;
    smdl->subst_type_names = g_new(char *, 1);
    smdl->subst_type_names[0] = NULL;
  }
  else if(strcmp(mdl_type, "KIMURA")==0) {
    smdl->id = SUBST_MDL_KIMURA;
    smdl->n_subst_type = 2;
    smdl->subst_type_ids = g_new(int, 2);
    smdl->subst_type_ids[0] = SUBST_TYPE_TRANSITION;
    smdl->subst_type_ids[1] = SUBST_TYPE_TRANSVERSION;
    smdl->subst_type_names = g_new(char *, 2);
    smdl->subst_type_names[0] = g_strdup("I");
    smdl->subst_type_names[1] = g_strdup("V");
  } else {
    g_error("get_mdl_type: SUBST_MODEL_TYPE should be JUKES_CANTOR "
	    "or KIMURA not '%s'", mdl_type);
  }

  return smdl;
}




/**
 * Reads mutation model type from configuration
 */
static int get_mu_type(Config *config) {
  char *mdl_type;

  mdl_type = config_get_str(config, "MU_MODEL_TYPE");
  
  if(strcmp(mdl_type, "SINGLE")==0) {
    return MU_TYPE_SINGLE;
  }
  if(strcmp(mdl_type, "CAT_QUAD")==0) {
    return MU_TYPE_CAT_QUAD;
  }
  if(strcmp(mdl_type, "CAT_LIN")==0) {
    return MU_TYPE_CAT_LIN;
  }

  g_error("get_mdl_type: mutation model type should be SINGLE, CAT_QUAD, "
	  "or CAT_LIN not '%s'", mdl_type);

  return MU_TYPE_SINGLE;
}



/**
 * Initializes configuration flags for model
 */
BkgdEvoMdlConfig *bkgd_evo_mdl_config_new(Config *config) {
  BkgdEvoMdlConfig *bem_conf;

  bem_conf = g_new(BkgdEvoMdlConfig,1);
  
  bem_conf->mdl_type = get_mdl_type(config);
  bem_conf->subst_mdl = subst_mdl_new(config);
  bem_conf->mu_type  = get_mu_type(config);
    
  if(config_has_key(config, "MIN_CAT")) {
    bem_conf->min_cat = config_get_double(config, "MIN_CAT");
    fprintf(stderr, "only considering bins with CAT >= %g\n", 
	    bem_conf->min_cat);
  } else {
    bem_conf->min_cat = 0.0;
  }

  if(config_has_key(config, "MAX_CAT")) {
    bem_conf->max_cat = config_get_double(config, "MAX_CAT");
    fprintf(stderr, "only considering bins with CAT <= %g\n", 
	    bem_conf->max_cat);
  } else {
    bem_conf->max_cat = 1.0;
  }

  bem_conf->use_double_substs = config_get_boolean(config, "USE_DOUBLE_SUBSTS");

  bem_conf->n_branch = 0;
  bem_conf->n_cltype = 0;
  bem_conf->cons_cltype = NULL;

  return bem_conf;
}



/**
 * Frees memory for model and associated data and configuration
 * information
 */
void bkgd_evo_mdl_free(Model *mdl) {
  BkgdEvoMdlData *data;

  g_free(mdl->config);

  data = mdl->data;
  if(data != NULL) {
    g_free(data->bin);
    g_free(data);
  }

  model_free(mdl);
}


void bkgd_evo_mdl_free_data(Model *mdl) {
  BkgdEvoMdlData *data;

  data = mdl->data;

  g_free(data->bin);
  g_free(data);

  mdl->data = NULL;
}



/**
 * Adds a parameter to the current model
 */
void add_param(Config *config, Model *mdl, const int allow_neg, 
	       const char *name) {
  char *uc_name;
  char *key;
  int is_locked;
  double val, scale;
  ModelParam *param;

  /* convert name to upper case for use as lookup key in config */
  uc_name = g_strdup(name);
  util_str_uc(uc_name);

  /* check config to see if parameter should be locked */
  key = g_strconcat("LOCK_", uc_name, NULL);
  if(config_has_key(config, key)) {     
    is_locked = config_get_boolean(config, key);
  } else {
    is_locked = FALSE;
  }
  g_free(key);

  /* get initial parameter value */
  key = g_strconcat("PARAM_", uc_name, NULL);
  val = config_get_double(config, key);
  g_free(key);

  scale = 1.0/val;

  /* add parameter to model */
  param = model_add_param(mdl, name, is_locked, allow_neg, scale);
  model_set_param_val(mdl, name, val);

  /* check if there is an upper bound on the parameter value */
  key = g_strconcat("UP_BOUND_", uc_name, NULL);
  if(config_has_key(config, key)) {
    param->up_bound = TRUE;
    param->max_val = config_get_double(config, key);

    fprintf(stderr, "using upper bound of %g for parameter '%s'\n", 
	    param->max_val, name);
  }
  g_free(key);


  g_free(uc_name);
}



char *get_clname(char *name, char *stype) {
  if(stype) {
    return g_strconcat(name, "_", stype, NULL);
  }
  return g_strdup(name);
}


void hc_init(Config *config, Model *mdl) {
  ColType *cltype;
  Branch *br;
  BkgdEvoMdlConfig *bem_conf;
  SubstModel *subst_mdl;
  char *clname;
  int i;

  bem_conf = mdl->config;
  subst_mdl = bem_conf->subst_mdl;

  if(config_get_boolean(config, "USE_DOUBLE_SUBSTS")) {
    g_error("hc_init: USE_DOUBLE_SUBSTS config option should "
            "be set to FALSE for HC model");
  }

  add_param(config, mdl, FALSE, "T_hc");
  add_param(config, mdl, FALSE, "N_hc");

  /* create H+C branch */
  br = add_branch(mdl->config, "H+C", &branch_set_len__hc_h_c,
		  &branch_set_dlen__hc_h_c);


  for(i = 0; i < subst_mdl->n_subst_type; i++) {
    /* create H+C column types */
    clname = get_clname("H+C", subst_mdl->subst_type_names[i]);
    cltype = cltype_new(clname, subst_mdl->subst_type_ids[i], 
			cltype_set_n__hc_h_c);
    cltype_set_br(cltype, br);
    add_cltype(mdl->config, cltype);
  }
  
  /* create conserved column type */
  cltype = cltype_new("CONS", SUBST_TYPE_CONSERVED, cltype_set_n__hc_cons);
  add_cltype(mdl->config, cltype);
}




void hcm_init(Config *config, Model *mdl) {
  ColType *cltype;
  Branch *h_br, *m_br;
  BkgdEvoMdlConfig *bem_conf;
  SubstModel *subst_mdl;
  char *clname;
  int i;

  bem_conf = mdl->config;
  subst_mdl = bem_conf->subst_mdl;

  add_param(config, mdl, FALSE, "T_hc");
  add_param(config, mdl, FALSE, "N_hc");
  add_param(config, mdl, FALSE, "T_hcm");
  add_param(config, mdl, FALSE, "N_hcm");

  /* create H+C branch */
  h_br = add_branch(mdl->config, "H+C", &branch_set_len__hcm_h_c,
		  &branch_set_dlen__hcm_h_c);

  /* create M branch */
  m_br = add_branch(mdl->config, "M", &branch_set_len__hcm_m,
		  &branch_set_dlen__hcm_m);


  for(i = 0; i < subst_mdl->n_subst_type; i++) {
    /* create H+C column type */
    clname = get_clname("H+C", subst_mdl->subst_type_names[i]);
    cltype = cltype_new(clname, subst_mdl->subst_type_ids[i], 
			cltype_set_n__hcm_h_c);
    cltype_set_br(cltype, h_br);
    add_cltype(mdl->config, cltype);
    g_free(clname);
  
    /* create M column type */
    clname = get_clname("M",  subst_mdl->subst_type_names[i]);
    cltype = cltype_new(clname, subst_mdl->subst_type_ids[i], 
			cltype_set_n__hcm_m);
    cltype_set_br(cltype, m_br);
    add_cltype(mdl->config, cltype);
  }
  
  /* create conserved column type */
  cltype = cltype_new("CONS", SUBST_TYPE_CONSERVED, cltype_set_n__hcm_cons);
  add_cltype(mdl->config, cltype);
}



void hcom_init(Config *config, Model *mdl) {
  ColType *cltype;
  Branch *h_br, *hc_br, *o_br, *m_br;
  BkgdEvoMdlConfig *bem_conf;
  SubstModel *subst_mdl;
  int i, stype;

  bem_conf = mdl->config;
  subst_mdl = bem_conf->subst_mdl;

  add_param(config, mdl, FALSE, "T_hc");
  add_param(config, mdl, FALSE, "N_hc");
  add_param(config, mdl, FALSE, "T_hco");
  add_param(config, mdl, FALSE, "N_hco");
  add_param(config, mdl, FALSE, "T_hcom");
  add_param(config, mdl, FALSE, "N_hcom");

  /* create H+C branch */
  h_br = add_branch(mdl->config, "H+C", &branch_set_len__hcom_h_c,
		  &branch_set_dlen__hcom_h_c);


  /* create HC branch */
  hc_br = add_branch(mdl->config, "HC", &branch_set_len__hcom_hc,
		  &branch_set_dlen__hcom_hc);


  /* create O branch */
  o_br = add_branch(mdl->config, "O", &branch_set_len__hcom_o,
		  &branch_set_dlen__hcom_o);


  /* create M branch */
  m_br = add_branch(mdl->config, "M", &branch_set_len__hcm_m,
		  &branch_set_dlen__hcom_m);


  for(i = 0; i < subst_mdl->n_subst_type; i++) {
    stype = subst_mdl->subst_type_ids[i];

    /* create H+C column type */
    cltype = cltype_new("H+C", stype, cltype_set_n__hcom_h_c);
    cltype_set_br(cltype, h_br);
    add_cltype(mdl->config, cltype);

    /* create HC column type */
    cltype = cltype_new("HC", stype, cltype_set_n__hcom_hc);
    cltype_set_br(cltype, hc_br);
    add_cltype(mdl->config, cltype);
    /* TODO: could add double subst branches to HC col type */

    /* create O column type */
    cltype = cltype_new("O", stype, cltype_set_n__hcom_o);
    cltype_set_br(cltype, o_br);
    add_cltype(mdl->config, cltype);
  
    /* create M column type */
    cltype = cltype_new("M", stype, cltype_set_n__hcom_m);
    cltype_set_br(cltype, m_br);
    add_cltype(mdl->config, cltype);
  }

  /* TODO: could add HO and CO col types */
  
  /* create conserved column type */
  cltype = cltype_new("CONS", SUBST_TYPE_CONSERVED, cltype_set_n__hcom_cons);
  add_cltype(mdl->config, cltype);
}





static void add_hcgo_branches(Config *config, Model *mdl) {
  ColType *cltype;
  Branch *h_br, *c_br, *g_br, *o_br, *hc_br, *hg_br, *cg_br, *hcg_br;

  BkgdEvoMdlConfig *bem_conf;
  SubstModel *subst_mdl;
  char *clname;
  int i, stype;

  bem_conf = mdl->config;
  subst_mdl = bem_conf->subst_mdl;

  /* add parameters to mdl */
  add_param(config, mdl, FALSE, "T_hc");
  add_param(config, mdl, FALSE, "T_hcg");
  add_param(config, mdl, FALSE, "T_hcgo");

  add_param(config, mdl, FALSE, "N_hc");
  add_param(config, mdl, FALSE, "N_hcg");
  add_param(config, mdl, FALSE, "N_hcgo");


  /* add branch definitions */
  fprintf(stderr, "  adding branches\n");
  h_br = add_branch(mdl->config, "H", branch_set_len__hcgom_h,
		    branch_set_dlen__hcgom_h);

  c_br = add_branch(mdl->config, "C", branch_set_len__hcgom_c,
		    branch_set_dlen__hcgom_c);

  g_br = add_branch(mdl->config, "G", branch_set_len__hcgom_g,
		    branch_set_dlen__hcgom_g);

  hc_br = add_branch(mdl->config, "HC", branch_set_len__hcgom_hc,
		     branch_set_dlen__hcgom_hc);

  hg_br = add_branch(mdl->config, "HG", branch_set_len__hcgom_hg,
			branch_set_dlen__hcgom_hg);

  cg_br = add_branch(mdl->config, "CG", branch_set_len__hcgom_cg,
			branch_set_dlen__hcgom_cg);

  hcg_br = add_branch(mdl->config, "HCG", branch_set_len__hcgom_hcg,
		      branch_set_dlen__hcgom_hcg);
  
  o_br = add_branch(mdl->config, "O", branch_set_len__hcgom_o,
		    branch_set_dlen__hcgom_o);

  /* add column types */

  fprintf(stderr, "  adding column types\n");
  for(i = 0; i < subst_mdl->n_subst_type; i++) {
    stype = subst_mdl->subst_type_ids[i];

    /*** H col ***/
    clname = get_clname("H", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_h);
    cltype_set_br(cltype, h_br);
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, hc_br, c_br);
      cltype_add_dbl_br(cltype, hg_br, g_br);
      cltype_add_dbl_br(cltype, hcg_br, cg_br);
    }
    add_cltype(mdl->config, cltype);
    g_free(clname);

    /*** C col ***/
    clname = get_clname("C", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_c);
    cltype_set_br(cltype, c_br);
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, hc_br, h_br);
      cltype_add_dbl_br(cltype, cg_br, g_br);  
      cltype_add_dbl_br(cltype, hcg_br, hg_br);
    }
    add_cltype(mdl->config, cltype);
    g_free(clname);

    /*** G col ***/
    clname = get_clname("G",subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_g);
    cltype_set_br(cltype, g_br);
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, hg_br, h_br);
      cltype_add_dbl_br(cltype, cg_br, c_br);
      cltype_add_dbl_br(cltype, hcg_br, hc_br);
    }
    add_cltype(mdl->config, cltype);
    g_free(clname);

    /*** HC col ***/    
    clname = get_clname("HC", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_hc);
    cltype_set_br(cltype, hc_br);  
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, h_br, c_br);
      cltype_add_dbl_br(cltype, hcg_br, g_br);  
    }
    add_cltype(mdl->config, cltype);
    g_free(clname);


    /*** HG col ***/
    clname = get_clname("HG", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_hg);
    cltype_set_br(cltype, hg_br);
    g_free(clname);
  
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, h_br, g_br);
      cltype_add_dbl_br(cltype, hcg_br, c_br);
    }
    add_cltype(mdl->config, cltype);

    /*** CG col ***/
    clname = get_clname("CG", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_cg);
    cltype_set_br(cltype, cg_br);
    g_free(clname);
  
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, c_br, g_br);
      cltype_add_dbl_br(cltype, hcg_br, h_br);
    }
    add_cltype(mdl->config, cltype);


    /*** HCG col ***/
    clname = get_clname("HCG", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_hcg);
    cltype_set_br(cltype, hcg_br);
    if(bem_conf->use_double_substs) {
      cltype_add_dbl_br(cltype, hc_br, g_br);
      cltype_add_dbl_br(cltype, cg_br, h_br);
      cltype_add_dbl_br(cltype, hg_br, c_br);
    }
    add_cltype(mdl->config, cltype);
    g_free(clname);

    /*** O col ***/
    clname = get_clname("O", subst_mdl->subst_type_names[i]);
    fprintf(stderr, "creating %s column\n", clname);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_o);
    cltype_set_br(cltype, o_br);
    add_cltype(mdl->config, cltype);
    g_free(clname);

    if(bem_conf->use_double_substs) {
      /*** HO and CO cols ***/
      clname = get_clname("HO+CO", subst_mdl->subst_type_names[i]);
      fprintf(stderr, "creating %s column\n", clname);
      cltype = cltype_new(clname, stype, cltype_set_n__hcgom_ho_co);
      cltype_add_dbl_br(cltype, h_br, o_br);
      cltype_add_dbl_br(cltype, c_br, o_br);
      add_cltype(mdl->config, cltype);
      g_free(clname);
    }
  }  
  /* TODO: could also have GO col type, etc. */


  fprintf(stderr, "done adding hcgo column types\n");
}



void hcgo_init(Config *config, Model *mdl) {
  ColType *cltype;

  add_hcgo_branches(config, mdl);

  /*** cons cols ***/
  cltype = cltype_new("CONS", SUBST_TYPE_CONSERVED, cltype_set_n__hcgo_cons);
  add_cltype(mdl->config, cltype);

}


/**
 *
 *
 */
void hcgom_init(Config *config, Model *mdl) {
  ColType *cltype, *o_cltype;
  Branch *m_br, *o_br, *hcg_br;
  BkgdEvoMdlConfig *bem_conf;
  SubstModel *subst_mdl;
  char *clname, *o_cltype_name;
  int i, stype, use_double_substs;

  use_double_substs = config_get_boolean(config, "USE_DOUBLE_SUBSTS");

  fprintf(stderr, "creating HCGOM model\n");
  bem_conf = mdl->config;
  subst_mdl = bem_conf->subst_mdl;


  add_hcgo_branches(config, mdl);

  /* add parameters to mdl */
  add_param(config, mdl, FALSE, "T_hcgom");
  add_param(config, mdl, FALSE, "N_hcgom");

  /* add M branch */
  m_br = add_branch(mdl->config, "M", branch_set_len__hcgom_m,
		    branch_set_dlen__hcgom_m);

  /* add column types */

  for(i = 0; i < subst_mdl->n_subst_type; i++) {
    /*** M col ***/
    stype = subst_mdl->subst_type_ids[i];
    clname = get_clname("M", subst_mdl->subst_type_names[i]);
    cltype = cltype_new(clname, stype, cltype_set_n__hcgom_m);
    cltype_set_br(cltype, m_br);
    if(use_double_substs) {      
      o_br = get_branch(bem_conf, "O");
      hcg_br = get_branch(bem_conf, "HCG");
      cltype_add_dbl_br(cltype, hcg_br, o_br);

      /* also add double subst for O branch */
      o_cltype_name = get_clname("O", subst_mdl->subst_type_names[i]);
      o_cltype = get_cltype(bem_conf, o_cltype_name);
      cltype_add_dbl_br(o_cltype, hcg_br, m_br);
      g_free(o_cltype_name);
    }
    add_cltype(mdl->config, cltype);
    g_free(clname);
  }


  /*** cons cols ***/
  cltype = cltype_new("CONS", SUBST_TYPE_CONSERVED, cltype_set_n__hcgom_cons);
  add_cltype(mdl->config, cltype);

  fprintf(stderr, "done creating HCGOM model\n");
}



void bkgd_evo_mdl_read_data(Model *mdl, Config *config) {
  BkgdEvoMdlConfig *bem_conf;
  BkgdEvoMdlData *bem_data;

  bem_conf = mdl->config;

  if(bem_conf->subst_mdl->id == SUBST_MDL_KIMURA) {
    /* read separate transition/transversion files */
    bem_data = bkgd_data_read_i_v_data(config);
  } else {
    bem_data = bkgd_data_read_data(config);
  }

  if(config_get_boolean(config, "JUKES_CANTOR_CORRECT_CATS")) {
    fprintf(stderr, "applying jukes-cantor correction to category values");
    bkgd_data_jukes_cantor_correct_cats(bem_data);
  }

  if(config_get_boolean(config, "IGNORE_B_NEX")) {
    fprintf(stderr, "collapsing B_nex values");
    bkgd_data_collapse_nex_bin(bem_data);
    model_set_param_val(mdl, "u_nex_scale", 1.0);
    model_set_param_val(mdl, "t_nex", 0.0);
    model_lock_param(mdl, "u_nex_scale");
    model_lock_param(mdl, "t_nex");
  }


  if(bem_conf->mu_type == MU_TYPE_SINGLE || 
     config_get_boolean(config, "USE_MDIV_CATS")) {
    /* combine category bins, since they will not be used */
    bkgd_data_combine_cat_bin(bem_data, bem_conf->min_cat, bem_conf->max_cat);
  }

  mdl->data = bem_data;

  bkgd_data_write_site_counts(stderr, mdl->data);
}



void bkgd_evo_mdl_bin_data(Model *mdl, Config *config) {
  long n_bin;
  BkgdEvoMdlData *bem_data;

  bem_data = mdl->data;
  n_bin = config_get_long(config, "N_DATA_BIN");
  bkgd_data_combine_bin(bem_data, n_bin);

  if(config_get_boolean(config, "USE_MDIV_CATS")) {
    /* use macaque divergence as category values ? */
    fprintf(stderr, "using macaque divergence as category values");
    bkgd_data_set_mdiv_cats(bem_data);
  }
}


/**
 * Creates and initializes a background selection evolutionary model
 * and associated parameters and data using information read from the
 * provided configuration.
 */
Model *bkgd_evo_mdl_new(Config *config) {
  Model *mdl;
  BkgdEvoMdlConfig *bem_conf;

  /* TODO: there should be a way to free the config when ML done */
  bem_conf = bkgd_evo_mdl_config_new(config);

  mdl = model_new();
  mdl->config = bem_conf;
  mdl->data = NULL;

  /* set selection coefficient parameters */
  add_param(config, mdl, FALSE, "t_ex");
  add_param(config, mdl, FALSE, "t_nex");


  /* set mutation rate parameter(s) */
  switch(bem_conf->mu_type) {
  case(MU_TYPE_SINGLE):
    if(bem_conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      /* separate mutation rates for transitions and transversions */
      add_param(config, mdl, FALSE, "mu_i");
      add_param(config, mdl, FALSE, "mu_v");
      if(bem_conf->use_double_substs) {
        add_param(config, mdl, FALSE, "lambda_i");
        add_param(config, mdl, FALSE, "lambda_v");
      }
    } else {
      add_param(config, mdl, FALSE, "mu");
      if(bem_conf->use_double_substs) {
        add_param(config, mdl, FALSE, "lambda");
      }
    }
    break;

  case(MU_TYPE_CAT_QUAD):
    if(bem_conf->subst_mdl->id == SUBST_MDL_KIMURA) {
      g_error("bkgd_evo_mdl_new: CAT_QUAD mutation model not implemented "
	      "for KIMURA substitution model");
    }

    /* three separate mutation rate parameters define quadratic
     * relationship with category values (e.g. GC content).  both mu_a
     * and mu_b coefficients are allowed to be negative.
     */
    add_param(config, mdl, TRUE, "mu_a");
    add_param(config, mdl, TRUE, "mu_b");
    add_param(config, mdl, TRUE, "mu_c");

    if(bem_conf->use_double_substs) {
      add_param(config, mdl, TRUE, "lambda");
    }
    break;

  case(MU_TYPE_CAT_LIN): 
    /* two mutation rate parameters define linear relationship with
     * category values (e.g. GC content).  
     */
    add_param(config, mdl, TRUE, "mu_a");
    add_param(config, mdl, TRUE, "mu_b");

    if(bem_conf->use_double_substs) {
      if(bem_conf->subst_mdl->id == SUBST_MDL_KIMURA) {
        /* separate mutation rates for transitions and transversions */
        add_param(config, mdl, FALSE, "lambda_i");
        add_param(config, mdl, FALSE, "lambda_v");
      } else {
        add_param(config, mdl, FALSE, "lambda");
      }
    }
    break;

  default:
    g_error("bkgd_evo_mdl_new: unknown mu model type\n");
  }

  /* set deleterious rate scaling parameters */
  add_param(config, mdl, FALSE, "u_ex_scale");
  add_param(config, mdl, FALSE, "u_nex_scale");

  /* set speciation time and ancestral effective pop size parameter(s) */
  switch(bem_conf->mdl_type) {
  case(MDL_TYPE_HC):
    hc_init(config, mdl);
    break;

  case(MDL_TYPE_HCM):
    hcm_init(config, mdl);
    break;

  case(MDL_TYPE_HCOM):
    hcom_init(config, mdl);
    break;

  case(MDL_TYPE_HCGOM):
    hcgom_init(config, mdl);
    break;

  case(MDL_TYPE_HCGO):
    hcgo_init(config, mdl);
    break;
    
  default:
    g_error("bkgd_evo_mdl_new: unknown model type");
  }

  model_write_param_ln(mdl, stderr);

  mdl->llhood = &bkgd_evo_mdl_calc_ll;

  return mdl;
}

