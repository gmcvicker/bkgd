
#include "branch.h"
#include "bkgd_evo_mdl.h"
#include <math.h>



void branch_set_prob(Branch *br, const BkgdEvoMdlParam *param,
		     const BkgdEvoMdlConfig *conf) {

  if(br->len < 0.0) {
    g_error("set_jc_prob: branch length (%g) for branch '%s' is "
	    " < 0.0", br->len, br->name);
  }

    
  /* calculate prob of observing substitution on branch using model of
   * DNA evolution
   */
  if(conf->subst_mdl->id == SUBST_MDL_JUKES_CANTOR) {
    /* jukes-cantor model: all substitution probs equal */
    br->jc = exp(-(4.0/3.0) * br->len * param->mu);
    br->prob = 0.75 * (1.0 - br->jc);
  }
  else if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
    /* kimura two parameter model: different probability of transition
     * and transversion substitutions
     */
    br->k_i = exp(-2.0 * (param->mu_i + param->mu_v) * br->len);
    br->k_v = exp(-4.0 * param->mu_v * br->len);

    br->prob_i = 0.25 - 0.5*br->k_i + 0.25*br->k_v;
    br->prob_v = 0.5 - 0.5*br->k_v;

    if(br->prob_i < BRANCH_LOW_PROB) {
     /*fprintf(stderr, "LOW prob_i for branch %s: %g\n", br->name, br->prob_i);*/
      br->prob_i = BRANCH_LOW_PROB;
    }
    if(br->prob_v < BRANCH_LOW_PROB) {
      /*fprintf(stderr, "LOW prob_v for branch %s: %g\n", br->name,br->prob_v);*/
      br->prob_v = BRANCH_LOW_PROB;
    }
				      

    br->prob = br->prob_i + br->prob_v;
  }
  else {
    g_error("branch_set_prob: unknown substitution model type");
  }

  if(br->prob < 0.0 || br->prob > 1.0 || isnan(br->prob) || isinf(br->prob)) {
    g_error("branch_set_prob: invalid probability (%g) for branch %s\n", 
	    br->prob, br->name);
  }

  /* don't allow branch probability to reach 0 (can because of
   * roundoff error on very short branches)
   */
  if(br->prob < BRANCH_LOW_PROB) {
    fprintf(stderr, "LOW prob for branch %s: %g\n", br->name,br->prob); 
    br->prob = BRANCH_LOW_PROB;
  }

}



/**
 * Sets partial derivatives of branch substituition probabilities for
 * jukes-cantor model.  Only sets dprob attributes, does not set
 * transition/transversion-specific dprob_i/dprob_v attributes.
 */
static void branch_set_dprob_jukes_cantor(Branch *br, const BkgdBin *bin, 
					  const BkgdEvoMdlParam *param, 
					  const BkgdEvoMdlConfig *conf) {
  double coef;

  /* calculate partial derivs of branch mutation probs (p)
   * using branch len partials
   */

  /* initialize prob deriv structure to 0 */
  param_set_zero(&br->dprob);

  coef = br->jc * param->mu;

  /* derivs w.r.t. pop size, speciation times and deleterious rates,
   * all depend on branch len derivs
   */
  br->dprob.N_hc    = coef * br->dlen.N_hc;
  br->dprob.N_hcg   = coef * br->dlen.N_hcg;
  br->dprob.N_hcgo  = coef * br->dlen.N_hcgo;
  br->dprob.N_hcgom = coef * br->dlen.N_hcgom;
  br->dprob.T_hc    = coef * br->dlen.T_hc;
  br->dprob.T_hcg   = coef * br->dlen.T_hcg;
  br->dprob.T_hcgo  = coef * br->dlen.T_hcgo;
  br->dprob.T_hcgom = coef * br->dlen.T_hcgom;

  /* HCM model params */
  br->dprob.N_hcm   = coef * br->dlen.N_hcm;
  br->dprob.T_hcm   = coef * br->dlen.T_hcm;

  /* HCOM model params */
  br->dprob.N_hco   = coef * br->dlen.N_hco;
  br->dprob.T_hco   = coef * br->dlen.T_hco;
  br->dprob.N_hcom  = coef * br->dlen.N_hcom;
  br->dprob.T_hcom  = coef * br->dlen.T_hcom;

  br->dprob.u_ex_scale  = coef * br->dlen.u_ex_scale;
  br->dprob.u_nex_scale = coef * br->dlen.u_nex_scale;

  /* derivs w.r.t. mu param do not depend on branch len derivs */
  if(conf->mu_type == MU_TYPE_SINGLE) {
    br->dprob.mu = br->len * br->jc;
    br->dprob.mu_a = 0.0;
    br->dprob.mu_b = 0.0;
    br->dprob.mu_c = 0.0;
  } 
  else if(conf->mu_type == MU_TYPE_CAT_LIN) {
    br->dprob.mu = 0.0;
    br->dprob.mu_c = 0.0;
    br->dprob.mu_b = br->len * br->jc;
    br->dprob.mu_a = br->dprob.mu_b * bin->cat;
  }
  else if(conf->mu_type == MU_TYPE_CAT_QUAD) {
    br->dprob.mu = 0.0;
    br->dprob.mu_c = br->len * br->jc;
    br->dprob.mu_b = br->dprob.mu_c * bin->cat;
    br->dprob.mu_a = br->dprob.mu_b * bin->cat;
  }
  else {
    g_error("branch_set_dprob_jukes_cantor: unknown mutation model type");
  }
}





/**
 * Sets partial derivatives of branch TRANSITION probabilities for
 * the kimura substition model.  Only sets dprob_i attributes, does
 * not set dprob or dprob_v attributes.
 */
static void branch_set_dprob_i_kimura(Branch *br, const BkgdBin *bin, 
					  const BkgdEvoMdlParam *param, 
					  const BkgdEvoMdlConfig *conf) {
  double coef;

  /* calculate partial derivs of branch mutation probs (p)
   * using branch len partials
   */

  /* initialize prob deriv structure to 0 */
  param_set_zero(&br->dprob_i);


  /* calculate partial derivs w.r.t. pop size, speciation times and
   * deleterious rates.  These all depend on branch len partial
   * derivs, and all have the same coefficient.
   */
  coef = br->k_i*(param->mu_i + param->mu_v) - br->k_v*param->mu_v;

  br->dprob_i.N_hc    = coef * br->dlen.N_hc;
  br->dprob_i.N_hcg   = coef * br->dlen.N_hcg;
  br->dprob_i.N_hcgo  = coef * br->dlen.N_hcgo;
  br->dprob_i.N_hcgom = coef * br->dlen.N_hcgom;
  br->dprob_i.T_hc    = coef * br->dlen.T_hc;
  br->dprob_i.T_hcg   = coef * br->dlen.T_hcg;
  br->dprob_i.T_hcgo  = coef * br->dlen.T_hcgo;
  br->dprob_i.T_hcgom = coef * br->dlen.T_hcgom;

  /* HCM model params */
  br->dprob_i.N_hcm   = coef * br->dlen.N_hcm;
  br->dprob_i.T_hcm   = coef * br->dlen.T_hcm;

  /* HCOM model params */
  br->dprob_i.N_hco   = coef * br->dlen.N_hco;
  br->dprob_i.T_hco   = coef * br->dlen.T_hco;
  br->dprob_i.N_hcom  = coef * br->dlen.N_hcom;
  br->dprob_i.T_hcom  = coef * br->dlen.T_hcom;

  br->dprob_i.u_ex_scale  = coef * br->dlen.u_ex_scale;
  br->dprob_i.u_nex_scale = coef * br->dlen.u_nex_scale;

  /* derivs w.r.t. mu param do not depend on branch len derivs */
  if(conf->mu_type == MU_TYPE_SINGLE) {
    br->dprob_i.mu_i = br->k_i * br->len;
    br->dprob_i.mu_v = br->k_i*br->len  - br->k_v*br->len;
  }
  else if(conf->mu_type == MU_TYPE_CAT_LIN) {
    br->dprob_i.mu_a = br->k_i * br->len * bin->cat;
    br->dprob_i.mu_b = (br->k_i*br->len  - br->k_v*br->len) * bin->cat;
  }
  else if(conf->mu_type == MU_TYPE_CAT_QUAD) {
    g_error("branch_set_dprob_i_kimura: CAT_QUAD mutation model not yet "
	    "implemented with kimura substitution model");
  }
  else {
    g_error("branch_set_dprob_i_kimura: unknown mutation model type");
  }
}




/**
 * Sets partial derivatives of branch TRANSVERSION probabilities for
 * the kimura substitution model. Only sets dprob_v, attributes, does
 * not set dprob or dprob_i attributes.
 */
static void branch_set_dprob_v_kimura(Branch *br, const BkgdBin *bin, 
				      const BkgdEvoMdlParam *param, 
				      const BkgdEvoMdlConfig *conf) {
  double coef;

  /* calculate partial derivs of branch mutation probs (p)
   * using branch len partials
   */

  /* initialize prob deriv structure to 0 */
  param_set_zero(&br->dprob_v);


  /* calculate partial derivs w.r.t. pop size, speciation times and
   * deleterious rates.  These all depend on branch len partial
   * derivs, and all have the same coefficient.
   */
  coef = 2.0 * br->k_v * param->mu_v;

  br->dprob_v.N_hc    = coef * br->dlen.N_hc;
  br->dprob_v.N_hcg   = coef * br->dlen.N_hcg;
  br->dprob_v.N_hcgo  = coef * br->dlen.N_hcgo;
  br->dprob_v.N_hcgom = coef * br->dlen.N_hcgom;
  br->dprob_v.T_hc    = coef * br->dlen.T_hc;
  br->dprob_v.T_hcg   = coef * br->dlen.T_hcg;
  br->dprob_v.T_hcgo  = coef * br->dlen.T_hcgo;
  br->dprob_v.T_hcgom = coef * br->dlen.T_hcgom;

  /* HCM model params */
  br->dprob_v.N_hcm   = coef * br->dlen.N_hcm;
  br->dprob_v.T_hcm   = coef * br->dlen.T_hcm;

  /* HCOM model params */
  br->dprob_v.N_hco   = coef * br->dlen.N_hco;
  br->dprob_v.T_hco   = coef * br->dlen.T_hco;
  br->dprob_v.N_hcom  = coef * br->dlen.N_hcom;
  br->dprob_v.T_hcom  = coef * br->dlen.T_hcom;

  br->dprob_v.u_ex_scale  = coef * br->dlen.u_ex_scale;
  br->dprob_v.u_nex_scale = coef * br->dlen.u_nex_scale;

  /* derivs w.r.t. mu param do not depend on branch len derivs */
  if(conf->mu_type == MU_TYPE_SINGLE) {
    br->dprob_v.mu_i = 0.0;
    br->dprob_v.mu_v = 2.0 * br->k_v * br->len;
  }
  else if(conf->mu_type == MU_TYPE_CAT_LIN) {
    br->dprob_v.mu_a = 0.0;
    br->dprob_v.mu_b = 2.0 * br->k_v * br->len * bin->cat;
  }
  else if(conf->mu_type == MU_TYPE_CAT_QUAD) {
    g_error("branch_set_dprob_v_kimura: CAT_QUAD mutation model not yet "
	    "implemented with kimura substitution model");
  }
  else {
    g_error("branch_set_dprob_v_kimura: unknown mutation model type");
  }
}




/**
 * Calculates partial derivatives of branch substitution probability
 * w.r.t. each of the model parameters. Assumes that branch len
 * partial derivatives have already been set through a call
 * to br->set_dlen()
 */
void branch_set_dprob(Branch *br, const BkgdBin *bin, 
		      const BkgdEvoMdlParam *param, 
		      const BkgdEvoMdlConfig *conf) {

  if(conf->subst_mdl->id == SUBST_MDL_JUKES_CANTOR) {
    branch_set_dprob_jukes_cantor(br, bin, param, conf);
  }

  else if(conf->subst_mdl->id == SUBST_MDL_KIMURA) {
    /* calculate transition and transversion prob partial derivs */
    branch_set_dprob_i_kimura(br, bin, param, conf);
    branch_set_dprob_v_kimura(br, bin, param, conf);

    /* subst prob partials are sums of transition + transversion partials */
    param_set_zero(&br->dprob);
    param_add(&br->dprob_i, 1.0, &br->dprob);
    param_add(&br->dprob_v, 1.0, &br->dprob);
  }
  else {
    g_error("branch_set_dprob: unknown substitution model\n");
  }
}				





/********************
 * HC model functions
 ********************/

void branch_set_len__hc_h_c(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = 2.0*param->T_hc + 4.0*param->N_hc*param->B;

  if(isnan(br->len)) {
    g_error("branch_set_len__hc_h_c: invalid len (%g)\n"
	    "  T_hc=%g, N_hc=%g, B=%g\n",
	    br->len, param->T_hc, param->N_hc, param->B);
  }
}


void branch_set_dlen__hc_h_c(Branch *br, const BkgdBin *bin ,
			     const BkgdEvoMdlParam *param) {
  br->dlen.T_hc = 2.0;
  br->dlen.N_hc = 4.0 * param->B;
  br->dlen.u_ex_scale  = 4.0 * param->N_hc * bin->lB_ex * param->B;
  br->dlen.u_nex_scale = 4.0 * param->N_hc * bin->lB_nex * param->B;
}



/********************
 * HCM model functions
 ********************/

void branch_set_len__hcm_h_c(Branch *br, const BkgdEvoMdlParam *param) {
  /* same as for HC model */
  branch_set_len__hc_h_c(br, param);
}



void branch_set_dlen__hcm_h_c(Branch *br, const BkgdBin *bin ,
			      const BkgdEvoMdlParam *param) {
  /* same as HC model */
  branch_set_dlen__hc_h_c(br, bin, param);
}



void branch_set_len__hcm_m(Branch *br, const BkgdEvoMdlParam *param) {
  /* calculate branch length */
  br->len = 4.0*param->N_hcm*param->B + param->T_hcm + 
    1.4*(param->T_hcm + param->T_hc) - 2.0*param->N_hc*param->B;
}


void branch_set_dlen__hcm_m(Branch *br, const BkgdBin *bin,
			    const BkgdEvoMdlParam *param) {
  br->dlen.N_hc = -2.0 * param->B;
  br->dlen.N_hcm = 4.0 * param->B;
  br->dlen.T_hc = 1.4;
  br->dlen.T_hcm = 2.4;
  br->dlen.u_ex_scale = 
    (4.0*param->N_hcm - 2.0*param->N_hc) * bin->lB_ex * param->B;
  br->dlen.u_nex_scale = 
    (4.0*param->N_hcm - 2.0*param->N_hc) * bin->lB_nex * param->B;
}





/********************
 * HCOM model functions
 ********************/

void branch_set_len__hcom_h_c(Branch *br, const BkgdEvoMdlParam *param) {
  /* calculate branch length */
  br->len = 2.0*param->T_hc + 4.0*param->N_hc*param->B;
}



void branch_set_dlen__hcom_h_c(Branch *br, const BkgdBin *bin ,
			       const BkgdEvoMdlParam *param) {
  br->dlen.T_hc = 2.0;
  br->dlen.N_hc = 4.0 * param->B;
  br->dlen.u_ex_scale  = 4.0 * param->N_hc * bin->lB_ex * param->B;
  br->dlen.u_nex_scale = 4.0 * param->N_hc * bin->lB_nex * param->B;
}


/* branch between orang speciation and human/chimp speciation */

void branch_set_len__hcom_hc(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hco + 2.0 * param->B * (param->N_hco - param->N_hc);

  if(br->len <= 0.0) {
    g_error("branch_set_len__hcom_hc: invalid branch len: %g\n"
	    "  B=%g, N_hco=%g, N_hc=%g, T_hco=%g", br->len, param->B, 
	    param->N_hco, param->N_hc, param->T_hco);
  }

}

void branch_set_dlen__hcom_hc(Branch *br, const BkgdBin *bin, 
			      const BkgdEvoMdlParam *param) {
  br->dlen.N_hc = -2.0 * param->B;
  br->dlen.N_hco = 2.0 * param->B;
  br->dlen.T_hco = 1.0;

  br->dlen.u_ex_scale = (2.0 * bin->lB_ex * param->B) * 
    (param->N_hco - param->N_hc);

  br->dlen.u_nex_scale = (2.0 * bin->lB_nex * param->B) * 
    (param->N_hco - param->N_hc);
}


void branch_set_len__hcom_o(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hc + param->T_hco + 2.0 * param->N_hco * param->B;
}


void branch_set_dlen__hcom_o(Branch *br, const BkgdBin *bin, 
			     const BkgdEvoMdlParam *param) {
  br->dlen.T_hc = 1.0;
  br->dlen.T_hco = 1.0;
  br->dlen.N_hco = 2.0 * param->B;
  br->dlen.u_ex_scale = 2.0 * param->N_hco * bin->lB_ex * param->B;
  br->dlen.u_nex_scale = 2.0 * param->N_hco * bin->lB_nex * param->B;
}



void branch_set_len__hcom_m(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hcom + (4.0 * param->N_hcom * param->B) - 
    (2.0 * param->N_hco * param->B) + 
    1.4 * (param->T_hcom + param->T_hco + param->T_hc);
}


void branch_set_dlen__hcom_m(Branch *br, const BkgdBin *bin ,
			     const BkgdEvoMdlParam *param) {
  br->dlen.T_hc = 1.4;
  br->dlen.T_hco = 1.4;
  br->dlen.T_hcom = 2.4;
  br->dlen.N_hco = -2.0 * param->B;
  br->dlen.N_hcom = 4.0 * param->B;
  br->dlen.u_ex_scale = (2.0 * bin->lB_ex * param->B) * 
    (2.0*param->N_hco - param->N_hc);

  br->dlen.u_nex_scale = (2.0 * bin->lB_nex * param->B) * 
    (2.0*param->N_hco - param->N_hc);

}



/********************************************************
 * HCGOM model functions
 ********************************************************/


/* helper function: calculates partial derivatives 
 * of kappa_hcg (the probability human chimp coalescent
 * is older than gorilla speciation event)
 *
 * caches results, to avoid re-doing same calculations
 */
static void set_k_hcg_deriv(const BkgdBin *bin,
			    const BkgdEvoMdlParam *param, 
			    BkgdEvoMdlParam *k_deriv) {

  static double last_N_hc = -1.0;
  static double last_T_hcg = -1.0;
  static double last_u_ex_scale = -1.0;
  static double last_u_nex_scale = -1.0;
  static double last_B_ex = -1.0;
  static double last_B_nex = -1.0;

  static double dk_dN_hc = 0.0;
  static double dk_dT_hcg = 0.0;
  static double dk_du_ex_scale = 0.0;
  static double dk_du_nex_scale = 0.0;

  double coef;

  param_set_zero(k_deriv);

  if((last_N_hc != param->N_hc) || 
     (last_T_hcg != param->T_hcg) ||
     (last_u_ex_scale != param->u_ex_scale) ||
     (last_u_nex_scale != param->u_nex_scale) || 
     (last_B_ex != bin->B_ex) ||
     (last_B_nex != bin->B_nex)) {

    /* params from last request were different, recompute
     * kappa partial derivatives and cache
     */
    coef = (param->k_hcg)/(2.0 * param->N_hc * param->B);
    dk_dN_hc = (coef * param->T_hcg) / param->N_hc;
    dk_dT_hcg = -coef;
    dk_du_ex_scale = coef * param->T_hcg * bin->lB_ex;
    dk_du_nex_scale = coef * param->T_hcg * bin->lB_nex;
      
    last_N_hc = param->N_hc;
    last_T_hcg = param->T_hcg;
    last_u_ex_scale = param->u_ex_scale;
    last_u_nex_scale = param->u_nex_scale;

    last_B_ex = bin->B_ex;
    last_B_nex = bin->B_nex;
  } 

  /* use cached values */
  k_deriv->N_hc = dk_dN_hc;
  k_deriv->T_hcg = dk_dT_hcg;
  k_deriv->u_ex_scale = dk_du_ex_scale;
  k_deriv->u_nex_scale = dk_du_nex_scale; 
}


/* H branch */
void branch_set_len__hcgom_h(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hc + (2.0 * param->B * param->N_hc) + 
    param->k_hcg * ((4.0/3.0)*param->B*param->N_hcg - 
		    2.0 * param->B*param->N_hc);
}

void branch_set_dlen__hcgom_h(Branch *br, const BkgdBin *bin,
			      const BkgdEvoMdlParam *param) {
  BkgdEvoMdlParam k_hcg_deriv;
  double c;

  set_k_hcg_deriv(bin, param, &k_hcg_deriv);

  c = 2.0 * param->B;

  br->dlen.N_hc = c * (1.0 - param->k_hcg + 
		    k_hcg_deriv.N_hc*((2.0/3.0)*param->N_hcg - param->N_hc));

  br->dlen.N_hcg = c * (2.0 / 3.0) * param->k_hcg;

  br->dlen.T_hc = 1.0;

  br->dlen.T_hcg = c * 
     k_hcg_deriv.T_hcg * (2.0/3.0*param->N_hcg - param->N_hc);
  
  br->dlen.u_ex_scale = 
    c * (bin->lB_ex * param->N_hc + 
	 (k_hcg_deriv.u_ex_scale + bin->lB_ex*param->k_hcg) *
	 ((2.0/3.0) * param->N_hcg - param->N_hc));
  
  br->dlen.u_nex_scale = 
    c * (bin->lB_nex * param->N_hc + 
	 (k_hcg_deriv.u_nex_scale + bin->lB_nex*param->k_hcg) *
	 ((2.0/3.0) * param->N_hcg - param->N_hc));
}




/* C branch */
void branch_set_len__hcgom_c(Branch *br, const BkgdEvoMdlParam *param) {
  /* same as H branch */
  branch_set_len__hcgom_h(br, param);
}

void branch_set_dlen__hcgom_c(Branch *br, const BkgdBin *bin, 
			      const BkgdEvoMdlParam *param) {
  /* same as H branch */
  branch_set_dlen__hcgom_h(br, bin, param);
}



/* G branch */
void branch_set_len__hcgom_g(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hc + param->T_hcg + 
    2.0*(1.0 - param->k_hcg/3.0) * param->B * param->N_hcg;
}

void branch_set_dlen__hcgom_g(Branch *br, const BkgdBin *bin,
			      const BkgdEvoMdlParam *param) {
  BkgdEvoMdlParam k_hcg_deriv;
  double c;

  set_k_hcg_deriv(bin, param, &k_hcg_deriv);

  c = 2.0*param->B;

  br->dlen.N_hc = -(c/3.0) * param->N_hcg * k_hcg_deriv.N_hc;

  br->dlen.N_hcg = c * (1.0 - param->k_hcg/3.0);
  
  br->dlen.T_hc = 1.0;

  br->dlen.T_hcg = 1.0 - (c/3.0)*param->N_hcg * k_hcg_deriv.T_hcg;


  br->dlen.u_ex_scale = 
    c * param->N_hcg * (bin->lB_ex * (1.0 - param->k_hcg / 3.0) -
			(1.0/3.0) * k_hcg_deriv.u_ex_scale);

  br->dlen.u_nex_scale = 
    c * param->N_hcg * (bin->lB_nex * (1.0 - param->k_hcg / 3.0) -
			(1.0/3.0) * k_hcg_deriv.u_nex_scale);
}



/* HC branch */
void branch_set_len__hcgom_hc(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hcg + 
    2.0*param->B*((param->N_hcg - param->N_hc)*(1.0-param->k_hcg) +
		  (param->k_hcg * param->N_hcg)/3.0);
}


/* HC branch partial derivatives */
void branch_set_dlen__hcgom_hc(Branch *br, const BkgdBin *bin,
			    const BkgdEvoMdlParam *param) {
  BkgdEvoMdlParam k_deriv;
  double c;

  set_k_hcg_deriv(bin, param, &k_deriv);

  c = 2.0 * param->B;

  br->dlen.N_hc = 
    c * (param->k_hcg - 1.0 - 
	 k_deriv.N_hc * ((2.0/3.0)*param->N_hcg - param->N_hc));

  br->dlen.N_hcg = c * (1.0 - (2.0/3.0)*param->k_hcg);

  br->dlen.T_hcg = 1.0 - 
    c * k_deriv.T_hcg * ((2.0/3.0)*param->N_hcg - param->N_hc);

  br->dlen.u_ex_scale = 
    c * bin->lB_ex * (param->N_hcg - param->N_hc - 
		      param->k_hcg * ((2.0/3.0)*param->N_hcg - param->N_hc)) +
    c * k_deriv.u_ex_scale * (param->N_hc - (2.0/3.0)*param->N_hcg);

  br->dlen.u_nex_scale = 
    c * bin->lB_nex * (param->N_hcg - param->N_hc - 
		       param->k_hcg * ((2.0/3.0)*param->N_hcg - param->N_hc)) +
    c * k_deriv.u_nex_scale * (param->N_hc - (2.0/3.0)*param->N_hcg);
}



/* HG+CG branches */
void branch_set_len__hcgom_hg_cg(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->k_hcg * (4.0/3.0)*param->B*param->N_hcg;
}


/* HG branch */
void branch_set_len__hcgom_hg(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->k_hcg * (2.0/3.0)*param->B*param->N_hcg;
}

/* CG branch */
void branch_set_len__hcgom_cg(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->k_hcg * (2.0/3.0)*param->B*param->N_hcg;
}



/* combined CG+HG branch partial derivatives */
void branch_set_dlen__hcgom_hg_cg(Branch *br, const BkgdBin *bin,
				  const BkgdEvoMdlParam *param) {
  BkgdEvoMdlParam k_hcg_deriv;
  double c;

  set_k_hcg_deriv(bin, param, &k_hcg_deriv);

  /* double HG deriv to count both HG and CG branches */
  c = (4.0 / 3.0) * param->B;

  br->dlen.N_hc = c * param->N_hcg * k_hcg_deriv.N_hc;

  br->dlen.N_hcg = c * param->k_hcg;

  br->dlen.T_hcg = c * k_hcg_deriv.T_hcg * param->N_hcg;

  br->dlen.u_ex_scale = 
    c * (bin->lB_ex * param->k_hcg + k_hcg_deriv.u_ex_scale) * param->N_hcg;

  br->dlen.u_nex_scale = 
    c * (bin->lB_nex * param->k_hcg + k_hcg_deriv.u_nex_scale) * param->N_hcg;

}


/* combined HG branch partial derivatives */
void branch_set_dlen__hcgom_hg(Branch *br, const BkgdBin *bin,
				  const BkgdEvoMdlParam *param) {
  BkgdEvoMdlParam k_hcg_deriv;
  double c;

  set_k_hcg_deriv(bin, param, &k_hcg_deriv);

  /* double HG deriv to count both HG and CG branches */
  c = (2.0 / 3.0) * param->B;

  br->dlen.N_hc = c * param->N_hcg * k_hcg_deriv.N_hc;

  br->dlen.N_hcg = c * param->k_hcg;

  br->dlen.T_hcg = c * k_hcg_deriv.T_hcg * param->N_hcg;

  br->dlen.u_ex_scale = 
    c * (bin->lB_ex * param->k_hcg + k_hcg_deriv.u_ex_scale) * param->N_hcg;

  br->dlen.u_nex_scale = 
    c * (bin->lB_nex * param->k_hcg + k_hcg_deriv.u_nex_scale) * param->N_hcg;

}


/* combined HG branch partial derivatives */
void branch_set_dlen__hcgom_cg(Branch *br, const BkgdBin *bin,
			       const BkgdEvoMdlParam *param) {
  /* same as HG branch */
  branch_set_dlen__hcgom_hg(br, bin, param);
}


/* HCG branch */
void branch_set_len__hcgom_hcg(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hcgo + 
    2.0*param->B * (param->N_hcgo - param->N_hcg*(1.0 + param->k_hcg/3.0));
}




/* HCG branch partial derivatives */
void branch_set_dlen__hcgom_hcg(Branch *br, const BkgdBin *bin,
			     const BkgdEvoMdlParam *param) {
  BkgdEvoMdlParam dk;
  double c;

  set_k_hcg_deriv(bin, param, &dk);

  c = 2.0 * param->B;

  br->dlen.N_hc = -(c/3.0) *  param->N_hcg * dk.N_hc;

  br->dlen.N_hcg = -c * (1.0 + param->k_hcg / 3.0);

  br->dlen.N_hcgo = c;

  br->dlen.T_hcg = -(c/3.0) * param->N_hcg * dk.T_hcg;

  br->dlen.T_hcgo = 1.0;

  
  br->dlen.u_ex_scale = c * (bin->lB_ex * param->N_hcgo -
			     bin->lB_ex * param->N_hcg -
			     (bin->lB_ex * param->N_hcg * param->k_hcg +
			      dk.u_ex_scale * param->N_hcg)/3.0);\

  br->dlen.u_nex_scale = c * (bin->lB_nex * param->N_hcgo -
			     bin->lB_nex * param->N_hcg -
			     (bin->lB_nex * param->N_hcg * param->k_hcg +
			      dk.u_nex_scale * param->N_hcg)/3.0);


}



/* 0 branch */
void branch_set_len__hcgom_o(Branch *br, const BkgdEvoMdlParam *param) {
  br->len = param->T_hc + param->T_hcg + param->T_hcgo + 
    2.0*param->B*param->N_hcgo;
}



/* O branch partial derivatives */
void branch_set_dlen__hcgom_o(Branch *br, const BkgdBin *bin,
			   const BkgdEvoMdlParam *param) {
  br->dlen.N_hcgo = 2.0 * param->B;
  br->dlen.T_hc = 1.0;
  br->dlen.T_hcg = 1.0;
  br->dlen.T_hcgo = 1.0;

  br->dlen.u_ex_scale = 2.0 * bin->lB_ex * param->B * param->N_hcgo;
  br->dlen.u_nex_scale = 2.0 * bin->lB_nex * param->B * param->N_hcgo;
}


/* M branch */
void branch_set_len__hcgom_m(Branch *br, const BkgdEvoMdlParam *param) {
  /* 1.4 is correction for accelerated mutation rate in old world monkeys */
  br->len = 1.4 * (param->T_hc+param->T_hcg+param->T_hcgo+param->T_hcgom) + 
    param->T_hcgom + 4.0*param->B*param->N_hcgom - 
    2.0 * param->B * param->N_hcgo;
}



/* M branch partial derivatives */
void branch_set_dlen__hcgom_m(Branch *br, const BkgdBin *bin,
			      const BkgdEvoMdlParam *param) {
  br->dlen.N_hcgo = -2.0*param->B;
  br->dlen.N_hcgom = 4.0 * param->B;
  br->dlen.T_hc = 1.4;
  br->dlen.T_hcg = 1.4;
  br->dlen.T_hcgo = 1.4;
  br->dlen.T_hcgom = 2.4;

  br->dlen.u_ex_scale = 
    2.0 * param->B * bin->lB_ex * (2.0*param->N_hcgom - param->N_hcgo); 

  br->dlen.u_nex_scale = 
    2.0 * param->B * bin->lB_nex * (2.0*param->N_hcgom - param->N_hcgo);
}

