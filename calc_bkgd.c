#include <math.h>
#include <string.h>

#include <glib.h>

#include "bio/seqcoord.h"
#include "bio/chr.h"
#include "util/config.h"
#include "util/util.h"
#include "rectab.h"
#include "bkgd.h"
#include "bkgd_param.h"
#include "bkgd_interp.h"


/**
 * Retrieves a table of recombination rates for current chromosome
 * from database
 */
static RecRateTable *get_rectab(Config *config, Chromosome *chr) {
  long n_sf, i, start;
  SeqFeature *sf;
  char *rec_rate_dir, *rec_rate_prefix, *rec_rate_postfix;
  char *rec_rate_filename;
  double rec_scale, rate;
  RecRateTable *rtab;
  FILE *f;
  char buf[1024];

  rec_rate_dir = config_get_str(config, "RECOMB_RATE_DIR");  
  rec_rate_prefix = config_get_str(config, "RECOMB_RATE_PREFIX");
  rec_rate_postfix = config_get_str(config, "RECOMB_RATE_POSTFIX");
  rec_scale = config_get_double(config, "RECOMB_RATE_SCALE");
  
  rec_rate_filename = util_str_concat(rec_rate_dir, "/", rec_rate_prefix,
				      chr->name, rec_rate_postfix, NULL);

  f = util_must_fopen(rec_rate_filename, "r");
  /* Subtract two from number of lines in file to determine number of 
   * features. First line is header, and there are 1 less recombination rate
   * interval features than there are markers in the map.
   */
  n_sf = util_fcount_lines(f) - 2;
  
  /* skip header line */
  if(fgets(buf, sizeof(buf), f) == NULL) {
    g_error("%s:%d: failed to read header line from file %s", 
	    __FILE__, __LINE__, rec_rate_filename);
  }

  sf = g_new(SeqFeature, n_sf);

  /* create a seq feature for each recombination rate block
   * set seq feature score equal to the recombination rate
   * Note that file actually gives markers, not blocks.
   */
  fprintf(stderr, "reading recombination rates from file '%s'\n", 
	  rec_rate_filename);
  for(i = 0; i < n_sf+1; i++) {
    if(fgets(buf, sizeof(buf), f) == NULL) {
      g_error("%s:%d: failed to read line %ld from file %s", 
	      __FILE__, __LINE__, i+1, rec_rate_filename);
    }

    /* read rate and start of this interval */
    sscanf(buf, "%*s %ld %lf %*s", &start, &rate);

    if(i < n_sf) {
      /* this is not the final line, so start a new feature */
      sf[i].c.chr = NULL;
      sf[i].c.seqname = util_str_dup(chr->name);

      sf[i].n_sub_feat = 0;
      sf[i].sub_feats = NULL;
      sf[i].name = NULL;
      sf[i].attrib = NULL;
      
      sf[i].c.start = start;
      sf[i].score = rate;

      if(sf[i].score < 0.0) {
	g_error("%s:%d: invalid recombination rate %g for position %ld\n",
		__FILE__, __LINE__, sf[i].score, sf[i].c.start);
		
      }
    }

    if(i > 0) {
      /* set end of previous feature */
      sf[i-1].c.end = start - 1;
    }
  }
  fprintf(stderr, "read %ld recombination regions\n", n_sf);

  /* convert these features to a 'recombination rate table' */
  rtab = rectab_from_feats(sf, n_sf, chr->len, rec_scale);

  /* free the features that we created */
  seqfeat_array_free(sf, n_sf);
  g_free(rec_rate_filename);

  return rtab;
}




static char *get_cons_filename(Config *config, Chromosome *chr) {
  char *filename;
  char *dir;
  char *postfix;

  dir = config_get_str(config, "CONS_DIR");
  postfix = config_get_str(config, "CONS_POSTFIX");
  
  filename = util_str_concat(dir, "/", chr->name, postfix, NULL);

  return filename;
}


/**
 * Converts provided recombination rates to recombination distances
 * along chromosome and retrieves a list of rec-dist coords
 * representing conserved blocks. Conserved blocks are split when
 * recombination rates change in order to make integral approximation
 * of background selection sum more efficient.
 */
static GList *get_cons_rec_dists(Config *config, Chromosome *chr,
				 RecRateTable *rtab) {
  long i, j, n_sf, ttl_left, ttl_right, ttl_cons, last_end;
  SeqFeature *sf;
  ConsBlock *cblk;
  GList *cons_list, *cur;
  char *cons_filename;

  /* read list of conserved features for this chromosome from a BED file */
  cons_filename = get_cons_filename(config, chr);
  fprintf(stderr, "reading conserved features from '%s'\n", cons_filename);
  sf = seqfeat_read_bed(cons_filename, &n_sf);
  g_free(cons_filename);

  /* order conserved elements */
  qsort(sf, n_sf, sizeof(SeqFeature), seqfeat_cmp_nostrand);

  /* Create list of conserved block coords splitting elements when
   * they span a change in recombination rate. The coordinates are
   * initially set to physical positions, and then updated to be
   * recombination distances.
   */
  cons_list = NULL;
  ttl_cons = 0;
  last_end = -1;
  for(i = 0; i < n_sf; i++) {
    cblk = g_new(ConsBlock, 1);
    cblk->start = sf[i].c.start;

    /* sanity check: we don't want overlapping blocks */
    if(cblk->start <= last_end) {
      g_error("get_cons_rec_dists: conserved blocks should not overlap: "
	      "cur_start=%ld, last_end=%ld", cblk->start, last_end);
    }
    
    ttl_cons += sf[i].c.end - sf[i].c.start + 1;

    cblk->r = rectab_rate(rtab, sf[i].c.start);
    for(j = sf[i].c.start; j < sf[i].c.end; j++) {
      if(rectab_rate(rtab,j+1) != cblk->r) {
	/* change in rate, end previous block */
	cblk->end = j;
	last_end = cblk->end;

	cons_list = g_list_append(cons_list, cblk);
	
/* 	fprintf(stderr, "splitting cons segment %ld-%ld into " */
/* 		"%ld-%ld and %ld-%ld (rec-rate change from %g to %g)\n", */
/* 		cblk->start, sf[i].c.end, */
/* 		cblk->start, cblk->end, */
/* 		j+1, sf[i].c.end, cblk->r, rectab_rate(rtab,j+1)); */

	/* start new block */
	cblk = g_new(ConsBlock, 1);
	cblk->start = j+1;
	cblk->r = rectab_rate(rtab, j+1);
      }
    }
    /* end current block */
    cblk->end = sf[i].c.end;
    last_end = cblk->end;

    cons_list = g_list_append(cons_list, cblk);
  }

  seqfeat_array_free(sf, n_sf);
  
  /* convert cons block phys coords to rec-dist coords */
  ttl_left = 0;
  ttl_right = ttl_cons;

  cur = cons_list;
  while(cur != NULL) {
    cblk = cur->data;
    cblk->r_start = rectab_rpos(rtab, cblk->start);
    cblk->r_end   = rectab_rpos(rtab, cblk->end);

    /* keep track of number of cons sites to left and right of each
     * blk (not including the block itself)
     */
    cblk->left_ttl = ttl_left;
    ttl_left  += cblk->end - cblk->start + 1;

    ttl_right -= cblk->end - cblk->start + 1;
    cblk->right_ttl = ttl_right;

    cur = g_list_next(cur);
  }

  return cons_list;
}



/* 
 * Calculates background selection strength at each position along a
 * chromosome using provided parameters, recombination map and list of
 * conserved elements.
 */
void calc_bkgd_chr(Chromosome *chr, RecRateTable *rtab, GList *cons_list,
		   BkgdParam *parm,  FILE *out_fh) {
  double b;
  long pos;
  GList *next_cons;
  int b_int, prev_b_int, b_len;
  BkgdInterp *bgi;

  next_cons = cons_list;

  /* b is background selection strength */
  prev_b_int = b_int = -1;

  b_len = 0;

  /* Create interpolator to estimate B values at positions along chr */
  bgi = bkgd_interp_new(rtab, chr->len, cons_list, parm);

  pos = 1;
  while(pos <= chr->len) {
    b = bkgd_interp_eval(bgi, pos);
    
    if((pos % 1000000)==0) {
      fprintf(stderr, ".");
    }

    /* fprintf(stderr, "pos=%ld, b=%g\n", pos, b); */

    /* truncate to 3 digit integer */
    b_int = (int)round(b* BKGD_SCALE);

    /* only print out value if rounded value is different from prev one */
    if(prev_b_int != b_int) {
      if(prev_b_int >= 0) {
	fprintf(out_fh, "%d %d\n", prev_b_int, b_len);
	/* fprintf(stderr, "%ld %d %d\n", pos, prev_b_int, b_len); */
      }
      prev_b_int = b_int;
      b_len = 0;
    }

    pos++;
    b_len++;
  }

  /* print out final value */
  if(b_len > 0) {
    fprintf(out_fh, "%d %d\n", b_int, b_len);
    /* fprintf(stderr, "%ld %d %d\n", pos, prev_b_int, b_len); */
  }

  fprintf(stderr, "\n");

  bkgd_interp_free(bgi);

}






int main(const int argc, const char **argv) {
  Config *config;
  GList *cons_list, *cur;
  char *out_dir, *out_path, *chr_filename;
  int n_chr, i;
  RecRateTable *rtab;
  FILE *out_fh;
  BkgdParam *parm;
  Chromosome *chrs;
  
  if(argc < 2) {
    fprintf(stderr, "usage: %s <config_file> [additional_config_args]\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");

  config = config_read_args(argc, argv);

  out_dir = config_get_str(config, "OUTPUT_DIR");
  parm = bkgd_param_new(config);

  /* get list of chromosome names to run on */
  chr_filename = config_get_str(config, "CHROMOSOME_FILE");
  chrs = chr_read_file(chr_filename, &n_chr);

  for(i = 0; i < n_chr; i++) {
    fprintf(stderr, "\nchromosome: %s\n", chrs[i].name);
    out_path = util_str_concat(out_dir, "/", chrs[i].name, ".bkgd", NULL);
    out_fh = util_must_fopen(out_path, "w");
    g_free(out_path);

    fprintf(stderr, "retrieving recombination rates\n");
    rtab = get_rectab(config, &chrs[i]);
    
    /* get recombination distances of conserved sites,
     * and convert rec rates to dists
     */
    fprintf(stderr, "calculating conserved rec dists\n");
    cons_list = get_cons_rec_dists(config, &chrs[i], rtab);

    fprintf(stderr, "total recomb dist for %s: %gM\n", 
	    chrs[i].name, rtab->chr_r_len);

    /* calculate strength of background selection at each position on chr */
    fprintf(stderr, "calculating b vals\n");
    calc_bkgd_chr(&chrs[i], rtab, cons_list, parm, out_fh);
    
    fprintf(stderr, "freeing recombination distances\n");
    rectab_free(rtab);
	    
    /* free conserved blocks */
    fprintf(stderr, "freeing conserved blocks\n");
    cur = cons_list;
    while(cur != NULL) {
      g_free(cur->data);
      cur = g_list_next(cur);
    }
    g_list_free(cons_list);

    fclose(out_fh);
  }

  fprintf(stderr, "freeing parameters and interpolation tables\n");
  bkgd_param_free(parm);

  fprintf(stderr, "freeing config\n");
  config_free(config);

  return 0;  
}

