
#include <glib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <bio/bkgd_reader.h>
#include <analysis.h>
#include <bio/seqmask.h>
#include <bio/seqcoord.h>
#include <util/config.h>

#define MAX_BKGD 1000
#define BKGD_TYPE_SINGLE 1
#define BKGD_TYPE_COMBINED 2


static BkgdReader *get_bkgd_reader(Analysis *an, char *dir) {
  char  *path;
  SeqCoord *region;
  BkgdReader *br;
  char *postfix;

  postfix = config_get_str(an->config, "BKGD_POSTFIX");

  region = analysis_get_region(an);
  path = g_strconcat(dir, region->chr->name, postfix, NULL);

  fprintf(stderr, "reading bkgd data from '%s'\n", path);

  br = bkgd_reader_new(path);

  g_free(path);

  return br;
}



static int combine_bkgd(int b1, double uscale1, int b2, double uscale2) {
  static int last_b1 = -1;
  static int last_b2 = -1;
  static int last_combined = -1;
  static double last_uscale1 = -1.0;
  static double last_uscale2 = -1.0;
  int b_combined_int;
  double b1_dbl, b2_dbl, b_combined;
  
  if((b1 == last_b1) && (b2 == last_b2) && 
     (last_uscale1 == uscale1) && (last_uscale2 == uscale2)) {
    /* don't recompute if values were the same as last time */
    return last_combined;
  }

  if((b1 > 1000) || (b1 < 0) || (b2 > 1000) || (b2 < 0)) {
    g_error("combine_bkgd: expected B integers (%d,%d)to be in "
	    "range 0-1000",b1,b2);
  }

  /* convert integers in range 0-1000 to floating points between 0-1.0 */
  b1_dbl = (double)b1 / 1000.0;
  b2_dbl = (double)b2 / 1000.0;
     
  if(uscale1 != 1.0 || uscale2 != 1.0) {
    b_combined = exp(log(b1_dbl)*uscale1 + log(b2_dbl)*uscale2);
  } else {
    b_combined = b1_dbl*b2_dbl;
  }

  /* convert back to integer in range 0-1000 */
  b_combined_int = round(b_combined*1000.0);
  if(b_combined_int < 0 || b_combined_int > 1000) {
    g_error("combine_bkgd: expected combined B integer (%d) to be in "
	    "range 0-1000", b_combined_int);
  }

  /* Remember the last parameters that were passed. Because we often
   * expect the same parameters to be passed multiple times in a row,
   * we can save repeating expensive computations.
   */
  last_b1 = b1;
  last_b2 = b2;
  last_uscale1 = uscale1;
  last_uscale2 = uscale2;
  last_combined = b_combined_int;
  
  return b_combined_int;
}




static void bkgd_hist(Analysis *an, long *hist) {
  char *bkgd_dir;
  int i, bkgd_int;
  long pos;
  BkgdReader *br;
  long seq_len;

  SeqCoord *region;

  bkgd_dir = config_get_str(an->config, "BKGD_DIR");
  region = analysis_get_region(an);

  seq_len = region->end - region->start + 1;

  fprintf(stderr, "%s, len=%ld\n", region->chr->name, seq_len);

  br = get_bkgd_reader(an, bkgd_dir);

  for(pos = region->start; pos <= region->end; pos++) {

    i = pos - region->start;

    /* skip masked sites */
    if(an->mask->mask[i] & an->mask_id) {
      continue;
    }
      
    /* get B value at position and add to histogram */
    bkgd_int = bkgd_reader_get_pos(br, pos);
    
    if(bkgd_int > MAX_BKGD) {
      g_error("%s:%d: B value is too large: %d", __FILE__, __LINE__, bkgd_int);
    }
    if(bkgd_int < 0) {
      continue;
    }
      
    hist[bkgd_int] += 1;
  }
  
  bkgd_reader_free(br);
}




static void bkgd_hist_combined(Analysis *an, long *hist) {
  char *bkgd1_dir, *bkgd2_dir;
  double uscale1, uscale2;
  int b_new, b1_int, b2_int;
  long len, i;
  BkgdReader *br1, *br2;
  SeqCoord *region;

  bkgd1_dir = config_get_str(an->config, "BKGD1_DIR");
  bkgd2_dir = config_get_str(an->config, "BKGD2_DIR");

  /* get deleterious mutation rate rescale */
  uscale1 = config_get_double(an->config, "BKGD1_USCALE");
  uscale2 = config_get_double(an->config, "BKGD2_USCALE");

  region = analysis_get_region(an);

  fprintf(stderr, "chromosome %s\n", region->chr->name);

  len = region->end - region->start + 1;

  br1 = get_bkgd_reader(an, bkgd1_dir);
  br2 = get_bkgd_reader(an, bkgd2_dir);
    
  for(i = 0; i < len; i++) {
    /* skip masked sites */
    if(an->mask->mask[i] & an->mask_id) {
      continue;
    }
      
    /* get B value at position and add to histogram */
    b1_int = bkgd_reader_get_pos(br1, i+1);
    b2_int = bkgd_reader_get_pos(br2, i+1);

    if((b1_int < 0) || (b2_int < 0)) {
      continue;
    }

    b_new = combine_bkgd(b1_int, uscale1, b2_int, uscale2);    
    if(b_new > MAX_BKGD) {
      g_error("%s:%d: B value is too large: %d", __FILE__, __LINE__, b_new);
    }
      
    hist[b_new] += 1;
  }

  bkgd_reader_free(br1);
  bkgd_reader_free(br2);

}



static int get_bkgd_type(Config *config) {
  char *str;

  str = config_get_str(config, "BKGD_TYPE");
  if(strcmp(str, "COMBINED")==0) {
    return BKGD_TYPE_COMBINED;
  }
  if(strcmp(str, "SINGLE")==0) {
    return BKGD_TYPE_SINGLE;
  }

  g_error("%s:%d: BKGD_TYPE should be 'COMBINED' or 'SINGLE' not '%s'",
	  __FILE__, __LINE__, str);

  return 0;
}


/*
 * Main program logic
 */
int main(const int argc, const char **argv) {  
  Config *config;
  Analysis *analysis;
  long *hist;
  int bkgd_type,i;

                
  if(argc < 2) {
    fprintf(stderr, "usage: %s <config_file> "
	    "[<CONF_KEY1>=<VAL1> [<CONF_KEY2>=<VAL2>] [...]]\n", argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");

  config   = config_read_args(argc, argv, CONFIG_MISSING_KEY_ERROR);
  analysis = analysis_new(config);

  bkgd_type = get_bkgd_type(config);
  hist = g_new0(long, MAX_BKGD+1);

  while(analysis_next_region(analysis)) {
    if(bkgd_type == BKGD_TYPE_SINGLE) {
      bkgd_hist(analysis, hist);
    } else {
      bkgd_hist_combined(analysis, hist);
    }
  }
  
  fprintf(analysis->output_fh[0], "B\tCOUNT\n");
  for(i = 0; i <= MAX_BKGD; i++) {
    fprintf(analysis->output_fh[0], "%d\t%lu\n", i, hist[i]);
  }

  g_free(hist);

  fprintf(stderr, "Freeing analysis\n");
  analysis_free(analysis);

  fprintf(stderr, "Freeing config\n");
  config_free(config);
  
  return 0;
}
