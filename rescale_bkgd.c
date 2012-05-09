
#include <glib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <util/config.h>

#include "../analysis/analysis.h"

#include "bkgd_reader.h"

#define B_UNDEF -1

BkgdReader *get_bkgd_reader(Analysis *an, char *dir) {
  char  *path;
  SeqCoord *region;
  BkgdReader *br;

  region = analysis_get_region(an);
  path = g_strconcat(dir, region->chr->name, ".bkgd.gz", NULL);

  br = bkgd_reader_new(path);

  g_free(path);

  return br;
}


FILE *get_output_fh(Analysis *an) {
  FILE *fh;
  SeqCoord *region;
  char *path, *dir;

  region = analysis_get_region(an);
  dir = config_get_str(an->config, "OUTPUT_DIR");

  path = g_strconcat(dir, region->chr->name, ".bkgd", NULL);

  fh = fopen(path, "w");
  if(fh == NULL) {
    g_error("get_output_fh: could not open file '%s'", path);
  }
  g_free(path);
  
  return fh;
}



static int rescale_bkgd(int b, double uscale) {
  static int last_b = -1;
  static int last_rescaled = -1;
  static double last_uscale = -1.0;
  int b_rescaled_int;
  double b_dbl, b_rescaled;
  
  if((b == last_b) && (uscale == last_uscale)) {
    /* don't recompute if values were the same as last time */
    return last_rescaled;
  }

  if((b > 1000) || (b < 0)) {
    g_error("combine_bkgd: expected B integer (%d)to be in "
	    "range 0-1000",b);
  }

  /* convert integer in range 0-1000 to floating points between 0-1.0 */
  b_dbl = (double)b / 1000.0;
     
  b_rescaled = exp(log(b_dbl)*uscale);

  /* convert back to integer in range 0-1000 */
  b_rescaled_int = round(b_rescaled*1000.0);
  if(b_rescaled_int < 0 || b_rescaled_int > 1000) {
    g_error("rescale_bkgd: expected combined B integer (%d) to be in "
	    "range 0-1000", b_rescaled_int);
  }

  /* Remember the last parameters that were passed. Because we often
   * expect the same parameters to be passed multiple times in a row,
   * we can save repeating expensive computations.
   */
  last_b = b;
  last_uscale = uscale;
  last_rescaled = b_rescaled_int;
  
  return b_rescaled_int;
}



/*
 * Main program logic
 */
int main(int argc, char **argv) {  
  Config *config;
  Analysis *analysis;
  char *bkgd_dir;
  SeqCoord *region;
  long len;
  int i, b_int, b_new, last_b, b_len;
  double uscale;
  BkgdReader *br;
  FILE *out_fh;
        
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");

  config   = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);
  analysis = analysis_new(config);

  bkgd_dir = config_get_str(config, "BKGD_DIR");

  /* get deleterious mutation rate rescale */
  uscale = config_get_double(config, "BKGD_USCALE");

  fprintf(stderr, "using uscale=%g\n", uscale);

  while(analysis_next_region(analysis)) {
    region = analysis_get_region(analysis);

    fprintf(stderr, "chromosome %s\n", region->chr->name);

    len = region->end - region->start + 1;

    out_fh = get_output_fh(analysis);

    br = get_bkgd_reader(analysis, bkgd_dir);

    last_b = B_UNDEF;
    b_len = 0;
    
    for(i = 0; i < len; i++) {
      b_int = bkgd_reader_get_pos(br, i+1);

      b_new = rescale_bkgd(b_int, uscale);

      if(b_new == last_b) {
	b_len++;
      } else {
	if(b_len > 0) {
	  fprintf(out_fh, "%d %d\n", b_new, b_len);
	}
	b_len = 1;
	last_b = b_new;
      }
    }

    bkgd_reader_free(br);
  }

  fprintf(stderr, "Freeing analysis\n");
  analysis_free(analysis);

  fprintf(stderr, "Freeing config\n");
  config_free(config);
  
  return 0;
}
