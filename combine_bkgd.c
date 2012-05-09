
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



/*
 * Main program logic
 */
int main(int argc, char **argv) {  
  Config *config;
  Analysis *analysis;
  char *bkgd1_dir, *bkgd2_dir;
  SeqCoord *region;
  long len;
  int i, b1_int, b2_int, b_new, last_b, b_len;
  double uscale1, uscale2;
  BkgdReader *br1, *br2;
  FILE *out_fh;
        
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>\n", 
	    argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");

  config   = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);
  analysis = analysis_new(config);

  bkgd1_dir = config_get_str(config, "BKGD1_DIR");
  bkgd2_dir = config_get_str(config, "BKGD2_DIR");

  /* get deleterious mutation rate rescale */
  uscale1 = config_get_double(config, "BKGD1_USCALE");
  uscale2 = config_get_double(config, "BKGD2_USCALE");

  fprintf(stderr, "using uscale1=%g, uscale2=%g\n", uscale1, uscale2);

  while(analysis_next_region(analysis)) {
    region = analysis_get_region(analysis);

    fprintf(stderr, "chromosome %s\n", region->chr->name);

    len = region->end - region->start + 1;

    out_fh = get_output_fh(analysis);

    br1 = get_bkgd_reader(analysis, bkgd1_dir);
    br2 = get_bkgd_reader(analysis, bkgd2_dir);

    last_b = B_UNDEF;
    b_len = 0;
    
    for(i = 0; i < len; i++) {
      b1_int = bkgd_reader_get_pos(br1, i+1);
      b2_int = bkgd_reader_get_pos(br2, i+1);

      b_new = combine_bkgd(b1_int, uscale1, b2_int, uscale2);

      if(b_new == last_b) {
	b_len++;
      } else {
	if(b_len > 0) {
	  fprintf(out_fh, "%d %d\n", last_b, b_len);
	}
	b_len = 1;
	last_b = b_new;
      }
    }

    /* write final region */
    if(b_len > 0) {
      fprintf(out_fh, "%d %d\n", last_b, b_len);
    }

    bkgd_reader_free(br1);
    bkgd_reader_free(br2);
  }

  fprintf(stderr, "Freeing analysis\n");
  analysis_free(analysis);

  fprintf(stderr, "Freeing config\n");
  config_free(config);
  
  return 0;
}
