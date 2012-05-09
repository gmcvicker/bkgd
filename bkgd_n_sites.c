
#include <stdio.h>
#include <stdlib.h>

#include <util/config.h>

#include "bkgd_data.h"


void write_data(FILE *fh, BkgdEvoMdlData *data) {
  long i, n;
  BkgdBin *bin;

  bin = data->bin;

  /* write header line */
  fprintf(fh, "B.EX\tB.NEX\tCAT\tN\t"
	  "H\tC\tG\tO\tM\t"
	  "HC\tHG\tCG\tHO\tCO\tHCG\t"
	  "CONS\n");
    
  for(i =0; i < data->n_bin; i++) {
    n = bin[i].h + bin[i].c + bin[i].g + bin[i].o + bin[i].m +
      bin[i].hc + bin[i].hg + bin[i].cg + bin[i].ho + bin[i].co +
      bin[i].hcg + bin[i].cons;

    fprintf(fh, "%g\t%g\t%g\t%ld\t"
	    "%ld\t%ld\t%ld\t%ld\t%ld\t"
	    "%ld\t%ld\t%ld\t%ld\t%ld\t"
	    "%ld\t%ld\n",
	    bin[i].B_ex, bin[i].B_nex, bin[i].cat, n,
	    bin[i].h, bin[i].c, bin[i].g, bin[i].o, bin[i].m, 
	    bin[i].hc, bin[i].hg, bin[i].cg, bin[i].ho, bin[i].co, 
	    bin[i].hcg, bin[i].cons);
  }

}



int main(int argc, char **argv) {
  BkgdEvoMdlData *data;
  Config *config;

  if(argc != 2) {
    fprintf(stderr, "usage: %s <config_file>", argv[0]);
    exit(2);
  }

  fprintf(stderr, "Reading config\n");
  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  data = bkgd_data_read_data(config);

  if(config_get_boolean(config, "COMBINE_CAT_BINS")) {
    double min_cat, max_cat;
    min_cat = config_get_double(config, "MIN_CAT");
    max_cat = config_get_double(config, "MAX_CAT");

    bkgd_data_combine_cat_bin(data, min_cat, max_cat);
  } 

  if(config_get_boolean(config, "COMBINE_DATA_BINS")) {
    long n_bin = config_get_long(config, "N_DATA_BIN");
    bkgd_data_combine_bin(data, n_bin);
  }

  if(config_get_boolean(config, "WRITE_COUNTS")) {
    bkgd_data_write_site_counts(stdout, data);
  }

  if(config_get_boolean(config, "WRITE_BINS")) {
     write_data(stdout, data);
  }

  g_free(data->bin);
  g_free(data);
  
  return 0;
}
