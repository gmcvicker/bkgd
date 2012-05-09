#include <glib.h>

#include <stdlib.h>

#include <util/config.h>

#include "bkgd_data.h"


/**
 * Main function
 */
int main(int argc, char **argv) {
  double mean_cat;
  BkgdEvoMdlData *data;
  Config *config;
    
  if(argc != 2) {
    fprintf(stderr, "usage: %s <config>\n", 
	    argv[0]);
    exit(2);
  }

  config = config_read_file(argv[1], CONFIG_MISSING_KEY_ERROR);

  data = bkgd_data_read_data(config);

  mean_cat = bkgd_data_mean_cat_val(data);

  fprintf(stderr, "mean category value: %g\n", mean_cat);

  free(data->bin);
  free(data);


  return 0;
}
