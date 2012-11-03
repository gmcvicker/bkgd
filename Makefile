CFLAGS = -c -Wall -O2 -std=c99
CFLAGS += `pkg-config glib-2.0 --cflags`
CFLAGS += `gsl-config --cflags`

LFLAGS=-lz -lm
LFLAGS += `pkg-config glib-2.0 --libs`
LFLAGS += `gsl-config --libs`
CC=gcc

objects=bkgd.o bkgd_interp.o interp_tab.o bkgd_intg.o bkgd_point.o bkgd_data.o util/util.o util/config.o

default: $(objects) calc_bkgd

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

calc_bkgd.o: calc_bkgd.c
	$(CC) -c $(CFLAGS) calc_bkgd.c -o calc_bkgd.o

calc_bkgd: $(objects) calc_bkgd.o
	$(CC) $(LFLAGS) -o calc_bkgd calc_bkgd.c $(objects)

clean:
	rm -f $(objects) calc_bkgd
