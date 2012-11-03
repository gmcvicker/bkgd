CFLAGS = -c -Wall -O2
CFLAGS += `pkg-config --libs`
CFLAGS += `gsl-config --cflags`

LFLAGS=-lz
CC=gcc

objects=bkgd.o bkgd_interp.o interp_tab.o bkgd_intg.o bkgd_point.o bkgd_data.o

default: $(objects) calc_bkgd

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

calc_bkgd: $(objects) calc_bkgd.c
	$(CC) $(LFLAGS) -o calc_bkgd calc_bkgd.c $(objects)

clean:
	rm -f $(objects) calc_bkgd
