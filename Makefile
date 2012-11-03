
GSL_CFLAGS=`gsl-config --cflags`
GLIB_CFLAGS=`pkg-config --libs`
CFLAGS"=-Wall -O4 $GSL_CFLAGS $GLIB_CFLAGS"

GSL_LFLAGS=`gsl-config --libs`
GLIB_LFLAGS=`pkg-config glib-2.0 --libs`
LFLAGS="-lz $GSL_LFLAGS $GLIB_LFLAGS"
CC=gcc

objects=bkgd.o bkgd_interp.o interp_tab.o bkgd_intg.o bkgd_point.o bkgd_data.o

default: $(objects) calc_bkgd

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

calc_bkgd: $(objects) calc_bkgd.c
	$(CC) $(LFLAGS) -o calc_bkgd calc_bkgd.c $(objects)

clean:
	rm -f $(objects) calc_bkgd
