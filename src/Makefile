###############################################################################
# Sloan Digital Sky Survey (SDSS)
# Support code for "template" product
# Assumes it will be built using "evilmake"
# To incorporate into another product, 
# add any necessary flags and change "SUBDIRS"
# to reflect the list of subdirectories.
###############################################################################

SHELL = /bin/sh

.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c

CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I../include

# List of subdirectories
SUBDIRS = gimg

# build each subdir 
all :
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

# clean each subdir 
clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done
