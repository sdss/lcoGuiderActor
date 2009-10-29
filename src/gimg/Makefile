
###############################################################################
# Sloan Digital Sky Survey (SDSS)
# Silly example for the template product
###############################################################################

SO_EXT = so
SHELL = /bin/sh
SDSS_CFLAGS =  -g -O2 -Wall

BIN = $(GUIDERACTOR_DIR)/bin
LIB = $(GUIDERACTOR_DIR)/lib
INC = $(GUIDERACTOR_DIR)/include

GIMG_TARGET = $(LIB)/libguide.$(SO_EXT)

all: $(GIMG_TARGET)

#
.c.o : 
	$(CC) -c $(CCCHK) $(CFLAGS) $(X_CFLAGS) $*.c

CFLAGS  = $(SDSS_CFLAGS) -fPIC -DCHECK_LEAKS -I$(INC)

OBJECTS = gutils.o shUtils.o ipGguide.o

GIMG_TARGET = $(LIB)/libguide.$(SO_EXT)
TARGETS = $(GIMG_TARGET)

all : $(TARGETS)

$(GIMG_TARGET): $(OBJECTS)
	$(LD) -shared -o $(GIMG_TARGET) $(OBJECTS) -lm

clean :
	- /bin/rm -f $(OBJECTS) $(TARGETS) *~ core *.o

ipGguide.o : $(INC)/shLegacy.h
shUtils.o : $(INC)/ipGguide.h $(INC)/shLegacy.h
