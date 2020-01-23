# Generic Makefile for PARI programs -- amd64 running linux (x86-64/GMP-6.1.2 kernel) 64-bit version
#
#  This file was created by Configure. Any change made to it will be
#  lost when Configure is run.
#
# make all will create
#  extgcd-dyn (linked dynamically with libpari)
#  extgcd-sta (linked statically)
#  libextgcd.so (to be used by "install" under GP)
#
# Under GP: install("extgcd", "GG&&", "gcdex", "./libextgcd.so") enables
# you to subsequently use gcdex to call extgcd (see the reference manual).
#

# change this TARGET to compile your own programs
SHELL  = /bin/sh

DBGFLAGS   = -g -Wall
CFLAGS     = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer
EXTRACFLAGS=
#CFLAGS    = $(DBGFLAGS)

# Various linkers use different flags to force static compilation. Choose
# the one which is relevant for your installation.
#
# Solaris ld (global)
#STATIC    = -dn

# Solaris ld (toggle: no shared object accepted until -B dynamic is seen
#STATIC    = -B static

# gcc
#STATIC    = -static

CC         = /usr/bin/gcc
CPPFLAGS   = -I. -I/home/nicolas/pari/GP/include
LD         = /usr/bin/gcc
LDFLAGS    = -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer    -Wl,--export-dynamic 
MODLD      = /usr/bin/gcc
MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared 
EXTRAMODLDFLAGS = -lc -lm -L/home/nicolas/pari/GP/lib -lpari
EXTRALIBS  =

RUNPTH     = -Wl,-rpath "/home/nicolas/pari/GP/lib"
DLCFLAGS   = -fPIC
LIBS       = -lm -L/home/nicolas/pari/GP/lib -lpari

RM = rm -f

all: libzn.so libmodjac.so

%.o: %.c
	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
clean:
	-$(RM) *.o $(ALL)


libzn.so: zn.o
	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) zn.o $(EXTRAMODLDFLAGS)

libmodjac.so: modjac.o
	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) modjac.o $(EXTRAMODLDFLAGS)
