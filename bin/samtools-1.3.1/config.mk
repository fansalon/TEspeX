#  Optional configure Makefile overrides for samtools.
#
#    Copyright (C) 2015 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# This is config.mk.  Generated from config.mk.in by configure.
#
# If you use configure, this file overrides variables and augments rules
# in the Makefile to reflect your configuration choices.  If you don't run
# configure, the main Makefile contains suitable conservative defaults.

prefix       = /home/fansalon/Software/TEspeX/bin/samtools-1.3.1
exec_prefix  = ${prefix}
bindir       = ${exec_prefix}/bin
datarootdir  = ${prefix}/share
mandir       = ${datarootdir}/man

CC       = gcc
CPPFLAGS = 
CFLAGS   = -I/home/fansalon/Software/Rpackages/include
LDFLAGS  = -L/home/fansalon/Software/Rpackages/lib
LIBS     = 

HTSDIR = htslib-1.3.1
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LIB = $(HTSLIB)
BGZIP = $(HTSDIR)/bgzip
HTSLIB_CPPFLAGS = -Ihtslib-1.3.1
#HTSLIB_LDFLAGS = -Lhtslib-1.3.1
#HTSLIB_LIB = -lhts

CURSES_LIB = 
