# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
# Copyright (c) 2002-2014 The ROMS/TOMS Group             Kate Hedstrom :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#  ROMS/TOMS Framework Master Makefile                                  :::
#                                                                       :::
#  This makefile is designed to work only with GNU Make version 3.80 or :::
#  higher. It can be used in any architecture provided that there is a  :::
#  machine/compiler rules file in the  "Compilers"  subdirectory.  You  :::
#  may need to modify the rules file to specify the  correct path  for  :::
#  the NetCDF and ARPACK libraries. The ARPACK library is only used in  :::
#  the Generalized Stability Theory analysis and Laczos algorithm.      :::
#                                                                       :::
#  If appropriate,  the USER needs to modify the  macro definitions in  :::
#  in user-defined section below.  To activate an option set the macro  :::
#  to "on". For example, if you want to compile with debugging options  :::
#  set:                                                                 :::
#                                                                       :::
#      USE_DEBUG := on                                                  :::
#                                                                       :::
#  Otherwise, leave macro definition blank.                             :::
#                                                                       :::
#  The USER needs to provide a value for the  macro FORT.  Choose  the  :::
#  appropriate value from the list below.                               :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

NEED_VERSION := 3.80 3.81 3.82 3.82.90 4.0
$(if $(filter $(MAKE_VERSION),$(NEED_VERSION)),,        \
 $(error This makefile requires one of GNU make version $(NEED_VERSION).))

#--------------------------------------------------------------------------
#  Initialize some things.
#--------------------------------------------------------------------------

  sources    :=
  c_sources  :=

#==========================================================================
#  Start of user-defined options. In some macro definitions below: "on" or
#  any other string means TRUE while blank (or spaces) is FALSE.
#==========================================================================
#
#  The CPP option defining a particular application is specified below.
#  See header file "ROMS/Include/cppdefs.h" for all available idealized
#  and realistic applications CPP flags. For example, to activate the
#  upwelling test case (UPWELLING) set:
#
#    ROMS_APPLICATION ?= UPWELLING
#
#  Notice that this makefile will include the associated application header
#  file, which is located either in the "ROMS/Include" or MY_HEADER_DIR
#  directory.  This makefile is designed to search in both directories.
#  The only constrain is that the application CPP option must be unique
#  and header file name is the lowercase value of ROMS_APPLICATION with
#  the .h extension. For example, the upwelling application includes the
#  "upwelling.h" header file.

ROMS_APPLICATION ?= UPWELLING

#  If application header files is not located in "ROMS/Include",
#  provide an alternate directory FULL PATH.

MY_HEADER_DIR ?=

#  If your application requires analytical expressions and they are
#  not located in "ROMS/Functionals", provide an alternate directory.
#  Notice that a set analytical expressions templates can be found in
#  "User/Functionals".
#
#  If applicable, also used this directory to place your customized
#  biology model header file (like fennel.h, nemuro.h, ecosim.h, etc).

MY_ANALYTICAL_DIR ?=

#  Sometimes it is desirable to activate one or more CPP options to
#  run different variants of the same application without modifying
#  its header file. If this is the case, specify such options here
#  using the -D syntax.  For example, to write time-averaged fields
#  set:
#
#    MY_CPP_FLAGS ?= -DAVERAGES
#

MY_CPP_FLAGS ?=

#  Activate debugging compiler options:

   USE_DEBUG ?=

#  If parallel applications, use at most one of these definitions
#  (leave both definitions blank in serial applications):

     USE_MPI ?=
  USE_OpenMP ?=

#  If distributed-memory, turn on compilation via the script "mpif90".
#  This is needed in some Linux operating systems. In some systems with
#  native MPI libraries the compilation does not require MPICH type
#  scripts. This macro is also convient when there are several fortran
#  compiliers (ifort, pgf90, pathf90) in the system that use mpif90.
#  In this, case the user need to select the desired compiler below and
#  turn on both USE_MPI and USE_MPIF90 macros.

  USE_MPIF90 ?=

#  If applicable, activate 64-bit compilation:

   USE_LARGE ?=

#  If applicable, link with NetCDF-4 library. Notice that the NetCDF-4
#  library needs both the HDF5 and MPI libraries.

 USE_NETCDF4 ?= on

#--------------------------------------------------------------------------
#  We are going to include a file with all the settings that depend on
#  the system and the compiler. We are going to build up the name of the
#  include file using information on both. Set your compiler here from
#  the following list:
#
#  Operating System        Compiler(s)
#
#     AIX:                    xlf
#     ALPHA:                  f90
#     CYGWIN:                 g95, df, ifort
#     Darwin:                 f90, xlf
#     IRIX:                   f90
#     Linux:                  ftn, ifc, ifort, pgi, path, g95, gfortran
#     SunOS:                  f95
#     UNICOS-mp:              ftn
#     SunOS/Linux:            ftn (Cray cross-compiler)
#
#  Feel free to send us additional rule files to include! Also, be sure
#  to check the appropriate file to make sure it has the right paths to
#  NetCDF and so on.
#--------------------------------------------------------------------------

        FORT ?= gfortran

#--------------------------------------------------------------------------
#  Set directory for executable.
#--------------------------------------------------------------------------

      BINDIR ?= .

#==========================================================================
#  End of user-defined options. See also the machine-dependent include
#  file being used above.
#==========================================================================

#--------------------------------------------------------------------------
#  Set directory for temporary objects.
#--------------------------------------------------------------------------

SCRATCH_DIR ?= Build
 clean_list := core *.ipo $(SCRATCH_DIR)

ifeq "$(strip $(SCRATCH_DIR))" "."
  clean_list := core *.o *.oo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif
ifeq "$(strip $(SCRATCH_DIR))" "./"
  clean_list := core *.o *.oo *.ipo *.mod *.f90 lib*.a *.bak
  clean_list += $(CURDIR)/*.ipo
endif

#--------------------------------------------------------------------------
#  Notice that the token "libraries" is initialize with the ROMS/Utility
#  library to account for calls to objects in other ROMS libraries or
#  cycling dependencies. These type of dependencies are problematic in
#  some compilers during linking. This library appears twice at linking
#  step (begining and almost the end of ROMS library list).
#--------------------------------------------------------------------------

  libraries  := $(SCRATCH_DIR)/libMODS.a

#--------------------------------------------------------------------------
#  Set Pattern rules.
#--------------------------------------------------------------------------

%.o: %.F

%.o: %.f90
	cd $(SCRATCH_DIR); $(FC) -c $(FFLAGS) $(notdir $<)

%.f90: %.F
	$(CPP) $(CPPFLAGS) $(MY_CPP_FLAGS) $< > $*.f90
	$(CLEAN) $*.f90

CLEAN := ROMS/Bin/cpp_clean

#--------------------------------------------------------------------------
#  Set C-preprocessing flags associated with ROMS application. They are
#  used in "ROMS/Include/cppdefs.h" to include the appropriate application
#  header file.
#--------------------------------------------------------------------------

ifdef ROMS_APPLICATION
        HEADER := $(addsuffix .h, \
			$(shell echo ${ROMS_APPLICATION} | tr [A-Z] [a-z]))
 ROMS_CPPFLAGS := -D$(ROMS_APPLICATION)
 ROMS_CPPFLAGS += -D'HEADER="$(HEADER)"'
 ifdef MY_HEADER_DIR
#   ROMS_CPPFLAGS += -D'ROMS_HEADER="$(MY_HEADER_DIR)/$(HEADER)"'
  ROMS_CPPFLAGS += -I$(MY_HEADER_DIR)
 endif
# else
  ROMS_CPPFLAGS += -D'ROMS_HEADER="$(HEADER)"'
# endif
 ifdef MY_CPP_FLAGS
  ROMS_CPPFLAGS += $(MY_CPP_FLAGS)
 endif
endif

#--------------------------------------------------------------------------
#  Internal macro definitions used to select the code to compile and
#  additional libraries to link. It uses the CPP activated in the
#  header file ROMS/Include/cppdefs.h to determine macro definitions.
#--------------------------------------------------------------------------

  COMPILERS ?= $(CURDIR)/Compilers

MAKE_MACROS := $(shell echo ${HOME} | sed 's| |\\ |g')/make_macros.mk

ifneq "$(MAKECMDGOALS)" "clean"
 MACROS := $(shell cpp -P $(ROMS_CPPFLAGS) Compilers/make_macros.h > \
		$(MAKE_MACROS); $(CLEAN) $(MAKE_MACROS))

 GET_MACROS := $(wildcard $(SCRATCH_DIR)/make_macros.*)

 ifdef GET_MACROS
  include $(SCRATCH_DIR)/make_macros.mk
  $(if ,, $(warning INCLUDING FILE $(SCRATCH_DIR)/make_macros.mk \
                    WHICH CONTAINS APPLICATION-DEPENDENT MAKE DEFINITIONS))
 else
  include $(MAKE_MACROS)
  $(if ,, $(warning INCLUDING FILE $(MAKE_MACROS) \
                   WHICH CONTAINS APPLICATION-DEPENDENT MAKE DEFINITIONS))
 endif
endif

clean_list += $(MAKE_MACROS)

#--------------------------------------------------------------------------
#  Make functions for putting the temporary files in $(SCRATCH_DIR)
#  DO NOT modify this section; spaces and blank lines are needed.
#--------------------------------------------------------------------------

# $(call source-dir-to-binary-dir, directory-list)
source-dir-to-binary-dir = $(addprefix $(SCRATCH_DIR)/, $(notdir $1))

# $(call source-to-object, source-file-list)
source-to-object = $(call source-dir-to-binary-dir,   \
                   $(subst .F,.o,$1))

# $(call source-to-object, source-file-list)
c-source-to-object = $(call source-dir-to-binary-dir,       \
                     $(subst .c,.o,$(filter %.c,$1))        \
                     $(subst .cc,.o,$(filter %.cc,$1)))

# $(call make-library, library-name, source-file-list)
define make-library
   libraries += $(SCRATCH_DIR)/$1
   sources   += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .F,.o,$2))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

# $(call make-c-library, library-name, source-file-list)
define make-c-library
   libraries += $(SCRATCH_DIR)/$1
   c_sources += $2

   $(SCRATCH_DIR)/$1: $(call source-dir-to-binary-dir,    \
                      $(subst .c,.o,$(filter %.c,$2))     \
                      $(subst .cc,.o,$(filter %.cc,$2)))
	$(AR) $(ARFLAGS) $$@ $$^
	$(RANLIB) $$@
endef

# $(call f90-source, source-file-list)
f90-source = $(call source-dir-to-binary-dir,     \
                   $(subst .F,.f90,$1))

# $(compile-rules)
define compile-rules
  $(foreach f, $(local_src),       \
    $(call one-compile-rule,$(call source-to-object,$f), \
    $(call f90-source,$f),$f))
endef

# $(c-compile-rules)
define c-compile-rules
  $(foreach f, $(local_c_src),       \
    $(call one-c-compile-rule,$(call c-source-to-object,$f), $f))
endef

# $(call one-compile-rule, binary-file, f90-file, source-file)
define one-compile-rule
  $1: $2 $3
	cd $$(SCRATCH_DIR); $$(FC) -c $$(FFLAGS) $(notdir $2)

  $2: $3
	$$(CPP) $$(CPPFLAGS) $$(MY_CPP_FLAGS) $$< > $$@
	$$(CLEAN) $$@

endef

# $(call one-c-compile-rule, binary-file, source-file)
define one-c-compile-rule
  $1: $2
	cd $$(SCRATCH_DIR); $$(CXX) -c $$(CXXFLAGS) $$<

endef

#--------------------------------------------------------------------------
#  Set ROMS/TOMS executable file name.
#--------------------------------------------------------------------------

BIN := $(BINDIR)/oceanS
ifdef USE_DEBUG
  BIN := $(BINDIR)/oceanG
else
 ifdef USE_MPI
   BIN := $(BINDIR)/oceanM
 endif
 ifdef USE_OpenMP
   BIN := $(BINDIR)/oceanO
 endif
endif

#--------------------------------------------------------------------------
#  Set name of module files for netCDF F90 interface. On some platforms
#  these will need to be overridden in the machine-dependent include file.
#--------------------------------------------------------------------------

   NETCDF_MODFILE := netcdf.mod
TYPESIZES_MODFILE := typesizes.mod

#--------------------------------------------------------------------------
#  "uname -s" should return the OS or kernel name and "uname -m" should
#  return the CPU or hardware name. In practice the results can be pretty
#  flaky. Run the results through sed to convert "/" and " " to "-",
#  then apply platform-specific conversions.
#--------------------------------------------------------------------------

OS := $(shell uname -s | sed 's/[\/ ]/-/g')
OS := $(patsubst CYGWIN_%,CYGWIN,$(OS))
OS := $(patsubst MINGW%,MINGW,$(OS))
OS := $(patsubst sn%,UNICOS-sn,$(OS))

CPU := $(shell uname -m | sed 's/[\/ ]/-/g')

SVNREV ?= $(shell svnversion -n .)

ROOTDIR := $(shell pwd)

ifndef FORT
  $(error Variable FORT not set)
endif

ifneq "$(MAKECMDGOALS)" "clean"
  include $(COMPILERS)/$(OS)-$(strip $(FORT)).mk
endif

ifdef USE_MPI
 ifdef USE_OpenMP
  $(error You cannot activate USE_MPI and USE_OpenMP at the same time!)
 endif
endif

#--------------------------------------------------------------------------
#  Pass the platform variables to the preprocessor as macros. Convert to
#  valid, upper-case identifiers. Attach ROMS application  CPP options.
#--------------------------------------------------------------------------

CPPFLAGS += -D$(shell echo ${OS} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${CPU} | tr "-" "_" | tr [a-z] [A-Z])
CPPFLAGS += -D$(shell echo ${FORT} | tr "-" "_" | tr [a-z] [A-Z])

CPPFLAGS += -D'ROOT_DIR="$(ROOTDIR)"'
ifdef ROMS_APPLICATION
  CPPFLAGS  += $(ROMS_CPPFLAGS)
  MDEPFLAGS += -DROMS_HEADER="$(HEADER)"
endif

ifndef MY_ANALYTICAL_DIR
  MY_ANALYTICAL_DIR := $(ROOTDIR)/ROMS/Functionals
endif
ifeq (,$(findstring ROMS/Functionals,$(MY_ANALYTICAL_DIR)))
  MY_ANALYTICAL := on
endif
CPPFLAGS += -D'ANALYTICAL_DIR="$(MY_ANALYTICAL_DIR)"'

ifdef MY_ANALYTICAL
  CPPFLAGS += -D'MY_ANALYTICAL="$(MY_ANALYTICAL)"'
endif

ifdef SVNREV
  CPPFLAGS += -D'SVN_REV="$(SVNREV)"'
else
  SVNREV := $(shell grep Revision ./ROMS/Version | sed 's/.* \([0-9]*\) .*/\1/')
  CPPFLAGS += -D'SVN_REV="$(SVNREV)"'
endif

#--------------------------------------------------------------------------
#  Build target directories.
#--------------------------------------------------------------------------

.PHONY: all

all: $(SCRATCH_DIR) $(SCRATCH_DIR)/MakeDepend $(BIN) rm_macros

 modules  :=
ifdef USE_ADJOINT
 modules  +=	ROMS/Adjoint \
		ROMS/Adjoint/Biology
endif
ifdef USE_REPRESENTER
 modules  +=	ROMS/Representer \
		ROMS/Representer/Biology
endif
ifdef USE_TANGENT
 modules  +=	ROMS/Tangent \
		ROMS/Tangent/Biology
endif
 modules  +=	ROMS/Nonlinear \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Sediment \
		ROMS/Functionals
ifdef USE_SEAICE
 modules  +=	ROMS/SeaIce
endif
 modules  +=	ROMS/Utility \
		ROMS/Modules

 includes :=	ROMS/Include
ifdef MY_ANALYTICAL
 includes +=	$(MY_ANALYTICAL_DIR)
endif
ifdef USE_ADJOINT
 includes +=	ROMS/Adjoint \
		ROMS/Adjoint/Biology
endif
ifdef USE_REPRESENTER
 includes +=	ROMS/Representer \
		ROMS/Representer/Biology
endif
ifdef USE_SEAICE
 includes +=	ROMS/SeaIce
endif
ifdef USE_TANGENT
 includes +=	ROMS/Tangent \
		ROMS/Tangent/Biology
endif
 includes +=	ROMS/Nonlinear \
		ROMS/Nonlinear/Biology \
		ROMS/Nonlinear/Sediment \
		ROMS/Utility \
		ROMS/Drivers \
                ROMS/Functionals
ifdef MY_HEADER_DIR
 includes +=	$(MY_HEADER_DIR)
endif

ifdef USE_SWAN
 modules  +=	Waves/SWAN/Src
 includes +=	Waves/SWAN/Src
endif

 modules  +=	Master
 includes +=	Master Compilers

vpath %.F $(modules)
vpath %.cc $(modules)
vpath %.h $(includes)
vpath %.f90 $(SCRATCH_DIR)
vpath %.o $(SCRATCH_DIR)

include $(addsuffix /Module.mk,$(modules))

MDEPFLAGS += $(patsubst %,-I %,$(includes)) --silent --moddir $(SCRATCH_DIR)

CPPFLAGS  += $(patsubst %,-I%,$(includes))

ifdef MY_HEADER_DIR
  CPPFLAGS += -D'HEADER_DIR="$(MY_HEADER_DIR)"'
else
  CPPFLAGS += -D'HEADER_DIR="$(ROOTDIR)/ROMS/Include"'
endif

$(SCRATCH_DIR):
	$(shell $(TEST) -d $(SCRATCH_DIR) || $(MKDIR) $(SCRATCH_DIR) )

#--------------------------------------------------------------------------
#  Add profiling.
#--------------------------------------------------------------------------

# FFLAGS += -check bounds                 # ifort
# FFLAGS += -C                            # pgi
# FFLAGS += -xpg                          # Sun
# FFLAGS += -pg                           # g95
# FFLAGS += -qp                           # ifort
# FFLAGS += -Mprof=func,lines             # pgi
# FFLAGS += -Mprof=mpi,lines              # pgi
# FFLAGS += -Mprof=mpi,hwcts              # pgi
# FFLAGS += -Mprof=func                   # pgi

#--------------------------------------------------------------------------
#  Special CPP macros for mod_strings.F
#--------------------------------------------------------------------------

$(SCRATCH_DIR)/mod_strings.f90: CPPFLAGS += -DMY_OS='"$(OS)"' \
              -DMY_CPU='"$(CPU)"' -DMY_FORT='"$(FORT)"' \
              -DMY_FC='"$(FC)"' -DMY_FFLAGS='"$(FFLAGS)"'

#--------------------------------------------------------------------------
#  ROMS/TOMS libraries.
#--------------------------------------------------------------------------

MYLIB := libocean.a

.PHONY: libraries

libraries: $(libraries)

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS dependecies.
#--------------------------------------------------------------------------

$(SCRATCH_DIR)/$(NETCDF_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(NETCDF_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/$(TYPESIZES_MODFILE): | $(SCRATCH_DIR)
	cp -f $(NETCDF_INCDIR)/$(TYPESIZES_MODFILE) $(SCRATCH_DIR)

$(SCRATCH_DIR)/MakeDepend: makefile \
                           $(SCRATCH_DIR)/$(NETCDF_MODFILE) \
                           $(SCRATCH_DIR)/$(TYPESIZES_MODFILE) \
                           | $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend
	cp -p $(MAKE_MACROS) $(SCRATCH_DIR)

.PHONY: depend

SFMAKEDEPEND := ./ROMS/Bin/sfmakedepend

depend: $(SCRATCH_DIR)
	$(SFMAKEDEPEND) $(MDEPFLAGS) $(sources) > $(SCRATCH_DIR)/MakeDepend

ifneq "$(MAKECMDGOALS)" "clean"
  -include $(SCRATCH_DIR)/MakeDepend
endif

#--------------------------------------------------------------------------
#  Target to create ROMS/TOMS tar file.
#--------------------------------------------------------------------------

.PHONY: tarfile

tarfile:
		tar --exclude=".svn" -cvf roms-3_0.tar *

.PHONY: zipfile

zipfile:
		zip -r roms-3_0.zip *

.PHONY: gzipfile

gzipfile:
		gzip -v roms-3_0.gzip *

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

.PHONY: clean

clean:
	$(RM) -r $(clean_list)

.PHONY: rm_macros

rm_macros:
	$(RM) -r $(MAKE_MACROS)

#--------------------------------------------------------------------------
#  A handy debugging target. This will allow to print the value of any
#  makefile defined macro (see http://tinyurl.com/8ax3j). For example,
#  to find the value of CPPFLAGS execute:
#
#        gmake print-CPPFLAGS
#  or
#        make print-CPPFLAGS
#--------------------------------------------------------------------------

.PHONY: print-%

print-%:
	@echo $* = $($*)
# DO NOT DELETE THIS LINE - used by make depend
esmf_roms.o: cppdefs.h globaldefs.h arctic.h
esmf_roms.o: /center/w/kate/Build_ARC/distribute.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_coupler.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_forces.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_grid.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_iounits.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_ncparam.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_ocean.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_parallel.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_param.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_scalars.o
esmf_roms.o: /center/w/kate/Build_ARC/mod_stepping.o
esmf_roms.o: /center/w/kate/Build_ARC/ocean_control.o
esmf_roms.o: /center/w/kate/Build_ARC/roms_import.o

master.o: esmf_coupler.h ocean.h mct_coupler.h cppdefs.h globaldefs.h arctic.h
master.o: /center/w/kate/Build_ARC/esmf_roms.o
master.o: /center/w/kate/Build_ARC/mod_coupler.o
master.o: /center/w/kate/Build_ARC/mod_iounits.o
master.o: /center/w/kate/Build_ARC/mod_parallel.o
master.o: /center/w/kate/Build_ARC/mod_param.o
master.o: /center/w/kate/Build_ARC/mod_scalars.o
master.o: /center/w/kate/Build_ARC/ocean_control.o
master.o: /center/w/kate/Build_ARC/ocean_coupler.o

ocean_control.o: tl_w4dpsas_ocean.h ad_ocean.h w4dpsas_ocean.h tl_ocean.h
ocean_control.o: correlation.h w4dvar_ocean.h hessian_so_ocean.h afte_ocean.h
ocean_control.o: fsv_ocean.h so_ocean.h obs_sen_w4dpsas.h so_semi_ocean.h
ocean_control.o: is4dvar_ocean.h cppdefs.h globaldefs.h arctic.h pert_ocean.h
ocean_control.o: nl_ocean.h optobs_ocean.h fte_ocean.h op_ocean.h
ocean_control.o: tlcheck_ocean.h obs_sen_w4dvar.h array_modes_w4dvar.h
ocean_control.o: obs_sen_is4dvar.h picard_ocean.h adsen_ocean.h
ocean_control.o: tl_w4dvar_ocean.h rp_ocean.h symmetry.h hessian_op_ocean.h
ocean_control.o: /center/w/kate/Build_ARC/analytical.o
ocean_control.o: /center/w/kate/Build_ARC/array_modes.o
ocean_control.o: /center/w/kate/Build_ARC/back_cost.o
ocean_control.o: /center/w/kate/Build_ARC/cgradient.o
ocean_control.o: /center/w/kate/Build_ARC/convolve.o
ocean_control.o: /center/w/kate/Build_ARC/cost_grad.o
ocean_control.o: /center/w/kate/Build_ARC/distribute.o
ocean_control.o: /center/w/kate/Build_ARC/dotproduct.o
ocean_control.o: /center/w/kate/Build_ARC/ini_adjust.o
ocean_control.o: /center/w/kate/Build_ARC/ini_fields.o
ocean_control.o: /center/w/kate/Build_ARC/mod_boundary.o
ocean_control.o: /center/w/kate/Build_ARC/mod_forces.o
ocean_control.o: /center/w/kate/Build_ARC/mod_fourdvar.o
ocean_control.o: /center/w/kate/Build_ARC/mod_iounits.o
ocean_control.o: /center/w/kate/Build_ARC/mod_ncparam.o
ocean_control.o: /center/w/kate/Build_ARC/mod_netcdf.o
ocean_control.o: /center/w/kate/Build_ARC/mod_ocean.o
ocean_control.o: /center/w/kate/Build_ARC/mod_parallel.o
ocean_control.o: /center/w/kate/Build_ARC/mod_param.o
ocean_control.o: /center/w/kate/Build_ARC/mod_scalars.o
ocean_control.o: /center/w/kate/Build_ARC/mod_stepping.o
ocean_control.o: /center/w/kate/Build_ARC/mod_storage.o
ocean_control.o: /center/w/kate/Build_ARC/normalization.o
ocean_control.o: /center/w/kate/Build_ARC/ocean_coupler.o
ocean_control.o: /center/w/kate/Build_ARC/packing.o
ocean_control.o: /center/w/kate/Build_ARC/posterior.o
ocean_control.o: /center/w/kate/Build_ARC/posterior_var.o
ocean_control.o: /center/w/kate/Build_ARC/propagator.o
ocean_control.o: /center/w/kate/Build_ARC/random_ic.o
ocean_control.o: /center/w/kate/Build_ARC/strings.o
ocean_control.o: /center/w/kate/Build_ARC/sum_grad.o
ocean_control.o: /center/w/kate/Build_ARC/zeta_balance.o

ocean_coupler.o: mct_roms_wrf.h tile.h mct_roms_swan.h cppdefs.h globaldefs.h
ocean_coupler.o: arctic.h
ocean_coupler.o: /center/w/kate/Build_ARC/distribute.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_coupler.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_forces.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_grid.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_iounits.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_kinds.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_ocean.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_parallel.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_param.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_scalars.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_sedbed.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_sediment.o
ocean_coupler.o: /center/w/kate/Build_ARC/mod_stepping.o
ocean_coupler.o: /center/w/kate/Build_ARC/roms_export.o
ocean_coupler.o: /center/w/kate/Build_ARC/roms_import.o

propagator.o: propagator_fte.h propagator_so_semi.h propagator_so.h
propagator.o: propagator_op.h propagator_hso.h propagator_afte.h cppdefs.h
propagator.o: globaldefs.h arctic.h propagator_fsv.h propagator_hop.h
propagator.o: /center/w/kate/Build_ARC/dotproduct.o
propagator.o: /center/w/kate/Build_ARC/ini_adjust.o
propagator.o: /center/w/kate/Build_ARC/inner2state.o
propagator.o: /center/w/kate/Build_ARC/mod_coupling.o
propagator.o: /center/w/kate/Build_ARC/mod_forces.o
propagator.o: /center/w/kate/Build_ARC/mod_iounits.o
propagator.o: /center/w/kate/Build_ARC/mod_kinds.o
propagator.o: /center/w/kate/Build_ARC/mod_netcdf.o
propagator.o: /center/w/kate/Build_ARC/mod_ocean.o
propagator.o: /center/w/kate/Build_ARC/mod_parallel.o
propagator.o: /center/w/kate/Build_ARC/mod_param.o
propagator.o: /center/w/kate/Build_ARC/mod_scalars.o
propagator.o: /center/w/kate/Build_ARC/mod_stepping.o
propagator.o: /center/w/kate/Build_ARC/mod_storage.o
propagator.o: /center/w/kate/Build_ARC/packing.o
propagator.o: /center/w/kate/Build_ARC/set_depth.o

roms_export.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
roms_export.o: /center/w/kate/Build_ARC/distribute.o
roms_export.o: /center/w/kate/Build_ARC/mod_kinds.o
roms_export.o: /center/w/kate/Build_ARC/mod_ncparam.o
roms_export.o: /center/w/kate/Build_ARC/mod_param.o

roms_import.o: cppdefs.h globaldefs.h arctic.h
roms_import.o: /center/w/kate/Build_ARC/distribute.o
roms_import.o: /center/w/kate/Build_ARC/exchange_2d.o
roms_import.o: /center/w/kate/Build_ARC/mod_kinds.o
roms_import.o: /center/w/kate/Build_ARC/mod_ncparam.o
roms_import.o: /center/w/kate/Build_ARC/mod_param.o
roms_import.o: /center/w/kate/Build_ARC/mod_scalars.o
roms_import.o: /center/w/kate/Build_ARC/mp_exchange.o

analytical.o: ana_ssh.h set_bounds.h tile.h ana_trc_psource.h ana_tobc.h
analytical.o: ana_ncep.h ana_fsobc.h ana_passive.h ana_drag.h ana_stflux.h
analytical.o: ana_spinning.h ana_m2clima.h ana_psource.h ana_rain.h ana_cloud.h
analytical.o: ana_specir.h ana_scope.h ana_hiobc.h ana_humid.h ana_initial.h
analytical.o: ana_sst.h ana_grid.h ana_hsnobc.h ana_nudgcoef.h ana_aiobc.h
analytical.o: ana_wwave.h ana_m3obc.h ana_sponge.h ana_btflux.h ana_ice.h
analytical.o: ana_sss.h ana_tclima.h ana_diag.h ana_mask.h ana_pair.h
analytical.o: ana_m2obc.h ana_winds.h ana_biology.h ana_m3clima.h ana_smflux.h
analytical.o: ana_perturb.h ana_vmix.h ana_dqdsst.h ana_tair.h cppdefs.h
analytical.o: globaldefs.h arctic.h ana_srflux.h ana_wtype.h ana_sediment.h
analytical.o: /center/w/kate/Build_ARC/distribute.o
analytical.o: /center/w/kate/Build_ARC/erf.o
analytical.o: /center/w/kate/Build_ARC/exchange_2d.o
analytical.o: /center/w/kate/Build_ARC/exchange_3d.o
analytical.o: /center/w/kate/Build_ARC/mod_biology.o
analytical.o: /center/w/kate/Build_ARC/mod_boundary.o
analytical.o: /center/w/kate/Build_ARC/mod_clima.o
analytical.o: /center/w/kate/Build_ARC/mod_eclight.o
analytical.o: /center/w/kate/Build_ARC/mod_forces.o
analytical.o: /center/w/kate/Build_ARC/mod_grid.o
analytical.o: /center/w/kate/Build_ARC/mod_ice.o
analytical.o: /center/w/kate/Build_ARC/mod_iounits.o
analytical.o: /center/w/kate/Build_ARC/mod_mixing.o
analytical.o: /center/w/kate/Build_ARC/mod_ncparam.o
analytical.o: /center/w/kate/Build_ARC/mod_ocean.o
analytical.o: /center/w/kate/Build_ARC/mod_parallel.o
analytical.o: /center/w/kate/Build_ARC/mod_param.o
analytical.o: /center/w/kate/Build_ARC/mod_scalars.o
analytical.o: /center/w/kate/Build_ARC/mod_sedbed.o
analytical.o: /center/w/kate/Build_ARC/mod_sediment.o
analytical.o: /center/w/kate/Build_ARC/mod_sources.o
analytical.o: /center/w/kate/Build_ARC/mod_stepping.o
analytical.o: /center/w/kate/Build_ARC/mod_trc_sources.o
analytical.o: /center/w/kate/Build_ARC/mp_exchange.o

mod_arrays.o: cppdefs.h globaldefs.h arctic.h
mod_arrays.o: /center/w/kate/Build_ARC/mod_average.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_average2.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_bbl.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_boundary.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_clima.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_coupling.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_diags.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_filter.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_forces.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_grid.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_ice.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_mixing.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_nesting.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_ocean.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_param.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_scalars.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_sedbed.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_sources.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_tides.o
mod_arrays.o: /center/w/kate/Build_ARC/mod_trc_sources.o

mod_average.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_average.o: /center/w/kate/Build_ARC/mod_biology.o
mod_average.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_average.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_average.o: /center/w/kate/Build_ARC/mod_param.o
mod_average.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_average2.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_average2.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_average2.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_average2.o: /center/w/kate/Build_ARC/mod_param.o
mod_average2.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_bbl.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_bbl.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_bbl.o: /center/w/kate/Build_ARC/mod_param.o

mod_behavior.o: oyster_floats_mod.h cppdefs.h globaldefs.h arctic.h
mod_behavior.o: /center/w/kate/Build_ARC/mod_param.o

mod_biology.o: goanpz_mod.h umaine_mod.h npzd_iron_mod.h cppdefs.h globaldefs.h
mod_biology.o: arctic.h npzd_Franks_mod.h npzd_Powell_mod.h nemuro_mod.h
mod_biology.o: ecosim_mod.h bestnpz_mod.h fennel_mod.h
mod_biology.o: /center/w/kate/Build_ARC/mod_param.o
mod_biology.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_boundary.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_boundary.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_boundary.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_boundary.o: /center/w/kate/Build_ARC/mod_param.o
mod_boundary.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_clima.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_clima.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_clima.o: /center/w/kate/Build_ARC/mod_param.o
mod_clima.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_coupler.o: cppdefs.h globaldefs.h arctic.h
mod_coupler.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_coupler.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_coupler.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_coupler.o: /center/w/kate/Build_ARC/mod_param.o
mod_coupler.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_coupling.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_coupling.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_coupling.o: /center/w/kate/Build_ARC/mod_param.o

mod_diags.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_diags.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_diags.o: /center/w/kate/Build_ARC/mod_param.o

mod_eclight.o: cppdefs.h globaldefs.h arctic.h
mod_eclight.o: /center/w/kate/Build_ARC/mod_biology.o
mod_eclight.o: /center/w/kate/Build_ARC/mod_param.o

mod_eoscoef.o: cppdefs.h globaldefs.h arctic.h
mod_eoscoef.o: /center/w/kate/Build_ARC/mod_kinds.o

mod_filter.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_filter.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_filter.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_filter.o: /center/w/kate/Build_ARC/mod_param.o
mod_filter.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_floats.o: cppdefs.h globaldefs.h arctic.h
mod_floats.o: /center/w/kate/Build_ARC/mod_param.o
mod_floats.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_forces.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_forces.o: /center/w/kate/Build_ARC/mod_biology.o
mod_forces.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_forces.o: /center/w/kate/Build_ARC/mod_param.o
mod_forces.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_fourdvar.o: cppdefs.h globaldefs.h arctic.h
mod_fourdvar.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_fourdvar.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_fourdvar.o: /center/w/kate/Build_ARC/mod_netcdf.o
mod_fourdvar.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_fourdvar.o: /center/w/kate/Build_ARC/mod_param.o
mod_fourdvar.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_grid.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_grid.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_grid.o: /center/w/kate/Build_ARC/mod_param.o
mod_grid.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_ice.o: tile.h cppdefs.h globaldefs.h arctic.h
mod_ice.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_ice.o: /center/w/kate/Build_ARC/mod_param.o

mod_iounits.o: cppdefs.h globaldefs.h arctic.h
mod_iounits.o: /center/w/kate/Build_ARC/mod_param.o


mod_mixing.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_mixing.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_mixing.o: /center/w/kate/Build_ARC/mod_param.o
mod_mixing.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_ncparam.o: umaine_var.h ecosim_var.h npzd_Franks_var.h npzd_iron_var.h
mod_ncparam.o: bestnpz_var.h nemuro_var.h goanpz_var.h cppdefs.h globaldefs.h
mod_ncparam.o: arctic.h npzd_Powell_var.h sediment_var.h fennel_var.h
mod_ncparam.o: /center/w/kate/Build_ARC/mod_biology.o
mod_ncparam.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_ncparam.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_ncparam.o: /center/w/kate/Build_ARC/mod_param.o
mod_ncparam.o: /center/w/kate/Build_ARC/mod_scalars.o
mod_ncparam.o: /center/w/kate/Build_ARC/mod_sediment.o

mod_nesting.o: cppdefs.h globaldefs.h arctic.h
mod_nesting.o: /center/w/kate/Build_ARC/mod_boundary.o
mod_nesting.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_nesting.o: /center/w/kate/Build_ARC/mod_param.o
mod_nesting.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_netcdf.o: cppdefs.h globaldefs.h arctic.h
mod_netcdf.o: /center/w/kate/Build_ARC/distribute.o
mod_netcdf.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_netcdf.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_netcdf.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_netcdf.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_netcdf.o: /center/w/kate/Build_ARC/mod_param.o
mod_netcdf.o: /center/w/kate/Build_ARC/mod_scalars.o
mod_netcdf.o: /center/w/kate/Build_ARC/strings.o

mod_ocean.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_ocean.o: /center/w/kate/Build_ARC/mod_biology.o
mod_ocean.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_ocean.o: /center/w/kate/Build_ARC/mod_param.o

mod_parallel.o: cppdefs.h globaldefs.h arctic.h
mod_parallel.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_parallel.o: /center/w/kate/Build_ARC/mod_param.o
mod_parallel.o: /center/w/kate/Build_ARC/mod_scalars.o
mod_parallel.o: /center/w/kate/Build_ARC/mod_strings.o

mod_param.o: cppdefs.h globaldefs.h arctic.h
mod_param.o: /center/w/kate/Build_ARC/mod_kinds.o

mod_scalars.o: cppdefs.h globaldefs.h arctic.h
mod_scalars.o: /center/w/kate/Build_ARC/mod_param.o

mod_sedbed.o: sedbed_mod.h set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_sedbed.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_sedbed.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_sedbed.o: /center/w/kate/Build_ARC/mod_param.o
mod_sedbed.o: /center/w/kate/Build_ARC/mod_sediment.o

mod_sediment.o: sediment_mod.h cppdefs.h globaldefs.h arctic.h
mod_sediment.o: /center/w/kate/Build_ARC/mod_param.o

mod_sources.o: cppdefs.h globaldefs.h arctic.h
mod_sources.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_sources.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_sources.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_sources.o: /center/w/kate/Build_ARC/mod_netcdf.o
mod_sources.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_sources.o: /center/w/kate/Build_ARC/mod_param.o
mod_sources.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_stepping.o: cppdefs.h globaldefs.h arctic.h
mod_stepping.o: /center/w/kate/Build_ARC/mod_param.o

mod_storage.o: cppdefs.h globaldefs.h arctic.h
mod_storage.o: /center/w/kate/Build_ARC/mod_param.o
mod_storage.o: /center/w/kate/Build_ARC/mod_scalars.o

mod_strings.o: cppdefs.h globaldefs.h arctic.h

mod_tides.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mod_tides.o: /center/w/kate/Build_ARC/mod_iounits.o
mod_tides.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_tides.o: /center/w/kate/Build_ARC/mod_ncparam.o
mod_tides.o: /center/w/kate/Build_ARC/mod_netcdf.o
mod_tides.o: /center/w/kate/Build_ARC/mod_parallel.o
mod_tides.o: /center/w/kate/Build_ARC/mod_param.o
mod_tides.o: /center/w/kate/Build_ARC/mod_scalars.o
mod_tides.o: /center/w/kate/Build_ARC/mod_stepping.o

mod_trc_sources.o: cppdefs.h globaldefs.h arctic.h
mod_trc_sources.o: /center/w/kate/Build_ARC/mod_kinds.o
mod_trc_sources.o: /center/w/kate/Build_ARC/mod_param.o

biology.o: umaine.h set_bounds.h tile.h fennel.h bestnpz.h cppdefs.h
biology.o: globaldefs.h arctic.h goanpz.h npzd_iron.h ecosim.h nemuro.h
biology.o: npzd_Franks.h npzd_Powell.h
biology.o: /center/w/kate/Build_ARC/mod_biology.o
biology.o: /center/w/kate/Build_ARC/mod_clima.o
biology.o: /center/w/kate/Build_ARC/mod_diags.o
biology.o: /center/w/kate/Build_ARC/mod_eclight.o
biology.o: /center/w/kate/Build_ARC/mod_forces.o
biology.o: /center/w/kate/Build_ARC/mod_grid.o
biology.o: /center/w/kate/Build_ARC/mod_ice.o
biology.o: /center/w/kate/Build_ARC/mod_iounits.o
biology.o: /center/w/kate/Build_ARC/mod_kinds.o
biology.o: /center/w/kate/Build_ARC/mod_ncparam.o
biology.o: /center/w/kate/Build_ARC/mod_ocean.o
biology.o: /center/w/kate/Build_ARC/mod_parallel.o
biology.o: /center/w/kate/Build_ARC/mod_param.o
biology.o: /center/w/kate/Build_ARC/mod_scalars.o
biology.o: /center/w/kate/Build_ARC/mod_stepping.o

biology_floats.o: oyster_floats.h cppdefs.h globaldefs.h arctic.h
biology_floats.o: /center/w/kate/Build_ARC/mod_behavior.o
biology_floats.o: /center/w/kate/Build_ARC/mod_biology.o
biology_floats.o: /center/w/kate/Build_ARC/mod_floats.o
biology_floats.o: /center/w/kate/Build_ARC/mod_grid.o
biology_floats.o: /center/w/kate/Build_ARC/mod_iounits.o
biology_floats.o: /center/w/kate/Build_ARC/mod_parallel.o
biology_floats.o: /center/w/kate/Build_ARC/mod_param.o
biology_floats.o: /center/w/kate/Build_ARC/mod_scalars.o
biology_floats.o: /center/w/kate/Build_ARC/mod_stepping.o

sed_bed.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
sed_bed.o: /center/w/kate/Build_ARC/bc_3d.o
sed_bed.o: /center/w/kate/Build_ARC/exchange_2d.o
sed_bed.o: /center/w/kate/Build_ARC/mod_bbl.o
sed_bed.o: /center/w/kate/Build_ARC/mod_forces.o
sed_bed.o: /center/w/kate/Build_ARC/mod_grid.o
sed_bed.o: /center/w/kate/Build_ARC/mod_ocean.o
sed_bed.o: /center/w/kate/Build_ARC/mod_param.o
sed_bed.o: /center/w/kate/Build_ARC/mod_scalars.o
sed_bed.o: /center/w/kate/Build_ARC/mod_sedbed.o
sed_bed.o: /center/w/kate/Build_ARC/mod_sediment.o
sed_bed.o: /center/w/kate/Build_ARC/mod_stepping.o
sed_bed.o: /center/w/kate/Build_ARC/mp_exchange.o

sed_bedload.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
sed_bedload.o: /center/w/kate/Build_ARC/bc_3d.o
sed_bedload.o: /center/w/kate/Build_ARC/exchange_2d.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_bbl.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_forces.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_grid.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_ncparam.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_ocean.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_param.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_scalars.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_sedbed.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_sediment.o
sed_bedload.o: /center/w/kate/Build_ARC/mod_stepping.o
sed_bedload.o: /center/w/kate/Build_ARC/mp_exchange.o

sed_fluxes.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
sed_fluxes.o: /center/w/kate/Build_ARC/mod_bbl.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_forces.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_grid.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_ocean.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_param.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_scalars.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_sedbed.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_sediment.o
sed_fluxes.o: /center/w/kate/Build_ARC/mod_stepping.o

sed_settling.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
sed_settling.o: /center/w/kate/Build_ARC/mod_forces.o
sed_settling.o: /center/w/kate/Build_ARC/mod_grid.o
sed_settling.o: /center/w/kate/Build_ARC/mod_ocean.o
sed_settling.o: /center/w/kate/Build_ARC/mod_param.o
sed_settling.o: /center/w/kate/Build_ARC/mod_scalars.o
sed_settling.o: /center/w/kate/Build_ARC/mod_sedbed.o
sed_settling.o: /center/w/kate/Build_ARC/mod_sediment.o
sed_settling.o: /center/w/kate/Build_ARC/mod_stepping.o

sed_surface.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
sed_surface.o: /center/w/kate/Build_ARC/bc_3d.o
sed_surface.o: /center/w/kate/Build_ARC/mod_ocean.o
sed_surface.o: /center/w/kate/Build_ARC/mod_param.o
sed_surface.o: /center/w/kate/Build_ARC/mod_scalars.o
sed_surface.o: /center/w/kate/Build_ARC/mod_sedbed.o
sed_surface.o: /center/w/kate/Build_ARC/mod_sediment.o
sed_surface.o: /center/w/kate/Build_ARC/mod_stepping.o
sed_surface.o: /center/w/kate/Build_ARC/mp_exchange.o

sediment.o: cppdefs.h globaldefs.h arctic.h
sediment.o: /center/w/kate/Build_ARC/sed_bed.o
sediment.o: /center/w/kate/Build_ARC/sed_bedload.o
sediment.o: /center/w/kate/Build_ARC/sed_fluxes.o
sediment.o: /center/w/kate/Build_ARC/sed_settling.o
sediment.o: /center/w/kate/Build_ARC/sed_surface.o

bbl.o: mb_bbl.h set_bounds.h tile.h sg_bbl.h ssw_bbl.h cppdefs.h globaldefs.h
bbl.o: arctic.h
bbl.o: /center/w/kate/Build_ARC/bc_2d.o /center/w/kate/Build_ARC/mod_bbl.o
bbl.o: /center/w/kate/Build_ARC/mod_forces.o
bbl.o: /center/w/kate/Build_ARC/mod_grid.o /center/w/kate/Build_ARC/mod_ocean.o
bbl.o: /center/w/kate/Build_ARC/mod_parallel.o
bbl.o: /center/w/kate/Build_ARC/mod_param.o
bbl.o: /center/w/kate/Build_ARC/mod_scalars.o
bbl.o: /center/w/kate/Build_ARC/mod_sedbed.o
bbl.o: /center/w/kate/Build_ARC/mod_sediment.o
bbl.o: /center/w/kate/Build_ARC/mod_stepping.o
bbl.o: /center/w/kate/Build_ARC/mp_exchange.o

bc_2d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
bc_2d.o: /center/w/kate/Build_ARC/exchange_2d.o
bc_2d.o: /center/w/kate/Build_ARC/mod_boundary.o
bc_2d.o: /center/w/kate/Build_ARC/mod_grid.o
bc_2d.o: /center/w/kate/Build_ARC/mod_ncparam.o
bc_2d.o: /center/w/kate/Build_ARC/mod_param.o
bc_2d.o: /center/w/kate/Build_ARC/mod_scalars.o

bc_3d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
bc_3d.o: /center/w/kate/Build_ARC/exchange_3d.o
bc_3d.o: /center/w/kate/Build_ARC/mod_boundary.o
bc_3d.o: /center/w/kate/Build_ARC/mod_grid.o
bc_3d.o: /center/w/kate/Build_ARC/mod_ncparam.o
bc_3d.o: /center/w/kate/Build_ARC/mod_param.o
bc_3d.o: /center/w/kate/Build_ARC/mod_scalars.o

bc_bry2d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
bc_bry2d.o: /center/w/kate/Build_ARC/mod_param.o
bc_bry2d.o: /center/w/kate/Build_ARC/mod_scalars.o

bc_bry3d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
bc_bry3d.o: /center/w/kate/Build_ARC/mod_param.o
bc_bry3d.o: /center/w/kate/Build_ARC/mod_scalars.o

bulk_flux.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
bulk_flux.o: /center/w/kate/Build_ARC/exchange_2d.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_clima.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_coupling.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_forces.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_grid.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_ice.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_kinds.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_mixing.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_ocean.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_param.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_scalars.o
bulk_flux.o: /center/w/kate/Build_ARC/mod_stepping.o
bulk_flux.o: /center/w/kate/Build_ARC/mp_exchange.o

bvf_mix.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
bvf_mix.o: /center/w/kate/Build_ARC/exchange_3d.o
bvf_mix.o: /center/w/kate/Build_ARC/mod_mixing.o
bvf_mix.o: /center/w/kate/Build_ARC/mod_param.o
bvf_mix.o: /center/w/kate/Build_ARC/mod_scalars.o
bvf_mix.o: /center/w/kate/Build_ARC/mp_exchange.o

ccsm_flux.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
ccsm_flux.o: /center/w/kate/Build_ARC/exchange_2d.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_biology.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_clima.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_coupling.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_forces.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_grid.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_ice.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_mixing.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_ocean.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_param.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_scalars.o
ccsm_flux.o: /center/w/kate/Build_ARC/mod_stepping.o
ccsm_flux.o: /center/w/kate/Build_ARC/mp_exchange.o

conv_2d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
conv_2d.o: /center/w/kate/Build_ARC/bc_2d.o
conv_2d.o: /center/w/kate/Build_ARC/mod_param.o
conv_2d.o: /center/w/kate/Build_ARC/mod_scalars.o
conv_2d.o: /center/w/kate/Build_ARC/mp_exchange.o

conv_3d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
conv_3d.o: /center/w/kate/Build_ARC/bc_3d.o
conv_3d.o: /center/w/kate/Build_ARC/mod_param.o
conv_3d.o: /center/w/kate/Build_ARC/mod_scalars.o
conv_3d.o: /center/w/kate/Build_ARC/mp_exchange.o

conv_bry2d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
conv_bry2d.o: /center/w/kate/Build_ARC/bc_bry2d.o
conv_bry2d.o: /center/w/kate/Build_ARC/mod_param.o
conv_bry2d.o: /center/w/kate/Build_ARC/mod_scalars.o
conv_bry2d.o: /center/w/kate/Build_ARC/mp_exchange.o

conv_bry3d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
conv_bry3d.o: /center/w/kate/Build_ARC/bc_bry3d.o
conv_bry3d.o: /center/w/kate/Build_ARC/mod_param.o
conv_bry3d.o: /center/w/kate/Build_ARC/mod_scalars.o
conv_bry3d.o: /center/w/kate/Build_ARC/mp_exchange.o

diag.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
diag.o: /center/w/kate/Build_ARC/analytical.o
diag.o: /center/w/kate/Build_ARC/distribute.o
diag.o: /center/w/kate/Build_ARC/mod_biology.o
diag.o: /center/w/kate/Build_ARC/mod_grid.o
diag.o: /center/w/kate/Build_ARC/mod_iounits.o
diag.o: /center/w/kate/Build_ARC/mod_ocean.o
diag.o: /center/w/kate/Build_ARC/mod_parallel.o
diag.o: /center/w/kate/Build_ARC/mod_param.o
diag.o: /center/w/kate/Build_ARC/mod_scalars.o
diag.o: /center/w/kate/Build_ARC/mod_stepping.o

exchange_2d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
exchange_2d.o: /center/w/kate/Build_ARC/mod_param.o
exchange_2d.o: /center/w/kate/Build_ARC/mod_scalars.o

exchange_3d.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
exchange_3d.o: /center/w/kate/Build_ARC/mod_param.o
exchange_3d.o: /center/w/kate/Build_ARC/mod_scalars.o

forcing.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
forcing.o: /center/w/kate/Build_ARC/mod_coupling.o
forcing.o: /center/w/kate/Build_ARC/mod_iounits.o
forcing.o: /center/w/kate/Build_ARC/mod_ocean.o
forcing.o: /center/w/kate/Build_ARC/mod_parallel.o
forcing.o: /center/w/kate/Build_ARC/mod_param.o
forcing.o: /center/w/kate/Build_ARC/mod_scalars.o

frc_adjust.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
frc_adjust.o: /center/w/kate/Build_ARC/mod_forces.o
frc_adjust.o: /center/w/kate/Build_ARC/mod_param.o
frc_adjust.o: /center/w/kate/Build_ARC/mod_scalars.o

get_data.o: cppdefs.h globaldefs.h arctic.h
get_data.o: /center/w/kate/Build_ARC/mod_biology.o
get_data.o: /center/w/kate/Build_ARC/mod_boundary.o
get_data.o: /center/w/kate/Build_ARC/mod_clima.o
get_data.o: /center/w/kate/Build_ARC/mod_forces.o
get_data.o: /center/w/kate/Build_ARC/mod_grid.o
get_data.o: /center/w/kate/Build_ARC/mod_iounits.o
get_data.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_data.o: /center/w/kate/Build_ARC/mod_ocean.o
get_data.o: /center/w/kate/Build_ARC/mod_param.o
get_data.o: /center/w/kate/Build_ARC/mod_scalars.o
get_data.o: /center/w/kate/Build_ARC/mod_sources.o
get_data.o: /center/w/kate/Build_ARC/mod_stepping.o

get_idata.o: cppdefs.h globaldefs.h arctic.h
get_idata.o: /center/w/kate/Build_ARC/mod_grid.o
get_idata.o: /center/w/kate/Build_ARC/mod_iounits.o
get_idata.o: /center/w/kate/Build_ARC/mod_mixing.o
get_idata.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_idata.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_idata.o: /center/w/kate/Build_ARC/mod_parallel.o
get_idata.o: /center/w/kate/Build_ARC/mod_param.o
get_idata.o: /center/w/kate/Build_ARC/mod_scalars.o
get_idata.o: /center/w/kate/Build_ARC/mod_sources.o
get_idata.o: /center/w/kate/Build_ARC/mod_stepping.o
get_idata.o: /center/w/kate/Build_ARC/mod_tides.o
get_idata.o: /center/w/kate/Build_ARC/nf_fread3d.o
get_idata.o: /center/w/kate/Build_ARC/nf_fread4d.o

gls_corstep.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
gls_corstep.o: /center/w/kate/Build_ARC/exchange_3d.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_forces.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_grid.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_mixing.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_ncparam.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_ocean.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_param.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_scalars.o
gls_corstep.o: /center/w/kate/Build_ARC/mod_stepping.o
gls_corstep.o: /center/w/kate/Build_ARC/mp_exchange.o
gls_corstep.o: /center/w/kate/Build_ARC/tkebc_im.o

gls_prestep.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
gls_prestep.o: /center/w/kate/Build_ARC/exchange_3d.o
gls_prestep.o: /center/w/kate/Build_ARC/mod_grid.o
gls_prestep.o: /center/w/kate/Build_ARC/mod_mixing.o
gls_prestep.o: /center/w/kate/Build_ARC/mod_ocean.o
gls_prestep.o: /center/w/kate/Build_ARC/mod_param.o
gls_prestep.o: /center/w/kate/Build_ARC/mod_scalars.o
gls_prestep.o: /center/w/kate/Build_ARC/mod_stepping.o
gls_prestep.o: /center/w/kate/Build_ARC/mp_exchange.o
gls_prestep.o: /center/w/kate/Build_ARC/tkebc_im.o

hmixing.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
hmixing.o: /center/w/kate/Build_ARC/exchange_3d.o
hmixing.o: /center/w/kate/Build_ARC/mod_grid.o
hmixing.o: /center/w/kate/Build_ARC/mod_mixing.o
hmixing.o: /center/w/kate/Build_ARC/mod_ncparam.o
hmixing.o: /center/w/kate/Build_ARC/mod_ocean.o
hmixing.o: /center/w/kate/Build_ARC/mod_param.o
hmixing.o: /center/w/kate/Build_ARC/mod_scalars.o
hmixing.o: /center/w/kate/Build_ARC/mod_stepping.o
hmixing.o: /center/w/kate/Build_ARC/mp_exchange.o

ini_fields.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
ini_fields.o: /center/w/kate/Build_ARC/exchange_2d.o
ini_fields.o: /center/w/kate/Build_ARC/exchange_3d.o
ini_fields.o: /center/w/kate/Build_ARC/mod_clima.o
ini_fields.o: /center/w/kate/Build_ARC/mod_coupling.o
ini_fields.o: /center/w/kate/Build_ARC/mod_grid.o
ini_fields.o: /center/w/kate/Build_ARC/mod_ice.o
ini_fields.o: /center/w/kate/Build_ARC/mod_ncparam.o
ini_fields.o: /center/w/kate/Build_ARC/mod_ocean.o
ini_fields.o: /center/w/kate/Build_ARC/mod_param.o
ini_fields.o: /center/w/kate/Build_ARC/mod_scalars.o
ini_fields.o: /center/w/kate/Build_ARC/mod_sedbed.o
ini_fields.o: /center/w/kate/Build_ARC/mod_sediment.o
ini_fields.o: /center/w/kate/Build_ARC/mod_stepping.o
ini_fields.o: /center/w/kate/Build_ARC/mp_exchange.o
ini_fields.o: /center/w/kate/Build_ARC/pt3dbc_im.o
ini_fields.o: /center/w/kate/Build_ARC/t3dbc_im.o
ini_fields.o: /center/w/kate/Build_ARC/u2dbc_im.o
ini_fields.o: /center/w/kate/Build_ARC/u3dbc_im.o
ini_fields.o: /center/w/kate/Build_ARC/v2dbc_im.o
ini_fields.o: /center/w/kate/Build_ARC/v3dbc_im.o
ini_fields.o: /center/w/kate/Build_ARC/zetabc.o

initial.o: cppdefs.h globaldefs.h arctic.h
initial.o: /center/w/kate/Build_ARC/analytical.o
initial.o: /center/w/kate/Build_ARC/distribute.o
initial.o: /center/w/kate/Build_ARC/ini_adjust.o
initial.o: /center/w/kate/Build_ARC/ini_hmixcoef.o
initial.o: /center/w/kate/Build_ARC/metrics.o
initial.o: /center/w/kate/Build_ARC/mod_bbl.o
initial.o: /center/w/kate/Build_ARC/mod_fourdvar.o
initial.o: /center/w/kate/Build_ARC/mod_grid.o
initial.o: /center/w/kate/Build_ARC/mod_iounits.o
initial.o: /center/w/kate/Build_ARC/mod_ncparam.o
initial.o: /center/w/kate/Build_ARC/mod_nesting.o
initial.o: /center/w/kate/Build_ARC/mod_ocean.o
initial.o: /center/w/kate/Build_ARC/mod_parallel.o
initial.o: /center/w/kate/Build_ARC/mod_param.o
initial.o: /center/w/kate/Build_ARC/mod_scalars.o
initial.o: /center/w/kate/Build_ARC/mod_stepping.o
initial.o: /center/w/kate/Build_ARC/nesting.o
initial.o: /center/w/kate/Build_ARC/ocean_coupler.o
initial.o: /center/w/kate/Build_ARC/omega.o /center/w/kate/Build_ARC/rho_eos.o
initial.o: /center/w/kate/Build_ARC/set_depth.o
initial.o: /center/w/kate/Build_ARC/set_masks.o
initial.o: /center/w/kate/Build_ARC/set_massflux.o
initial.o: /center/w/kate/Build_ARC/stiffness.o
initial.o: /center/w/kate/Build_ARC/wetdry.o /center/w/kate/Build_ARC/wpoints.o

interp_floats.o: cppdefs.h globaldefs.h arctic.h
interp_floats.o: /center/w/kate/Build_ARC/mod_floats.o
interp_floats.o: /center/w/kate/Build_ARC/mod_ncparam.o
interp_floats.o: /center/w/kate/Build_ARC/mod_param.o
interp_floats.o: /center/w/kate/Build_ARC/mod_scalars.o

interp_floats_diapW.o: cppdefs.h globaldefs.h arctic.h
interp_floats_diapW.o: /center/w/kate/Build_ARC/mod_floats.o
interp_floats_diapW.o: /center/w/kate/Build_ARC/mod_ncparam.o
interp_floats_diapW.o: /center/w/kate/Build_ARC/mod_param.o
interp_floats_diapW.o: /center/w/kate/Build_ARC/mod_scalars.o

lmd_bkpp.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
lmd_bkpp.o: /center/w/kate/Build_ARC/bc_2d.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_forces.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_grid.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_mixing.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_ocean.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_param.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_scalars.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mod_stepping.o
lmd_bkpp.o: /center/w/kate/Build_ARC/mp_exchange.o
lmd_bkpp.o: /center/w/kate/Build_ARC/shapiro.o

lmd_skpp.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
lmd_skpp.o: /center/w/kate/Build_ARC/bc_2d.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_clima.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_forces.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_grid.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_mixing.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_ocean.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_param.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_scalars.o
lmd_skpp.o: /center/w/kate/Build_ARC/mod_stepping.o
lmd_skpp.o: /center/w/kate/Build_ARC/mp_exchange.o
lmd_skpp.o: /center/w/kate/Build_ARC/shapiro.o

lmd_swfrac.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
lmd_swfrac.o: /center/w/kate/Build_ARC/mod_mixing.o
lmd_swfrac.o: /center/w/kate/Build_ARC/mod_param.o
lmd_swfrac.o: /center/w/kate/Build_ARC/mod_scalars.o

lmd_vmix.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
lmd_vmix.o: /center/w/kate/Build_ARC/bc_3d.o
lmd_vmix.o: /center/w/kate/Build_ARC/lmd_bkpp.o
lmd_vmix.o: /center/w/kate/Build_ARC/lmd_skpp.o
lmd_vmix.o: /center/w/kate/Build_ARC/mod_grid.o
lmd_vmix.o: /center/w/kate/Build_ARC/mod_mixing.o
lmd_vmix.o: /center/w/kate/Build_ARC/mod_ocean.o
lmd_vmix.o: /center/w/kate/Build_ARC/mod_param.o
lmd_vmix.o: /center/w/kate/Build_ARC/mod_scalars.o
lmd_vmix.o: /center/w/kate/Build_ARC/mod_stepping.o
lmd_vmix.o: /center/w/kate/Build_ARC/mp_exchange.o

main2d.o: cppdefs.h globaldefs.h arctic.h
main2d.o: /center/w/kate/Build_ARC/bulk_flux.o
main2d.o: /center/w/kate/Build_ARC/ccsm_flux.o /center/w/kate/Build_ARC/diag.o
main2d.o: /center/w/kate/Build_ARC/dotproduct.o
main2d.o: /center/w/kate/Build_ARC/forcing.o
main2d.o: /center/w/kate/Build_ARC/frc_adjust.o
main2d.o: /center/w/kate/Build_ARC/ini_fields.o
main2d.o: /center/w/kate/Build_ARC/mod_coupler.o
main2d.o: /center/w/kate/Build_ARC/mod_iounits.o
main2d.o: /center/w/kate/Build_ARC/mod_nesting.o
main2d.o: /center/w/kate/Build_ARC/mod_parallel.o
main2d.o: /center/w/kate/Build_ARC/mod_param.o
main2d.o: /center/w/kate/Build_ARC/mod_scalars.o
main2d.o: /center/w/kate/Build_ARC/mod_stepping.o
main2d.o: /center/w/kate/Build_ARC/nesting.o
main2d.o: /center/w/kate/Build_ARC/obc_adjust.o
main2d.o: /center/w/kate/Build_ARC/ocean_coupler.o
main2d.o: /center/w/kate/Build_ARC/radiation_stress.o
main2d.o: /center/w/kate/Build_ARC/set_avg.o
main2d.o: /center/w/kate/Build_ARC/set_tides.o
main2d.o: /center/w/kate/Build_ARC/set_vbc.o /center/w/kate/Build_ARC/step2d.o
main2d.o: /center/w/kate/Build_ARC/step_floats.o

main3d.o: cppdefs.h globaldefs.h arctic.h
main3d.o: /center/w/kate/Build_ARC/analytical.o /center/w/kate/Build_ARC/bbl.o
main3d.o: /center/w/kate/Build_ARC/biology.o
main3d.o: /center/w/kate/Build_ARC/bulk_flux.o
main3d.o: /center/w/kate/Build_ARC/bvf_mix.o
main3d.o: /center/w/kate/Build_ARC/cawdir_eval.o
main3d.o: /center/w/kate/Build_ARC/ccsm_flux.o /center/w/kate/Build_ARC/diag.o
main3d.o: /center/w/kate/Build_ARC/dotproduct.o
main3d.o: /center/w/kate/Build_ARC/forcing.o
main3d.o: /center/w/kate/Build_ARC/frc_adjust.o
main3d.o: /center/w/kate/Build_ARC/gls_corstep.o
main3d.o: /center/w/kate/Build_ARC/gls_prestep.o
main3d.o: /center/w/kate/Build_ARC/hmixing.o
main3d.o: /center/w/kate/Build_ARC/ini_fields.o
main3d.o: /center/w/kate/Build_ARC/lmd_vmix.o
main3d.o: /center/w/kate/Build_ARC/mod_coupler.o
main3d.o: /center/w/kate/Build_ARC/mod_iounits.o
main3d.o: /center/w/kate/Build_ARC/mod_nesting.o
main3d.o: /center/w/kate/Build_ARC/mod_parallel.o
main3d.o: /center/w/kate/Build_ARC/mod_param.o
main3d.o: /center/w/kate/Build_ARC/mod_scalars.o
main3d.o: /center/w/kate/Build_ARC/mod_stepping.o
main3d.o: /center/w/kate/Build_ARC/my25_corstep.o
main3d.o: /center/w/kate/Build_ARC/my25_prestep.o
main3d.o: /center/w/kate/Build_ARC/nesting.o
main3d.o: /center/w/kate/Build_ARC/obc_adjust.o
main3d.o: /center/w/kate/Build_ARC/ocean_coupler.o
main3d.o: /center/w/kate/Build_ARC/omega.o
main3d.o: /center/w/kate/Build_ARC/radiation_stress.o
main3d.o: /center/w/kate/Build_ARC/rho_eos.o /center/w/kate/Build_ARC/rhs3d.o
main3d.o: /center/w/kate/Build_ARC/sediment.o
main3d.o: /center/w/kate/Build_ARC/set_avg.o
main3d.o: /center/w/kate/Build_ARC/set_avg2.o
main3d.o: /center/w/kate/Build_ARC/set_depth.o
main3d.o: /center/w/kate/Build_ARC/set_massflux.o
main3d.o: /center/w/kate/Build_ARC/set_tides.o
main3d.o: /center/w/kate/Build_ARC/set_vbc.o
main3d.o: /center/w/kate/Build_ARC/set_zeta.o /center/w/kate/Build_ARC/step2d.o
main3d.o: /center/w/kate/Build_ARC/step3d_t.o
main3d.o: /center/w/kate/Build_ARC/step3d_uv.o
main3d.o: /center/w/kate/Build_ARC/step_floats.o
main3d.o: /center/w/kate/Build_ARC/wvelocity.o

main3d_offline.o: cppdefs.h globaldefs.h arctic.h
main3d_offline.o: /center/w/kate/Build_ARC/analytical.o
main3d_offline.o: /center/w/kate/Build_ARC/bbl.o
main3d_offline.o: /center/w/kate/Build_ARC/biology.o
main3d_offline.o: /center/w/kate/Build_ARC/bulk_flux.o
main3d_offline.o: /center/w/kate/Build_ARC/bvf_mix.o
main3d_offline.o: /center/w/kate/Build_ARC/cawdir_eval.o
main3d_offline.o: /center/w/kate/Build_ARC/ccsm_flux.o
main3d_offline.o: /center/w/kate/Build_ARC/diag.o
main3d_offline.o: /center/w/kate/Build_ARC/dotproduct.o
main3d_offline.o: /center/w/kate/Build_ARC/gls_corstep.o
main3d_offline.o: /center/w/kate/Build_ARC/gls_prestep.o
main3d_offline.o: /center/w/kate/Build_ARC/hmixing.o
main3d_offline.o: /center/w/kate/Build_ARC/ini_fields.o
main3d_offline.o: /center/w/kate/Build_ARC/lmd_vmix.o
main3d_offline.o: /center/w/kate/Build_ARC/mod_coupler.o
main3d_offline.o: /center/w/kate/Build_ARC/mod_iounits.o
main3d_offline.o: /center/w/kate/Build_ARC/mod_parallel.o
main3d_offline.o: /center/w/kate/Build_ARC/mod_param.o
main3d_offline.o: /center/w/kate/Build_ARC/mod_scalars.o
main3d_offline.o: /center/w/kate/Build_ARC/mod_stepping.o
main3d_offline.o: /center/w/kate/Build_ARC/my25_corstep.o
main3d_offline.o: /center/w/kate/Build_ARC/my25_prestep.o
main3d_offline.o: /center/w/kate/Build_ARC/nesting.o
main3d_offline.o: /center/w/kate/Build_ARC/ocean_coupler.o
main3d_offline.o: /center/w/kate/Build_ARC/omega.o
main3d_offline.o: /center/w/kate/Build_ARC/radiation_stress.o
main3d_offline.o: /center/w/kate/Build_ARC/rho_eos.o
main3d_offline.o: /center/w/kate/Build_ARC/rhs3d.o
main3d_offline.o: /center/w/kate/Build_ARC/sediment.o
main3d_offline.o: /center/w/kate/Build_ARC/set_avg.o
main3d_offline.o: /center/w/kate/Build_ARC/set_avg2.o
main3d_offline.o: /center/w/kate/Build_ARC/set_massflux.o
main3d_offline.o: /center/w/kate/Build_ARC/set_tides.o
main3d_offline.o: /center/w/kate/Build_ARC/set_vbc.o
main3d_offline.o: /center/w/kate/Build_ARC/set_zeta.o
main3d_offline.o: /center/w/kate/Build_ARC/step2d.o
main3d_offline.o: /center/w/kate/Build_ARC/step3d_t.o
main3d_offline.o: /center/w/kate/Build_ARC/step3d_uv.o
main3d_offline.o: /center/w/kate/Build_ARC/step_floats.o
main3d_offline.o: /center/w/kate/Build_ARC/wvelocity.o

mpdata_adiff.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mpdata_adiff.o: /center/w/kate/Build_ARC/mod_ncparam.o
mpdata_adiff.o: /center/w/kate/Build_ARC/mod_param.o
mpdata_adiff.o: /center/w/kate/Build_ARC/mod_scalars.o

my25_corstep.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
my25_corstep.o: /center/w/kate/Build_ARC/exchange_3d.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_forces.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_grid.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_mixing.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_ncparam.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_ocean.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_param.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_scalars.o
my25_corstep.o: /center/w/kate/Build_ARC/mod_stepping.o
my25_corstep.o: /center/w/kate/Build_ARC/mp_exchange.o
my25_corstep.o: /center/w/kate/Build_ARC/tkebc_im.o

my25_prestep.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
my25_prestep.o: /center/w/kate/Build_ARC/exchange_3d.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_grid.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_mixing.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_ncparam.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_ocean.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_param.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_scalars.o
my25_prestep.o: /center/w/kate/Build_ARC/mod_stepping.o
my25_prestep.o: /center/w/kate/Build_ARC/mp_exchange.o
my25_prestep.o: /center/w/kate/Build_ARC/tkebc_im.o

nesting.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
nesting.o: /center/w/kate/Build_ARC/distribute.o
nesting.o: /center/w/kate/Build_ARC/exchange_2d.o
nesting.o: /center/w/kate/Build_ARC/exchange_3d.o
nesting.o: /center/w/kate/Build_ARC/mod_clima.o
nesting.o: /center/w/kate/Build_ARC/mod_coupling.o
nesting.o: /center/w/kate/Build_ARC/mod_forces.o
nesting.o: /center/w/kate/Build_ARC/mod_grid.o
nesting.o: /center/w/kate/Build_ARC/mod_ncparam.o
nesting.o: /center/w/kate/Build_ARC/mod_nesting.o
nesting.o: /center/w/kate/Build_ARC/mod_ocean.o
nesting.o: /center/w/kate/Build_ARC/mod_parallel.o
nesting.o: /center/w/kate/Build_ARC/mod_param.o
nesting.o: /center/w/kate/Build_ARC/mod_scalars.o
nesting.o: /center/w/kate/Build_ARC/mod_stepping.o
nesting.o: /center/w/kate/Build_ARC/mp_exchange.o
nesting.o: /center/w/kate/Build_ARC/set_depth.o

obc_adjust.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
obc_adjust.o: /center/w/kate/Build_ARC/mod_boundary.o
obc_adjust.o: /center/w/kate/Build_ARC/mod_ncparam.o
obc_adjust.o: /center/w/kate/Build_ARC/mod_param.o
obc_adjust.o: /center/w/kate/Build_ARC/mod_scalars.o

obc_volcons.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
obc_volcons.o: /center/w/kate/Build_ARC/distribute.o
obc_volcons.o: /center/w/kate/Build_ARC/mod_grid.o
obc_volcons.o: /center/w/kate/Build_ARC/mod_ocean.o
obc_volcons.o: /center/w/kate/Build_ARC/mod_parallel.o
obc_volcons.o: /center/w/kate/Build_ARC/mod_param.o
obc_volcons.o: /center/w/kate/Build_ARC/mod_scalars.o
obc_volcons.o: /center/w/kate/Build_ARC/mp_exchange.o

omega.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
omega.o: /center/w/kate/Build_ARC/bc_3d.o
omega.o: /center/w/kate/Build_ARC/exchange_3d.o
omega.o: /center/w/kate/Build_ARC/mod_grid.o
omega.o: /center/w/kate/Build_ARC/mod_ncparam.o
omega.o: /center/w/kate/Build_ARC/mod_ocean.o
omega.o: /center/w/kate/Build_ARC/mod_param.o
omega.o: /center/w/kate/Build_ARC/mod_scalars.o
omega.o: /center/w/kate/Build_ARC/mod_sedbed.o
omega.o: /center/w/kate/Build_ARC/mod_sources.o
omega.o: /center/w/kate/Build_ARC/mod_stepping.o
omega.o: /center/w/kate/Build_ARC/mp_exchange.o

output.o: cppdefs.h globaldefs.h arctic.h
output.o: /center/w/kate/Build_ARC/distribute.o
output.o: /center/w/kate/Build_ARC/mod_filter.o
output.o: /center/w/kate/Build_ARC/mod_floats.o
output.o: /center/w/kate/Build_ARC/mod_fourdvar.o
output.o: /center/w/kate/Build_ARC/mod_iounits.o
output.o: /center/w/kate/Build_ARC/mod_ncparam.o
output.o: /center/w/kate/Build_ARC/mod_netcdf.o
output.o: /center/w/kate/Build_ARC/mod_parallel.o
output.o: /center/w/kate/Build_ARC/mod_param.o
output.o: /center/w/kate/Build_ARC/mod_scalars.o

pre_step3d.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
pre_step3d.o: /center/w/kate/Build_ARC/exchange_3d.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_biology.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_clima.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_diags.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_forces.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_grid.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_mixing.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_ocean.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_param.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_scalars.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_sources.o
pre_step3d.o: /center/w/kate/Build_ARC/mod_stepping.o
pre_step3d.o: /center/w/kate/Build_ARC/mp_exchange.o
pre_step3d.o: /center/w/kate/Build_ARC/pt3dbc_im.o
pre_step3d.o: /center/w/kate/Build_ARC/t3dbc_im.o

prsgrd.o: prsgrd32.h set_bounds.h tile.h prsgrd31.h prsgrd42.h prsgrd40.h
prsgrd.o: prsgrd44.h cppdefs.h globaldefs.h arctic.h
prsgrd.o: /center/w/kate/Build_ARC/mod_diags.o
prsgrd.o: /center/w/kate/Build_ARC/mod_forces.o
prsgrd.o: /center/w/kate/Build_ARC/mod_grid.o
prsgrd.o: /center/w/kate/Build_ARC/mod_ocean.o
prsgrd.o: /center/w/kate/Build_ARC/mod_param.o
prsgrd.o: /center/w/kate/Build_ARC/mod_scalars.o
prsgrd.o: /center/w/kate/Build_ARC/mod_stepping.o
prsgrd.o: /center/w/kate/Build_ARC/mod_tides.o

pt3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
pt3dbc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
pt3dbc_im.o: /center/w/kate/Build_ARC/mod_grid.o
pt3dbc_im.o: /center/w/kate/Build_ARC/mod_ocean.o
pt3dbc_im.o: /center/w/kate/Build_ARC/mod_param.o
pt3dbc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
pt3dbc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

radiation_stress.o: nearshore_mellor08.h set_bounds.h tile.h
radiation_stress.o: nearshore_mellor05.h cppdefs.h globaldefs.h arctic.h
radiation_stress.o: /center/w/kate/Build_ARC/bc_2d.o
radiation_stress.o: /center/w/kate/Build_ARC/bc_3d.o
radiation_stress.o: /center/w/kate/Build_ARC/exchange_2d.o
radiation_stress.o: /center/w/kate/Build_ARC/exchange_3d.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_coupling.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_diags.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_forces.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_grid.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_mixing.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_ocean.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_param.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_scalars.o
radiation_stress.o: /center/w/kate/Build_ARC/mod_stepping.o
radiation_stress.o: /center/w/kate/Build_ARC/mp_exchange.o

rho_eos.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
rho_eos.o: /center/w/kate/Build_ARC/exchange_2d.o
rho_eos.o: /center/w/kate/Build_ARC/exchange_3d.o
rho_eos.o: /center/w/kate/Build_ARC/mod_coupling.o
rho_eos.o: /center/w/kate/Build_ARC/mod_eoscoef.o
rho_eos.o: /center/w/kate/Build_ARC/mod_grid.o
rho_eos.o: /center/w/kate/Build_ARC/mod_mixing.o
rho_eos.o: /center/w/kate/Build_ARC/mod_ocean.o
rho_eos.o: /center/w/kate/Build_ARC/mod_param.o
rho_eos.o: /center/w/kate/Build_ARC/mod_scalars.o
rho_eos.o: /center/w/kate/Build_ARC/mod_sediment.o
rho_eos.o: /center/w/kate/Build_ARC/mod_stepping.o
rho_eos.o: /center/w/kate/Build_ARC/mp_exchange.o

rhs3d.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
rhs3d.o: /center/w/kate/Build_ARC/mod_clima.o
rhs3d.o: /center/w/kate/Build_ARC/mod_coupling.o
rhs3d.o: /center/w/kate/Build_ARC/mod_diags.o
rhs3d.o: /center/w/kate/Build_ARC/mod_forces.o
rhs3d.o: /center/w/kate/Build_ARC/mod_grid.o
rhs3d.o: /center/w/kate/Build_ARC/mod_mixing.o
rhs3d.o: /center/w/kate/Build_ARC/mod_ocean.o
rhs3d.o: /center/w/kate/Build_ARC/mod_param.o
rhs3d.o: /center/w/kate/Build_ARC/mod_scalars.o
rhs3d.o: /center/w/kate/Build_ARC/mod_stepping.o
rhs3d.o: /center/w/kate/Build_ARC/pre_step3d.o
rhs3d.o: /center/w/kate/Build_ARC/prsgrd.o /center/w/kate/Build_ARC/t3dmix.o
rhs3d.o: /center/w/kate/Build_ARC/uv3dmix.o

set_avg.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_avg.o: /center/w/kate/Build_ARC/exchange_2d.o
set_avg.o: /center/w/kate/Build_ARC/exchange_3d.o
set_avg.o: /center/w/kate/Build_ARC/mod_average.o
set_avg.o: /center/w/kate/Build_ARC/mod_biology.o
set_avg.o: /center/w/kate/Build_ARC/mod_coupling.o
set_avg.o: /center/w/kate/Build_ARC/mod_forces.o
set_avg.o: /center/w/kate/Build_ARC/mod_grid.o
set_avg.o: /center/w/kate/Build_ARC/mod_ice.o
set_avg.o: /center/w/kate/Build_ARC/mod_mixing.o
set_avg.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_avg.o: /center/w/kate/Build_ARC/mod_ocean.o
set_avg.o: /center/w/kate/Build_ARC/mod_param.o
set_avg.o: /center/w/kate/Build_ARC/mod_scalars.o
set_avg.o: /center/w/kate/Build_ARC/mod_sedbed.o
set_avg.o: /center/w/kate/Build_ARC/mod_sediment.o
set_avg.o: /center/w/kate/Build_ARC/mod_stepping.o
set_avg.o: /center/w/kate/Build_ARC/mod_tides.o
set_avg.o: /center/w/kate/Build_ARC/mp_exchange.o
set_avg.o: /center/w/kate/Build_ARC/set_masks.o
set_avg.o: /center/w/kate/Build_ARC/uv_rotate.o
set_avg.o: /center/w/kate/Build_ARC/vorticity.o

set_avg2.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_avg2.o: /center/w/kate/Build_ARC/mod_average2.o
set_avg2.o: /center/w/kate/Build_ARC/mod_coupling.o
set_avg2.o: /center/w/kate/Build_ARC/mod_forces.o
set_avg2.o: /center/w/kate/Build_ARC/mod_grid.o
set_avg2.o: /center/w/kate/Build_ARC/mod_ice.o
set_avg2.o: /center/w/kate/Build_ARC/mod_mixing.o
set_avg2.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_avg2.o: /center/w/kate/Build_ARC/mod_ocean.o
set_avg2.o: /center/w/kate/Build_ARC/mod_param.o
set_avg2.o: /center/w/kate/Build_ARC/mod_scalars.o
set_avg2.o: /center/w/kate/Build_ARC/mod_stepping.o
set_avg2.o: /center/w/kate/Build_ARC/set_masks.o
set_avg2.o: /center/w/kate/Build_ARC/uv_rotate.o

set_data.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_data.o: /center/w/kate/Build_ARC/analytical.o
set_data.o: /center/w/kate/Build_ARC/distribute.o
set_data.o: /center/w/kate/Build_ARC/exchange_2d.o
set_data.o: /center/w/kate/Build_ARC/exchange_3d.o
set_data.o: /center/w/kate/Build_ARC/mod_biology.o
set_data.o: /center/w/kate/Build_ARC/mod_boundary.o
set_data.o: /center/w/kate/Build_ARC/mod_clima.o
set_data.o: /center/w/kate/Build_ARC/mod_forces.o
set_data.o: /center/w/kate/Build_ARC/mod_grid.o
set_data.o: /center/w/kate/Build_ARC/mod_mixing.o
set_data.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_data.o: /center/w/kate/Build_ARC/mod_ocean.o
set_data.o: /center/w/kate/Build_ARC/mod_param.o
set_data.o: /center/w/kate/Build_ARC/mod_scalars.o
set_data.o: /center/w/kate/Build_ARC/mod_sources.o
set_data.o: /center/w/kate/Build_ARC/mod_stepping.o
set_data.o: /center/w/kate/Build_ARC/mod_trc_sources.o
set_data.o: /center/w/kate/Build_ARC/mp_exchange.o
set_data.o: /center/w/kate/Build_ARC/set_2dfld.o
set_data.o: /center/w/kate/Build_ARC/set_3dfld.o

set_depth.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_depth.o: /center/w/kate/Build_ARC/exchange_2d.o
set_depth.o: /center/w/kate/Build_ARC/exchange_3d.o
set_depth.o: /center/w/kate/Build_ARC/mod_boundary.o
set_depth.o: /center/w/kate/Build_ARC/mod_coupling.o
set_depth.o: /center/w/kate/Build_ARC/mod_grid.o
set_depth.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_depth.o: /center/w/kate/Build_ARC/mod_ocean.o
set_depth.o: /center/w/kate/Build_ARC/mod_param.o
set_depth.o: /center/w/kate/Build_ARC/mod_scalars.o
set_depth.o: /center/w/kate/Build_ARC/mod_sedbed.o
set_depth.o: /center/w/kate/Build_ARC/mod_stepping.o
set_depth.o: /center/w/kate/Build_ARC/mp_exchange.o

set_massflux.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_massflux.o: /center/w/kate/Build_ARC/exchange_3d.o
set_massflux.o: /center/w/kate/Build_ARC/mod_coupling.o
set_massflux.o: /center/w/kate/Build_ARC/mod_grid.o
set_massflux.o: /center/w/kate/Build_ARC/mod_ocean.o
set_massflux.o: /center/w/kate/Build_ARC/mod_param.o
set_massflux.o: /center/w/kate/Build_ARC/mod_scalars.o
set_massflux.o: /center/w/kate/Build_ARC/mod_stepping.o
set_massflux.o: /center/w/kate/Build_ARC/mp_exchange.o

set_tides.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_tides.o: /center/w/kate/Build_ARC/distribute.o
set_tides.o: /center/w/kate/Build_ARC/exchange_2d.o
set_tides.o: /center/w/kate/Build_ARC/mod_boundary.o
set_tides.o: /center/w/kate/Build_ARC/mod_clima.o
set_tides.o: /center/w/kate/Build_ARC/mod_grid.o
set_tides.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_tides.o: /center/w/kate/Build_ARC/mod_ocean.o
set_tides.o: /center/w/kate/Build_ARC/mod_param.o
set_tides.o: /center/w/kate/Build_ARC/mod_scalars.o
set_tides.o: /center/w/kate/Build_ARC/mod_stepping.o
set_tides.o: /center/w/kate/Build_ARC/mod_tides.o
set_tides.o: /center/w/kate/Build_ARC/mp_exchange.o

set_vbc.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_vbc.o: /center/w/kate/Build_ARC/bc_2d.o
set_vbc.o: /center/w/kate/Build_ARC/mod_forces.o
set_vbc.o: /center/w/kate/Build_ARC/mod_grid.o
set_vbc.o: /center/w/kate/Build_ARC/mod_ice.o
set_vbc.o: /center/w/kate/Build_ARC/mod_ocean.o
set_vbc.o: /center/w/kate/Build_ARC/mod_param.o
set_vbc.o: /center/w/kate/Build_ARC/mod_scalars.o
set_vbc.o: /center/w/kate/Build_ARC/mod_stepping.o
set_vbc.o: /center/w/kate/Build_ARC/mp_exchange.o

set_zeta.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_zeta.o: /center/w/kate/Build_ARC/exchange_2d.o
set_zeta.o: /center/w/kate/Build_ARC/mod_coupling.o
set_zeta.o: /center/w/kate/Build_ARC/mod_forces.o
set_zeta.o: /center/w/kate/Build_ARC/mod_grid.o
set_zeta.o: /center/w/kate/Build_ARC/mod_ocean.o
set_zeta.o: /center/w/kate/Build_ARC/mod_param.o
set_zeta.o: /center/w/kate/Build_ARC/mod_scalars.o
set_zeta.o: /center/w/kate/Build_ARC/mp_exchange.o

step2d.o: step2d_LF_AM3.h set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
step2d.o: /center/w/kate/Build_ARC/exchange_2d.o
step2d.o: /center/w/kate/Build_ARC/mod_clima.o
step2d.o: /center/w/kate/Build_ARC/mod_coupling.o
step2d.o: /center/w/kate/Build_ARC/mod_diags.o
step2d.o: /center/w/kate/Build_ARC/mod_forces.o
step2d.o: /center/w/kate/Build_ARC/mod_grid.o
step2d.o: /center/w/kate/Build_ARC/mod_mixing.o
step2d.o: /center/w/kate/Build_ARC/mod_ncparam.o
step2d.o: /center/w/kate/Build_ARC/mod_ocean.o
step2d.o: /center/w/kate/Build_ARC/mod_param.o
step2d.o: /center/w/kate/Build_ARC/mod_scalars.o
step2d.o: /center/w/kate/Build_ARC/mod_sedbed.o
step2d.o: /center/w/kate/Build_ARC/mod_sediment.o
step2d.o: /center/w/kate/Build_ARC/mod_sources.o
step2d.o: /center/w/kate/Build_ARC/mod_stepping.o
step2d.o: /center/w/kate/Build_ARC/mp_exchange.o
step2d.o: /center/w/kate/Build_ARC/obc_volcons.o
step2d.o: /center/w/kate/Build_ARC/u2dbc_im.o
step2d.o: /center/w/kate/Build_ARC/v2dbc_im.o /center/w/kate/Build_ARC/wetdry.o
step2d.o: /center/w/kate/Build_ARC/zetabc.o

step3d_t.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
step3d_t.o: /center/w/kate/Build_ARC/exchange_3d.o
step3d_t.o: /center/w/kate/Build_ARC/mod_biology.o
step3d_t.o: /center/w/kate/Build_ARC/mod_clima.o
step3d_t.o: /center/w/kate/Build_ARC/mod_diags.o
step3d_t.o: /center/w/kate/Build_ARC/mod_grid.o
step3d_t.o: /center/w/kate/Build_ARC/mod_mixing.o
step3d_t.o: /center/w/kate/Build_ARC/mod_ncparam.o
step3d_t.o: /center/w/kate/Build_ARC/mod_nesting.o
step3d_t.o: /center/w/kate/Build_ARC/mod_ocean.o
step3d_t.o: /center/w/kate/Build_ARC/mod_param.o
step3d_t.o: /center/w/kate/Build_ARC/mod_scalars.o
step3d_t.o: /center/w/kate/Build_ARC/mod_sources.o
step3d_t.o: /center/w/kate/Build_ARC/mod_stepping.o
step3d_t.o: /center/w/kate/Build_ARC/mod_trc_sources.o
step3d_t.o: /center/w/kate/Build_ARC/mp_exchange.o
step3d_t.o: /center/w/kate/Build_ARC/mpdata_adiff.o
step3d_t.o: /center/w/kate/Build_ARC/nesting.o
step3d_t.o: /center/w/kate/Build_ARC/pt3dbc_im.o
step3d_t.o: /center/w/kate/Build_ARC/t3dbc_im.o

step3d_uv.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
step3d_uv.o: /center/w/kate/Build_ARC/exchange_2d.o
step3d_uv.o: /center/w/kate/Build_ARC/exchange_3d.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_coupling.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_diags.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_forces.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_grid.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_mixing.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_ncparam.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_ocean.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_param.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_scalars.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_sources.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_stepping.o
step3d_uv.o: /center/w/kate/Build_ARC/mod_tides.o
step3d_uv.o: /center/w/kate/Build_ARC/mp_exchange.o
step3d_uv.o: /center/w/kate/Build_ARC/u3dbc_im.o
step3d_uv.o: /center/w/kate/Build_ARC/v3dbc_im.o

step_floats.o: cppdefs.h globaldefs.h arctic.h
step_floats.o: /center/w/kate/Build_ARC/biology_floats.o
step_floats.o: /center/w/kate/Build_ARC/distribute.o
step_floats.o: /center/w/kate/Build_ARC/interp_floats.o
step_floats.o: /center/w/kate/Build_ARC/interp_floats_diapW.o
step_floats.o: /center/w/kate/Build_ARC/mod_floats.o
step_floats.o: /center/w/kate/Build_ARC/mod_grid.o
step_floats.o: /center/w/kate/Build_ARC/mod_iounits.o
step_floats.o: /center/w/kate/Build_ARC/mod_ncparam.o
step_floats.o: /center/w/kate/Build_ARC/mod_ocean.o
step_floats.o: /center/w/kate/Build_ARC/mod_parallel.o
step_floats.o: /center/w/kate/Build_ARC/mod_param.o
step_floats.o: /center/w/kate/Build_ARC/mod_scalars.o
step_floats.o: /center/w/kate/Build_ARC/mod_stepping.o
step_floats.o: /center/w/kate/Build_ARC/vwalk_floats.o

t3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
t3dbc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_clima.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_grid.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_ncparam.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_ocean.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_param.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
t3dbc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

t3dmix.o: t3dmix4_geo.h set_bounds.h tile.h t3dmix2_geo.h t3dmix2_s.h
t3dmix.o: t3dmix4_iso.h t3dmix4_s.h t3dmix2_iso.h cppdefs.h globaldefs.h
t3dmix.o: arctic.h
t3dmix.o: /center/w/kate/Build_ARC/mod_biology.o
t3dmix.o: /center/w/kate/Build_ARC/mod_clima.o
t3dmix.o: /center/w/kate/Build_ARC/mod_diags.o
t3dmix.o: /center/w/kate/Build_ARC/mod_grid.o
t3dmix.o: /center/w/kate/Build_ARC/mod_mixing.o
t3dmix.o: /center/w/kate/Build_ARC/mod_ncparam.o
t3dmix.o: /center/w/kate/Build_ARC/mod_ocean.o
t3dmix.o: /center/w/kate/Build_ARC/mod_param.o
t3dmix.o: /center/w/kate/Build_ARC/mod_scalars.o
t3dmix.o: /center/w/kate/Build_ARC/mod_stepping.o

tkebc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
tkebc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
tkebc_im.o: /center/w/kate/Build_ARC/mod_grid.o
tkebc_im.o: /center/w/kate/Build_ARC/mod_mixing.o
tkebc_im.o: /center/w/kate/Build_ARC/mod_ncparam.o
tkebc_im.o: /center/w/kate/Build_ARC/mod_param.o
tkebc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
tkebc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

u2dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
u2dbc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_clima.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_forces.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_grid.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_ncparam.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_ocean.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_param.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
u2dbc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

u3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
u3dbc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_clima.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_grid.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_ncparam.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_ocean.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_param.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
u3dbc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

uv3dmix.o: uv3dmix2_s.h set_bounds.h tile.h uv3dmix2_geo.h uv3dmix4_s.h
uv3dmix.o: uv3dmix4_geo.h cppdefs.h globaldefs.h arctic.h
uv3dmix.o: /center/w/kate/Build_ARC/mod_coupling.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_diags.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_grid.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_mixing.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_ncparam.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_ocean.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_param.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_scalars.o
uv3dmix.o: /center/w/kate/Build_ARC/mod_stepping.o

v2dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
v2dbc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_clima.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_forces.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_grid.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_ncparam.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_ocean.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_param.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
v2dbc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

v3dbc_im.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
v3dbc_im.o: /center/w/kate/Build_ARC/mod_boundary.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_clima.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_grid.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_ncparam.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_ocean.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_param.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_scalars.o
v3dbc_im.o: /center/w/kate/Build_ARC/mod_stepping.o

vwalk_floats.o: cppdefs.h globaldefs.h arctic.h
vwalk_floats.o: /center/w/kate/Build_ARC/distribute.o
vwalk_floats.o: /center/w/kate/Build_ARC/interp_floats.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_floats.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_grid.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_mixing.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_ncparam.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_ocean.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_parallel.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_param.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_scalars.o
vwalk_floats.o: /center/w/kate/Build_ARC/mod_stepping.o
vwalk_floats.o: /center/w/kate/Build_ARC/nrutil.o

wetdry.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
wetdry.o: /center/w/kate/Build_ARC/exchange_2d.o
wetdry.o: /center/w/kate/Build_ARC/mod_coupling.o
wetdry.o: /center/w/kate/Build_ARC/mod_grid.o
wetdry.o: /center/w/kate/Build_ARC/mod_ncparam.o
wetdry.o: /center/w/kate/Build_ARC/mod_ocean.o
wetdry.o: /center/w/kate/Build_ARC/mod_param.o
wetdry.o: /center/w/kate/Build_ARC/mod_scalars.o
wetdry.o: /center/w/kate/Build_ARC/mod_sources.o
wetdry.o: /center/w/kate/Build_ARC/mp_exchange.o

wvelocity.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
wvelocity.o: /center/w/kate/Build_ARC/bc_3d.o
wvelocity.o: /center/w/kate/Build_ARC/exchange_2d.o
wvelocity.o: /center/w/kate/Build_ARC/mod_coupling.o
wvelocity.o: /center/w/kate/Build_ARC/mod_grid.o
wvelocity.o: /center/w/kate/Build_ARC/mod_ncparam.o
wvelocity.o: /center/w/kate/Build_ARC/mod_ocean.o
wvelocity.o: /center/w/kate/Build_ARC/mod_param.o
wvelocity.o: /center/w/kate/Build_ARC/mod_scalars.o
wvelocity.o: /center/w/kate/Build_ARC/mod_stepping.o
wvelocity.o: /center/w/kate/Build_ARC/mp_exchange.o

zetabc.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
zetabc.o: /center/w/kate/Build_ARC/mod_boundary.o
zetabc.o: /center/w/kate/Build_ARC/mod_grid.o
zetabc.o: /center/w/kate/Build_ARC/mod_ncparam.o
zetabc.o: /center/w/kate/Build_ARC/mod_ocean.o
zetabc.o: /center/w/kate/Build_ARC/mod_param.o
zetabc.o: /center/w/kate/Build_ARC/mod_scalars.o
zetabc.o: /center/w/kate/Build_ARC/mod_stepping.o

abort.o: cppdefs.h globaldefs.h arctic.h
abort.o: /center/w/kate/Build_ARC/ocean_control.o

array_modes.o: cppdefs.h globaldefs.h arctic.h
array_modes.o: /center/w/kate/Build_ARC/distribute.o
array_modes.o: /center/w/kate/Build_ARC/mod_fourdvar.o
array_modes.o: /center/w/kate/Build_ARC/mod_iounits.o
array_modes.o: /center/w/kate/Build_ARC/mod_parallel.o
array_modes.o: /center/w/kate/Build_ARC/mod_param.o
array_modes.o: /center/w/kate/Build_ARC/mod_scalars.o

back_cost.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
back_cost.o: /center/w/kate/Build_ARC/distribute.o
back_cost.o: /center/w/kate/Build_ARC/mod_boundary.o
back_cost.o: /center/w/kate/Build_ARC/mod_forces.o
back_cost.o: /center/w/kate/Build_ARC/mod_fourdvar.o
back_cost.o: /center/w/kate/Build_ARC/mod_grid.o
back_cost.o: /center/w/kate/Build_ARC/mod_ncparam.o
back_cost.o: /center/w/kate/Build_ARC/mod_ocean.o
back_cost.o: /center/w/kate/Build_ARC/mod_parallel.o
back_cost.o: /center/w/kate/Build_ARC/mod_param.o
back_cost.o: /center/w/kate/Build_ARC/mod_scalars.o

cawdir_eval.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
cawdir_eval.o: /center/w/kate/Build_ARC/bc_2d.o
cawdir_eval.o: /center/w/kate/Build_ARC/mod_forces.o
cawdir_eval.o: /center/w/kate/Build_ARC/mod_grid.o
cawdir_eval.o: /center/w/kate/Build_ARC/mod_ice.o
cawdir_eval.o: /center/w/kate/Build_ARC/mod_param.o
cawdir_eval.o: /center/w/kate/Build_ARC/mod_scalars.o
cawdir_eval.o: /center/w/kate/Build_ARC/mod_stepping.o
cawdir_eval.o: /center/w/kate/Build_ARC/mp_exchange.o

cgradient.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
cgradient.o: /center/w/kate/Build_ARC/distribute.o
cgradient.o: /center/w/kate/Build_ARC/mod_boundary.o
cgradient.o: /center/w/kate/Build_ARC/mod_coupling.o
cgradient.o: /center/w/kate/Build_ARC/mod_forces.o
cgradient.o: /center/w/kate/Build_ARC/mod_fourdvar.o
cgradient.o: /center/w/kate/Build_ARC/mod_grid.o
cgradient.o: /center/w/kate/Build_ARC/mod_iounits.o
cgradient.o: /center/w/kate/Build_ARC/mod_ncparam.o
cgradient.o: /center/w/kate/Build_ARC/mod_netcdf.o
cgradient.o: /center/w/kate/Build_ARC/mod_ocean.o
cgradient.o: /center/w/kate/Build_ARC/mod_parallel.o
cgradient.o: /center/w/kate/Build_ARC/mod_param.o
cgradient.o: /center/w/kate/Build_ARC/mod_scalars.o
cgradient.o: /center/w/kate/Build_ARC/mod_stepping.o
cgradient.o: /center/w/kate/Build_ARC/nf_fread2d.o
cgradient.o: /center/w/kate/Build_ARC/nf_fread2d_bry.o
cgradient.o: /center/w/kate/Build_ARC/nf_fread3d.o
cgradient.o: /center/w/kate/Build_ARC/nf_fread3d_bry.o
cgradient.o: /center/w/kate/Build_ARC/state_addition.o
cgradient.o: /center/w/kate/Build_ARC/state_copy.o
cgradient.o: /center/w/kate/Build_ARC/state_dotprod.o
cgradient.o: /center/w/kate/Build_ARC/state_initialize.o
cgradient.o: /center/w/kate/Build_ARC/state_scale.o

check_multifile.o: cppdefs.h globaldefs.h arctic.h
check_multifile.o: /center/w/kate/Build_ARC/mod_iounits.o
check_multifile.o: /center/w/kate/Build_ARC/mod_netcdf.o
check_multifile.o: /center/w/kate/Build_ARC/mod_parallel.o
check_multifile.o: /center/w/kate/Build_ARC/mod_param.o
check_multifile.o: /center/w/kate/Build_ARC/mod_scalars.o

checkadj.o: cppdefs.h globaldefs.h arctic.h
checkadj.o: /center/w/kate/Build_ARC/mod_iounits.o
checkadj.o: /center/w/kate/Build_ARC/mod_ncparam.o
checkadj.o: /center/w/kate/Build_ARC/mod_parallel.o
checkadj.o: /center/w/kate/Build_ARC/mod_param.o
checkadj.o: /center/w/kate/Build_ARC/mod_scalars.o
checkadj.o: /center/w/kate/Build_ARC/mod_strings.o
checkadj.o: /center/w/kate/Build_ARC/strings.o

checkdefs.o: cppdefs.h globaldefs.h arctic.h
checkdefs.o: /center/w/kate/Build_ARC/mod_iounits.o
checkdefs.o: /center/w/kate/Build_ARC/mod_parallel.o
checkdefs.o: /center/w/kate/Build_ARC/mod_param.o
checkdefs.o: /center/w/kate/Build_ARC/mod_scalars.o
checkdefs.o: /center/w/kate/Build_ARC/mod_strings.o
checkdefs.o: /center/w/kate/Build_ARC/strings.o

checkerror.o: cppdefs.h globaldefs.h arctic.h
checkerror.o: /center/w/kate/Build_ARC/mod_iounits.o
checkerror.o: /center/w/kate/Build_ARC/mod_parallel.o
checkerror.o: /center/w/kate/Build_ARC/mod_param.o
checkerror.o: /center/w/kate/Build_ARC/mod_scalars.o

checkvars.o: cppdefs.h globaldefs.h arctic.h
checkvars.o: /center/w/kate/Build_ARC/mod_biology.o
checkvars.o: /center/w/kate/Build_ARC/mod_iounits.o
checkvars.o: /center/w/kate/Build_ARC/mod_ncparam.o
checkvars.o: /center/w/kate/Build_ARC/mod_netcdf.o
checkvars.o: /center/w/kate/Build_ARC/mod_parallel.o
checkvars.o: /center/w/kate/Build_ARC/mod_param.o
checkvars.o: /center/w/kate/Build_ARC/mod_scalars.o
checkvars.o: /center/w/kate/Build_ARC/mod_sediment.o

close_io.o: cppdefs.h globaldefs.h arctic.h
close_io.o: /center/w/kate/Build_ARC/mod_biology.o
close_io.o: /center/w/kate/Build_ARC/mod_iounits.o
close_io.o: /center/w/kate/Build_ARC/mod_ncparam.o
close_io.o: /center/w/kate/Build_ARC/mod_netcdf.o
close_io.o: /center/w/kate/Build_ARC/mod_parallel.o
close_io.o: /center/w/kate/Build_ARC/mod_param.o
close_io.o: /center/w/kate/Build_ARC/mod_scalars.o

congrad.o: cppdefs.h globaldefs.h arctic.h
congrad.o: /center/w/kate/Build_ARC/distribute.o
congrad.o: /center/w/kate/Build_ARC/mod_fourdvar.o
congrad.o: /center/w/kate/Build_ARC/mod_iounits.o
congrad.o: /center/w/kate/Build_ARC/mod_ncparam.o
congrad.o: /center/w/kate/Build_ARC/mod_netcdf.o
congrad.o: /center/w/kate/Build_ARC/mod_parallel.o
congrad.o: /center/w/kate/Build_ARC/mod_param.o
congrad.o: /center/w/kate/Build_ARC/mod_scalars.o

convolve.o: cppdefs.h globaldefs.h arctic.h
convolve.o: /center/w/kate/Build_ARC/ini_adjust.o
convolve.o: /center/w/kate/Build_ARC/mod_forces.o
convolve.o: /center/w/kate/Build_ARC/mod_iounits.o
convolve.o: /center/w/kate/Build_ARC/mod_ocean.o
convolve.o: /center/w/kate/Build_ARC/mod_parallel.o
convolve.o: /center/w/kate/Build_ARC/mod_param.o
convolve.o: /center/w/kate/Build_ARC/mod_scalars.o
convolve.o: /center/w/kate/Build_ARC/mod_stepping.o

cost_grad.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
cost_grad.o: /center/w/kate/Build_ARC/mod_boundary.o
cost_grad.o: /center/w/kate/Build_ARC/mod_forces.o
cost_grad.o: /center/w/kate/Build_ARC/mod_ncparam.o
cost_grad.o: /center/w/kate/Build_ARC/mod_ocean.o
cost_grad.o: /center/w/kate/Build_ARC/mod_param.o
cost_grad.o: /center/w/kate/Build_ARC/mod_scalars.o

def_avg.o: cppdefs.h globaldefs.h arctic.h
def_avg.o: /center/w/kate/Build_ARC/def_var.o
def_avg.o: /center/w/kate/Build_ARC/mod_biology.o
def_avg.o: /center/w/kate/Build_ARC/mod_filter.o
def_avg.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_avg.o: /center/w/kate/Build_ARC/mod_iounits.o
def_avg.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_avg.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_avg.o: /center/w/kate/Build_ARC/mod_parallel.o
def_avg.o: /center/w/kate/Build_ARC/mod_param.o
def_avg.o: /center/w/kate/Build_ARC/mod_scalars.o
def_avg.o: /center/w/kate/Build_ARC/mod_sediment.o

def_avg2.o: cppdefs.h globaldefs.h arctic.h
def_avg2.o: /center/w/kate/Build_ARC/def_var.o
def_avg2.o: /center/w/kate/Build_ARC/distribute.o
def_avg2.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_avg2.o: /center/w/kate/Build_ARC/mod_iounits.o
def_avg2.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_avg2.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_avg2.o: /center/w/kate/Build_ARC/mod_parallel.o
def_avg2.o: /center/w/kate/Build_ARC/mod_param.o
def_avg2.o: /center/w/kate/Build_ARC/mod_scalars.o
def_avg2.o: /center/w/kate/Build_ARC/mod_sediment.o

def_diags.o: cppdefs.h globaldefs.h arctic.h
def_diags.o: /center/w/kate/Build_ARC/def_var.o
def_diags.o: /center/w/kate/Build_ARC/mod_biology.o
def_diags.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_diags.o: /center/w/kate/Build_ARC/mod_iounits.o
def_diags.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_diags.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_diags.o: /center/w/kate/Build_ARC/mod_parallel.o
def_diags.o: /center/w/kate/Build_ARC/mod_param.o
def_diags.o: /center/w/kate/Build_ARC/mod_scalars.o
def_diags.o: /center/w/kate/Build_ARC/mod_sediment.o

def_dim.o: cppdefs.h globaldefs.h arctic.h
def_dim.o: /center/w/kate/Build_ARC/distribute.o
def_dim.o: /center/w/kate/Build_ARC/mod_iounits.o
def_dim.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_dim.o: /center/w/kate/Build_ARC/mod_parallel.o
def_dim.o: /center/w/kate/Build_ARC/mod_param.o
def_dim.o: /center/w/kate/Build_ARC/mod_scalars.o

def_error.o: cppdefs.h globaldefs.h arctic.h
def_error.o: /center/w/kate/Build_ARC/def_var.o
def_error.o: /center/w/kate/Build_ARC/mod_biology.o
def_error.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_error.o: /center/w/kate/Build_ARC/mod_iounits.o
def_error.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_error.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_error.o: /center/w/kate/Build_ARC/mod_parallel.o
def_error.o: /center/w/kate/Build_ARC/mod_param.o
def_error.o: /center/w/kate/Build_ARC/mod_scalars.o
def_error.o: /center/w/kate/Build_ARC/mod_sediment.o

def_filt.o: cppdefs.h globaldefs.h arctic.h
def_filt.o: /center/w/kate/Build_ARC/def_var.o
def_filt.o: /center/w/kate/Build_ARC/distribute.o
def_filt.o: /center/w/kate/Build_ARC/mod_filter.o
def_filt.o: /center/w/kate/Build_ARC/mod_iounits.o
def_filt.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_filt.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_filt.o: /center/w/kate/Build_ARC/mod_parallel.o
def_filt.o: /center/w/kate/Build_ARC/mod_param.o
def_filt.o: /center/w/kate/Build_ARC/mod_scalars.o

def_floats.o: cppdefs.h globaldefs.h arctic.h
def_floats.o: /center/w/kate/Build_ARC/def_var.o
def_floats.o: /center/w/kate/Build_ARC/mod_biology.o
def_floats.o: /center/w/kate/Build_ARC/mod_floats.o
def_floats.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_floats.o: /center/w/kate/Build_ARC/mod_grid.o
def_floats.o: /center/w/kate/Build_ARC/mod_iounits.o
def_floats.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_floats.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_floats.o: /center/w/kate/Build_ARC/mod_parallel.o
def_floats.o: /center/w/kate/Build_ARC/mod_param.o
def_floats.o: /center/w/kate/Build_ARC/mod_scalars.o
def_floats.o: /center/w/kate/Build_ARC/mod_sediment.o

def_gst.o: cppdefs.h globaldefs.h arctic.h
def_gst.o: /center/w/kate/Build_ARC/def_var.o
def_gst.o: /center/w/kate/Build_ARC/mod_iounits.o
def_gst.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_gst.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_gst.o: /center/w/kate/Build_ARC/mod_parallel.o
def_gst.o: /center/w/kate/Build_ARC/mod_param.o
def_gst.o: /center/w/kate/Build_ARC/mod_scalars.o
def_gst.o: /center/w/kate/Build_ARC/mod_storage.o

def_hessian.o: cppdefs.h globaldefs.h arctic.h
def_hessian.o: /center/w/kate/Build_ARC/def_var.o
def_hessian.o: /center/w/kate/Build_ARC/mod_biology.o
def_hessian.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_hessian.o: /center/w/kate/Build_ARC/mod_iounits.o
def_hessian.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_hessian.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_hessian.o: /center/w/kate/Build_ARC/mod_parallel.o
def_hessian.o: /center/w/kate/Build_ARC/mod_param.o
def_hessian.o: /center/w/kate/Build_ARC/mod_scalars.o
def_hessian.o: /center/w/kate/Build_ARC/mod_sediment.o

def_his.o: cppdefs.h globaldefs.h arctic.h
def_his.o: /center/w/kate/Build_ARC/def_var.o
def_his.o: /center/w/kate/Build_ARC/mod_biology.o
def_his.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_his.o: /center/w/kate/Build_ARC/mod_iounits.o
def_his.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_his.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_his.o: /center/w/kate/Build_ARC/mod_parallel.o
def_his.o: /center/w/kate/Build_ARC/mod_param.o
def_his.o: /center/w/kate/Build_ARC/mod_scalars.o
def_his.o: /center/w/kate/Build_ARC/mod_sediment.o

def_his2.o: cppdefs.h globaldefs.h arctic.h
def_his2.o: /center/w/kate/Build_ARC/def_var.o
def_his2.o: /center/w/kate/Build_ARC/distribute.o
def_his2.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_his2.o: /center/w/kate/Build_ARC/mod_iounits.o
def_his2.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_his2.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_his2.o: /center/w/kate/Build_ARC/mod_parallel.o
def_his2.o: /center/w/kate/Build_ARC/mod_param.o
def_his2.o: /center/w/kate/Build_ARC/mod_scalars.o
def_his2.o: /center/w/kate/Build_ARC/mod_sediment.o

def_impulse.o: cppdefs.h globaldefs.h arctic.h
def_impulse.o: /center/w/kate/Build_ARC/def_var.o
def_impulse.o: /center/w/kate/Build_ARC/mod_biology.o
def_impulse.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_impulse.o: /center/w/kate/Build_ARC/mod_iounits.o
def_impulse.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_impulse.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_impulse.o: /center/w/kate/Build_ARC/mod_parallel.o
def_impulse.o: /center/w/kate/Build_ARC/mod_param.o
def_impulse.o: /center/w/kate/Build_ARC/mod_scalars.o
def_impulse.o: /center/w/kate/Build_ARC/mod_sediment.o

def_info.o: sediment_def.h npzd_Franks_def.h npzd_iron_def.h fennel_def.h
def_info.o: cppdefs.h globaldefs.h arctic.h oyster_floats_def.h ecosim_def.h
def_info.o: npzd_Powell_def.h umaine_def.h nemuro_def.h
def_info.o: /center/w/kate/Build_ARC/def_var.o
def_info.o: /center/w/kate/Build_ARC/distribute.o
def_info.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_info.o: /center/w/kate/Build_ARC/mod_grid.o
def_info.o: /center/w/kate/Build_ARC/mod_iounits.o
def_info.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_info.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_info.o: /center/w/kate/Build_ARC/mod_parallel.o
def_info.o: /center/w/kate/Build_ARC/mod_param.o
def_info.o: /center/w/kate/Build_ARC/mod_scalars.o
def_info.o: /center/w/kate/Build_ARC/mod_strings.o
def_info.o: /center/w/kate/Build_ARC/strings.o

def_ini.o: cppdefs.h globaldefs.h arctic.h
def_ini.o: /center/w/kate/Build_ARC/def_var.o
def_ini.o: /center/w/kate/Build_ARC/mod_iounits.o
def_ini.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_ini.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_ini.o: /center/w/kate/Build_ARC/mod_parallel.o
def_ini.o: /center/w/kate/Build_ARC/mod_param.o
def_ini.o: /center/w/kate/Build_ARC/mod_scalars.o

def_lanczos.o: cppdefs.h globaldefs.h arctic.h
def_lanczos.o: /center/w/kate/Build_ARC/def_var.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_biology.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_iounits.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_parallel.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_param.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_scalars.o
def_lanczos.o: /center/w/kate/Build_ARC/mod_sediment.o

def_mod.o: cppdefs.h globaldefs.h arctic.h
def_mod.o: /center/w/kate/Build_ARC/def_var.o
def_mod.o: /center/w/kate/Build_ARC/distribute.o
def_mod.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_mod.o: /center/w/kate/Build_ARC/mod_iounits.o
def_mod.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_mod.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_mod.o: /center/w/kate/Build_ARC/mod_parallel.o
def_mod.o: /center/w/kate/Build_ARC/mod_param.o
def_mod.o: /center/w/kate/Build_ARC/mod_scalars.o
def_mod.o: /center/w/kate/Build_ARC/mod_strings.o

def_norm.o: cppdefs.h globaldefs.h arctic.h
def_norm.o: /center/w/kate/Build_ARC/def_var.o
def_norm.o: /center/w/kate/Build_ARC/mod_biology.o
def_norm.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_norm.o: /center/w/kate/Build_ARC/mod_iounits.o
def_norm.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_norm.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_norm.o: /center/w/kate/Build_ARC/mod_parallel.o
def_norm.o: /center/w/kate/Build_ARC/mod_param.o
def_norm.o: /center/w/kate/Build_ARC/mod_scalars.o
def_norm.o: /center/w/kate/Build_ARC/mod_sediment.o

def_rst.o: cppdefs.h globaldefs.h arctic.h
def_rst.o: /center/w/kate/Build_ARC/def_var.o
def_rst.o: /center/w/kate/Build_ARC/mod_biology.o
def_rst.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_rst.o: /center/w/kate/Build_ARC/mod_iounits.o
def_rst.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_rst.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_rst.o: /center/w/kate/Build_ARC/mod_parallel.o
def_rst.o: /center/w/kate/Build_ARC/mod_param.o
def_rst.o: /center/w/kate/Build_ARC/mod_scalars.o
def_rst.o: /center/w/kate/Build_ARC/mod_sediment.o

def_station.o: cppdefs.h globaldefs.h arctic.h
def_station.o: /center/w/kate/Build_ARC/def_var.o
def_station.o: /center/w/kate/Build_ARC/mod_biology.o
def_station.o: /center/w/kate/Build_ARC/mod_fourdvar.o
def_station.o: /center/w/kate/Build_ARC/mod_iounits.o
def_station.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_station.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_station.o: /center/w/kate/Build_ARC/mod_parallel.o
def_station.o: /center/w/kate/Build_ARC/mod_param.o
def_station.o: /center/w/kate/Build_ARC/mod_scalars.o
def_station.o: /center/w/kate/Build_ARC/mod_sediment.o

def_tides.o: cppdefs.h globaldefs.h arctic.h
def_tides.o: /center/w/kate/Build_ARC/def_var.o
def_tides.o: /center/w/kate/Build_ARC/mod_iounits.o
def_tides.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_tides.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_tides.o: /center/w/kate/Build_ARC/mod_parallel.o
def_tides.o: /center/w/kate/Build_ARC/mod_param.o
def_tides.o: /center/w/kate/Build_ARC/mod_scalars.o
def_tides.o: /center/w/kate/Build_ARC/mod_stepping.o
def_tides.o: /center/w/kate/Build_ARC/mod_tides.o

def_var.o: cppdefs.h globaldefs.h arctic.h
def_var.o: /center/w/kate/Build_ARC/distribute.o
def_var.o: /center/w/kate/Build_ARC/mod_iounits.o
def_var.o: /center/w/kate/Build_ARC/mod_ncparam.o
def_var.o: /center/w/kate/Build_ARC/mod_netcdf.o
def_var.o: /center/w/kate/Build_ARC/mod_parallel.o
def_var.o: /center/w/kate/Build_ARC/mod_param.o
def_var.o: /center/w/kate/Build_ARC/mod_scalars.o

distribute.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
distribute.o: /center/w/kate/Build_ARC/mod_iounits.o
distribute.o: /center/w/kate/Build_ARC/mod_ncparam.o
distribute.o: /center/w/kate/Build_ARC/mod_parallel.o
distribute.o: /center/w/kate/Build_ARC/mod_param.o
distribute.o: /center/w/kate/Build_ARC/mod_scalars.o

dotproduct.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
dotproduct.o: /center/w/kate/Build_ARC/distribute.o
dotproduct.o: /center/w/kate/Build_ARC/mod_fourdvar.o
dotproduct.o: /center/w/kate/Build_ARC/mod_grid.o
dotproduct.o: /center/w/kate/Build_ARC/mod_iounits.o
dotproduct.o: /center/w/kate/Build_ARC/mod_ncparam.o
dotproduct.o: /center/w/kate/Build_ARC/mod_ocean.o
dotproduct.o: /center/w/kate/Build_ARC/mod_parallel.o
dotproduct.o: /center/w/kate/Build_ARC/mod_param.o
dotproduct.o: /center/w/kate/Build_ARC/mod_scalars.o
dotproduct.o: /center/w/kate/Build_ARC/mod_stepping.o

erf.o: cppdefs.h globaldefs.h arctic.h
erf.o: /center/w/kate/Build_ARC/mod_iounits.o
erf.o: /center/w/kate/Build_ARC/mod_kinds.o
erf.o: /center/w/kate/Build_ARC/mod_param.o
erf.o: /center/w/kate/Build_ARC/mod_scalars.o

extract_obs.o: cppdefs.h globaldefs.h arctic.h
extract_obs.o: /center/w/kate/Build_ARC/mod_kinds.o
extract_obs.o: /center/w/kate/Build_ARC/mod_param.o

extract_sta.o: cppdefs.h globaldefs.h arctic.h
extract_sta.o: /center/w/kate/Build_ARC/distribute.o
extract_sta.o: /center/w/kate/Build_ARC/mod_grid.o
extract_sta.o: /center/w/kate/Build_ARC/mod_ncparam.o
extract_sta.o: /center/w/kate/Build_ARC/mod_parallel.o
extract_sta.o: /center/w/kate/Build_ARC/mod_param.o
extract_sta.o: /center/w/kate/Build_ARC/mod_scalars.o

frc_weak.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
frc_weak.o: /center/w/kate/Build_ARC/mod_forces.o
frc_weak.o: /center/w/kate/Build_ARC/mod_fourdvar.o
frc_weak.o: /center/w/kate/Build_ARC/mod_ocean.o
frc_weak.o: /center/w/kate/Build_ARC/mod_param.o
frc_weak.o: /center/w/kate/Build_ARC/mod_scalars.o
frc_weak.o: /center/w/kate/Build_ARC/mod_stepping.o

gasdev.o: cppdefs.h globaldefs.h arctic.h
gasdev.o: /center/w/kate/Build_ARC/mod_kinds.o
gasdev.o: /center/w/kate/Build_ARC/nrutil.o

get_2dfld.o: cppdefs.h globaldefs.h arctic.h
get_2dfld.o: /center/w/kate/Build_ARC/mod_iounits.o
get_2dfld.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_2dfld.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_2dfld.o: /center/w/kate/Build_ARC/mod_parallel.o
get_2dfld.o: /center/w/kate/Build_ARC/mod_param.o
get_2dfld.o: /center/w/kate/Build_ARC/mod_scalars.o
get_2dfld.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_2dfld.o: /center/w/kate/Build_ARC/nf_fread3d.o

get_2dfldr.o: cppdefs.h globaldefs.h arctic.h
get_2dfldr.o: /center/w/kate/Build_ARC/mod_iounits.o
get_2dfldr.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_2dfldr.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_2dfldr.o: /center/w/kate/Build_ARC/mod_parallel.o
get_2dfldr.o: /center/w/kate/Build_ARC/mod_param.o
get_2dfldr.o: /center/w/kate/Build_ARC/mod_scalars.o
get_2dfldr.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_2dfldr.o: /center/w/kate/Build_ARC/nf_fread3d.o

get_3dfld.o: cppdefs.h globaldefs.h arctic.h
get_3dfld.o: /center/w/kate/Build_ARC/mod_iounits.o
get_3dfld.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_3dfld.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_3dfld.o: /center/w/kate/Build_ARC/mod_parallel.o
get_3dfld.o: /center/w/kate/Build_ARC/mod_param.o
get_3dfld.o: /center/w/kate/Build_ARC/mod_scalars.o
get_3dfld.o: /center/w/kate/Build_ARC/nf_fread3d.o

get_3dfldr.o: cppdefs.h globaldefs.h arctic.h
get_3dfldr.o: /center/w/kate/Build_ARC/mod_iounits.o
get_3dfldr.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_3dfldr.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_3dfldr.o: /center/w/kate/Build_ARC/mod_parallel.o
get_3dfldr.o: /center/w/kate/Build_ARC/mod_param.o
get_3dfldr.o: /center/w/kate/Build_ARC/mod_scalars.o
get_3dfldr.o: /center/w/kate/Build_ARC/nf_fread3d.o

get_bounds.o: cppdefs.h globaldefs.h arctic.h
get_bounds.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_bounds.o: /center/w/kate/Build_ARC/mod_nesting.o
get_bounds.o: /center/w/kate/Build_ARC/mod_parallel.o
get_bounds.o: /center/w/kate/Build_ARC/mod_param.o
get_bounds.o: /center/w/kate/Build_ARC/mod_scalars.o

get_cycle.o: cppdefs.h globaldefs.h arctic.h
get_cycle.o: /center/w/kate/Build_ARC/mod_iounits.o
get_cycle.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_cycle.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_cycle.o: /center/w/kate/Build_ARC/mod_parallel.o
get_cycle.o: /center/w/kate/Build_ARC/mod_param.o
get_cycle.o: /center/w/kate/Build_ARC/mod_scalars.o

get_date.o: cppdefs.h globaldefs.h arctic.h
get_date.o: /center/w/kate/Build_ARC/mod_kinds.o

get_filter.o: cppdefs.h globaldefs.h arctic.h
get_filter.o: /center/w/kate/Build_ARC/mod_filter.o
get_filter.o: /center/w/kate/Build_ARC/mod_grid.o
get_filter.o: /center/w/kate/Build_ARC/mod_iounits.o
get_filter.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_filter.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_filter.o: /center/w/kate/Build_ARC/mod_ocean.o
get_filter.o: /center/w/kate/Build_ARC/mod_parallel.o
get_filter.o: /center/w/kate/Build_ARC/mod_param.o
get_filter.o: /center/w/kate/Build_ARC/mod_scalars.o
get_filter.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_filter.o: /center/w/kate/Build_ARC/nf_fread3d.o

get_grid.o: cppdefs.h globaldefs.h arctic.h
get_grid.o: /center/w/kate/Build_ARC/exchange_2d.o
get_grid.o: /center/w/kate/Build_ARC/mod_grid.o
get_grid.o: /center/w/kate/Build_ARC/mod_iounits.o
get_grid.o: /center/w/kate/Build_ARC/mod_mixing.o
get_grid.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_grid.o: /center/w/kate/Build_ARC/mod_nesting.o
get_grid.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_grid.o: /center/w/kate/Build_ARC/mod_parallel.o
get_grid.o: /center/w/kate/Build_ARC/mod_param.o
get_grid.o: /center/w/kate/Build_ARC/mod_scalars.o
get_grid.o: /center/w/kate/Build_ARC/mp_exchange.o
get_grid.o: /center/w/kate/Build_ARC/nesting.o
get_grid.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_grid.o: /center/w/kate/Build_ARC/strings.o

get_gst.o: cppdefs.h globaldefs.h arctic.h
get_gst.o: /center/w/kate/Build_ARC/distribute.o
get_gst.o: /center/w/kate/Build_ARC/mod_iounits.o
get_gst.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_gst.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_gst.o: /center/w/kate/Build_ARC/mod_parallel.o
get_gst.o: /center/w/kate/Build_ARC/mod_param.o
get_gst.o: /center/w/kate/Build_ARC/mod_scalars.o
get_gst.o: /center/w/kate/Build_ARC/mod_storage.o

get_ngfld.o: cppdefs.h globaldefs.h arctic.h
get_ngfld.o: /center/w/kate/Build_ARC/mod_iounits.o
get_ngfld.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_ngfld.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_ngfld.o: /center/w/kate/Build_ARC/mod_parallel.o
get_ngfld.o: /center/w/kate/Build_ARC/mod_param.o
get_ngfld.o: /center/w/kate/Build_ARC/mod_scalars.o

get_ngfldr.o: cppdefs.h globaldefs.h arctic.h
get_ngfldr.o: /center/w/kate/Build_ARC/mod_iounits.o
get_ngfldr.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_ngfldr.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_ngfldr.o: /center/w/kate/Build_ARC/mod_parallel.o
get_ngfldr.o: /center/w/kate/Build_ARC/mod_param.o
get_ngfldr.o: /center/w/kate/Build_ARC/mod_scalars.o

get_nudgcoef.o: cppdefs.h globaldefs.h arctic.h
get_nudgcoef.o: /center/w/kate/Build_ARC/exchange_2d.o
get_nudgcoef.o: /center/w/kate/Build_ARC/exchange_3d.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_clima.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_grid.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_iounits.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_parallel.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_param.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mod_scalars.o
get_nudgcoef.o: /center/w/kate/Build_ARC/mp_exchange.o
get_nudgcoef.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_nudgcoef.o: /center/w/kate/Build_ARC/nf_fread3d.o
get_nudgcoef.o: /center/w/kate/Build_ARC/strings.o

get_state.o: cppdefs.h globaldefs.h arctic.h
get_state.o: /center/w/kate/Build_ARC/mod_biology.o
get_state.o: /center/w/kate/Build_ARC/mod_boundary.o
get_state.o: /center/w/kate/Build_ARC/mod_filter.o
get_state.o: /center/w/kate/Build_ARC/mod_forces.o
get_state.o: /center/w/kate/Build_ARC/mod_fourdvar.o
get_state.o: /center/w/kate/Build_ARC/mod_grid.o
get_state.o: /center/w/kate/Build_ARC/mod_ice.o
get_state.o: /center/w/kate/Build_ARC/mod_iounits.o
get_state.o: /center/w/kate/Build_ARC/mod_mixing.o
get_state.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_state.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_state.o: /center/w/kate/Build_ARC/mod_ocean.o
get_state.o: /center/w/kate/Build_ARC/mod_parallel.o
get_state.o: /center/w/kate/Build_ARC/mod_param.o
get_state.o: /center/w/kate/Build_ARC/mod_scalars.o
get_state.o: /center/w/kate/Build_ARC/mod_sedbed.o
get_state.o: /center/w/kate/Build_ARC/mod_sediment.o
get_state.o: /center/w/kate/Build_ARC/mod_stepping.o
get_state.o: /center/w/kate/Build_ARC/mod_strings.o
get_state.o: /center/w/kate/Build_ARC/mp_exchange.o
get_state.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_state.o: /center/w/kate/Build_ARC/nf_fread2d_bry.o
get_state.o: /center/w/kate/Build_ARC/nf_fread3d.o
get_state.o: /center/w/kate/Build_ARC/nf_fread3d_bry.o
get_state.o: /center/w/kate/Build_ARC/nf_fread4d.o
get_state.o: /center/w/kate/Build_ARC/strings.o

get_varcoords.o: cppdefs.h globaldefs.h arctic.h
get_varcoords.o: /center/w/kate/Build_ARC/mod_grid.o
get_varcoords.o: /center/w/kate/Build_ARC/mod_iounits.o
get_varcoords.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_varcoords.o: /center/w/kate/Build_ARC/mod_parallel.o
get_varcoords.o: /center/w/kate/Build_ARC/mod_param.o
get_varcoords.o: /center/w/kate/Build_ARC/mod_scalars.o

get_wetdry.o: cppdefs.h globaldefs.h arctic.h
get_wetdry.o: /center/w/kate/Build_ARC/exchange_2d.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_grid.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_iounits.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_ncparam.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_netcdf.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_parallel.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_param.o
get_wetdry.o: /center/w/kate/Build_ARC/mod_scalars.o
get_wetdry.o: /center/w/kate/Build_ARC/mp_exchange.o
get_wetdry.o: /center/w/kate/Build_ARC/nf_fread2d.o
get_wetdry.o: /center/w/kate/Build_ARC/strings.o

grid_coords.o: cppdefs.h globaldefs.h arctic.h
grid_coords.o: /center/w/kate/Build_ARC/distribute.o
grid_coords.o: /center/w/kate/Build_ARC/interpolate.o
grid_coords.o: /center/w/kate/Build_ARC/mod_floats.o
grid_coords.o: /center/w/kate/Build_ARC/mod_grid.o
grid_coords.o: /center/w/kate/Build_ARC/mod_parallel.o
grid_coords.o: /center/w/kate/Build_ARC/mod_param.o
grid_coords.o: /center/w/kate/Build_ARC/mod_scalars.o

ini_adjust.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
ini_adjust.o: /center/w/kate/Build_ARC/exchange_2d.o
ini_adjust.o: /center/w/kate/Build_ARC/exchange_3d.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_boundary.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_coupling.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_forces.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_fourdvar.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_grid.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_iounits.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_ncparam.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_ocean.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_parallel.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_param.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_scalars.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_sedbed.o
ini_adjust.o: /center/w/kate/Build_ARC/mod_stepping.o
ini_adjust.o: /center/w/kate/Build_ARC/mp_exchange.o
ini_adjust.o: /center/w/kate/Build_ARC/set_depth.o
ini_adjust.o: /center/w/kate/Build_ARC/state_addition.o
ini_adjust.o: /center/w/kate/Build_ARC/state_copy.o
ini_adjust.o: /center/w/kate/Build_ARC/t3dbc_im.o
ini_adjust.o: /center/w/kate/Build_ARC/u2dbc_im.o
ini_adjust.o: /center/w/kate/Build_ARC/u3dbc_im.o
ini_adjust.o: /center/w/kate/Build_ARC/v2dbc_im.o
ini_adjust.o: /center/w/kate/Build_ARC/v3dbc_im.o
ini_adjust.o: /center/w/kate/Build_ARC/zetabc.o

ini_hmixcoef.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
ini_hmixcoef.o: /center/w/kate/Build_ARC/exchange_2d.o
ini_hmixcoef.o: /center/w/kate/Build_ARC/mod_grid.o
ini_hmixcoef.o: /center/w/kate/Build_ARC/mod_mixing.o
ini_hmixcoef.o: /center/w/kate/Build_ARC/mod_ncparam.o
ini_hmixcoef.o: /center/w/kate/Build_ARC/mod_param.o
ini_hmixcoef.o: /center/w/kate/Build_ARC/mod_scalars.o
ini_hmixcoef.o: /center/w/kate/Build_ARC/mp_exchange.o

ini_lanczos.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
ini_lanczos.o: /center/w/kate/Build_ARC/distribute.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_boundary.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_forces.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_fourdvar.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_grid.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_iounits.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_ncparam.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_netcdf.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_ocean.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_parallel.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_param.o
ini_lanczos.o: /center/w/kate/Build_ARC/mod_scalars.o
ini_lanczos.o: /center/w/kate/Build_ARC/nf_fread2d.o
ini_lanczos.o: /center/w/kate/Build_ARC/nf_fread2d_bry.o
ini_lanczos.o: /center/w/kate/Build_ARC/nf_fread3d.o
ini_lanczos.o: /center/w/kate/Build_ARC/nf_fread3d_bry.o
ini_lanczos.o: /center/w/kate/Build_ARC/state_addition.o
ini_lanczos.o: /center/w/kate/Build_ARC/state_dotprod.o
ini_lanczos.o: /center/w/kate/Build_ARC/state_initialize.o
ini_lanczos.o: /center/w/kate/Build_ARC/state_scale.o

inner2state.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
inner2state.o: /center/w/kate/Build_ARC/distribute.o
inner2state.o: /center/w/kate/Build_ARC/mod_boundary.o
inner2state.o: /center/w/kate/Build_ARC/mod_forces.o
inner2state.o: /center/w/kate/Build_ARC/mod_fourdvar.o
inner2state.o: /center/w/kate/Build_ARC/mod_grid.o
inner2state.o: /center/w/kate/Build_ARC/mod_iounits.o
inner2state.o: /center/w/kate/Build_ARC/mod_ncparam.o
inner2state.o: /center/w/kate/Build_ARC/mod_netcdf.o
inner2state.o: /center/w/kate/Build_ARC/mod_ocean.o
inner2state.o: /center/w/kate/Build_ARC/mod_parallel.o
inner2state.o: /center/w/kate/Build_ARC/mod_param.o
inner2state.o: /center/w/kate/Build_ARC/mod_scalars.o
inner2state.o: /center/w/kate/Build_ARC/mod_stepping.o
inner2state.o: /center/w/kate/Build_ARC/nf_fread2d.o
inner2state.o: /center/w/kate/Build_ARC/nf_fread2d_bry.o
inner2state.o: /center/w/kate/Build_ARC/nf_fread3d.o
inner2state.o: /center/w/kate/Build_ARC/nf_fread3d_bry.o
inner2state.o: /center/w/kate/Build_ARC/state_addition.o
inner2state.o: /center/w/kate/Build_ARC/state_dotprod.o
inner2state.o: /center/w/kate/Build_ARC/state_initialize.o
inner2state.o: /center/w/kate/Build_ARC/state_scale.o

inp_par.o: cppdefs.h globaldefs.h arctic.h
inp_par.o: /center/w/kate/Build_ARC/distribute.o
inp_par.o: /center/w/kate/Build_ARC/mod_filter.o
inp_par.o: /center/w/kate/Build_ARC/mod_iounits.o
inp_par.o: /center/w/kate/Build_ARC/mod_kinds.o
inp_par.o: /center/w/kate/Build_ARC/mod_ncparam.o
inp_par.o: /center/w/kate/Build_ARC/mod_netcdf.o
inp_par.o: /center/w/kate/Build_ARC/mod_parallel.o
inp_par.o: /center/w/kate/Build_ARC/mod_param.o
inp_par.o: /center/w/kate/Build_ARC/mod_scalars.o
inp_par.o: /center/w/kate/Build_ARC/mod_strings.o
inp_par.o: /center/w/kate/Build_ARC/mod_tides.o
inp_par.o: /center/w/kate/Build_ARC/ran_state.o
inp_par.o: /center/w/kate/Build_ARC/strings.o

inquire.o: cppdefs.h globaldefs.h arctic.h
inquire.o: /center/w/kate/Build_ARC/mod_iounits.o
inquire.o: /center/w/kate/Build_ARC/mod_ncparam.o
inquire.o: /center/w/kate/Build_ARC/mod_netcdf.o
inquire.o: /center/w/kate/Build_ARC/mod_parallel.o
inquire.o: /center/w/kate/Build_ARC/mod_param.o
inquire.o: /center/w/kate/Build_ARC/mod_scalars.o

interpolate.o: cppdefs.h globaldefs.h arctic.h
interpolate.o: /center/w/kate/Build_ARC/mod_kinds.o
interpolate.o: /center/w/kate/Build_ARC/mod_param.o
interpolate.o: /center/w/kate/Build_ARC/mod_scalars.o

lbc.o: cppdefs.h globaldefs.h arctic.h
lbc.o: /center/w/kate/Build_ARC/distribute.o
lbc.o: /center/w/kate/Build_ARC/mod_iounits.o
lbc.o: /center/w/kate/Build_ARC/mod_ncparam.o
lbc.o: /center/w/kate/Build_ARC/mod_netcdf.o
lbc.o: /center/w/kate/Build_ARC/mod_parallel.o
lbc.o: /center/w/kate/Build_ARC/mod_param.o
lbc.o: /center/w/kate/Build_ARC/mod_scalars.o

lubksb.o: cppdefs.h globaldefs.h arctic.h
lubksb.o: /center/w/kate/Build_ARC/mod_kinds.o

ludcmp.o: cppdefs.h globaldefs.h arctic.h
ludcmp.o: /center/w/kate/Build_ARC/mod_kinds.o

metrics.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
metrics.o: /center/w/kate/Build_ARC/distribute.o
metrics.o: /center/w/kate/Build_ARC/exchange_2d.o
metrics.o: /center/w/kate/Build_ARC/mod_fourdvar.o
metrics.o: /center/w/kate/Build_ARC/mod_grid.o
metrics.o: /center/w/kate/Build_ARC/mod_iounits.o
metrics.o: /center/w/kate/Build_ARC/mod_mixing.o
metrics.o: /center/w/kate/Build_ARC/mod_ncparam.o
metrics.o: /center/w/kate/Build_ARC/mod_nesting.o
metrics.o: /center/w/kate/Build_ARC/mod_ocean.o
metrics.o: /center/w/kate/Build_ARC/mod_parallel.o
metrics.o: /center/w/kate/Build_ARC/mod_param.o
metrics.o: /center/w/kate/Build_ARC/mod_scalars.o
metrics.o: /center/w/kate/Build_ARC/mod_sedbed.o
metrics.o: /center/w/kate/Build_ARC/mod_stepping.o
metrics.o: /center/w/kate/Build_ARC/mp_exchange.o
metrics.o: /center/w/kate/Build_ARC/set_depth.o

mp_exchange.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
mp_exchange.o: /center/w/kate/Build_ARC/mod_iounits.o
mp_exchange.o: /center/w/kate/Build_ARC/mod_parallel.o
mp_exchange.o: /center/w/kate/Build_ARC/mod_param.o
mp_exchange.o: /center/w/kate/Build_ARC/mod_scalars.o

mp_routines.o: cppdefs.h globaldefs.h arctic.h
mp_routines.o: /center/w/kate/Build_ARC/mod_kinds.o

nf_fread2d.o: cppdefs.h globaldefs.h arctic.h
nf_fread2d.o: /center/w/kate/Build_ARC/distribute.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_grid.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_iounits.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_param.o
nf_fread2d.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fread2d_bry.o: cppdefs.h globaldefs.h arctic.h
nf_fread2d_bry.o: /center/w/kate/Build_ARC/distribute.o
nf_fread2d_bry.o: /center/w/kate/Build_ARC/mod_iounits.o
nf_fread2d_bry.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fread2d_bry.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fread2d_bry.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fread2d_bry.o: /center/w/kate/Build_ARC/mod_param.o
nf_fread2d_bry.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fread3d.o: cppdefs.h globaldefs.h arctic.h
nf_fread3d.o: /center/w/kate/Build_ARC/distribute.o
nf_fread3d.o: /center/w/kate/Build_ARC/mod_iounits.o
nf_fread3d.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fread3d.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fread3d.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fread3d.o: /center/w/kate/Build_ARC/mod_param.o
nf_fread3d.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fread3d_bry.o: cppdefs.h globaldefs.h arctic.h
nf_fread3d_bry.o: /center/w/kate/Build_ARC/distribute.o
nf_fread3d_bry.o: /center/w/kate/Build_ARC/mod_iounits.o
nf_fread3d_bry.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fread3d_bry.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fread3d_bry.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fread3d_bry.o: /center/w/kate/Build_ARC/mod_param.o
nf_fread3d_bry.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fread4d.o: cppdefs.h globaldefs.h arctic.h
nf_fread4d.o: /center/w/kate/Build_ARC/distribute.o
nf_fread4d.o: /center/w/kate/Build_ARC/mod_iounits.o
nf_fread4d.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fread4d.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fread4d.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fread4d.o: /center/w/kate/Build_ARC/mod_param.o
nf_fread4d.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fwrite2d.o: cppdefs.h globaldefs.h arctic.h
nf_fwrite2d.o: /center/w/kate/Build_ARC/distribute.o
nf_fwrite2d.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fwrite2d.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fwrite2d.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fwrite2d.o: /center/w/kate/Build_ARC/mod_param.o
nf_fwrite2d.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fwrite2d_bry.o: cppdefs.h globaldefs.h arctic.h
nf_fwrite2d_bry.o: /center/w/kate/Build_ARC/distribute.o
nf_fwrite2d_bry.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fwrite2d_bry.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fwrite2d_bry.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fwrite2d_bry.o: /center/w/kate/Build_ARC/mod_param.o
nf_fwrite2d_bry.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fwrite3d.o: cppdefs.h globaldefs.h arctic.h
nf_fwrite3d.o: /center/w/kate/Build_ARC/distribute.o
nf_fwrite3d.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fwrite3d.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fwrite3d.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fwrite3d.o: /center/w/kate/Build_ARC/mod_param.o
nf_fwrite3d.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fwrite3d_bry.o: cppdefs.h globaldefs.h arctic.h
nf_fwrite3d_bry.o: /center/w/kate/Build_ARC/distribute.o
nf_fwrite3d_bry.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fwrite3d_bry.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fwrite3d_bry.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fwrite3d_bry.o: /center/w/kate/Build_ARC/mod_param.o
nf_fwrite3d_bry.o: /center/w/kate/Build_ARC/mod_scalars.o

nf_fwrite4d.o: cppdefs.h globaldefs.h arctic.h
nf_fwrite4d.o: /center/w/kate/Build_ARC/distribute.o
nf_fwrite4d.o: /center/w/kate/Build_ARC/mod_ncparam.o
nf_fwrite4d.o: /center/w/kate/Build_ARC/mod_netcdf.o
nf_fwrite4d.o: /center/w/kate/Build_ARC/mod_parallel.o
nf_fwrite4d.o: /center/w/kate/Build_ARC/mod_param.o
nf_fwrite4d.o: /center/w/kate/Build_ARC/mod_scalars.o

normalization.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
normalization.o: /center/w/kate/Build_ARC/bc_2d.o
normalization.o: /center/w/kate/Build_ARC/bc_3d.o
normalization.o: /center/w/kate/Build_ARC/bc_bry2d.o
normalization.o: /center/w/kate/Build_ARC/bc_bry3d.o
normalization.o: /center/w/kate/Build_ARC/distribute.o
normalization.o: /center/w/kate/Build_ARC/mod_boundary.o
normalization.o: /center/w/kate/Build_ARC/mod_forces.o
normalization.o: /center/w/kate/Build_ARC/mod_fourdvar.o
normalization.o: /center/w/kate/Build_ARC/mod_grid.o
normalization.o: /center/w/kate/Build_ARC/mod_iounits.o
normalization.o: /center/w/kate/Build_ARC/mod_kinds.o
normalization.o: /center/w/kate/Build_ARC/mod_mixing.o
normalization.o: /center/w/kate/Build_ARC/mod_ncparam.o
normalization.o: /center/w/kate/Build_ARC/mod_netcdf.o
normalization.o: /center/w/kate/Build_ARC/mod_ocean.o
normalization.o: /center/w/kate/Build_ARC/mod_parallel.o
normalization.o: /center/w/kate/Build_ARC/mod_param.o
normalization.o: /center/w/kate/Build_ARC/mod_scalars.o
normalization.o: /center/w/kate/Build_ARC/mod_sedbed.o
normalization.o: /center/w/kate/Build_ARC/mod_stepping.o
normalization.o: /center/w/kate/Build_ARC/mp_exchange.o
normalization.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
normalization.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
normalization.o: /center/w/kate/Build_ARC/set_depth.o
normalization.o: /center/w/kate/Build_ARC/white_noise.o

nrutil.o: /center/w/kate/Build_ARC/mod_kinds.o

ntimestep.o: cppdefs.h globaldefs.h arctic.h
ntimestep.o: /center/w/kate/Build_ARC/mod_iounits.o
ntimestep.o: /center/w/kate/Build_ARC/mod_parallel.o
ntimestep.o: /center/w/kate/Build_ARC/mod_param.o
ntimestep.o: /center/w/kate/Build_ARC/mod_scalars.o

obs_cost.o: cppdefs.h globaldefs.h arctic.h
obs_cost.o: /center/w/kate/Build_ARC/mod_fourdvar.o
obs_cost.o: /center/w/kate/Build_ARC/mod_parallel.o
obs_cost.o: /center/w/kate/Build_ARC/mod_param.o
obs_cost.o: /center/w/kate/Build_ARC/mod_scalars.o

obs_depth.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
obs_depth.o: /center/w/kate/Build_ARC/mod_fourdvar.o
obs_depth.o: /center/w/kate/Build_ARC/mod_ncparam.o
obs_depth.o: /center/w/kate/Build_ARC/mod_param.o
obs_depth.o: /center/w/kate/Build_ARC/mod_scalars.o

obs_initial.o: cppdefs.h globaldefs.h arctic.h
obs_initial.o: /center/w/kate/Build_ARC/mod_fourdvar.o
obs_initial.o: /center/w/kate/Build_ARC/mod_iounits.o
obs_initial.o: /center/w/kate/Build_ARC/mod_ncparam.o
obs_initial.o: /center/w/kate/Build_ARC/mod_netcdf.o
obs_initial.o: /center/w/kate/Build_ARC/mod_parallel.o
obs_initial.o: /center/w/kate/Build_ARC/mod_param.o
obs_initial.o: /center/w/kate/Build_ARC/mod_scalars.o

obs_read.o: cppdefs.h globaldefs.h arctic.h
obs_read.o: /center/w/kate/Build_ARC/mod_fourdvar.o
obs_read.o: /center/w/kate/Build_ARC/mod_iounits.o
obs_read.o: /center/w/kate/Build_ARC/mod_ncparam.o
obs_read.o: /center/w/kate/Build_ARC/mod_netcdf.o
obs_read.o: /center/w/kate/Build_ARC/mod_parallel.o
obs_read.o: /center/w/kate/Build_ARC/mod_param.o
obs_read.o: /center/w/kate/Build_ARC/mod_scalars.o

obs_write.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
obs_write.o: /center/w/kate/Build_ARC/distribute.o
obs_write.o: /center/w/kate/Build_ARC/extract_obs.o
obs_write.o: /center/w/kate/Build_ARC/mod_fourdvar.o
obs_write.o: /center/w/kate/Build_ARC/mod_grid.o
obs_write.o: /center/w/kate/Build_ARC/mod_iounits.o
obs_write.o: /center/w/kate/Build_ARC/mod_ncparam.o
obs_write.o: /center/w/kate/Build_ARC/mod_netcdf.o
obs_write.o: /center/w/kate/Build_ARC/mod_ocean.o
obs_write.o: /center/w/kate/Build_ARC/mod_parallel.o
obs_write.o: /center/w/kate/Build_ARC/mod_param.o
obs_write.o: /center/w/kate/Build_ARC/mod_scalars.o
obs_write.o: /center/w/kate/Build_ARC/mod_stepping.o
obs_write.o: /center/w/kate/Build_ARC/nf_fwrite2d.o

packing.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
packing.o: /center/w/kate/Build_ARC/distribute.o
packing.o: /center/w/kate/Build_ARC/mod_forces.o
packing.o: /center/w/kate/Build_ARC/mod_grid.o
packing.o: /center/w/kate/Build_ARC/mod_iounits.o
packing.o: /center/w/kate/Build_ARC/mod_ncparam.o
packing.o: /center/w/kate/Build_ARC/mod_netcdf.o
packing.o: /center/w/kate/Build_ARC/mod_ocean.o
packing.o: /center/w/kate/Build_ARC/mod_parallel.o
packing.o: /center/w/kate/Build_ARC/mod_param.o
packing.o: /center/w/kate/Build_ARC/mod_scalars.o
packing.o: /center/w/kate/Build_ARC/mod_stepping.o
packing.o: /center/w/kate/Build_ARC/mod_storage.o
packing.o: /center/w/kate/Build_ARC/nf_fread2d.o
packing.o: /center/w/kate/Build_ARC/nf_fread3d.o

posterior.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
posterior.o: /center/w/kate/Build_ARC/distribute.o
posterior.o: /center/w/kate/Build_ARC/mod_boundary.o
posterior.o: /center/w/kate/Build_ARC/mod_coupling.o
posterior.o: /center/w/kate/Build_ARC/mod_forces.o
posterior.o: /center/w/kate/Build_ARC/mod_fourdvar.o
posterior.o: /center/w/kate/Build_ARC/mod_grid.o
posterior.o: /center/w/kate/Build_ARC/mod_iounits.o
posterior.o: /center/w/kate/Build_ARC/mod_ncparam.o
posterior.o: /center/w/kate/Build_ARC/mod_netcdf.o
posterior.o: /center/w/kate/Build_ARC/mod_ocean.o
posterior.o: /center/w/kate/Build_ARC/mod_parallel.o
posterior.o: /center/w/kate/Build_ARC/mod_param.o
posterior.o: /center/w/kate/Build_ARC/mod_scalars.o
posterior.o: /center/w/kate/Build_ARC/mod_stepping.o
posterior.o: /center/w/kate/Build_ARC/nf_fread2d.o
posterior.o: /center/w/kate/Build_ARC/nf_fread2d_bry.o
posterior.o: /center/w/kate/Build_ARC/nf_fread3d.o
posterior.o: /center/w/kate/Build_ARC/nf_fread3d_bry.o
posterior.o: /center/w/kate/Build_ARC/state_addition.o
posterior.o: /center/w/kate/Build_ARC/state_copy.o
posterior.o: /center/w/kate/Build_ARC/state_dotprod.o
posterior.o: /center/w/kate/Build_ARC/state_initialize.o
posterior.o: /center/w/kate/Build_ARC/state_scale.o

posterior_var.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
posterior_var.o: /center/w/kate/Build_ARC/distribute.o
posterior_var.o: /center/w/kate/Build_ARC/mod_boundary.o
posterior_var.o: /center/w/kate/Build_ARC/mod_coupling.o
posterior_var.o: /center/w/kate/Build_ARC/mod_forces.o
posterior_var.o: /center/w/kate/Build_ARC/mod_fourdvar.o
posterior_var.o: /center/w/kate/Build_ARC/mod_grid.o
posterior_var.o: /center/w/kate/Build_ARC/mod_iounits.o
posterior_var.o: /center/w/kate/Build_ARC/mod_ncparam.o
posterior_var.o: /center/w/kate/Build_ARC/mod_ocean.o
posterior_var.o: /center/w/kate/Build_ARC/mod_parallel.o
posterior_var.o: /center/w/kate/Build_ARC/mod_param.o
posterior_var.o: /center/w/kate/Build_ARC/mod_scalars.o
posterior_var.o: /center/w/kate/Build_ARC/mod_stepping.o
posterior_var.o: /center/w/kate/Build_ARC/posterior.o
posterior_var.o: /center/w/kate/Build_ARC/state_addition.o
posterior_var.o: /center/w/kate/Build_ARC/state_copy.o
posterior_var.o: /center/w/kate/Build_ARC/state_initialize.o
posterior_var.o: /center/w/kate/Build_ARC/state_product.o

ran1.o: cppdefs.h globaldefs.h arctic.h
ran1.o: /center/w/kate/Build_ARC/mod_kinds.o
ran1.o: /center/w/kate/Build_ARC/ran_state.o

ran_state.o: /center/w/kate/Build_ARC/mod_kinds.o
ran_state.o: /center/w/kate/Build_ARC/nrutil.o

random_ic.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
random_ic.o: /center/w/kate/Build_ARC/distribute.o
random_ic.o: /center/w/kate/Build_ARC/mod_boundary.o
random_ic.o: /center/w/kate/Build_ARC/mod_coupling.o
random_ic.o: /center/w/kate/Build_ARC/mod_forces.o
random_ic.o: /center/w/kate/Build_ARC/mod_fourdvar.o
random_ic.o: /center/w/kate/Build_ARC/mod_grid.o
random_ic.o: /center/w/kate/Build_ARC/mod_iounits.o
random_ic.o: /center/w/kate/Build_ARC/mod_ncparam.o
random_ic.o: /center/w/kate/Build_ARC/mod_netcdf.o
random_ic.o: /center/w/kate/Build_ARC/mod_ocean.o
random_ic.o: /center/w/kate/Build_ARC/mod_parallel.o
random_ic.o: /center/w/kate/Build_ARC/mod_param.o
random_ic.o: /center/w/kate/Build_ARC/mod_scalars.o
random_ic.o: /center/w/kate/Build_ARC/mod_stepping.o
random_ic.o: /center/w/kate/Build_ARC/mp_exchange.o
random_ic.o: /center/w/kate/Build_ARC/white_noise.o

read_asspar.o: cppdefs.h globaldefs.h arctic.h
read_asspar.o: /center/w/kate/Build_ARC/mod_fourdvar.o
read_asspar.o: /center/w/kate/Build_ARC/mod_iounits.o
read_asspar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_asspar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_asspar.o: /center/w/kate/Build_ARC/mod_param.o
read_asspar.o: /center/w/kate/Build_ARC/mod_scalars.o

read_biopar.o: npzd_Powell_inp.h fennel_inp.h npzd_iron_inp.h nemuro_inp.h
read_biopar.o: cppdefs.h globaldefs.h arctic.h ecosim_inp.h bestnpz_inp.h
read_biopar.o: umaine_inp.h npzd_Franks_inp.h goanpz_inp.h
read_biopar.o: /center/w/kate/Build_ARC/mod_biology.o
read_biopar.o: /center/w/kate/Build_ARC/mod_eclight.o
read_biopar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_biopar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_biopar.o: /center/w/kate/Build_ARC/mod_param.o
read_biopar.o: /center/w/kate/Build_ARC/mod_scalars.o

read_couplepar.o: cppdefs.h globaldefs.h arctic.h
read_couplepar.o: /center/w/kate/Build_ARC/distribute.o
read_couplepar.o: /center/w/kate/Build_ARC/mod_coupler.o
read_couplepar.o: /center/w/kate/Build_ARC/mod_iounits.o
read_couplepar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_couplepar.o: /center/w/kate/Build_ARC/mod_param.o
read_couplepar.o: /center/w/kate/Build_ARC/mod_scalars.o

read_fltbiopar.o: oyster_floats_inp.h cppdefs.h globaldefs.h arctic.h
read_fltbiopar.o: /center/w/kate/Build_ARC/mod_behavior.o
read_fltbiopar.o: /center/w/kate/Build_ARC/mod_iounits.o
read_fltbiopar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_fltbiopar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_fltbiopar.o: /center/w/kate/Build_ARC/mod_param.o
read_fltbiopar.o: /center/w/kate/Build_ARC/mod_scalars.o

read_fltpar.o: cppdefs.h globaldefs.h arctic.h
read_fltpar.o: /center/w/kate/Build_ARC/mod_floats.o
read_fltpar.o: /center/w/kate/Build_ARC/mod_iounits.o
read_fltpar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_fltpar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_fltpar.o: /center/w/kate/Build_ARC/mod_param.o
read_fltpar.o: /center/w/kate/Build_ARC/mod_scalars.o

read_icepar.o: cppdefs.h globaldefs.h arctic.h
read_icepar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_icepar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_icepar.o: /center/w/kate/Build_ARC/mod_param.o
read_icepar.o: /center/w/kate/Build_ARC/mod_scalars.o

read_phypar.o: cppdefs.h globaldefs.h arctic.h
read_phypar.o: /center/w/kate/Build_ARC/mod_biology.o
read_phypar.o: /center/w/kate/Build_ARC/mod_coupler.o
read_phypar.o: /center/w/kate/Build_ARC/mod_fourdvar.o
read_phypar.o: /center/w/kate/Build_ARC/mod_iounits.o
read_phypar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_phypar.o: /center/w/kate/Build_ARC/mod_nesting.o
read_phypar.o: /center/w/kate/Build_ARC/mod_netcdf.o
read_phypar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_phypar.o: /center/w/kate/Build_ARC/mod_param.o
read_phypar.o: /center/w/kate/Build_ARC/mod_scalars.o
read_phypar.o: /center/w/kate/Build_ARC/mod_sediment.o
read_phypar.o: /center/w/kate/Build_ARC/mod_stepping.o
read_phypar.o: /center/w/kate/Build_ARC/mod_storage.o
read_phypar.o: /center/w/kate/Build_ARC/mod_strings.o

read_sedpar.o: sediment_inp.h cppdefs.h globaldefs.h arctic.h
read_sedpar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_sedpar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_sedpar.o: /center/w/kate/Build_ARC/mod_param.o
read_sedpar.o: /center/w/kate/Build_ARC/mod_scalars.o
read_sedpar.o: /center/w/kate/Build_ARC/mod_sediment.o

read_stapar.o: cppdefs.h globaldefs.h arctic.h
read_stapar.o: /center/w/kate/Build_ARC/mod_iounits.o
read_stapar.o: /center/w/kate/Build_ARC/mod_ncparam.o
read_stapar.o: /center/w/kate/Build_ARC/mod_parallel.o
read_stapar.o: /center/w/kate/Build_ARC/mod_param.o
read_stapar.o: /center/w/kate/Build_ARC/mod_scalars.o
read_stapar.o: /center/w/kate/Build_ARC/mod_sediment.o

regrid.o: cppdefs.h globaldefs.h arctic.h
regrid.o: /center/w/kate/Build_ARC/distribute.o
regrid.o: /center/w/kate/Build_ARC/interpolate.o
regrid.o: /center/w/kate/Build_ARC/mod_iounits.o
regrid.o: /center/w/kate/Build_ARC/mod_ncparam.o
regrid.o: /center/w/kate/Build_ARC/mod_parallel.o
regrid.o: /center/w/kate/Build_ARC/mod_param.o
regrid.o: /center/w/kate/Build_ARC/mod_scalars.o

rep_matrix.o: cppdefs.h globaldefs.h arctic.h
rep_matrix.o: /center/w/kate/Build_ARC/distribute.o
rep_matrix.o: /center/w/kate/Build_ARC/mod_fourdvar.o
rep_matrix.o: /center/w/kate/Build_ARC/mod_iounits.o
rep_matrix.o: /center/w/kate/Build_ARC/mod_parallel.o
rep_matrix.o: /center/w/kate/Build_ARC/mod_param.o
rep_matrix.o: /center/w/kate/Build_ARC/mod_scalars.o

set_2dfld.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
set_2dfld.o: /center/w/kate/Build_ARC/exchange_2d.o
set_2dfld.o: /center/w/kate/Build_ARC/mod_iounits.o
set_2dfld.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_2dfld.o: /center/w/kate/Build_ARC/mod_parallel.o
set_2dfld.o: /center/w/kate/Build_ARC/mod_param.o
set_2dfld.o: /center/w/kate/Build_ARC/mod_scalars.o
set_2dfld.o: /center/w/kate/Build_ARC/mp_exchange.o

set_2dfldr.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
set_2dfldr.o: /center/w/kate/Build_ARC/exchange_2d.o
set_2dfldr.o: /center/w/kate/Build_ARC/mod_iounits.o
set_2dfldr.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_2dfldr.o: /center/w/kate/Build_ARC/mod_parallel.o
set_2dfldr.o: /center/w/kate/Build_ARC/mod_param.o
set_2dfldr.o: /center/w/kate/Build_ARC/mod_scalars.o
set_2dfldr.o: /center/w/kate/Build_ARC/mp_exchange.o

set_3dfld.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
set_3dfld.o: /center/w/kate/Build_ARC/exchange_3d.o
set_3dfld.o: /center/w/kate/Build_ARC/mod_iounits.o
set_3dfld.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_3dfld.o: /center/w/kate/Build_ARC/mod_parallel.o
set_3dfld.o: /center/w/kate/Build_ARC/mod_param.o
set_3dfld.o: /center/w/kate/Build_ARC/mod_scalars.o
set_3dfld.o: /center/w/kate/Build_ARC/mp_exchange.o

set_3dfldr.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
set_3dfldr.o: /center/w/kate/Build_ARC/exchange_3d.o
set_3dfldr.o: /center/w/kate/Build_ARC/mod_iounits.o
set_3dfldr.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_3dfldr.o: /center/w/kate/Build_ARC/mod_parallel.o
set_3dfldr.o: /center/w/kate/Build_ARC/mod_param.o
set_3dfldr.o: /center/w/kate/Build_ARC/mod_scalars.o
set_3dfldr.o: /center/w/kate/Build_ARC/mp_exchange.o

set_contact.o: cppdefs.h globaldefs.h arctic.h
set_contact.o: /center/w/kate/Build_ARC/mod_iounits.o
set_contact.o: /center/w/kate/Build_ARC/mod_nesting.o
set_contact.o: /center/w/kate/Build_ARC/mod_netcdf.o
set_contact.o: /center/w/kate/Build_ARC/mod_parallel.o
set_contact.o: /center/w/kate/Build_ARC/mod_param.o
set_contact.o: /center/w/kate/Build_ARC/mod_scalars.o
set_contact.o: /center/w/kate/Build_ARC/strings.o

set_diags.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_diags.o: /center/w/kate/Build_ARC/bc_2d.o /center/w/kate/Build_ARC/bc_3d.o
set_diags.o: /center/w/kate/Build_ARC/exchange_2d.o
set_diags.o: /center/w/kate/Build_ARC/mod_biology.o
set_diags.o: /center/w/kate/Build_ARC/mod_diags.o
set_diags.o: /center/w/kate/Build_ARC/mod_grid.o
set_diags.o: /center/w/kate/Build_ARC/mod_ocean.o
set_diags.o: /center/w/kate/Build_ARC/mod_param.o
set_diags.o: /center/w/kate/Build_ARC/mod_scalars.o
set_diags.o: /center/w/kate/Build_ARC/mod_stepping.o
set_diags.o: /center/w/kate/Build_ARC/mp_exchange.o

set_filter.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_filter.o: /center/w/kate/Build_ARC/mod_average.o
set_filter.o: /center/w/kate/Build_ARC/mod_filter.o
set_filter.o: /center/w/kate/Build_ARC/mod_forces.o
set_filter.o: /center/w/kate/Build_ARC/mod_grid.o
set_filter.o: /center/w/kate/Build_ARC/mod_ice.o
set_filter.o: /center/w/kate/Build_ARC/mod_mixing.o
set_filter.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_filter.o: /center/w/kate/Build_ARC/mod_ocean.o
set_filter.o: /center/w/kate/Build_ARC/mod_parallel.o
set_filter.o: /center/w/kate/Build_ARC/mod_param.o
set_filter.o: /center/w/kate/Build_ARC/mod_scalars.o
set_filter.o: /center/w/kate/Build_ARC/mod_stepping.o
set_filter.o: /center/w/kate/Build_ARC/set_masks.o
set_filter.o: /center/w/kate/Build_ARC/uv_rotate.o

set_masks.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
set_masks.o: /center/w/kate/Build_ARC/exchange_2d.o
set_masks.o: /center/w/kate/Build_ARC/mod_grid.o
set_masks.o: /center/w/kate/Build_ARC/mod_param.o
set_masks.o: /center/w/kate/Build_ARC/mod_scalars.o
set_masks.o: /center/w/kate/Build_ARC/mod_sources.o
set_masks.o: /center/w/kate/Build_ARC/mp_exchange.o

set_ngfld.o: cppdefs.h globaldefs.h arctic.h
set_ngfld.o: /center/w/kate/Build_ARC/mod_iounits.o
set_ngfld.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_ngfld.o: /center/w/kate/Build_ARC/mod_parallel.o
set_ngfld.o: /center/w/kate/Build_ARC/mod_param.o
set_ngfld.o: /center/w/kate/Build_ARC/mod_scalars.o

set_ngfldr.o: cppdefs.h globaldefs.h arctic.h
set_ngfldr.o: /center/w/kate/Build_ARC/mod_iounits.o
set_ngfldr.o: /center/w/kate/Build_ARC/mod_ncparam.o
set_ngfldr.o: /center/w/kate/Build_ARC/mod_parallel.o
set_ngfldr.o: /center/w/kate/Build_ARC/mod_param.o
set_ngfldr.o: /center/w/kate/Build_ARC/mod_scalars.o

set_scoord.o: cppdefs.h globaldefs.h arctic.h
set_scoord.o: /center/w/kate/Build_ARC/mod_grid.o
set_scoord.o: /center/w/kate/Build_ARC/mod_iounits.o
set_scoord.o: /center/w/kate/Build_ARC/mod_parallel.o
set_scoord.o: /center/w/kate/Build_ARC/mod_param.o
set_scoord.o: /center/w/kate/Build_ARC/mod_scalars.o

set_weights.o: cppdefs.h globaldefs.h arctic.h
set_weights.o: /center/w/kate/Build_ARC/mod_iounits.o
set_weights.o: /center/w/kate/Build_ARC/mod_parallel.o
set_weights.o: /center/w/kate/Build_ARC/mod_param.o
set_weights.o: /center/w/kate/Build_ARC/mod_scalars.o

shapiro.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
shapiro.o: /center/w/kate/Build_ARC/mod_param.o

sqlq.o: cppdefs.h globaldefs.h arctic.h
sqlq.o: /center/w/kate/Build_ARC/mod_kinds.o

state_addition.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
state_addition.o: /center/w/kate/Build_ARC/mod_ncparam.o
state_addition.o: /center/w/kate/Build_ARC/mod_param.o
state_addition.o: /center/w/kate/Build_ARC/mod_scalars.o

state_copy.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
state_copy.o: /center/w/kate/Build_ARC/mod_ncparam.o
state_copy.o: /center/w/kate/Build_ARC/mod_param.o
state_copy.o: /center/w/kate/Build_ARC/mod_scalars.o

state_dotprod.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
state_dotprod.o: /center/w/kate/Build_ARC/distribute.o
state_dotprod.o: /center/w/kate/Build_ARC/mod_ncparam.o
state_dotprod.o: /center/w/kate/Build_ARC/mod_parallel.o
state_dotprod.o: /center/w/kate/Build_ARC/mod_param.o
state_dotprod.o: /center/w/kate/Build_ARC/mod_scalars.o

state_initialize.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
state_initialize.o: /center/w/kate/Build_ARC/mod_ncparam.o
state_initialize.o: /center/w/kate/Build_ARC/mod_param.o
state_initialize.o: /center/w/kate/Build_ARC/mod_scalars.o

state_product.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
state_product.o: /center/w/kate/Build_ARC/distribute.o
state_product.o: /center/w/kate/Build_ARC/mod_ncparam.o
state_product.o: /center/w/kate/Build_ARC/mod_parallel.o
state_product.o: /center/w/kate/Build_ARC/mod_param.o
state_product.o: /center/w/kate/Build_ARC/mod_scalars.o

state_scale.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
state_scale.o: /center/w/kate/Build_ARC/mod_ncparam.o
state_scale.o: /center/w/kate/Build_ARC/mod_param.o
state_scale.o: /center/w/kate/Build_ARC/mod_scalars.o

stats_modobs.o: cppdefs.h globaldefs.h arctic.h
stats_modobs.o: /center/w/kate/Build_ARC/mod_fourdvar.o
stats_modobs.o: /center/w/kate/Build_ARC/mod_iounits.o
stats_modobs.o: /center/w/kate/Build_ARC/mod_ncparam.o
stats_modobs.o: /center/w/kate/Build_ARC/mod_netcdf.o
stats_modobs.o: /center/w/kate/Build_ARC/mod_parallel.o
stats_modobs.o: /center/w/kate/Build_ARC/mod_param.o
stats_modobs.o: /center/w/kate/Build_ARC/mod_scalars.o

stiffness.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
stiffness.o: /center/w/kate/Build_ARC/distribute.o
stiffness.o: /center/w/kate/Build_ARC/mod_grid.o
stiffness.o: /center/w/kate/Build_ARC/mod_iounits.o
stiffness.o: /center/w/kate/Build_ARC/mod_ocean.o
stiffness.o: /center/w/kate/Build_ARC/mod_parallel.o
stiffness.o: /center/w/kate/Build_ARC/mod_param.o
stiffness.o: /center/w/kate/Build_ARC/mod_scalars.o

strings.o: cppdefs.h globaldefs.h arctic.h

sum_grad.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
sum_grad.o: /center/w/kate/Build_ARC/mod_boundary.o
sum_grad.o: /center/w/kate/Build_ARC/mod_forces.o
sum_grad.o: /center/w/kate/Build_ARC/mod_ncparam.o
sum_grad.o: /center/w/kate/Build_ARC/mod_ocean.o
sum_grad.o: /center/w/kate/Build_ARC/mod_param.o
sum_grad.o: /center/w/kate/Build_ARC/mod_scalars.o

timers.o: cppdefs.h globaldefs.h arctic.h
timers.o: /center/w/kate/Build_ARC/distribute.o
timers.o: /center/w/kate/Build_ARC/mod_iounits.o
timers.o: /center/w/kate/Build_ARC/mod_parallel.o
timers.o: /center/w/kate/Build_ARC/mod_param.o
timers.o: /center/w/kate/Build_ARC/mod_strings.o

uv_rotate.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
uv_rotate.o: /center/w/kate/Build_ARC/exchange_2d.o
uv_rotate.o: /center/w/kate/Build_ARC/exchange_3d.o
uv_rotate.o: /center/w/kate/Build_ARC/mod_param.o
uv_rotate.o: /center/w/kate/Build_ARC/mod_scalars.o
uv_rotate.o: /center/w/kate/Build_ARC/mp_exchange.o

vorticity.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
vorticity.o: /center/w/kate/Build_ARC/exchange_2d.o
vorticity.o: /center/w/kate/Build_ARC/exchange_3d.o
vorticity.o: /center/w/kate/Build_ARC/mod_average.o
vorticity.o: /center/w/kate/Build_ARC/mod_grid.o
vorticity.o: /center/w/kate/Build_ARC/mod_ncparam.o
vorticity.o: /center/w/kate/Build_ARC/mod_ocean.o
vorticity.o: /center/w/kate/Build_ARC/mod_param.o
vorticity.o: /center/w/kate/Build_ARC/mod_scalars.o
vorticity.o: /center/w/kate/Build_ARC/mod_stepping.o
vorticity.o: /center/w/kate/Build_ARC/mp_exchange.o

white_noise.o: cppdefs.h globaldefs.h arctic.h
white_noise.o: /center/w/kate/Build_ARC/distribute.o
white_noise.o: /center/w/kate/Build_ARC/mod_kinds.o
white_noise.o: /center/w/kate/Build_ARC/mod_parallel.o
white_noise.o: /center/w/kate/Build_ARC/mod_param.o
white_noise.o: /center/w/kate/Build_ARC/mod_scalars.o
white_noise.o: /center/w/kate/Build_ARC/mp_exchange.o
white_noise.o: /center/w/kate/Build_ARC/nrutil.o

wpoints.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
wpoints.o: /center/w/kate/Build_ARC/distribute.o
wpoints.o: /center/w/kate/Build_ARC/mod_grid.o
wpoints.o: /center/w/kate/Build_ARC/mod_iounits.o
wpoints.o: /center/w/kate/Build_ARC/mod_ncparam.o
wpoints.o: /center/w/kate/Build_ARC/mod_parallel.o
wpoints.o: /center/w/kate/Build_ARC/mod_param.o
wpoints.o: /center/w/kate/Build_ARC/mod_scalars.o
wpoints.o: /center/w/kate/Build_ARC/mod_storage.o

wrt_avg.o: cppdefs.h globaldefs.h arctic.h
wrt_avg.o: /center/w/kate/Build_ARC/mod_average.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_biology.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_ice.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_param.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_sedbed.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_avg.o: /center/w/kate/Build_ARC/mod_tides.o
wrt_avg.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_avg.o: /center/w/kate/Build_ARC/nf_fwrite3d.o

wrt_avg2.o: cppdefs.h globaldefs.h arctic.h
wrt_avg2.o: /center/w/kate/Build_ARC/distribute.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_average2.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_ice.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_param.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_avg2.o: /center/w/kate/Build_ARC/mod_tides.o
wrt_avg2.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_avg2.o: /center/w/kate/Build_ARC/nf_fwrite3d.o

wrt_diags.o: cppdefs.h globaldefs.h arctic.h
wrt_diags.o: /center/w/kate/Build_ARC/mod_biology.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_diags.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_param.o
wrt_diags.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_diags.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_diags.o: /center/w/kate/Build_ARC/nf_fwrite3d.o

wrt_error.o: cppdefs.h globaldefs.h arctic.h
wrt_error.o: /center/w/kate/Build_ARC/mod_boundary.o
wrt_error.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_error.o: /center/w/kate/Build_ARC/mod_fourdvar.o
wrt_error.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_error.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_error.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_error.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_error.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_error.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_error.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_error.o: /center/w/kate/Build_ARC/mod_param.o
wrt_error.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_error.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_error.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_error.o: /center/w/kate/Build_ARC/nf_fwrite2d_bry.o
wrt_error.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_error.o: /center/w/kate/Build_ARC/nf_fwrite3d_bry.o

wrt_filt.o: cppdefs.h globaldefs.h arctic.h
wrt_filt.o: /center/w/kate/Build_ARC/mod_filter.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_param.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_filt.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_filt.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_filt.o: /center/w/kate/Build_ARC/nf_fwrite3d.o

wrt_floats.o: cppdefs.h globaldefs.h arctic.h
wrt_floats.o: /center/w/kate/Build_ARC/mod_floats.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_param.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_floats.o: /center/w/kate/Build_ARC/mod_stepping.o

wrt_gst.o: cppdefs.h globaldefs.h arctic.h
wrt_gst.o: /center/w/kate/Build_ARC/distribute.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_param.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_gst.o: /center/w/kate/Build_ARC/mod_storage.o

wrt_hessian.o: cppdefs.h globaldefs.h arctic.h
wrt_hessian.o: /center/w/kate/Build_ARC/mod_boundary.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_param.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_hessian.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_hessian.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_hessian.o: /center/w/kate/Build_ARC/nf_fwrite2d_bry.o
wrt_hessian.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_hessian.o: /center/w/kate/Build_ARC/nf_fwrite3d_bry.o

wrt_his.o: cppdefs.h globaldefs.h arctic.h
wrt_his.o: /center/w/kate/Build_ARC/mod_bbl.o
wrt_his.o: /center/w/kate/Build_ARC/mod_biology.o
wrt_his.o: /center/w/kate/Build_ARC/mod_boundary.o
wrt_his.o: /center/w/kate/Build_ARC/mod_coupling.o
wrt_his.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_his.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_his.o: /center/w/kate/Build_ARC/mod_ice.o
wrt_his.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_his.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_his.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_his.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_his.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_his.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_his.o: /center/w/kate/Build_ARC/mod_param.o
wrt_his.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_his.o: /center/w/kate/Build_ARC/mod_sedbed.o
wrt_his.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_his.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_his.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_his.o: /center/w/kate/Build_ARC/nf_fwrite2d_bry.o
wrt_his.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_his.o: /center/w/kate/Build_ARC/nf_fwrite3d_bry.o
wrt_his.o: /center/w/kate/Build_ARC/omega.o
wrt_his.o: /center/w/kate/Build_ARC/uv_rotate.o

wrt_his2.o: cppdefs.h globaldefs.h arctic.h
wrt_his2.o: /center/w/kate/Build_ARC/distribute.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_ice.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_param.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_his2.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_his2.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_his2.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_his2.o: /center/w/kate/Build_ARC/uv_rotate.o

wrt_impulse.o: set_bounds.h cppdefs.h globaldefs.h arctic.h
wrt_impulse.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_param.o
wrt_impulse.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_impulse.o: /center/w/kate/Build_ARC/nf_fread2d.o
wrt_impulse.o: /center/w/kate/Build_ARC/nf_fread3d.o
wrt_impulse.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_impulse.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_impulse.o: /center/w/kate/Build_ARC/strings.o

wrt_info.o: sediment_wrt.h npzd_Powell_wrt.h fennel_wrt.h umaine_wrt.h
wrt_info.o: npzd_iron_wrt.h oyster_floats_wrt.h npzd_Franks_wrt.h ecosim_wrt.h
wrt_info.o: cppdefs.h globaldefs.h arctic.h nemuro_wrt.h
wrt_info.o: /center/w/kate/Build_ARC/distribute.o
wrt_info.o: /center/w/kate/Build_ARC/extract_sta.o
wrt_info.o: /center/w/kate/Build_ARC/mod_behavior.o
wrt_info.o: /center/w/kate/Build_ARC/mod_biology.o
wrt_info.o: /center/w/kate/Build_ARC/mod_fourdvar.o
wrt_info.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_info.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_info.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_info.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_info.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_info.o: /center/w/kate/Build_ARC/mod_param.o
wrt_info.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_info.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_info.o: /center/w/kate/Build_ARC/mod_sources.o
wrt_info.o: /center/w/kate/Build_ARC/mod_storage.o
wrt_info.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_info.o: /center/w/kate/Build_ARC/strings.o

wrt_ini.o: cppdefs.h globaldefs.h arctic.h
wrt_ini.o: /center/w/kate/Build_ARC/mod_boundary.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_fourdvar.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_param.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_ini.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_ini.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_ini.o: /center/w/kate/Build_ARC/nf_fwrite2d_bry.o
wrt_ini.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_ini.o: /center/w/kate/Build_ARC/nf_fwrite3d_bry.o

wrt_rst.o: cppdefs.h globaldefs.h arctic.h
wrt_rst.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_ice.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_param.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_sedbed.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_rst.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_rst.o: /center/w/kate/Build_ARC/nf_fwrite2d.o
wrt_rst.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_rst.o: /center/w/kate/Build_ARC/nf_fwrite4d.o

wrt_station.o: cppdefs.h globaldefs.h arctic.h
wrt_station.o: /center/w/kate/Build_ARC/extract_sta.o
wrt_station.o: /center/w/kate/Build_ARC/mod_bbl.o
wrt_station.o: /center/w/kate/Build_ARC/mod_forces.o
wrt_station.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_station.o: /center/w/kate/Build_ARC/mod_ice.o
wrt_station.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_station.o: /center/w/kate/Build_ARC/mod_mixing.o
wrt_station.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_station.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_station.o: /center/w/kate/Build_ARC/mod_ocean.o
wrt_station.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_station.o: /center/w/kate/Build_ARC/mod_param.o
wrt_station.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_station.o: /center/w/kate/Build_ARC/mod_sedbed.o
wrt_station.o: /center/w/kate/Build_ARC/mod_sediment.o
wrt_station.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_station.o: /center/w/kate/Build_ARC/uv_rotate.o

wrt_tides.o: cppdefs.h globaldefs.h arctic.h
wrt_tides.o: /center/w/kate/Build_ARC/mod_grid.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_iounits.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_ncparam.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_netcdf.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_parallel.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_param.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_scalars.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_stepping.o
wrt_tides.o: /center/w/kate/Build_ARC/mod_tides.o
wrt_tides.o: /center/w/kate/Build_ARC/nf_fwrite3d.o
wrt_tides.o: /center/w/kate/Build_ARC/nf_fwrite4d.o

zeta_balance.o: set_bounds.h tile.h cppdefs.h globaldefs.h arctic.h
zeta_balance.o: /center/w/kate/Build_ARC/distribute.o
zeta_balance.o: /center/w/kate/Build_ARC/exchange_2d.o
zeta_balance.o: /center/w/kate/Build_ARC/exchange_3d.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_coupling.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_fourdvar.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_grid.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_iounits.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_mixing.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_ncparam.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_ocean.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_parallel.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_param.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_scalars.o
zeta_balance.o: /center/w/kate/Build_ARC/mod_stepping.o
zeta_balance.o: /center/w/kate/Build_ARC/mp_exchange.o
zeta_balance.o: /center/w/kate/Build_ARC/rho_eos.o
zeta_balance.o: /center/w/kate/Build_ARC/set_depth.o

/center/w/kate/Build_ARC/analytical_mod.mod: analytical.o
/center/w/kate/Build_ARC/array_modes_mod.mod: array_modes.o
/center/w/kate/Build_ARC/back_cost_mod.mod: back_cost.o
/center/w/kate/Build_ARC/bbl_mod.mod: bbl.o
/center/w/kate/Build_ARC/bc_2d_mod.mod: bc_2d.o
/center/w/kate/Build_ARC/bc_3d_mod.mod: bc_3d.o
/center/w/kate/Build_ARC/bc_bry2d_mod.mod: bc_bry2d.o
/center/w/kate/Build_ARC/bc_bry3d_mod.mod: bc_bry3d.o
/center/w/kate/Build_ARC/biology_floats_mod.mod: biology_floats.o
/center/w/kate/Build_ARC/biology_mod.mod: biology.o
/center/w/kate/Build_ARC/bulk_flux_mod.mod: bulk_flux.o
/center/w/kate/Build_ARC/bvf_mix_mod.mod: bvf_mix.o
/center/w/kate/Build_ARC/cawdir_eval_mod.mod: cawdir_eval.o
/center/w/kate/Build_ARC/ccsm_flux_mod.mod: ccsm_flux.o
/center/w/kate/Build_ARC/cgradient_mod.mod: cgradient.o
/center/w/kate/Build_ARC/conv_2d_mod.mod: conv_2d.o
/center/w/kate/Build_ARC/conv_3d_bry_mod.mod: conv_bry3d.o
/center/w/kate/Build_ARC/conv_3d_mod.mod: conv_3d.o
/center/w/kate/Build_ARC/conv_bry2d_mod.mod: conv_bry2d.o
/center/w/kate/Build_ARC/convolve_mod.mod: convolve.o
/center/w/kate/Build_ARC/cost_grad_mod.mod: cost_grad.o
/center/w/kate/Build_ARC/def_var_mod.mod: def_var.o
/center/w/kate/Build_ARC/diag_mod.mod: diag.o
/center/w/kate/Build_ARC/distribute_mod.mod: distribute.o
/center/w/kate/Build_ARC/dotproduct_mod.mod: dotproduct.o
/center/w/kate/Build_ARC/erf_mod.mod: erf.o
/center/w/kate/Build_ARC/exchange_2d_mod.mod: exchange_2d.o
/center/w/kate/Build_ARC/exchange_3d_mod.mod: exchange_3d.o
/center/w/kate/Build_ARC/extract_obs_mod.mod: extract_obs.o
/center/w/kate/Build_ARC/extract_sta_mod.mod: extract_sta.o
/center/w/kate/Build_ARC/forcing_mod.mod: forcing.o
/center/w/kate/Build_ARC/frc_adjust_mod.mod: frc_adjust.o
/center/w/kate/Build_ARC/frc_weak_mod.mod: frc_weak.o
/center/w/kate/Build_ARC/gls_corstep_mod.mod: gls_corstep.o
/center/w/kate/Build_ARC/gls_prestep_mod.mod: gls_prestep.o
/center/w/kate/Build_ARC/hmixing_mod.mod: hmixing.o
/center/w/kate/Build_ARC/ini_adjust_mod.mod: ini_adjust.o
/center/w/kate/Build_ARC/ini_fields_mod.mod: ini_fields.o
/center/w/kate/Build_ARC/ini_hmixcoef_mod.mod: ini_hmixcoef.o
/center/w/kate/Build_ARC/ini_lanczos_mod.mod: ini_lanczos.o
/center/w/kate/Build_ARC/inner2state_mod.mod: inner2state.o
/center/w/kate/Build_ARC/interp_floats_diapw_mod.mod: interp_floats_diapW.o
/center/w/kate/Build_ARC/interp_floats_mod.mod: interp_floats.o
/center/w/kate/Build_ARC/interpolate_mod.mod: interpolate.o
/center/w/kate/Build_ARC/lmd_bkpp_mod.mod: lmd_bkpp.o
/center/w/kate/Build_ARC/lmd_skpp_mod.mod: lmd_skpp.o
/center/w/kate/Build_ARC/lmd_vmix_mod.mod: lmd_vmix.o
/center/w/kate/Build_ARC/metrics_mod.mod: metrics.o
/center/w/kate/Build_ARC/mp_exchange_mod.mod: mp_exchange.o
/center/w/kate/Build_ARC/mpdata_adiff_mod.mod: mpdata_adiff.o
/center/w/kate/Build_ARC/my25_corstep_mod.mod: my25_corstep.o
/center/w/kate/Build_ARC/my25_prestep_mod.mod: my25_prestep.o
/center/w/kate/Build_ARC/nesting_mod.mod: nesting.o
/center/w/kate/Build_ARC/nf_fread2d_bry_mod.mod: nf_fread2d_bry.o
/center/w/kate/Build_ARC/nf_fread2d_mod.mod: nf_fread2d.o
/center/w/kate/Build_ARC/nf_fread3d_bry_mod.mod: nf_fread3d_bry.o
/center/w/kate/Build_ARC/nf_fread3d_mod.mod: nf_fread3d.o
/center/w/kate/Build_ARC/nf_fread4d_mod.mod: nf_fread4d.o
/center/w/kate/Build_ARC/nf_fwrite2d_bry_mod.mod: nf_fwrite2d_bry.o
/center/w/kate/Build_ARC/nf_fwrite2d_mod.mod: nf_fwrite2d.o
/center/w/kate/Build_ARC/nf_fwrite3d_bry_mod.mod: nf_fwrite3d_bry.o
/center/w/kate/Build_ARC/nf_fwrite3d_mod.mod: nf_fwrite3d.o
/center/w/kate/Build_ARC/nf_fwrite4d_mod.mod: nf_fwrite4d.o
/center/w/kate/Build_ARC/normalization_mod.mod: normalization.o
/center/w/kate/Build_ARC/obc_adjust_mod.mod: obc_adjust.o
/center/w/kate/Build_ARC/obc_volcons_mod.mod: obc_volcons.o
/center/w/kate/Build_ARC/ocean_control_mod.mod: ocean_control.o
/center/w/kate/Build_ARC/ocean_coupler_mod.mod: ocean_coupler.o
/center/w/kate/Build_ARC/omega_mod.mod: omega.o
/center/w/kate/Build_ARC/packing_mod.mod: packing.o
/center/w/kate/Build_ARC/posterior_mod.mod: posterior.o
/center/w/kate/Build_ARC/posterior_var_mod.mod: posterior_var.o
/center/w/kate/Build_ARC/pre_step3d_mod.mod: pre_step3d.o
/center/w/kate/Build_ARC/propagator_mod.mod: propagator.o
/center/w/kate/Build_ARC/prsgrd_mod.mod: prsgrd.o
/center/w/kate/Build_ARC/pt3dbc_mod.mod: pt3dbc_im.o
/center/w/kate/Build_ARC/radiation_stress_mod.mod: radiation_stress.o
/center/w/kate/Build_ARC/random_ic_mod.mod: random_ic.o
/center/w/kate/Build_ARC/rho_eos_mod.mod: rho_eos.o
/center/w/kate/Build_ARC/rhs3d_mod.mod: rhs3d.o
/center/w/kate/Build_ARC/roms_export_mod.mod: roms_export.o
/center/w/kate/Build_ARC/roms_import_mod.mod: roms_import.o
/center/w/kate/Build_ARC/sed_bed_mod.mod: sed_bed.o
/center/w/kate/Build_ARC/sed_bedload_mod.mod: sed_bedload.o
/center/w/kate/Build_ARC/sed_fluxes_mod.mod: sed_fluxes.o
/center/w/kate/Build_ARC/sed_settling_mod.mod: sed_settling.o
/center/w/kate/Build_ARC/sed_surface_mod.mod: sed_surface.o
/center/w/kate/Build_ARC/sediment_mod.mod: sediment.o
/center/w/kate/Build_ARC/set_2dfld_mod.mod: set_2dfld.o
/center/w/kate/Build_ARC/set_2dfldr_mod.mod: set_2dfldr.o
/center/w/kate/Build_ARC/set_3dfld_mod.mod: set_3dfld.o
/center/w/kate/Build_ARC/set_3dfldr_mod.mod: set_3dfldr.o
/center/w/kate/Build_ARC/set_avg2_mod.mod: set_avg2.o
/center/w/kate/Build_ARC/set_avg_mod.mod: set_avg.o
/center/w/kate/Build_ARC/set_depth_mod.mod: set_depth.o
/center/w/kate/Build_ARC/set_masks_mod.mod: set_masks.o
/center/w/kate/Build_ARC/set_massflux_mod.mod: set_massflux.o
/center/w/kate/Build_ARC/set_tides_mod.mod: set_tides.o
/center/w/kate/Build_ARC/set_vbc_mod.mod: set_vbc.o
/center/w/kate/Build_ARC/set_zeta_mod.mod: set_zeta.o
/center/w/kate/Build_ARC/shapiro_mod.mod: shapiro.o
/center/w/kate/Build_ARC/state_addition_mod.mod: state_addition.o
/center/w/kate/Build_ARC/state_copy_mod.mod: state_copy.o
/center/w/kate/Build_ARC/state_dotprod_mod.mod: state_dotprod.o
/center/w/kate/Build_ARC/state_initialize_mod.mod: state_initialize.o
/center/w/kate/Build_ARC/state_product_mod.mod: state_product.o
/center/w/kate/Build_ARC/state_scale_mod.mod: state_scale.o
/center/w/kate/Build_ARC/step2d_mod.mod: step2d.o
/center/w/kate/Build_ARC/step3d_t_mod.mod: step3d_t.o
/center/w/kate/Build_ARC/step3d_uv_mod.mod: step3d_uv.o
/center/w/kate/Build_ARC/step_floats_mod.mod: step_floats.o
/center/w/kate/Build_ARC/stiffness_mod.mod: stiffness.o
/center/w/kate/Build_ARC/strings_mod.mod: strings.o
/center/w/kate/Build_ARC/sum_grad_mod.mod: sum_grad.o
/center/w/kate/Build_ARC/t3dbc_mod.mod: t3dbc_im.o
/center/w/kate/Build_ARC/t3dmix_mod.mod: t3dmix.o
/center/w/kate/Build_ARC/tkebc_mod.mod: tkebc_im.o
/center/w/kate/Build_ARC/u2dbc_mod.mod: u2dbc_im.o
/center/w/kate/Build_ARC/u3dbc_mod.mod: u3dbc_im.o
/center/w/kate/Build_ARC/uv3dmix_mod.mod: uv3dmix.o
/center/w/kate/Build_ARC/uv_rotate_mod.mod: uv_rotate.o
/center/w/kate/Build_ARC/v2dbc_mod.mod: v2dbc_im.o
/center/w/kate/Build_ARC/v3dbc_mod.mod: v3dbc_im.o
/center/w/kate/Build_ARC/vorticity_mod.mod: vorticity.o
/center/w/kate/Build_ARC/vwalk_floats_mod.mod: vwalk_floats.o
/center/w/kate/Build_ARC/wetdry_mod.mod: wetdry.o
/center/w/kate/Build_ARC/white_noise_mod.mod: white_noise.o
/center/w/kate/Build_ARC/wpoints_mod.mod: wpoints.o
/center/w/kate/Build_ARC/wvelocity_mod.mod: wvelocity.o
/center/w/kate/Build_ARC/zeta_balance_mod.mod: zeta_balance.o
/center/w/kate/Build_ARC/zetabc_mod.mod: zetabc.o
