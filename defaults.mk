DEBUG ?= 0 # Build in release mode by default

# To use the default compiler flags, set $(COMPILER_FAMILY)
# If you do not want to do this, be sure to set OUTPUT_MODULE_COMMAND to the correct
# one for your compiler in Makefile.inc.

# build in directories to not clutter the repository
MODDIR := .mod
OBJDIR := .obj
DEPDIR := .dep
$(shell mkdir -p $(MODDIR) $(OBJDIR) $(DEPDIR) >/dev/null)
INCLUDES += -I$(MODDIR)

# Detect the compiler vendors (sort to remove duplicates)
F_COMPILER_FAMILY :=$(sort $(shell $(FC) --version | grep -oim 1 'intel\|gcc\|gnu' | tr A-Z a-z | sed 's/gcc/gnu/'))
C_COMPILER_FAMILY :=$(sort $(shell $(CC) --version | grep -oim 1 'intel\|gcc\|gnu' | tr A-Z a-z | sed 's/gcc/gnu/'))
CXX_COMPILER_FAMILY :=$(sort $(shell $(CXX) --version | grep -oim 1 'intel\|gcc\|gnu\|g[+][+]' | tr A-Z a-z | sed -e 's/gcc/gnu/' -e 's/g[+][+]/gnu/'))
ifneq ($(F_COMPILER_FAMILY),$(C_COMPILER_FAMILY))
  $(error "Fortran compiler ($(F_COMPILER_FAMILY)) must be same as C compiler ($(C_COMPILER_FAMILY))")
endif
ifneq ($(CXX_COMPILER_FAMILY),$(C_COMPILER_FAMILY))
  $(error "C++ compiler ($(CXX_COMPILER_FAMILY)) must be same as C compiler ($(C_COMPILER_FAMILY))")
endif
ifeq ($(COMPILER_FAMILY),)
  COMPILER_FAMILY = $(F_COMPILER_FAMILY)
endif
ifneq ($(F_COMPILER_FAMILY),$(COMPILER_FAMILY))
  $(warning "*******************************************************************************************")
  $(warning "Set compiler family ($(COMPILER_FAMILY)) is different from $(FC) ($(F_COMPILER_FAMILY))")
  $(warning "*******************************************************************************************")
endif

# Clear some variables (makefile bug?)
F77FLAGS=
F90FLAGS=

# Default (GNU) preprocessor
GCPP ?= cpp

# Default flags for GCC
ifeq ($(COMPILER_FAMILY), gnu)
  FLAGS += -cpp -fopenmp
  FLAGS += -Wall -Wextra
  FLAGS += -Wno-unused-variable
  FFLAGS += -Wintrinsics-std
  FFLAGS += -Wcharacter-truncation
  FFLAGS += -Wsurprising -Wno-tabs
  FFLAGS += -ffree-line-length-none
  F77FLAGS += -fdefault-real-8 -fdefault-double-8
  ifeq ($(DEBUG), 1)
    FLAGS  += -g -Og -ggdb -fno-lto
    FLAGS  += -Wunused-variable
    FFLAGS += -fcheck=all
    FLAGS  += -ffpe-trap=invalid,zero,overflow,denormal
    FFLAGS += -ftrapv
    FFLAGS += -finit-real=snan -finit-integer=12345678
    FFLAGS += -Wconversion
    F90FLAGS += -fimplicit-none
  endif
  ifeq ($(DEBUG), 2)
    FFLAGS += -Wimplicit-interface -Wimplicit-procedure
  endif

  FFLAGS +=-J$(MODDIR)
endif

# Default flags for intel
ifeq ($(COMPILER_FAMILY), intel)
  COMPILER_MAJOR_VERSION=$(shell $(FC) -V 2>&1 | grep -o "Version [0-9]*" | cut -d' ' -f 2)
  FFLAGS += -align
  ifeq ($(shell test $(COMPILER_MAJOR_VERSION) -ge 15; echo $$?),0)
    FLAGS += -qopenmp
  else
    FLAGS += -openmp
  endif
  FFLAGS += -warn all
  FFLAGS += -warn nointerfaces
  FFLAGS += -warn nounused
  FFLAGS += -fpp
  FFLAGS += -r8
  F77FLAGS += -warn nodeclarations
  ifeq ($(DEBUG), 1)
    # Debug flags for ifort, see http://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
    FLAGS += -O0 -g -traceback
    FLAGS += -ftrapuv
    FLAGS += -debug all -debug-parameters
    FLAGS += -fstack-security-check
    FLAGS += -fpe0
    FFLAGS += -check all,noarg_temp_created
    FFLAGS += -check bounds
    FFLAGS += -check uninit
    FFLAGS += -init=snan -init=zero
    FFLAGS += -gen-interfaces -warn-interfaces
    F90FLAGS += -implicitnone
  endif

  FFLAGS +=-module $(MODDIR)
endif

#TODO identify good default flags for XLF
# Correct preprocessor-defines for IBM XLF Compiler
IBM_DEFINES = `echo $(DEFINES) | sed -e 's/^/-WF,/' -e 's/  */,/g'`
ifdef IBMFC
  DEFINES  := $(IBM_DEFINES) -DIBM_MACHINE
endif
# TODO set the option to output module files to a specific directory for XLF

# Make rules for specific files
# This is needed because the file stems must match and we do not really want to recreate the directory structure in $(OBJDIR) and $(DEPDIR)
# Also, it will make everything slightly nicer later
# Template for generating object files from source files
# Touch the .mod file again if it exists (because it is not written if there is no change to the interfaces, and this messes with the make rules)
# Note that gfortran outputs module files in lowercase of the module name, while intel remains case-sensitive. Best to use lowercase module names
define O_TEMPLATE
$(OBJDIR)/%.o $(MODDIR)/%.mod:: $(1)%.f90
	$$(FC) $$(FLAGS) $$(FFLAGS) $$(F90FLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o
	@test -e $(MODDIR)/$$*.mod && touch $(MODDIR)/$$*.mod || true

$(OBJDIR)/%.o:: $(1)%.f
	$$(FC) $$(FLAGS) $$(FFLAGS) $$(F77FLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o

$(OBJDIR)/%.o:: $(1)%.c
	$$(CC) $$(FLAGS) $$(CFLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o

$(OBJDIR)/%.o:: $(1)%.cpp
	$$(CXX) $$(FLAGS) $$(CXXFLAGS) $$(DEFINES) $$(INCLUDES) $$(EXTRA_FLAGS) -c $$< -o $(OBJDIR)/$$*.o
endef
# Template for generating dependencies from source file
define F90_D_TEMPLATE
$(DEPDIR)/%.d: $(1)%.f90 | $(MODDIR)/version.h
	@echo "Generating dependencies for $$<"
	@$(GCPP) -traditional-cpp -dI $$(DEFINES) $$(INCLUDES) $$< | util/makedepend $$< - $(DIRS) > $(DEPDIR)/$$(*F).d
endef

# This template defines a program $(file_stem)
# which has prerequisites $(OBJDIR)/$(file_stem).o and as determined by the output of obj_deps.sh
# Since make will try to build the .d files first and re-exec every time
# these will be up-to-date
file_stem=$(notdir $(basename $(1)))
define PROGRAM_TEMPLATE
$(notdir $(basename $(1))): $(OBJDIR)/$(file_stem).o $(shell ./util/obj_deps $(DEPDIR)/$(file_stem).d)
	$$(FC) $$(FLAGS) $$(EXTRA_FLAGS) $$(DEFINES) $$(INCLUDES) -o $(file_stem) $$^ $$(LIBS)
endef



LIBS += $(LIBLAPACK) $(LIBBLAS) $(OPENMPLIB)
DEFINES += -DJOREK_MODEL=$(MODEL_NUMBER) -DUSE_MPI

# Full-MHD models flags
ifeq (model710, $(MODEL))
  DEFINES  := $(DEFINES) -Dfullmhd
endif
ifeq (model711, $(MODEL))
  DEFINES  := $(DEFINES) -Dfullmhd
endif
ifeq (model712, $(MODEL))
  DEFINES  := $(DEFINES) -Dfullmhd
endif

# Generic models flags
ifeq (model303, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Vpar
endif
ifeq (model305, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Vpar
endif
ifeq (model306, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Vpar
endif
ifeq (model333, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Vpar
endif
ifeq (model401, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_TiTe
endif
ifeq (model500, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Neutrals
endif
ifeq (model501, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Impurities
endif
ifeq (model502, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_Impurities -DWITH_TiTe
endif
ifeq (model711, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_TiTe
endif
ifeq (model712, $(MODEL))
  DEFINES  := $(DEFINES) -DWITH_TiTe -DWITH_Neutrals
endif


ifeq (1, $(USE_FFTW))
  LIBS     := $(LIBS) $(LIBFFTW)
  DEFINES  := $(DEFINES) -DUSE_FFTW
  INCLUDES := $(INCLUDES) $(INC_FFTW)
endif

ifeq (1, $(USE_PASTIX_MURGE))
  LIBS     := $(LIBS) $(LIB_PASTIX_MURGE) $(LIB_PASTIX_BLAS)
  DEFINES  := $(DEFINES) -DUSE_MURGE
  INCLUDES := $(INCLUDES) $(INC_PASTIX)
endif

ifeq (1, $(USE_PASTIX))
  DEFINES  := $(DEFINES) -DUSE_PASTIX
  ifeq (0, $(USE_PASTIX_MURGE))
    LIBS     := $(LIBS) $(LIB_PASTIX) $(LIB_PASTIX_BLAS)
    INCLUDES := $(INCLUDES) $(INC_PASTIX)
  endif
  PASTIX_MEMORY_USAGE?=1
  ifeq (1, $(PASTIX_MEMORY_USAGE))
    DEFINES := $(DEFINES) -DMEMORY_USAGE
  endif
else
  # This is a hack to remove the linking problems that otherwise arise
  DEFINES += -Dpastix_fortran=fake_pastix_fortran
endif

ifeq (1, $(USE_PASTIX6))
  DEFINES  := $(DEFINES) -DUSE_PASTIX6
  LIBS     := $(LIBS) $(LIB_PASTIX6)  $(LIB_PASTIX6_BLAS)
  INCLUDES := $(INCLUDES) $(INC_PASTIX6)
endif

ifeq (1, $(USE_WSMP))
  DEFINES  := $(DEFINES) -DUSE_WSMP
  LIBS     := $(LIBS) $(LIB_WSMP)
endif

ifeq (1, $(USE_HIPS))
  LIBS := $(LIBS) $(LIBHIPS)
  INCLUDES := $(INCLUDES) $(INCHIPS)
  DEFINES  := $(DEFINES) -DUSE_HIPS
endif

ifeq (1, $(USE_MUMPS))
  LIBS := $(LIBS) $(LIB_MUMPS) $(ORDLIB) $(SCALAP) $(BLACS) $(LIBLAPACK) $(LIBBLAS) $(PPPLIB) $(OPENMP_LIB) $(LIBFFTW)
  INCLUDES := $(INCLUDES) -I$(INC_MUMPS) $(INC_MUMPS_EXTRA)
  DEFINES := $(DEFINES) -DUSE_MUMPS
endif

ifeq (1, $(USE_HDF5))
  LIBS     := $(LIBS) $(HDF5LIB)
  INCLUDES := $(INCLUDES) -I$(HDF5INCLUDE)
  DEFINES  := $(DEFINES) -DUSE_HDF5
else
  $(warning "USE_HDF5=1 is recommended for input/output")
endif

ifeq (1, $(USE_BLOCK))
  DEFINES := $(DEFINES) -DUSE_BLOCK
endif

ifeq (1, $(USE_DIRECT_CONSTRUCTION))
  DEFINES  := $(DEFINES) -DDIRECT_CONSTRUCTION
endif

ifeq (1, $(USE_COMPLEX_PRECOND))
  DEFINES  := $(DEFINES) -DUSE_COMPLEX_PRECOND
endif

ifeq (1, $(USE_INTSIZE64))
  DEFINES  := $(DEFINES) -DINTSIZE64
endif

ifeq (1, $(USE_STRUMPACK))
  DEFINES  := $(DEFINES) -DUSE_STRUMPACK
  LIBS     := $(LIBS) $(STRUMPACKLIB)
  INCLUDES := $(INCLUDES) $(STRUMPACKINC)
  EXTRA_FLAGS := $(EXTRA_FLAGS) -lstdc++ -std=c++14
endif


# Do not check to make these files to speed up and clean -d output
Makefile: ;
Makefile.inc: ;
%.mk: ;
%.f90: ;

# Try to create .mod/version.h, but only overwrite it if the contents have changed
$(MODDIR)/version.h:
	@echo "Generate .mod/version.h"
	@rm -f $@.tmp
	@echo "#define RCS_VERSION '`git describe --always --dirty --abbrev 2> /dev/null`'" >> $@.tmp
	@echo "#define RCS_LABEL '`git log -1 --format="%s (%D)" 2> /dev/null | sed -e "s/'/''/g" `'" >> $@.tmp
	@echo "#define RCS_TIME '`git log -1 --format="%ad" 2> /dev/null`'" >> $@.tmp
	@echo "#define compile_command '$(FC)'" >> $@.tmp
	@echo "#define compile_flags '$(FLAGS) $(FFLAGS) $(F90FLAGS) $(EXTRA_FLAGS)'" >> $@.tmp
	@echo "#define compile_includes '$(INCLUDES)'" >> $@.tmp
	@echo "#define compile_defines '$(DEFINES)'" >> $@.tmp
	@echo "#define compile_libs '$(LIBS)'" >> $@.tmp
	-@echo "#define compile_dir '`pwd`'" >> $@.tmp
	-@echo "#define compile_user '`whoami`'" >> $@.tmp
	-@echo "#define compile_machine '`hostname`'" >> $@.tmp
	@echo "#define compile_modules '$(LOADEDMODULES)'" >> $@.tmp
	@[ -f $@ ] && cmp --silent $@ $@.tmp || mv $@.tmp $@
	-@rm -f $@.tmp
