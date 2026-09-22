---
title: "Example SCOTCH Configuration"
nav_order: 8
parent: "Compilation and Libraries"
grand_parent: "Compiling and Running"
layout: default
render_with_liquid: false
---

# Example configuration file for Scotch 5.1.12b

```make
EXE		=
LIB		= .a
OBJ		= .o

MAKE		= make
AR		= ar
ARFLAGS		= -ruv
CAT		= cat
CCS		= icc
CCP		= mpicc -cc=icc
CCD		= icc
CFLAGS		= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER -DSCOTCH_PTHREAD -restrict -DIDXSIZE64 -DINTSIZE32
CLIBFLAGS	=
LDFLAGS		= -g -lz -lm -lrt -I${I_MPI_ROOT}/intel64/include
CP		= cp
LEX		= flex -Pscotchyy -olex.yy.c
LN		= ln
MKDIR		= mkdir
MV		= mv
RANLIB		= ranlib
YACC		= bison -pscotchyy -y -b y
```

Note: it may be necessary to alter `LDFLAGS` to include `-lpthread`, see Scotch bugtracker [#15048](https://gforge.inria.fr/tracker/?func=detail&atid=1079&aid=15048&group_id=248).

Set `prefix` to install in a non-standard path here.
