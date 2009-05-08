#include "cppdefs.h"
/*
** svn $Id: master.F 975 2009-05-05 22:51:13Z kate $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2009 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Master program to run ROMS/TOMS as single ocean model or coupled  **
**  to other models using the MCT or ESMF libraries.                  **
**                                                                    **
************************************************************************
*/

#if defined MODEL_COUPLING
# if defined MCT_LIB
#  include "mct_coupler.h"
# elif defined ESMF_LIB
#  include "esmf_coupler.h"
# endif
#else
# include "ocean.h"
#endif
