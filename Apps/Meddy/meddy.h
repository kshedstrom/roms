/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2016 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for Mediterranean Eddy Example.
*/

#undef STATIONS
#define UV_ADV
#define UV_COR
#undef TS_A4HADVECTION
#undef TS_A4VADVECTION
#define  TS_U3HADVECTION
#define DJ_GRADPS
#define UV_VIS2
#define TS_DIF2
#define MIX_S_UV
#define MIX_S_TS
#define SOLVE3D
#define FLOATS
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_VMIX
#undef MY25_MIXING
#define  UV_LDRAG
#define  DIAGNOSTICS_UV
#define  DIAGNOSTICS_TS

#undef BIOLOGY
#ifdef BIOLOGY
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define SHORTWAVE
# define ANA_SRFLUX
# define ANA_CLOUD
# define ANA_HUMIDITY
# define ANA_TAIR
#endif

#undef OFFLINE_BIOLOGY
#undef OFFLINE_FLOATS
