/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Franks et al. (1986) ecosystem   **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

              CASE ('idTvar(iNO3)')
                idTvar(iNO3)=varid
              CASE ('idTvar(iNH4)')
                idTvar(iNH4)=varid
              CASE ('idTvar(iPhS)')
                idTvar(iPhS)=varid
              CASE ('idTvar(iPhL)')
                idTvar(iPhL)=varid
              CASE ('idTvar(iMZS)')
                idTvar(iMZS)=varid
              CASE ('idTvar(iMZL)')
                idTvar(iMZL)=varid
              CASE ('idTvar(iCop)')
                idTvar(iCop)=varid
              CASE ('idTvar(iNCa)')
                idTvar(iNCa)=varid
              CASE ('idTvar(iEup)')
                idTvar(iEup)=varid
              CASE ('idTvar(iDet)')
                idTvar(iDet)=varid
#ifdef IRON_LIMIT
              CASE ('idTvar(iFe)')
                idTvar(iFe)=varid
#endif

             !-----------------------------
             !Stationary production tracers
             !-----------------------------
              CASE ('idTSvar(iPhSprd)')
                idTSvar(iPhSprd) = varid
              CASE ('idTSvar(iPhLprd)')
                idTSvar(iPhLprd) = varid
              CASE ('idTSvar(iMZSprd)')
                idTSvar(iMZSprd) = varid
              CASE ('idTSvar(iMZLprd)')
                idTSvar(iMZLprd) = varid
              CASE ('idTSvar(iCopPrd)')
                idTSvar(iCopPrd) = varid
              CASE ('idTSvar(iNCaPrd)')
                idTSvar(iNCaPrd) = varid
              CASE ('idTSvar(iEupPrd)')
                idTSvar(iEupPrd) = varid

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iNO3)')
                idTbry(iwest,iNO3)=varid
              CASE ('idTbry(ieast,iNO3)')
                idTbry(ieast,iNO3)=varid
              CASE ('idTbry(isouth,iNO3)')
                idTbry(isouth,iNO3)=varid
              CASE ('idTbry(inorth,iNO3)')
                idTbry(inorth,iNO3)=varid

              CASE ('idTbry(iwest,iNH4)')
                idTbry(iwest,iNH4)=varid
              CASE ('idTbry(ieast,iNH4)')
                idTbry(ieast,iNH4)=varid
              CASE ('idTbry(isouth,iNH4)')
                idTbry(isouth,iNH4)=varid
              CASE ('idTbry(inorth,iNH4)')
                idTbry(inorth,iNH4)=varid

              CASE ('idTbry(iwest,iPhS)')
                idTbry(iwest,iPhS)=varid
              CASE ('idTbry(ieast,iPhS)')
                idTbry(ieast,iPhS)=varid
              CASE ('idTbry(isouth,iPhS)')
                idTbry(isouth,iPhS)=varid
              CASE ('idTbry(inorth,iPhS)')
                idTbry(inorth,iPhS)=varid

              CASE ('idTbry(iwest,iPhL)')
                idTbry(iwest,iPhL)=varid
              CASE ('idTbry(ieast,iPhL)')
                idTbry(ieast,iPhL)=varid
              CASE ('idTbry(isouth,iPhL)')
                idTbry(isouth,iPhL)=varid
              CASE ('idTbry(inorth,iPhL)')
                idTbry(inorth,iPhL)=varid

              CASE ('idTbry(iwest,iMZS)')
                idTbry(iwest,iMZS)=varid
              CASE ('idTbry(ieast,iMZS)')
                idTbry(ieast,iMZS)=varid
              CASE ('idTbry(isouth,iMZS)')
                idTbry(isouth,iMZS)=varid
              CASE ('idTbry(inorth,iMZS)')
                idTbry(inorth,iMZS)=varid

              CASE ('idTbry(iwest,iMZL)')
                idTbry(iwest,iMZL)=varid
              CASE ('idTbry(ieast,iMZL)')
                idTbry(ieast,iMZL)=varid
              CASE ('idTbry(isouth,iMZL)')
                idTbry(isouth,iMZL)=varid
              CASE ('idTbry(inorth,iMZL)')
                idTbry(inorth,iMZL)=varid

              CASE ('idTbry(iwest,iCop)')
                idTbry(iwest,iCop)=varid
              CASE ('idTbry(ieast,iCop)')
                idTbry(ieast,iCop)=varid
              CASE ('idTbry(isouth,iCop)')
                idTbry(isouth,iCop)=varid
              CASE ('idTbry(inorth,iCop)')
                idTbry(inorth,iCop)=varid

              CASE ('idTbry(iwest,iNCa)')
                idTbry(iwest,iNCa)=varid
              CASE ('idTbry(ieast,iNCa)')
                idTbry(ieast,iNCa)=varid
              CASE ('idTbry(isouth,iNCa)')
                idTbry(isouth,iNCa)=varid
              CASE ('idTbry(inorth,iNCa)')
                idTbry(inorth,iNCa)=varid

              CASE ('idTbry(iwest,iEup)')
                idTbry(iwest,iEup)=varid
              CASE ('idTbry(ieast,iEup)')
                idTbry(ieast,iEup)=varid
              CASE ('idTbry(isouth,iEup)')
                idTbry(isouth,iEup)=varid
              CASE ('idTbry(inorth,iEup)')
                idTbry(inorth,iEup)=varid

              CASE ('idTbry(iwest,iDet)')
                idTbry(iwest,iDet)=varid
              CASE ('idTbry(ieast,iDet)')
                idTbry(ieast,iDet)=varid
              CASE ('idTbry(isouth,iDet)')
                idTbry(isouth,iDet)=varid
              CASE ('idTbry(inorth,iDet)')
                idTbry(inorth,iDet)=varid

#ifdef IRON_LIMIT
              CASE ('idTbry(iwest,iFe)')
                idTbry(iwest,iFe)=varid
              CASE ('idTbry(ieast,iFe)')
                idTbry(ieast,iFe)=varid
              CASE ('idTbry(isouth,iFe)')
                idTbry(isouth,iFe)=varid
              CASE ('idTbry(inorth,iFe)')
                idTbry(inorth,iFe)=varid
#endif

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iNO3)')
                idRtrc(iNO3)=varid
              CASE ('idRtrc(iNH4)')
                idRtrc(iNH4)=varid
#endif
