/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2011 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices to the Nemuro ecosystem model            **
**  variables that are used in input and output NetCDF files.         **
**  The metadata nformation is read from "varinfo.dat".               **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

!
!  Model state biological tracers.
!

              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iSiOH)')
                idTvar(iSiOH)=varid
              CASE ('idTvar(iSphy)')
                idTvar(iSphy)=varid
              CASE ('idTvar(iLphy)')
                idTvar(iLphy)=varid
              CASE ('idTvar(iSzoo)')
                idTvar(iSzoo)=varid
              CASE ('idTvar(iLzoo)')
                idTvar(iLzoo)=varid
              CASE ('idTvar(iSDet)')
                idTvar(iSDet)=varid
              CASE ('idTvar(iopal)')
                idTvar(iopal)=varid
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
#ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
#endif

! Do it again for climatologies
              CASE ('idTclm(iNO3_)')
                idTclm(iNO3_)=varid
              CASE ('idTclm(iNH4_)')
                idTclm(iNH4_)=varid
              CASE ('idTclm(iSiOH)')
                idTclm(iSiOH)=varid
              CASE ('idTclm(iSphy)')
                idTclm(iSphy)=varid
              CASE ('idTclm(iLphy)')
                idTclm(iLphy)=varid
              CASE ('idTclm(iSzoo)')
                idTclm(iSzoo)=varid
              CASE ('idTclm(iLzoo)')
                idTclm(iLzoo)=varid
              CASE ('idTclm(iSDet)')
                idTclm(iSDet)=varid
              CASE ('idTclm(iopal)')
                idTclm(iopal)=varid
              CASE ('idTclm(iPO4_)')
                idTclm(iPO4_)=varid
#ifdef OXYGEN
              CASE ('idTclm(iOxyg)')
                idTclm(iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTclm(iTIC_)')
                idTclm(iTIC_)=varid
              CASE ('idTclm(iTAlk)')
                idTclm(iTAlk)=varid
#endif

!
!  Biological tracers open boundary conditions.
!

              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid
              CASE ('idTbry(iwest,iSiOH)')
                idTbry(iwest,iSiOH)=varid
              CASE ('idTbry(ieast,iSiOH)')
                idTbry(ieast,iSiOH)=varid
              CASE ('idTbry(isouth,iSiOH)')
                idTbry(isouth,iSiOH)=varid
              CASE ('idTbry(inorth,iSiOH)')
                idTbry(inorth,iSiOH)=varid
              CASE ('idTbry(iwest,iSphy)')
                idTbry(iwest,iSphy)=varid
              CASE ('idTbry(ieast,iSphy)')
                idTbry(ieast,iSphy)=varid
              CASE ('idTbry(isouth,iSphy)')
                idTbry(isouth,iSphy)=varid
              CASE ('idTbry(inorth,iSphy)')
                idTbry(inorth,iSphy)=varid
              CASE ('idTbry(iwest,iLphy)')
                idTbry(iwest,iLphy)=varid
              CASE ('idTbry(ieast,iLphy)')
                idTbry(ieast,iLphy)=varid
              CASE ('idTbry(isouth,iLphy)')
                idTbry(isouth,iLphy)=varid
              CASE ('idTbry(inorth,iLphy)')
                idTbry(inorth,iLphy)=varid
              CASE ('idTbry(iwest,iSzoo)')
                idTbry(iwest,iSzoo)=varid
              CASE ('idTbry(ieast,iSzoo)')
                idTbry(ieast,iSzoo)=varid
              CASE ('idTbry(isouth,iSzoo)')
                idTbry(isouth,iSzoo)=varid
              CASE ('idTbry(inorth,iSzoo)')
                idTbry(inorth,iSzoo)=varid
              CASE ('idTbry(iwest,iLzoo)')
                idTbry(iwest,iLzoo)=varid
              CASE ('idTbry(ieast,iLzoo)')
                idTbry(ieast,iLzoo)=varid
              CASE ('idTbry(isouth,iLzoo)')
                idTbry(isouth,iLzoo)=varid
              CASE ('idTbry(inorth,iLzoo)')
                idTbry(inorth,iLzoo)=varid
              CASE ('idTbry(iwest,iSDet)')
                idTbry(iwest,iSDet)=varid
              CASE ('idTbry(ieast,iSDet)')
                idTbry(ieast,iSDet)=varid
              CASE ('idTbry(isouth,iSDet)')
                idTbry(isouth,iSDet)=varid
              CASE ('idTbry(inorth,iSDet)')
                idTbry(inorth,iSDet)=varid
              CASE ('idTbry(iwest,iopal)')
                idTbry(iwest,iopal)=varid
              CASE ('idTbry(ieast,iopal)')
                idTbry(ieast,iopal)=varid
              CASE ('idTbry(isouth,iopal)')
                idTbry(isouth,iopal)=varid
              CASE ('idTbry(inorth,iopal)')
                idTbry(inorth,iopal)=varid
              CASE ('idTbry(iwest,iPO4_)')
                idTbry(iwest,iPO4_)=varid
              CASE ('idTbry(ieast,iPO4_)')
                idTbry(ieast,iPO4_)=varid
              CASE ('idTbry(isouth,iPO4_)')
                idTbry(isouth,iPO4_)=varid
              CASE ('idTbry(inorth,iPO4_)')
                idTbry(inorth,iPO4_)=varid
#ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid
              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid
#endif

#ifdef TS_PSOURCE

!
!  Biological tracers point Source/Sinks (river runoff).
!

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iSiOH)')
                idRtrc(iSiOH)=varid
              CASE ('idRtrc(iSphy)')
                idRtrc(iSphy)=varid
              CASE ('idRtrc(iLphy)')
                idRtrc(iLphy)=varid
              CASE ('idRtrc(iSzoo)')
                idRtrc(iSzoo)=varid
              CASE ('idRtrc(iLzoo)')
                idRtrc(iLzoo)=varid
              CASE ('idRtrc(iSDet)')
                idRtrc(iSDet)=varid
              CASE ('idRtrc(iopal)')
                idRtrc(iopal)=varid
              CASE ('idRtrc(iPO4_)')
                idRtrc(iPO4_)=varid
# ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
# endif
# ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
# endif
#endif
