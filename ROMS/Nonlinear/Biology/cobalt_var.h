/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2013 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices to the Cobalt ecosystem model            **
**  variables that are used in input and output NetCDF files.         **
**  The metadata nformation is read from "varinfo.dat".               **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

!TODO : what about diag variables
!
!  Model state biological tracers.
!

               ! Nitrogen Dynamics
               CASE ('idTvar(insm)')
                 idTvar(insm)=varid
               CASE ('idTvar(inlg)')
                 idTvar(inlg)=varid
               CASE ('idTvar(indi)')
                 idTvar(indi)=varid
               CASE ('idTvar(insmz)')
                 idTvar(insmz)=varid
               CASE ('idTvar(inmdz)')
                 idTvar(inmdz)=varid
               CASE ('idTvar(inlgz)')
                 idTvar(inlgz)=varid
               CASE ('idTvar(ildon)')
                 idTvar(ildon)=varid
               CASE ('idTvar(isldon)')
                 idTvar(isldon)=varid
               CASE ('idTvar(isrdon)')
                 idTvar(isrdon)=varid
               CASE ('idTvar(inbact)')
                 idTvar(inbact)=varid
               CASE ('idTvar(inh4)')
                 idTvar(inh4)=varid
               CASE ('idTvar(ino3)')
                 idTvar(ino3)=varid
               CASE ('idTvar(indet)')
                 idTvar(indet)=varid
#ifdef COBALT_MINERALS
               ! Biogenic Minerals and Lithogenic Materials
               CASE ('idTvar(isio4)')
                 idTvar(isio4)=varid
               CASE ('idTvar(isilg)')
                 idTvar(isilg)=varid
               CASE ('idTvar(isidet)')
                 idTvar(isidet)=varid
               CASE ('idTvar(icadet_calc)')
                 idTvar(icadet_calc)=varid
               CASE ('idTvar(icadet_arag)')
                 idTvar(icadet_arag)=varid
               CASE ('idTvar(ilith)')
                 idTvar(ilith)=varid
               CASE ('idTvar(ilithdet)')
                 idTvar(ilithdet)=varid
#endif
#ifdef COBALT_PHOSPHORUS
               ! Phosporus Dynamics
               CASE ('idTvar(ildop)')
                 idTvar(ildop)=varid
               CASE ('idTvar(isldop)')
                 idTvar(isldop)=varid
               CASE ('idTvar(isrdop)')
                 idTvar(isrdop)=varid
               CASE ('idTvar(ipo4)')
                 idTvar(ipo4)=varid
               CASE ('idTvar(ipdet)')
                 idTvar(ipdet)=varid
#else
               CASE ('idTvar(ipo4)')
                 idTvar(ipo4)=varid
#endif
#ifdef COBALT_IRON
               ! Iron Dynamics
               CASE ('idTvar(ifesm)')
                 idTvar(ifesm)=varid
               CASE ('idTvar(ifedi)')
                 idTvar(ifedi)=varid
               CASE ('idTvar(ifelg)')
                 idTvar(ifelg)=varid
               CASE ('idTvar(ifed)')
                 idTvar(ifed)=varid
               CASE ('idTvar(ifedet)')
                 idTvar(ifedet)=varid
#endif
#ifdef COBALT_CARBON
               ! Oxygen, Carbon and Alkalinity
               CASE ('idTvar(io2)')
                 idTvar(io2)=varid
               CASE ('idTvar(idic)')
                 idTvar(idic)=varid
               CASE ('idTvar(ialk)')
                 idTvar(ialk)=varid
#endif
#ifdef COASTDIAT
               CASE ('idTvar(inmd)')
                 idTvar(inmd)=varid
               CASE ('idTvar(isimd)')
                 idTvar(isimd)=varid
               CASE ('idTvar(ifemd)')
                 idTvar(ifemd)=varid
#endif

! ----- Other BGC variables
               CASE ('iDobgc(iochl)')
                 iDobgc(iochl)=varid
               CASE ('iDobgc(ioirr_mem)')
                 iDobgc(ioirr_mem)=varid
               CASE ('iDobgc(iohtotal)')
                 iDobgc(iohtotal)=varid
               CASE ('iDobgc(ioco3_ion)')
                 iDobgc(ioco3_ion)=varid
               CASE ('iDobgc(iomu_mem_sm)')
                 iDobgc(iomu_mem_sm)=varid
               CASE ('iDobgc(iomu_mem_di)')
                 iDobgc(iomu_mem_di)=varid
               CASE ('iDobgc(iomu_mem_lg)')
                 iDobgc(iomu_mem_lg)=varid
#ifdef COASTDIAT
               CASE ('iDobgc(iomu_mem_md)')
                 iDobgc(iomu_mem_md)=varid
#endif

!RD : diag variables
!               CASE ('idTvar(icased)')
!                 idTvar(icased)=varid
!               CASE ('idTvar(ichl)')
!                 idTvar(ichl)=varid
!               CASE ('idTvar(ico3_ion)')
!                 idTvar(ico3_ion)=varid
!               CASE ('idTvar(icadet_arag_btf)')
!                 idTvar(icadet_arag_btf)=varid
!               CASE ('idTvar(icadet_calc_btf)')
!                 idTvar(icadet_calc_btf)=varid
!               CASE ('idTvar(indet_btf)')
!                 idTvar(indet_btf)=varid
!               CASE ('idTvar(ipdet_btf)')
!                 idTvar(ipdet_btf)=varid
!               CASE ('idTvar(isidet_btf)')
!                 idTvar(isidet_btf)=varid
!               CASE ('idTvar(ihtotal)')
!                 idTvar(ihtotal)=varid
!               CASE ('idTvar(iirr_mem)')
!                 idTvar(iirr_mem)=varid

! Do it again for climatologies

               ! Nitrogen Dynamics
               CASE ('idTclm(insm)')
                 idTclm(insm)=varid
               CASE ('idTclm(inlg)')
                 idTclm(inlg)=varid
               CASE ('idTclm(indi)')
                 idTclm(indi)=varid
               CASE ('idTclm(insmz)')
                 idTclm(insmz)=varid
               CASE ('idTclm(inmdz)')
                 idTclm(inmdz)=varid
               CASE ('idTclm(inlgz)')
                 idTclm(inlgz)=varid
               CASE ('idTclm(ildon)')
                 idTclm(ildon)=varid
               CASE ('idTclm(isldon)')
                 idTclm(isldon)=varid
               CASE ('idTclm(isrdon)')
                 idTclm(isrdon)=varid
               CASE ('idTclm(inbact)')
                 idTclm(inbact)=varid
               CASE ('idTclm(inh4)')
                 idTclm(inh4)=varid
               CASE ('idTclm(ino3)')
                 idTclm(ino3)=varid
               CASE ('idTclm(indet)')
                 idTclm(indet)=varid
#ifdef COBALT_MINERALS
               ! Biogenic Minerals and Lithogenic Materials
               CASE ('idTclm(isio4)')
                 idTclm(isio4)=varid
               CASE ('idTclm(isilg)')
                 idTclm(isilg)=varid
               CASE ('idTclm(isidet)')
                 idTclm(isidet)=varid
               CASE ('idTclm(icadet_calc)')
                 idTclm(icadet_calc)=varid
               CASE ('idTclm(icadet_arag)')
                 idTclm(icadet_arag)=varid
               CASE ('idTclm(ilith)')
                 idTclm(ilith)=varid
               CASE ('idTclm(ilithdet)')
                 idTclm(ilithdet)=varid
#endif
#ifdef COBALT_PHOSPHORUS
               ! Phosporus Dynamics
               CASE ('idTclm(ildop)')
                 idTclm(ildop)=varid
               CASE ('idTclm(isldop)')
                 idTclm(isldop)=varid
               CASE ('idTclm(isrdop)')
                 idTclm(isrdop)=varid
               CASE ('idTclm(ipo4)')
                 idTclm(ipo4)=varid
               CASE ('idTclm(ipdet)')
                 idTclm(ipdet)=varid
#endif
#ifdef COBALT_IRON
               ! Iron Dynamics
               CASE ('idTclm(ifesm)')
                 idTclm(ifesm)=varid
               CASE ('idTclm(ifedi)')
                 idTclm(ifedi)=varid
               CASE ('idTclm(ifelg)')
                 idTclm(ifelg)=varid
               CASE ('idTclm(ifed)')
                 idTclm(ifed)=varid
               CASE ('idTclm(ifedet)')
                 idTclm(ifedet)=varid
#endif
#ifdef COBALT_CARBON
               ! Oxygen, Carbon and Alkalinity
               CASE ('idTclm(io2)')
                 idTclm(io2)=varid
               CASE ('idTclm(idic)')
                 idTclm(idic)=varid
               CASE ('idTclm(ialk)')
                 idTclm(ialk)=varid
#endif
#ifdef COASTDIAT
               CASE ('idTclm(inmd)')
                 idTclm(inmd)=varid
               CASE ('idTclm(isimd)')
                 idTclm(isimd)=varid
               CASE ('idTclm(ifemd)')
                 idTclm(ifemd)=varid
#endif

! ----- Other BGC variables


!RD : diag variables
!               CASE ('idTclm(icased)')
!                 idTclm(icased)=varid
!               CASE ('idTclm(ichl)')
!                 idTclm(ichl)=varid
!               CASE ('idTclm(ico3_ion)')
!                 idTclm(ico3_ion)=varid
!               CASE ('idTclm(icadet_arag_btf)')
!                 idTclm(icadet_arag_btf)=varid
!               CASE ('idTclm(icadet_calc_btf)')
!                 idTclm(icadet_calc_btf)=varid
!               CASE ('idTclm(indet_btf)')
!                 idTclm(indet_btf)=varid
!               CASE ('idTclm(ipdet_btf)')
!                 idTclm(ipdet_btf)=varid
!               CASE ('idTclm(isidet_btf)')
!                 idTclm(isidet_btf)=varid
!               CASE ('idTclm(ihtotal)')
!                 idTclm(ihtotal)=varid
!               CASE ('idTclm(iirr_mem)')
!                 idTclm(iirr_mem)=varid

!
!  Biological tracers open boundary conditions.
!

               ! Nitrogen Dynamics
               CASE ('idTbry(iwest,insm)')
                 idTbry(iwest,insm)=varid
               CASE ('idTbry(ieast,insm)')
                 idTbry(ieast,insm)=varid
               CASE ('idTbry(isouth,insm)')
                 idTbry(isouth,insm)=varid
               CASE ('idTbry(inorth,insm)')
                 idTbry(inorth,insm)=varid
               CASE ('idTbry(iwest,inlg)')
                 idTbry(iwest,inlg)=varid
               CASE ('idTbry(ieast,inlg)')
                 idTbry(ieast,inlg)=varid
               CASE ('idTbry(isouth,inlg)')
                 idTbry(isouth,inlg)=varid
               CASE ('idTbry(inorth,inlg)')
                 idTbry(inorth,inlg)=varid
               CASE ('idTbry(iwest,indi)')
                 idTbry(iwest,indi)=varid
               CASE ('idTbry(ieast,indi)')
                 idTbry(ieast,indi)=varid
               CASE ('idTbry(isouth,indi)')
                 idTbry(isouth,indi)=varid
               CASE ('idTbry(inorth,indi)')
                 idTbry(inorth,indi)=varid
               CASE ('idTbry(iwest,insmz)')
                 idTbry(iwest,insmz)=varid
               CASE ('idTbry(ieast,insmz)')
                 idTbry(ieast,insmz)=varid
               CASE ('idTbry(isouth,insmz)')
                 idTbry(isouth,insmz)=varid
               CASE ('idTbry(inorth,insmz)')
                 idTbry(inorth,insmz)=varid
               CASE ('idTbry(iwest,inmdz)')
                 idTbry(iwest,inmdz)=varid
               CASE ('idTbry(ieast,inmdz)')
                 idTbry(ieast,inmdz)=varid
               CASE ('idTbry(isouth,inmdz)')
                 idTbry(isouth,inmdz)=varid
               CASE ('idTbry(inorth,inmdz)')
                 idTbry(inorth,inmdz)=varid
               CASE ('idTbry(iwest,inlgz)')
                 idTbry(iwest,inlgz)=varid
               CASE ('idTbry(ieast,inlgz)')
                 idTbry(ieast,inlgz)=varid
               CASE ('idTbry(isouth,inlgz)')
                 idTbry(isouth,inlgz)=varid
               CASE ('idTbry(inorth,inlgz)')
                 idTbry(inorth,inlgz)=varid
               CASE ('idTbry(iwest,ildon)')
                 idTbry(iwest,ildon)=varid
               CASE ('idTbry(ieast,ildon)')
                 idTbry(ieast,ildon)=varid
               CASE ('idTbry(isouth,ildon)')
                 idTbry(isouth,ildon)=varid
               CASE ('idTbry(inorth,ildon)')
                 idTbry(inorth,ildon)=varid
               CASE ('idTbry(iwest,isldon)')
                 idTbry(iwest,isldon)=varid
               CASE ('idTbry(ieast,isldon)')
                 idTbry(ieast,isldon)=varid
               CASE ('idTbry(isouth,isldon)')
                 idTbry(isouth,isldon)=varid
               CASE ('idTbry(inorth,isldon)')
                 idTbry(inorth,isldon)=varid
               CASE ('idTbry(iwest,isrdon)')
                 idTbry(iwest,isrdon)=varid
               CASE ('idTbry(ieast,isrdon)')
                 idTbry(ieast,isrdon)=varid
               CASE ('idTbry(isouth,isrdon)')
                 idTbry(isouth,isrdon)=varid
               CASE ('idTbry(inorth,isrdon)')
                 idTbry(inorth,isrdon)=varid
               CASE ('idTbry(iwest,inbact)')
                 idTbry(iwest,inbact)=varid
               CASE ('idTbry(ieast,inbact)')
                 idTbry(ieast,inbact)=varid
               CASE ('idTbry(isouth,inbact)')
                 idTbry(isouth,inbact)=varid
               CASE ('idTbry(inorth,inbact)')
                 idTbry(inorth,inbact)=varid
               CASE ('idTbry(iwest,inh4)')
                 idTbry(iwest,inh4)=varid
               CASE ('idTbry(ieast,inh4)')
                 idTbry(ieast,inh4)=varid
               CASE ('idTbry(isouth,inh4)')
                 idTbry(isouth,inh4)=varid
               CASE ('idTbry(inorth,inh4)')
                 idTbry(inorth,inh4)=varid
               CASE ('idTbry(iwest,ino3)')
                 idTbry(iwest,ino3)=varid
               CASE ('idTbry(ieast,ino3)')
                 idTbry(ieast,ino3)=varid
               CASE ('idTbry(isouth,ino3)')
                 idTbry(isouth,ino3)=varid
               CASE ('idTbry(inorth,ino3)')
                 idTbry(inorth,ino3)=varid
               CASE ('idTbry(iwest,indet)')
                 idTbry(iwest,indet)=varid
               CASE ('idTbry(ieast,indet)')
                 idTbry(ieast,indet)=varid
               CASE ('idTbry(isouth,indet)')
                 idTbry(isouth,indet)=varid
               CASE ('idTbry(inorth,indet)')
                 idTbry(inorth,indet)=varid
#ifdef COBALT_MINERALS
               ! Biogenic Minerals and Lithogenic Materials
               CASE ('idTbry(iwest,isio4)')
                 idTbry(iwest,isio4)=varid
               CASE ('idTbry(ieast,isio4)')
                 idTbry(ieast,isio4)=varid
               CASE ('idTbry(isouth,isio4)')
                 idTbry(isouth,isio4)=varid
               CASE ('idTbry(inorth,isio4)')
                 idTbry(inorth,isio4)=varid
               CASE ('idTbry(iwest,isilg)')
                 idTbry(iwest,isilg)=varid
               CASE ('idTbry(ieast,isilg)')
                 idTbry(ieast,isilg)=varid
               CASE ('idTbry(isouth,isilg)')
                 idTbry(isouth,isilg)=varid
               CASE ('idTbry(inorth,isilg)')
                 idTbry(inorth,isilg)=varid
               CASE ('idTbry(iwest,isidet)')
                 idTbry(iwest,isidet)=varid
               CASE ('idTbry(ieast,isidet)')
                 idTbry(ieast,isidet)=varid
               CASE ('idTbry(isouth,isidet)')
                 idTbry(isouth,isidet)=varid
               CASE ('idTbry(inorth,isidet)')
                 idTbry(inorth,isidet)=varid
               CASE ('idTbry(iwest,icadet_calc)')
                 idTbry(iwest,icadet_calc)=varid
               CASE ('idTbry(ieast,icadet_calc)')
                 idTbry(ieast,icadet_calc)=varid
               CASE ('idTbry(isouth,icadet_calc)')
                 idTbry(isouth,icadet_calc)=varid
               CASE ('idTbry(inorth,icadet_calc)')
                 idTbry(inorth,icadet_calc)=varid
               CASE ('idTbry(iwest,icadet_arag)')
                 idTbry(iwest,icadet_arag)=varid
               CASE ('idTbry(ieast,icadet_arag)')
                 idTbry(ieast,icadet_arag)=varid
               CASE ('idTbry(isouth,icadet_arag)')
                 idTbry(isouth,icadet_arag)=varid
               CASE ('idTbry(inorth,icadet_arag)')
                 idTbry(inorth,icadet_arag)=varid
               CASE ('idTbry(iwest,ilith)')
                 idTbry(iwest,ilith)=varid
               CASE ('idTbry(ieast,ilith)')
                 idTbry(ieast,ilith)=varid
               CASE ('idTbry(isouth,ilith)')
                 idTbry(isouth,ilith)=varid
               CASE ('idTbry(inorth,ilith)')
                 idTbry(inorth,ilith)=varid
               CASE ('idTbry(iwest,ilithdet)')
                 idTbry(iwest,ilithdet)=varid
               CASE ('idTbry(ieast,ilithdet)')
                 idTbry(ieast,ilithdet)=varid
               CASE ('idTbry(isouth,ilithdet)')
                 idTbry(isouth,ilithdet)=varid
               CASE ('idTbry(inorth,ilithdet)')
                 idTbry(inorth,ilithdet)=varid
#endif
#ifdef COBALT_PHOSPHORUS
               ! Phosporus Dynamics
               CASE ('idTbry(iwest,ildop)')
                 idTbry(iwest,ildop)=varid
               CASE ('idTbry(ieast,ildop)')
                 idTbry(ieast,ildop)=varid
               CASE ('idTbry(isouth,ildop)')
                 idTbry(isouth,ildop)=varid
               CASE ('idTbry(inorth,ildop)')
                 idTbry(inorth,ildop)=varid
               CASE ('idTbry(iwest,isldop)')
                 idTbry(iwest,isldop)=varid
               CASE ('idTbry(ieast,isldop)')
                 idTbry(ieast,isldop)=varid
               CASE ('idTbry(isouth,isldop)')
                 idTbry(isouth,isldop)=varid
               CASE ('idTbry(inorth,isldop)')
                 idTbry(inorth,isldop)=varid
               CASE ('idTbry(iwest,isrdop)')
                 idTbry(iwest,isrdop)=varid
               CASE ('idTbry(ieast,isrdop)')
                 idTbry(ieast,isrdop)=varid
               CASE ('idTbry(isouth,isrdop)')
                 idTbry(isouth,isrdop)=varid
               CASE ('idTbry(inorth,isrdop)')
                 idTbry(inorth,isrdop)=varid
               CASE ('idTbry(iwest,ipo4)')
                 idTbry(iwest,ipo4)=varid
               CASE ('idTbry(ieast,ipo4)')
                 idTbry(ieast,ipo4)=varid
               CASE ('idTbry(isouth,ipo4)')
                 idTbry(isouth,ipo4)=varid
               CASE ('idTbry(inorth,ipo4)')
                 idTbry(inorth,ipo4)=varid
               CASE ('idTbry(iwest,ipdet)')
                 idTbry(iwest,ipdet)=varid
               CASE ('idTbry(ieast,ipdet)')
                 idTbry(ieast,ipdet)=varid
               CASE ('idTbry(isouth,ipdet)')
                 idTbry(isouth,ipdet)=varid
               CASE ('idTbry(inorth,ipdet)')
                 idTbry(inorth,ipdet)=varid
#else
               CASE ('idTbry(iwest,ipo4)')
                 idTbry(iwest,ipo4)=varid
               CASE ('idTbry(ieast,ipo4)')
                 idTbry(ieast,ipo4)=varid
               CASE ('idTbry(isouth,ipo4)')
                 idTbry(isouth,ipo4)=varid
               CASE ('idTbry(inorth,ipo4)')
                 idTbry(inorth,ipo4)=varid
#endif
#ifdef COBALT_IRON
               ! Iron Dynamics
               CASE ('idTbry(iwest,ifesm)')
                 idTbry(iwest,ifesm)=varid
               CASE ('idTbry(ieast,ifesm)')
                 idTbry(ieast,ifesm)=varid
               CASE ('idTbry(isouth,ifesm)')
                 idTbry(isouth,ifesm)=varid
               CASE ('idTbry(inorth,ifesm)')
                 idTbry(inorth,ifesm)=varid
               CASE ('idTbry(iwest,ifedi)')
                 idTbry(iwest,ifedi)=varid
               CASE ('idTbry(ieast,ifedi)')
                 idTbry(ieast,ifedi)=varid
               CASE ('idTbry(isouth,ifedi)')
                 idTbry(isouth,ifedi)=varid
               CASE ('idTbry(inorth,ifedi)')
                 idTbry(inorth,ifedi)=varid
               CASE ('idTbry(iwest,ifelg)')
                 idTbry(iwest,ifelg)=varid
               CASE ('idTbry(ieast,ifelg)')
                 idTbry(ieast,ifelg)=varid
               CASE ('idTbry(isouth,ifelg)')
                 idTbry(isouth,ifelg)=varid
               CASE ('idTbry(inorth,ifelg)')
                 idTbry(inorth,ifelg)=varid
               CASE ('idTbry(iwest,ifed)')
                 idTbry(iwest,ifed)=varid
               CASE ('idTbry(ieast,ifed)')
                 idTbry(ieast,ifed)=varid
               CASE ('idTbry(isouth,ifed)')
                 idTbry(isouth,ifed)=varid
               CASE ('idTbry(inorth,ifed)')
                 idTbry(inorth,ifed)=varid
               CASE ('idTbry(iwest,ifedet)')
                 idTbry(iwest,ifedet)=varid
               CASE ('idTbry(ieast,ifedet)')
                 idTbry(ieast,ifedet)=varid
               CASE ('idTbry(isouth,ifedet)')
                 idTbry(isouth,ifedet)=varid
               CASE ('idTbry(inorth,ifedet)')
                 idTbry(inorth,ifedet)=varid
#endif
#ifdef COBALT_CARBON
               ! Oxygen, Carbon and Alkalinity
               CASE ('idTbry(iwest,io2)')
                 idTbry(iwest,io2)=varid
               CASE ('idTbry(ieast,io2)')
                 idTbry(ieast,io2)=varid
               CASE ('idTbry(isouth,io2)')
                 idTbry(isouth,io2)=varid
               CASE ('idTbry(inorth,io2)')
                 idTbry(inorth,io2)=varid
               CASE ('idTbry(iwest,idic)')
                 idTbry(iwest,idic)=varid
               CASE ('idTbry(ieast,idic)')
                 idTbry(ieast,idic)=varid
               CASE ('idTbry(isouth,idic)')
                 idTbry(isouth,idic)=varid
               CASE ('idTbry(inorth,idic)')
                 idTbry(inorth,idic)=varid
               CASE ('idTbry(iwest,ialk)')
                 idTbry(iwest,ialk)=varid
               CASE ('idTbry(ieast,ialk)')
                 idTbry(ieast,ialk)=varid
               CASE ('idTbry(isouth,ialk)')
                 idTbry(isouth,ialk)=varid
               CASE ('idTbry(inorth,ialk)')
                 idTbry(inorth,ialk)=varid
#endif
#ifdef COASTDIAT
               CASE ('idTbry(iwest,inmd)')
                 idTbry(iwest,inmd)=varid
               CASE ('idTbry(ieast,inmd)')
                 idTbry(ieast,inmd)=varid
               CASE ('idTbry(isouth,inmd)')
                 idTbry(isouth,inmd)=varid
               CASE ('idTbry(inorth,inmd)')
                 idTbry(inorth,inmd)=varid
               CASE ('idTbry(iwest,isimd)')
                 idTbry(iwest,isimd)=varid
               CASE ('idTbry(ieast,isimd)')
                 idTbry(ieast,isimd)=varid
               CASE ('idTbry(isouth,isimd)')
                 idTbry(isouth,isimd)=varid
               CASE ('idTbry(inorth,isimd)')
                 idTbry(inorth,isimd)=varid
               CASE ('idTbry(iwest,ifemd)')
                 idTbry(iwest,ifemd)=varid
               CASE ('idTbry(ieast,ifemd)')
                 idTbry(ieast,ifemd)=varid
               CASE ('idTbry(isouth,ifemd)')
                 idTbry(isouth,ifemd)=varid
               CASE ('idTbry(inorth,ifemd)')
                 idTbry(inorth,ifemd)=varid
#endif

! ----- Other BGC variables


!               CASE ('idTbry(iwest,icased)')
!                idTbry(iwest,icased)=varid
!              CASE ('idTbry(ieast,icased)')
!                idTbry(ieast,icased)=varid
!              CASE ('idTbry(isouth,icased)')
!                idTbry(isouth,icased)=varid
!              CASE ('idTbry(inorth,icased)')
!                idTbry(inorth,icased)=varid
!              CASE ('idTbry(iwest,ichl)')
!                idTbry(iwest,ichl)=varid
!              CASE ('idTbry(ieast,ichl)')
!                idTbry(ieast,ichl)=varid
!              CASE ('idTbry(isouth,ichl)')
!                idTbry(isouth,ichl)=varid
!              CASE ('idTbry(inorth,ichl)')
!                idTbry(inorth,ichl)=varid
!              CASE ('idTbry(iwest,ico3_ion)')
!                idTbry(iwest,ico3_ion)=varid
!              CASE ('idTbry(ieast,ico3_ion)')
!                idTbry(ieast,ico3_ion)=varid
!              CASE ('idTbry(isouth,ico3_ion)')
!                idTbry(isouth,ico3_ion)=varid
!              CASE ('idTbry(inorth,ico3_ion)')
!                idTbry(inorth,ico3_ion)=varid
!              CASE ('idTbry(iwest,icadet_arag_btf)')
!                idTbry(iwest,icadet_arag_btf)=varid
!              CASE ('idTbry(ieast,icadet_arag_btf)')
!                idTbry(ieast,icadet_arag_btf)=varid
!              CASE ('idTbry(isouth,icadet_arag_btf)')
!                idTbry(isouth,icadet_arag_btf)=varid
!              CASE ('idTbry(inorth,icadet_arag_btf)')
!                idTbry(inorth,icadet_arag_btf)=varid
!              CASE ('idTbry(iwest,icadet_calc_btf)')
!                idTbry(iwest,icadet_calc_btf)=varid
!              CASE ('idTbry(ieast,icadet_calc_btf)')
!                idTbry(ieast,icadet_calc_btf)=varid
!              CASE ('idTbry(isouth,icadet_calc_btf)')
!                idTbry(isouth,icadet_calc_btf)=varid
!              CASE ('idTbry(inorth,icadet_calc_btf)')
!                idTbry(inorth,icadet_calc_btf)=varid
!              CASE ('idTbry(iwest,indet_btf)')
!                idTbry(iwest,indet_btf)=varid
!              CASE ('idTbry(ieast,indet_btf)')
!                idTbry(ieast,indet_btf)=varid
!              CASE ('idTbry(isouth,indet_btf)')
!                idTbry(isouth,indet_btf)=varid
!              CASE ('idTbry(inorth,indet_btf)')
!                idTbry(inorth,indet_btf)=varid
!              CASE ('idTbry(iwest,ipdet_btf)')
!                idTbry(iwest,ipdet_btf)=varid
!              CASE ('idTbry(ieast,ipdet_btf)')
!                idTbry(ieast,ipdet_btf)=varid
!              CASE ('idTbry(isouth,ipdet_btf)')
!                idTbry(isouth,ipdet_btf)=varid
!              CASE ('idTbry(inorth,ipdet_btf)')
!                idTbry(inorth,ipdet_btf)=varid
!              CASE ('idTbry(iwest,isidet_btf)')
!                idTbry(iwest,isidet_btf)=varid
!              CASE ('idTbry(ieast,isidet_btf)')
!                idTbry(ieast,isidet_btf)=varid
!              CASE ('idTbry(isouth,isidet_btf)')
!                idTbry(isouth,isidet_btf)=varid
!              CASE ('idTbry(inorth,isidet_btf)')
!                idTbry(inorth,isidet_btf)=varid
!              CASE ('idTbry(iwest,ihtotal)')
!                idTbry(iwest,ihtotal)=varid
!              CASE ('idTbry(ieast,ihtotal)')
!                idTbry(ieast,ihtotal)=varid
!              CASE ('idTbry(isouth,ihtotal)')
!                idTbry(isouth,ihtotal)=varid
!              CASE ('idTbry(inorth,ihtotal)')
!                idTbry(inorth,ihtotal)=varid
!              CASE ('idTbry(iwest,iirr_mem)')
!                idTbry(iwest,iirr_mem)=varid
!              CASE ('idTbry(ieast,iirr_mem)')
!                idTbry(ieast,iirr_mem)=varid
!              CASE ('idTbry(isouth,iirr_mem)')
!                idTbry(isouth,iirr_mem)=varid
!              CASE ('idTbry(inorth,iirr_mem)')
!                idTbry(inorth,iirr_mem)=varid

!
!  Biological tracers point Source/Sinks (river runoff).
!

               ! Nitrogen Dynamics
               CASE ('idRtrc(insm)')
                 idRtrc(insm)=varid
               CASE ('idRtrc(inlg)')
                 idRtrc(inlg)=varid
               CASE ('idRtrc(indi)')
                 idRtrc(indi)=varid
               CASE ('idRtrc(insmz)')
                 idRtrc(insmz)=varid
               CASE ('idRtrc(inmdz)')
                 idRtrc(inmdz)=varid
               CASE ('idRtrc(inlgz)')
                 idRtrc(inlgz)=varid
               CASE ('idRtrc(ildon)')
                 idRtrc(ildon)=varid
               CASE ('idRtrc(isldon)')
                 idRtrc(isldon)=varid
               CASE ('idRtrc(isrdon)')
                 idRtrc(isrdon)=varid
               CASE ('idRtrc(inbact)')
                 idRtrc(inbact)=varid
               CASE ('idRtrc(inh4)')
                 idRtrc(inh4)=varid
               CASE ('idRtrc(ino3)')
                 idRtrc(ino3)=varid
               CASE ('idRtrc(indet)')
                 idRtrc(indet)=varid
# ifdef COBALT_MINERALS
               ! Biogenic Minerals and Lithogenic Materials
               CASE ('idRtrc(isio4)')
                 idRtrc(isio4)=varid
               CASE ('idRtrc(isilg)')
                 idRtrc(isilg)=varid
               CASE ('idRtrc(isidet)')
                 idRtrc(isidet)=varid
               CASE ('idRtrc(icadet_calc)')
                 idRtrc(icadet_calc)=varid
               CASE ('idRtrc(icadet_arag)')
                 idRtrc(icadet_arag)=varid
               CASE ('idRtrc(ilith)')
                 idRtrc(ilith)=varid
               CASE ('idRtrc(ilithdet)')
                 idRtrc(ilithdet)=varid
# endif
# ifdef COBALT_PHOSPHORUS
               ! Phosporus Dynamics
               CASE ('idRtrc(ildop)')
                 idRtrc(ildop)=varid
               CASE ('idRtrc(isldop)')
                 idRtrc(isldop)=varid
               CASE ('idRtrc(isrdop)')
                 idRtrc(isrdop)=varid
               CASE ('idRtrc(ipo4)')
                 idRtrc(ipo4)=varid
               CASE ('idRtrc(ipdet)')
                 idRtrc(ipdet)=varid
# else
               CASE ('idRtrc(ipo4)')
                 idRtrc(ipo4)=varid
# endif
# ifdef COBALT_IRON
               ! Iron Dynamics
               CASE ('idRtrc(ifesm)')
                 idRtrc(ifesm)=varid
               CASE ('idRtrc(ifedi)')
                 idRtrc(ifedi)=varid
               CASE ('idRtrc(ifelg)')
                 idRtrc(ifelg)=varid
               CASE ('idRtrc(ifed)')
                 idRtrc(ifed)=varid
               CASE ('idRtrc(ifedet)')
                 idRtrc(ifedet)=varid
# endif
# ifdef COBALT_CARBON
               ! Oxygen, Carbon and Alkalinity
               CASE ('idRtrc(io2)')
                 idRtrc(io2)=varid
               CASE ('idRtrc(idic)')
                 idRtrc(idic)=varid
               CASE ('idRtrc(ialk)')
                 idRtrc(ialk)=varid
# endif
# ifdef COASTDIAT
               CASE ('idRtrc(inmd)')
                 idRtrc(inmd)=varid
               CASE ('idRtrc(isimd)')
                 idRtrc(isimd)=varid
               CASE ('idRtrc(ifemd)')
                 idRtrc(ifemd)=varid
# endif

! river as external forcing
               CASE ('idriver_no3')
                 idriver_no3=varid
               CASE ('idriver_ldon')
                 idriver_ldon=varid
               CASE ('idriver_sldon')
                 idriver_sldon=varid
               CASE ('idriver_srdon')
                 idriver_srdon=varid
               CASE ('idriver_ndet')
                 idriver_ndet=varid
               CASE ('idriver_po4')
                 idriver_po4=varid
               CASE ('idriver_ldop')
                 idriver_ldop=varid
               CASE ('idriver_sldop')
                 idriver_sldop=varid
               CASE ('idriver_srdop')
                 idriver_srdop=varid
#ifdef COBALT_IRON
               CASE ('idriver_fed')
                 idriver_fed=varid
#endif
#ifdef DIAGNOSTICS_BIO
!------------------------------------------------------
!---  Biological tracers term diagnostics.
!------------------------------------------------------
              CASE ('iDbio2(icased)')
                iDbio2(icased)=varid
              CASE ('iDbio2(icadet_arag_btf)')
                iDbio2(icadet_arag_btf)=varid
              CASE ('iDbio2(icadet_calc_btf)')
                iDbio2(icadet_calc_btf)=varid
              CASE ('iDbio2(indet_btf)')
                iDbio2(indet_btf)=varid
              CASE ('iDbio2(ipdet_btf)')
                iDbio2(ipdet_btf)=varid
              CASE ('iDbio2(isidet_btf)')
                iDbio2(isidet_btf)=varid
              CASE ('iDbio2(ialk_btf)')
                iDbio2(ialk_btf)=varid
              CASE ('iDbio2(idic_btf)')
                iDbio2(idic_btf)=varid
              CASE ('iDbio2(ifed_btf)')
                iDbio2(ifed_btf)=varid
              CASE ('iDbio2(inh4_btf)')
                iDbio2(inh4_btf)=varid
              CASE ('iDbio2(ino3_btf)')
                iDbio2(ino3_btf)=varid
              CASE ('iDbio2(io2_btf)')
                iDbio2(io2_btf)=varid
              CASE ('iDbio2(ipo4_btf)')
                iDbio2(ipo4_btf)=varid
              CASE ('iDbio2(isio4_btf)')
                iDbio2(isio4_btf)=varid
              CASE ('iDbio2(imxl_depth)')
                iDbio2(imxl_depth)=varid
              CASE ('iDbio2(imxl_level)')
                iDbio2(imxl_level)=varid

!              CASE ('iDbio2(ijprod_ndet_100)')
!                iDbio2(ijprod_ndet_100)=varid
!              CASE ('iDbio2(ijremin_ndet_100)')
!                iDbio2(ijremin_ndet_100)=varid

              CASE ('iDbio2(ialpha)')
                iDbio2(ialpha)=varid
              CASE ('iDbio2(ico2star)')
                iDbio2(ico2star)=varid
              CASE ('iDbio2(ipco2surf)')
                iDbio2(ipco2surf)=varid
              CASE ('iDbio2(ico2_flx)')
                iDbio2(ico2_flx)=varid
              CASE ('iDbio2(io2_flx)')
                iDbio2(io2_flx)=varid
              CASE ('iDbio2(iironsed_flx)')
                iDbio2(iironsed_flx)=varid
              CASE ('iDbio2(inpp_100)')
                iDbio2(inpp_100)=varid
              CASE ('iDbio2(imesozoo_200)')
                iDbio2(imesozoo_200)=varid

              CASE ('iDbio2(iprod_n_100_sm)')
                iDbio2(iprod_n_100_sm)=varid
              CASE ('iDbio2(iaggloss_n_100_sm)')
                iDbio2(iaggloss_n_100_sm)=varid
              CASE ('iDbio2(izloss_n_100_sm)')
                iDbio2(izloss_n_100_sm)=varid
              CASE ('iDbio2(iprod_n_100_lg)')
                iDbio2(iprod_n_100_lg)=varid
              CASE ('iDbio2(iaggloss_n_100_lg)')
                iDbio2(iaggloss_n_100_lg)=varid
              CASE ('iDbio2(izloss_n_100_lg)')
                iDbio2(izloss_n_100_lg)=varid
              CASE ('iDbio2(iprod_n_100_di)')
                iDbio2(iprod_n_100_di)=varid
              CASE ('iDbio2(iaggloss_n_100_di)')
                iDbio2(iaggloss_n_100_di)=varid
              CASE ('iDbio2(izloss_n_100_di)')
                iDbio2(izloss_n_100_di)=varid
              CASE ('iDbio2(iprod_n_100_smz)')
                iDbio2(iprod_n_100_smz)=varid
              CASE ('iDbio2(iingest_n_100_smz)')
                iDbio2(iingest_n_100_smz)=varid
              CASE ('iDbio2(izloss_n_100_smz)')
                iDbio2(izloss_n_100_smz)=varid
              CASE ('iDbio2(ihploss_n_100_smz)')
                iDbio2(ihploss_n_100_smz)=varid
              CASE ('iDbio2(iprod_ndet_100_smz)')
                iDbio2(iprod_ndet_100_smz)=varid
              CASE ('iDbio2(iprod_n_100_mdz)')
                iDbio2(iprod_n_100_mdz)=varid
              CASE ('iDbio2(iingest_n_100_mdz)')
                iDbio2(iingest_n_100_mdz)=varid
              CASE ('iDbio2(izloss_n_100_mdz)')
                iDbio2(izloss_n_100_mdz)=varid
              CASE ('iDbio2(ihploss_n_100_mdz)')
                iDbio2(ihploss_n_100_mdz)=varid
              CASE ('iDbio2(iprod_ndet_100_mdz)')
                iDbio2(iprod_ndet_100_mdz)=varid
              CASE ('iDbio2(iprod_n_100_lgz)')
                iDbio2(iprod_n_100_lgz)=varid
              CASE ('iDbio2(iingest_n_100_lgz)')
                iDbio2(iingest_n_100_lgz)=varid
              CASE ('iDbio2(izloss_n_100_lgz)')
                iDbio2(izloss_n_100_lgz)=varid
              CASE ('iDbio2(ihploss_n_100_lgz)')
                iDbio2(ihploss_n_100_lgz)=varid
              CASE ('iDbio2(iprod_ndet_100_lgz)')
                iDbio2(iprod_ndet_100_lgz)=varid
              CASE ('iDbio2(iprod_n_100_bact)')
                iDbio2(iprod_n_100_bact)=varid
              CASE ('iDbio2(izloss_n_100_bact)')
                iDbio2(izloss_n_100_bact)=varid
!RDIAGS
              CASE ('iDbio2(imesozooprod_200)')
                iDbio2(imesozooprod_200)=varid
              CASE ('iDbio2(iuptake_din_100)')
                iDbio2(iuptake_din_100)=varid
              CASE ('iDbio2(iuptake_no3_n2_100)')
                iDbio2(iuptake_no3_n2_100)=varid
              CASE ('iDbio2(iprod_mesozoo_100)')
                iDbio2(iprod_mesozoo_100)=varid
              CASE ('iDbio2(iz_ratio_100)')
                iDbio2(iz_ratio_100)=varid
              CASE ('iDbio2(ipe_ratio_100)')
                iDbio2(ipe_ratio_100)=varid
              CASE ('iDbio2(if_ratio_100)')
                iDbio2(if_ratio_100)=varid
              CASE ('iDbio2(iprod_don_100_smz)')
                iDbio2(iprod_don_100_smz)=varid
              CASE ('iDbio2(iprod_don_100_mdz)')
                iDbio2(iprod_don_100_mdz)=varid
              CASE ('iDbio2(iprod_don_100_lgz)')
                iDbio2(iprod_don_100_lgz)=varid
              CASE ('iDbio2(ijno3denit_wc_vint)')
                iDbio2(ijno3denit_wc_vint)=varid

#ifdef COASTDIAT
              CASE ('iDbio2(iprod_n_100_md)')
                iDbio2(iprod_n_100_md)=varid
              CASE ('iDbio2(iaggloss_n_100_md)')
                iDbio2(iaggloss_n_100_md)=varid
              CASE ('iDbio2(izloss_n_100_md)')
                iDbio2(izloss_n_100_md)=varid
#endif
              CASE ('iDbio3(ichl)')
                iDbio3(ichl)=varid
              CASE ('iDbio3(ico3_ion)')
                iDbio3(ico3_ion)=varid
              CASE ('iDbio3(ihtotal)')
                iDbio3(ihtotal)=varid
              CASE ('iDbio3(iirr_mem)')
                iDbio3(iirr_mem)=varid
              CASE ('iDbio3(iirr_mix)')
                iDbio3(iirr_mix)=varid
              CASE ('iDbio3(iirr_inst)')
                iDbio3(iirr_inst)=varid

!              CASE ('iDbio3(ijprod_cadet_arag)')
!                iDbio3(ijprod_cadet_arag)=varid
!              CASE ('iDbio3(ijprod_cadet_calc)')
!                iDbio3(ijprod_cadet_calc)=varid
!              CASE ('iDbio3(ijprod_fedet)')
!                iDbio3(ijprod_fedet)=varid
!              CASE ('iDbio3(ijprod_lithdet)')
!                iDbio3(ijprod_lithdet)=varid
!              CASE ('iDbio3(ijprod_ndet)')
!                iDbio3(ijprod_ndet)=varid
!              CASE ('iDbio3(ijprod_pdet)')
!                iDbio3(ijprod_pdet)=varid
!              CASE ('iDbio3(ijprod_sidet)')
!                iDbio3(ijprod_sidet)=varid
!              CASE ('iDbio3(ijdiss_cadet_arag)')
!                iDbio3(ijdiss_cadet_arag)=varid
!              CASE ('iDbio3(ijdiss_cadet_calc)')
!                iDbio3(ijdiss_cadet_calc)=varid
!              CASE ('iDbio3(ijremin_ndet)')
!                iDbio3(ijremin_ndet)=varid
!              CASE ('iDbio3(ijremin_pdet)')
!                iDbio3(ijremin_pdet)=varid
!              CASE ('iDbio3(ijremin_fedet)')
!                iDbio3(ijremin_fedet)=varid
!              CASE ('iDbio3(idet_jzloss_n)')
!                iDbio3(idet_jzloss_n)=varid
!              CASE ('iDbio3(idet_jzloss_p)')
!                iDbio3(idet_jzloss_p)=varid
!              CASE ('iDbio3(idet_jzloss_fe)')
!                iDbio3(idet_jzloss_fe)=varid
!              CASE ('iDbio3(idet_jhploss_n)')
!                iDbio3(idet_jhploss_n)=varid
!              CASE ('iDbio3(idet_jhploss_p)')
!                iDbio3(idet_jhploss_p)=varid
!              CASE ('iDbio3(idet_jhploss_fe)')
!                iDbio3(idet_jhploss_fe)=varid
              CASE ('iDbio3(ico3_sol_calc)')
                iDbio3(ico3_sol_calc)=varid
              CASE ('iDbio3(ico3_sol_arag)')
                iDbio3(ico3_sol_arag)=varid
              CASE ('iDbio3(ife_bulk_flx)')
                iDbio3(ife_bulk_flx)=varid


!              CASE ('iDbio3(indet_b4sink)')
!                iDbio3(indet_b4sink)=varid
!              CASE ('iDbio3(ipdet_b4sink)')
!                iDbio3(ipdet_b4sink)=varid
!              CASE ('iDbio3(ifedet_b4sink)')
!                iDbio3(ifedet_b4sink)=varid
!              CASE ('iDbio3(indet_afsink)')
!                iDbio3(indet_afsink)=varid
!              CASE ('iDbio3(ipdet_afsink)')
!                iDbio3(ipdet_afsink)=varid
!              CASE ('iDbio3(ifedet_afsink)')
!                iDbio3(ifedet_afsink)=varid
!              CASE ('iDbio3(indet_flx)')
!                iDbio3(indet_flx)=varid

              CASE ('iDbio3(iomega_cadet_calc)')
                iDbio3(iomega_cadet_calc)=varid
              CASE ('iDbio3(iomega_cadet_arag)')
                iDbio3(iomega_cadet_arag)=varid
              CASE ('iDbio3(iswdk)')
                iDbio3(iswdk)=varid
              CASE ('iDbio3(imu_mem_sm)')
                iDbio3(imu_mem_sm)=varid
              CASE ('iDbio3(imu_mem_di)')
                iDbio3(imu_mem_di)=varid
              CASE ('iDbio3(imu_mem_lg)')
                iDbio3(imu_mem_lg)=varid
              CASE ('iDbio3(iagg_lim_sm)')
                iDbio3(iagg_lim_sm)=varid
              CASE ('iDbio3(iagg_lim_di)')
                iDbio3(iagg_lim_di)=varid
              CASE ('iDbio3(iagg_lim_lg)')
                iDbio3(iagg_lim_lg)=varid

              CASE ('iDbio3(iaggloss_di)')
                iDbio3(iaggloss_di)=varid
              CASE ('iDbio3(iaggloss_sm)')
                iDbio3(iaggloss_sm)=varid
              CASE ('iDbio3(iaggloss_lg)')
                iDbio3(iaggloss_lg)=varid
              CASE ('iDbio3(ivirloss_di)')
                iDbio3(ivirloss_di)=varid
              CASE ('iDbio3(ivirloss_sm)')
                iDbio3(ivirloss_sm)=varid
              CASE ('iDbio3(ivirloss_lg)')
                iDbio3(ivirloss_lg)=varid
              CASE ('iDbio3(izloss_di)')
                iDbio3(izloss_di)=varid
              CASE ('iDbio3(izloss_sm)')
                iDbio3(izloss_sm)=varid
              CASE ('iDbio3(izloss_lg)')
                iDbio3(izloss_lg)=varid

              CASE ('iDbio3(idef_fe_sm)')
                iDbio3(idef_fe_sm)=varid
              CASE ('iDbio3(idef_fe_di)')
                iDbio3(idef_fe_di)=varid
              CASE ('iDbio3(idef_fe_lg)')
                iDbio3(idef_fe_lg)=varid
              CASE ('iDbio3(ifelim_sm)')
                iDbio3(ifelim_sm)=varid
              CASE ('iDbio3(ifelim_di)')
                iDbio3(ifelim_di)=varid
              CASE ('iDbio3(ifelim_lg)')
                iDbio3(ifelim_lg)=varid
              CASE ('iDbio3(ino3lim_sm)')
                iDbio3(ino3lim_sm)=varid
              CASE ('iDbio3(ino3lim_di)')
                iDbio3(ino3lim_di)=varid
              CASE ('iDbio3(ino3lim_lg)')
                iDbio3(ino3lim_lg)=varid
              CASE ('iDbio3(inh4lim_sm)')
                iDbio3(inh4lim_sm)=varid
              CASE ('iDbio3(inh4lim_di)')
                iDbio3(inh4lim_di)=varid
              CASE ('iDbio3(inh4lim_lg)')
                iDbio3(inh4lim_lg)=varid
              CASE ('iDbio3(ipo4lim_sm)')
                iDbio3(ipo4lim_sm)=varid
              CASE ('iDbio3(ipo4lim_di)')
                iDbio3(ipo4lim_di)=varid
              CASE ('iDbio3(ipo4lim_lg)')
                iDbio3(ipo4lim_lg)=varid
              CASE ('iDbio3(ichl_di)')
                iDbio3(ichl_di)=varid
              CASE ('iDbio3(iC_2_chl_di)')
                iDbio3(iC_2_chl_di)=varid
              CASE ('iDbio3(ichl_sm)')
                iDbio3(ichl_sm)=varid
              CASE ('iDbio3(iC_2_chl_sm)')
                iDbio3(iC_2_chl_sm)=varid
              CASE ('iDbio3(ichl_lg)')
                iDbio3(ichl_lg)=varid
              CASE ('iDbio3(iC_2_chl_lg)')
                iDbio3(iC_2_chl_lg)=varid
#ifdef COASTDIAT
              CASE ('iDbio3(imu_mem_md)')
                iDbio3(imu_mem_md)=varid
              CASE ('iDbio3(iagg_lim_md)')
                iDbio3(iagg_lim_md)=varid
              CASE ('iDbio3(iaggloss_md)')
                iDbio3(iaggloss_md)=varid
              CASE ('iDbio3(ivirloss_md)')
                iDbio3(ivirloss_md)=varid
              CASE ('iDbio3(izloss_md)')
                iDbio3(izloss_md)=varid
              CASE ('iDbio3(idef_fe_md)')
                iDbio3(idef_fe_md)=varid
              CASE ('iDbio3(ifelim_md)')
                iDbio3(ifelim_md)=varid
              CASE ('iDbio3(ino3lim_md)')
                iDbio3(ino3lim_md)=varid
              CASE ('iDbio3(inh4lim_md)')
                iDbio3(inh4lim_md)=varid
              CASE ('iDbio3(ipo4lim_md)')
                iDbio3(ipo4lim_md)=varid
              CASE ('iDbio3(ichl_md)')
                iDbio3(ichl_md)=varid
              CASE ('iDbio3(iC_2_chl_md)')
                iDbio3(iC_2_chl_md)=varid
#endif


#endif
#ifdef BENTHIC
              CASE ('idBeTvar(icased)')
               idBeTvar(icased)=varid
!               PRINT *, 'RD in mod_ncparam, icased=', icased
!               PRINT *, 'RD in mod_ncparam, idBeTvar(icased)=', idBeTvar(icased)
              CASE ('idBeTvar(icadet_arag_btf)')
                idBeTvar(icadet_arag_btf)=varid
!               PRINT *, 'RD in mod_ncparam, icadet_arag_btf=', icadet_arag_btf
!               PRINT *, 'RD in mod_ncparam, idBeTvar(icadet_arag_btf)=', idBeTvar(icadet_arag_btf)
              CASE ('idBeTvar(icadet_calc_btf)')
                idBeTvar(icadet_calc_btf)=varid
!               PRINT *, 'RD in mod_ncparam, icadet_calc_btf=', icadet_calc_btf
!               PRINT *, 'RD in mod_ncparam, idBeTvar(icadet_calc_btf)=', idBeTvar(icadet_calc_btf)
              CASE ('idBeTvar(indet_btf)')
               idBeTvar(indet_btf)=varid
!               PRINT *, 'RD in mod_ncparam, indet_btf=', indet_btf
!               PRINT *, 'RD in mod_ncparam, idBeTvar(indet_btf)=', idBeTvar(indet_btf)
              CASE ('idBeTvar(ipdet_btf)')
               idBeTvar(ipdet_btf)=varid
!               PRINT *, 'RD in mod_ncparam, ipdet_btf=', ipdet_btf
!               PRINT *, 'RD in mod_ncparam, idBeTvar(ipdet_btf)=', idBeTvar(ipdet_btf)
              CASE ('idBeTvar(isidet_btf)')
               idBeTvar(isidet_btf)=varid
!               PRINT *, 'RD in mod_ncparam, isidet_btf=', isidet_btf
!               PRINT *, 'RD in mod_ncparam, idBeTvar(isidet_btf)=', idBeTvar(isidet_btf)
#endif
