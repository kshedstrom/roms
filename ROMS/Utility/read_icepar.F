#include "cppdefs.h"
#ifdef ICE_MODEL
      SUBROUTINE read_IcePar (inp, out, Lwrite)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads and reports ice model input parameters.          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Lstr, Lval, Npts, Nval, i, ng, itrc, status

      integer :: decode_line, lenstr, load_i, load_l, load_r

      real(r8), dimension(200) :: Rval

      character (len=40) :: KeyWord
      character (len=80) :: line
      character (len=80), dimension(200) :: Cval

!
!-----------------------------------------------------------------------
!  Read in ice model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lice') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lice)
          ELSE IF (TRIM(KeyWord).eq.'DTICE') THEN
            Npts=load_r(Nval, Rval, Ngrids, dtice)
          ELSE IF (TRIM(KeyWord).eq.'DTICE_EQ') THEN
            Npts=load_r(Nval, Rval, Ngrids, dtice_eq)
          ELSE IF (TRIM(KeyWord).eq.'nstrs') THEN
            Npts=load_i(Nval, Rval, Ngrids, nstrs)
          ELSE IF (TRIM(KeyWord).eq.'nevp') THEN
            Npts=load_i(Nval, Rval, Ngrids, nevp)
          ELSE IF (TRIM(KeyWord).eq.'rhoice') THEN
            Npts=load_r(Nval, Rval, Ngrids, rhoice)
          ELSE IF (TRIM(KeyWord).eq.'cdiw') THEN
            Npts=load_r(Nval, Rval, Ngrids, cdiw)
          ELSE IF (TRIM(KeyWord).eq.'cdai') THEN
            Npts=load_r(Nval, Rval, Ngrids, cdai)
          ELSE IF (TRIM(KeyWord).eq.'rho_air') THEN
            Npts=load_r(Nval, Rval, Ngrids, rho_air)
          ELSE IF (TRIM(KeyWord).eq.'rhosnow_dry') THEN
            Npts=load_r(Nval, Rval, Ngrids, rhosnow_dry)
          ELSE IF (TRIM(KeyWord).eq.'rhosnow_wet') THEN
            Npts=load_r(Nval, Rval, Ngrids, rhosnow_wet)
          ELSE IF (TRIM(KeyWord).eq.'pstar') THEN
            Npts=load_r(Nval, Rval, Ngrids, pstar)
          ELSE IF (TRIM(KeyWord).eq.'astren') THEN
            Npts=load_r(Nval, Rval, Ngrids, astren)
          ELSE IF (TRIM(KeyWord).eq.'zetamax') THEN
            Npts=load_r(Nval, Rval, Ngrids, zetamax)
          ELSE IF (TRIM(KeyWord).eq.'zetamin') THEN
            Npts=load_r(Nval, Rval, Ngrids, zetamin)
          ELSE IF (TRIM(KeyWord).eq.'ellip_sq') THEN
            Npts=load_r(Nval, Rval, Ngrids, ellip_sq)
          ELSE IF (TRIM(KeyWord).eq.'alphai') THEN
            Npts=load_r(Nval, Rval, Ngrids, alphai)
            do ng=1,Ngrids
               alphai(ng) = alphai(ng)*deg2rad
            enddo
          ELSE IF (TRIM(KeyWord).eq.'tol') THEN
            Npts=load_r(Nval, Rval, 1, tol)
          ELSE IF (TRIM(KeyWord).eq.'min_h') THEN
            Npts=load_r(Nval, Rval, Ngrids, min_h)
          ELSE IF (TRIM(KeyWord).eq.'min_a') THEN
            Npts=load_r(Nval, Rval, Ngrids, min_a)
          ELSE IF (TRIM(KeyWord).eq.'max_a') THEN
            Npts=load_r(Nval, Rval, Ngrids, max_a)
          ELSE IF (TRIM(KeyWord).eq.'sfwat_max') THEN
            Npts=load_r(Nval, Rval, Ngrids, sfwat_max)
          ELSE IF (TRIM(KeyWord).eq.'stressang') THEN
            Npts=load_r(Nval, Rval, Ngrids, stressang)
            do ng=1,Ngrids
               stressang(ng) = stressang(ng)*deg2rad
            enddo
          ELSE IF (TRIM(KeyWord).eq.'ice_emiss') THEN
            Npts=load_r(Nval, Rval, 1, ice_emiss)
          ELSE IF (TRIM(KeyWord).eq.'spec_heat_air') THEN
            Npts=load_r(Nval, Rval, 1, spec_heat_air)
          ELSE IF (TRIM(KeyWord).eq.'trans_coeff') THEN
            Npts=load_r(Nval, Rval, 1, trans_coeff)
          ELSE IF (TRIM(KeyWord).eq.'sublim_latent_heat') THEN
            Npts=load_r(Nval, Rval, 1, sublim_latent_heat)
          ELSE IF (TRIM(KeyWord).eq.'t0deg') THEN
            Npts=load_r(Nval, Rval, 1, t0deg)
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUice)') THEN
            IF (idUice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUice,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVice)') THEN
            IF (idVice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVice,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUiceE)') THEN
            IF (idUiceE.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUiceE'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUiceE,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idViceN)') THEN
            IF (idViceN.eq.0) THEN
              IF (Master) WRITE (out,280) 'idViceN'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idViceN,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idAice)') THEN
            IF (idAice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idAice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idAice,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHice)') THEN
            IF (idHice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHice,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHsno)') THEN
            IF (idHsno.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHsno'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHsno,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTice)') THEN
            IF (idTice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTice,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTimid)') THEN
            IF (idTimid.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTimid'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTimid,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTauiw)') THEN
            IF (idTauiw.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTauiw'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTauiw,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idChuiw)') THEN
            IF (idChuiw.eq.0) THEN
              IF (Master) WRITE (out,280) 'idChuiw'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idChuiw,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSfwat)') THEN
            IF (idSfwat.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSfwat'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSfwat,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idAgeice)') THEN
            IF (idAgeice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idAgeice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idAgeice,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idIomflx)') THEN
            IF (idIomflx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idIomflx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idIomflx,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSig11)') THEN
            IF (idSig11.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSig11'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSig11,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSig12)') THEN
            IF (idSig12.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSig12'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSig12,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSig22)') THEN
            IF (idSig22.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSig22'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSig22,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idT0mk)') THEN
            IF (idT0mk.eq.0) THEN
              IF (Master) WRITE (out,280) 'idT0mk'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idT0mk,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idS0mk)') THEN
            IF (idS0mk.eq.0) THEN
              IF (Master) WRITE (out,280) 'idS0mk'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idS0mk,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWfr)') THEN
            IF (idWfr.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWfr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWfr,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWai)') THEN
            IF (idWai.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWai'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWai,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWao)') THEN
            IF (idWao.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWao'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWao,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWio)') THEN
            IF (idWio.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWio'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWio,:))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWro)') THEN
            IF (idWro.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWro'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWro,:))
#ifdef AVERAGES
          ELSE IF (TRIM(KeyWord).eq.'Aout(idUice)') THEN
            IF (idUice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idUice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idVice)') THEN
            IF (idVice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idVice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idVice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idUiceE)') THEN
            IF (idUiceE.eq.0) THEN
              IF (Master) WRITE (out,280) 'idUiceE'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idUiceE,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idViceN)') THEN
            IF (idViceN.eq.0) THEN
              IF (Master) WRITE (out,280) 'idViceN'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idViceN,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idAice)') THEN
            IF (idAice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idAice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idAice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idHice)') THEN
            IF (idHice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idHice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idHsno)') THEN
            IF (idHsno.eq.0) THEN
              IF (Master) WRITE (out,280) 'idHsno'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idHsno,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idTice)') THEN
            IF (idTice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idTice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idTimid)') THEN
            IF (idTimid.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTimid'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idTimid,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idTauiw)') THEN
            IF (idTauiw.eq.0) THEN
              IF (Master) WRITE (out,280) 'idTauiw'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idTauiw,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idChuiw)') THEN
            IF (idChuiw.eq.0) THEN
              IF (Master) WRITE (out,280) 'idChuiw'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idChuiw,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idSfwat)') THEN
            IF (idSfwat.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSfwat'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idSfwat,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idAgeice)') THEN
            IF (idAgeice.eq.0) THEN
              IF (Master) WRITE (out,280) 'idAgeice'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idAgeice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idIomflx)') THEN
            IF (idIomflx.eq.0) THEN
              IF (Master) WRITE (out,280) 'idIomflx'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idIomflx,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idSig11)') THEN
            IF (idSig11.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSig11'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idSig11,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idSig12)') THEN
            IF (idSig12.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSig12'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idSig12,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idSig22)') THEN
            IF (idSig22.eq.0) THEN
              IF (Master) WRITE (out,280) 'idSig22'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idSig22,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idT0mk)') THEN
            IF (idT0mk.eq.0) THEN
              IF (Master) WRITE (out,280) 'idT0mk'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idT0mk,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idS0mk)') THEN
            IF (idS0mk.eq.0) THEN
              IF (Master) WRITE (out,280) 'idS0mk'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idS0mk,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idWfr)') THEN
            IF (idWfr.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWfr'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idWfr,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idWai)') THEN
            IF (idWai.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWai'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idWai,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idWao)') THEN
            IF (idWao.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWao'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idWao,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idWio)') THEN
            IF (idWio.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWio'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idWio,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout(idWro)') THEN
            IF (idWro.eq.0) THEN
              IF (Master) WRITE (out,280) 'idWro'
              exit_flag=5
              RETURN
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Aout(idWro,:))
#endif
#ifdef AVERAGES2
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idUice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idUice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idVice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idVice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idUiceE)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idUiceE,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idViceN)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idViceN,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idAice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idAice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idHice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idHice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idHsno)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idHsno,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idTice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idTice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idTimid)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idTimid,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idTauiw)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idTauiw,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idChuiw)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idChuiw,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idSfwat)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idSfwat,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idAgeice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idAgeice,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idIomflx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idIomflx,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idSig11)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idSig11,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idSig12)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idSig12,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idSig22)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idSig22,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idT0mk)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idT0mk,:))
          ELSE IF (TRIM(KeyWord).eq.'Aout2(idS0mk)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Aout2(idS0mk,:))
#endif
#ifdef STATIONS
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUice)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUice,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVice)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVice,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUiceE)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUiceE,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idViceN)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idViceN,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idAice)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idAice,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHice)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHice,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTice)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTice,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTimid)') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTimid,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSfwat)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSfwat,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idAgeice)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idAgeice,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idIomflx)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idIomflx,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSig11)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSig11,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSig12)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSig12,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSig22)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSig22,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWfr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWfr,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTauiw)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTauiw,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idChuiw)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idChuiw,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idT0mk)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idT0mk,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idS0mk)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idS0mk,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWfr)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWfr,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWai)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWai,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWao)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWao,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWio)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWio,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWro)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWro,1))
#endif
          END IF
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CLOSE (inp)

! Set ice time step to ocean time step

      DO ng = 1,Ngrids
        dtice(ng) = dt(ng)
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (.not. Lice(ng)) RETURN
          WRITE (out,40) ng
          WRITE(out,*) 'Ice time step = ocean time step'
          WRITE (out,100) dtice(ng), 'DTICE',                           &
     &          'Ice model time step (s).'
          WRITE (out,100) dtice_eq(ng), 'DTICE_EQ',                     &
     &          'Ice drift update (equilibrium) time step (s).'
          WRITE (out,50) nstrs(ng), 'nstrs',                            &
     &          'Number of iterations for nonlinear ice dynamics.'
          WRITE (out,50) nevp(ng), 'nevp',                              &
     &          'Number of elastic steps per plastic step in EVP.'
          WRITE (out,100) rhoice(ng), 'rhoice',                         &
     &          'Density of sea ice (kg/m3).'
          WRITE (out,100) cdiw(ng), 'cdiw',                             &
     &          'Ice-water drag coefficient (nondimensional).'
          WRITE (out,100) cdai(ng), 'cdai',                             &
     &          'Air-ice drag coefficient (nondimensional).'
          WRITE (out,100) rho_air(ng), 'rho_air',                       &
     &          'Air density (kg/m3).'
          WRITE (out,100) rhosnow_dry(ng), 'rhosnow_dry',               &
     &          'Dry snow density (kg/m3).'
          WRITE (out,100) rhosnow_wet(ng), 'rhosnow_wet',               &
     &          'Wet snow density (kg/m3).'
          WRITE (out,100) alphai(ng)*rad2deg, 'alphai',                 &
     &          'Mohr-Coulomb stress angle (degrees).'
          WRITE (out,100) min_h(ng), 'min_h',                           &
     &          'Minimum average ice thickness (m).'
          WRITE (out,100) min_a(ng), 'min_a',                           &
     &          'Minimum ice concentration (nondimensional).'
          WRITE (out,100) max_a(ng), 'max_a',                           &
     &          'Maximum ice concentration (nondimensional).'
          WRITE (out,100) sfwat_max(ng)*rad2deg, 'sfwat_max',           &
     &          'Maximum surface fresh water (nondimensional).'
          WRITE (out,100) stressang(ng)*rad2deg, 'stressang',           &
     &          'Turning angle for ice-water drag (degrees).'
          WRITE (out,100) tol, 'tol',                                   &
     &          'Numerical tolerance in rheology calculations .'
          WRITE (out,100) ice_emiss, 'ice_emiss',                       &
     &          'Ice emissivity.'
          WRITE (out,100) spec_heat_air, 'spec_heat_air',               &
     &          'Specific heat of air.'
          WRITE (out,100) trans_coeff, 'trans_coeff',                   &
     &          'Transfer coefficient.'
          WRITE (out,100) sublim_latent_heat, 'sublim_latent_heat',     &
     &          'Latent_heat of sublimation.'
          WRITE (out,100) t0deg, 't0deg',                               &
     &          'Zero degrees Celsius in degrees Kelvin.'
          IF (Hout(idUice,ng)) WRITE (out,170) Hout(idUice,ng),         &
     &       'Hout(idUice)',                                            &
     &       'Write out U-component ice velocity.'
          IF (Hout(idVice,ng)) WRITE (out,170) Hout(idVice,ng),         &
     &       'Hout(idVice)',                                            &
     &       'Write out V-component ice velocity.'
          IF (Hout(idUiceE,ng)) WRITE (out,170) Hout(idUiceE,ng),       &
     &       'Hout(idUiceE)',                                           &
     &       'Write out East component ice velocity.'
          IF (Hout(idViceN,ng)) WRITE (out,170) Hout(idViceN,ng),       &
     &       'Hout(idViceN)',                                           &
     &       'Write out North component ice velocity.'
          IF (Hout(idAice,ng)) WRITE (out,170) Hout(idAice,ng),         &
     &       'Hout(idAice)',                                            &
     &       'Write out ice concentration.'
          IF (Hout(idHice,ng)) WRITE (out,170) Hout(idHice,ng),         &
     &       'Hout(idHice)',                                            &
     &       'Write out average ice thickness.'
          IF (Hout(idHsno,ng)) WRITE (out,170) Hout(idHsno,ng),         &
     &       'Hout(idHsno)',                                            &
     &       'Write out snow thickness.'
          IF (Hout(idTice,ng)) WRITE (out,170) Hout(idTice,ng),         &
     &       'Hout(idTice)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Hout(idTimid,ng)) WRITE (out,170) Hout(idTimid,ng),       &
     &       'Hout(idTimid)',                                           &
     &       'Write out interior ice temperature.'
          IF (Hout(idSfwat,ng)) WRITE (out,170) Hout(idSfwat,ng),       &
     &       'Hout(idSfwat)',                                           &
     &       'Write out surface water (on ice) thickness.'
          IF (Hout(idAgeice,ng)) WRITE (out,170) Hout(idAgeice,ng),     &
     &       'Hout(idAgeice)',                                          &
     &       'Write out ice age.'
          IF (Hout(idIomflx,ng)) WRITE (out,170) Hout(idIomflx,ng),     &
     &       'Hout(idIomflx)',                                          &
     &       'Write out ice-ocean mass flux'
          IF (Hout(idSig11,ng)) WRITE (out,170) Hout(idSig11,ng),       &
     &       'Hout(idSig11)',                                           &
     &       'Write out internal ice stress component 11.'
          IF (Hout(idSig12,ng)) WRITE (out,170) Hout(idSig12,ng),       &
     &       'Hout(idSig12)',                                           &
     &       'Write out internal ice stress component 12.'
          IF (Hout(idSig22,ng)) WRITE (out,170) Hout(idSig22,ng),       &
     &       'Hout(idSig22)',                                           &
     &       'Write out internal ice stress component 22.'
          IF (Hout(idTauiw,ng)) WRITE (out,170) Hout(idTauiw,ng),       &
     &       'Hout(idTauiw)',                                           &
     &       'Write out ice-water friction velocity.'
          IF (Hout(idChuiw,ng)) WRITE (out,170) Hout(idChuiw,ng),       &
     &       'Hout(idChuiw)',                                           &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Hout(idT0mk,ng)) WRITE (out,170) Hout(idT0mk,ng),         &
     &       'Hout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Hout(idS0mk,ng)) WRITE (out,170) Hout(idS0mk,ng),         &
     &       'Hout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
          IF (Hout(idWfr,ng)) WRITE (out,170) Hout(idWfr,ng),           &
     &       'Hout(idWfr)',                                             &
     &       'Write out frazil ice growth rate.'
          IF (Hout(idWai,ng)) WRITE (out,170) Hout(idWai,ng),           &
     &       'Hout(idWai)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Hout(idWao,ng)) WRITE (out,170) Hout(idWao,ng),           &
     &       'Hout(idWao)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Hout(idWio,ng)) WRITE (out,170) Hout(idWio,ng),           &
     &       'Hout(idWio)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Hout(idWro,ng)) WRITE (out,170) Hout(idWro,ng),           &
     &       'Hout(idWro)',                                             &
     &       'Write out ice melt runoff rate.'
#ifdef AVERAGES
          IF (Aout(idUice,ng)) WRITE (out,170) Aout(idUice,ng),         &
     &       'Aout(idUice)',                                            &
     &       'Write out time-averaged U-component ice velocity.'
          IF (Aout(idVice,ng)) WRITE (out,170) Aout(idVice,ng),         &
     &       'Aout(idVice)',                                            &
     &       'Write out time-averaged V-component ice velocity.'
          IF (Aout(idUiceE,ng)) WRITE (out,170) Aout(idUiceE,ng),       &
     &       'Aout(idUiceE)',                                           &
     &       'Write out time-averaged East component ice velocity.'
          IF (Aout(idViceN,ng)) WRITE (out,170) Aout(idViceN,ng),       &
     &       'Aout(idViceN)',                                           &
     &       'Write out time-averaged North component ice velocity.'
          IF (Aout(idAice,ng)) WRITE (out,170) Aout(idAice,ng),         &
     &       'Aout(idAice)',                                            &
     &       'Write out time-averaged ice concentration.'
          IF (Aout(idHice,ng)) WRITE (out,170) Aout(idHice,ng),         &
     &       'Aout(idHice)',                                            &
     &       'Write out time-averaged average ice thickness.'
          IF (Aout(idHsno,ng)) WRITE (out,170) Aout(idHsno,ng),         &
     &       'Aout(idHsno)',                                            &
     &       'Write out time-averaged snow thickness.'
          IF (Aout(idTice,ng)) WRITE (out,170) Aout(idTice,ng),         &
     &       'Aout(idTice)',                                            &
     &       'Write out time-averaged ice/snow surface temperature.'
          IF (Aout(idTimid,ng)) WRITE (out,170) Aout(idTimid,ng),       &
     &       'Aout(idTimid)',                                           &
     &       'Write out time-averaged interior ice temperature.'
          IF (Aout(idSfwat,ng)) WRITE (out,170) Aout(idSfwat,ng),       &
     &       'Aout(idSfwat)',                                           &
     &       'Write out time-averaged surface water (on ice) thickness.'
          IF (Aout(idAgeice,ng)) WRITE (out,170) Aout(idAgeice,ng),     &
     &       'Aout(idAgeice)',                                          &
     &       'Write out time-averaged ice age.'
          IF (Aout(idIomflx,ng)) WRITE (out,170) Aout(idIomflx,ng),     &
     &       'Aout(idIomflx)',                                          &
     &       'Write out time-averaged ice-ocean mass flux'
          IF (Aout(idSig11,ng)) WRITE (out,170) Aout(idSig11,ng),       &
     &       'Aout(idSig11)',                                           &
     &       'Write out time-averaged internal ice stress component 11.'
          IF (Aout(idSig12,ng)) WRITE (out,170) Aout(idSig12,ng),       &
     &       'Aout(idSig12)',                                           &
     &       'Write out time-averaged internal ice stress component 12.'
          IF (Aout(idSig22,ng)) WRITE (out,170) Aout(idSig22,ng),       &
     &       'Aout(idSig22)',                                           &
     &       'Write out time-averaged internal ice stress component 22.'
          IF (Aout(idTauiw,ng)) WRITE (out,170) Aout(idTauiw,ng),       &
     &       'Aout(idTauiw)',                                           &
     &       'Write out time-averaged ice-water friction velocity.'
          IF (Aout(idChuiw,ng)) WRITE (out,170) Aout(idChuiw,ng),       &
     &       'Aout(idChuiw)',                                           &
     &       'Write out time-averaged ice-water transfer coefficient.'
          IF (Aout(idT0mk,ng)) WRITE (out,170) Aout(idT0mk,ng),         &
     &       'Aout(idT0mk)',                                            &
     &       'Write out time-averaged under ice temperature.'
          IF (Aout(idS0mk,ng)) WRITE (out,170) Aout(idS0mk,ng),         &
     &       'Aout(idS0mk)',                                            &
     &       'Write out time-averaged under ice salinity.'
          IF (Aout(idWfr,ng)) WRITE (out,170) Aout(idWfr,ng),           &
     &       'Aout(idWfr)',                                             &
     &       'Write out time-averaged frazil ice growth rate.'
          IF (Aout(idWai,ng)) WRITE (out,170) Aout(idWai,ng),           &
     &       'Aout(idWai)',                                             &
     &       'Write out time-averaged ice growth/melt rate.'
          IF (Aout(idWao,ng)) WRITE (out,170) Aout(idWao,ng),           &
     &       'Aout(idWao)',                                             &
     &       'Write out time-averaged ice growth/melt rate.'
          IF (Aout(idWio,ng)) WRITE (out,170) Aout(idWio,ng),           &
     &       'Aout(idWio)',                                             &
     &       'Write out time-averaged ice growth/melt rate.'
          IF (Aout(idWro,ng)) WRITE (out,170) Aout(idWro,ng),           &
     &       'Aout(idWro)',                                             &
     &       'Write out time-averaged ice melt runoff rate.'
#endif
#ifdef AVERAGES2
          IF (Aout2(idUice,ng)) WRITE (out,170) Aout2(idUice,ng),       &
     &       'Aout2(idUice)',                                           &
     &       'Write out U-component ice velocity.'
          IF (Aout2(idVice,ng)) WRITE (out,170) Aout2(idVice,ng),       &
     &       'Aout2(idVice)',                                           &
     &       'Write out V-component ice velocity.'
          IF (Aout2(idUiceE,ng)) WRITE (out,170) Aout2(idUiceE,ng),     &
     &       'Aout2(idUiceE)',                                          &
     &       'Write out East component ice velocity.'
          IF (Aout2(idViceN,ng)) WRITE (out,170) Aout2(idViceN,ng),     &
     &       'Aout2(idViceN)',                                          &
     &       'Write out North component ice velocity.'
          IF (Aout2(idAice,ng)) WRITE (out,170) Aout2(idAice,ng),       &
     &       'Aout2(idAice)',                                           &
     &       'Write out ice concentration.'
          IF (Aout2(idHice,ng)) WRITE (out,170) Aout2(idHice,ng),       &
     &       'Aout2(idHice)',                                           &
     &       'Write out average ice thickness.'
          IF (Aout2(idHsno,ng)) WRITE (out,170) Aout2(idHsno,ng),       &
     &       'Aout2(idHsno)',                                           &
     &       'Write out snow thickness.'
          IF (Aout2(idTice,ng)) WRITE (out,170) Aout2(idTice,ng),       &
     &       'Aout2(idTice)',                                           &
     &       'Write out ice/snow surface temperature.'
          IF (Aout2(idTimid,ng)) WRITE (out,170) Aout2(idTimid,ng),     &
     &       'Aout2(idTimid)',                                          &
     &       'Write out interior ice temperature.'
          IF (Aout2(idSfwat,ng)) WRITE (out,170) Aout2(idSfwat,ng),     &
     &       'Aout2(idSfwat)',                                          &
     &       'Write out surface water (on ice) thickness.'
          IF (Aout2(idAgeice,ng)) WRITE (out,170) Aout2(idAgeice,ng),   &
     &       'Aout2(idAgeice)',                                         &
     &       'Write out ice age.'
          IF (Aout2(idIomflx,ng)) WRITE (out,170) Aout2(idIomflx,ng),   &
     &       'Aout2(idIomflx)',                                         &
     &       'Write out ice-ocean mass flux'
          IF (Aout2(idSig11,ng)) WRITE (out,170) Aout2(idSig11,ng),     &
     &       'Aout2(idSig11)',                                          &
     &       'Write out internal ice stress component 11.'
          IF (Aout2(idSig12,ng)) WRITE (out,170) Aout2(idSig12,ng),     &
     &       'Aout2(idSig12)',                                          &
     &       'Write out internal ice stress component 12.'
          IF (Aout2(idSig22,ng)) WRITE (out,170) Aout2(idSig22,ng),     &
     &       'Aout2(idSig22)',                                          &
     &       'Write out internal ice stress component 22.'
          IF (Aout2(idTauiw,ng)) WRITE (out,170) Aout2(idTauiw,ng),     &
     &       'Aout2(idTauiw)',                                          &
     &       'Write out ice-water friction velocity.'
          IF (Aout2(idChuiw,ng)) WRITE (out,170) Aout2(idChuiw,ng),     &
     &       'Aout2(idChuiw)',                                          &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Aout2(idT0mk,ng)) WRITE (out,170) Aout2(idT0mk,ng),       &
     &       'Aout2(idT0mk)',                                           &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Aout2(idS0mk,ng)) WRITE (out,170) Aout2(idS0mk,ng),       &
     &       'Aout2(idS0mk)',                                           &
     &       'Write out salinity of molecular sublayer under ice.'
#endif
#ifdef STATIONS
          IF (Sout(idUice,ng)) WRITE (out,170) Sout(idUice,ng),         &
     &       'Sout(idUice)',                                            &
     &       'Write out U-component ice velocity.'
          IF (Sout(idVice,ng)) WRITE (out,170) Sout(idVice,ng),         &
     &       'Sout(idVice)',                                            &
     &       'Write out V-component ice velocity.'
          IF (Sout(idUiceE,ng)) WRITE (out,170) Sout(idUiceE,ng),       &
     &       'Sout(idUiceE)',                                           &
     &       'Write out East-component ice velocity.'
          IF (Sout(idViceN,ng)) WRITE (out,170) Sout(idViceN,ng),       &
     &       'Sout(idViceN)',                                           &
     &       'Write out North-component ice velocity.'
          IF (Sout(idAice,ng)) WRITE (out,170) Sout(idAice,ng),         &
     &       'Sout(idAice)',                                            &
     &       'Write out ice concentration.'
          IF (Sout(idHice,ng)) WRITE (out,170) Sout(idHice,ng),         &
     &       'Sout(idHice)',                                            &
     &       'Write out average ice thickness.'
          IF (Sout(idHsno,ng)) WRITE (out,170) Sout(idHsno,ng),         &
     &       'Sout(idHsno)',                                            &
     &       'Write out snow thickness.'
          IF (Sout(idTice,ng)) WRITE (out,170) Sout(idTice,ng),         &
     &       'Sout(idTice)',                                            &
     &       'Write out ice/snow surface temperature.'
          IF (Sout(idTimid,ng)) WRITE (out,170) Sout(idTimid,ng),       &
     &       'Sout(idTimid)',                                           &
     &       'Write out interior ice temperature.'
          IF (Sout(idSfwat,ng)) WRITE (out,170) Sout(idSfwat,ng),       &
     &       'Sout(idSfwat)',                                           &
     &       'Write out surface water (on ice) thickness.'
          IF (Sout(idAgeice,ng)) WRITE (out,170) Sout(idAgeice,ng),     &
     &       'Sout(idAgeice)',                                          &
     &       'Write out surface water (on ice) thickness.'
          IF (Sout(idIomflx,ng)) WRITE (out,170) Sout(idIomflx,ng),     &
     &       'Sout(idIomflx)',                                          &
     &       'Write out ice-ocean mass flux.'
          IF (Sout(idSig11,ng)) WRITE (out,170) Sout(idSig11,ng),       &
     &       'Sout(idSig11)',                                           &
     &       'Write out internal ice stress component 11.'
          IF (Sout(idSig12,ng)) WRITE (out,170) Sout(idSig12,ng),       &
     &       'Sout(idSig12)',                                           &
     &       'Write out internal ice stress component 12.'
          IF (Sout(idSig22,ng)) WRITE (out,170) Sout(idSig22,ng),       &
     &       'Sout(idSig22)',                                           &
     &       'Write out internal ice stress component 22.'
          IF (Sout(idTauiw,ng)) WRITE (out,170) Sout(idTauiw,ng),       &
     &       'Sout(idTauiw)',                                           &
     &       'Write out ice-water friction velocity.'
          IF (Sout(idChuiw,ng)) WRITE (out,170) Sout(idChuiw,ng),       &
     &       'Hout(idChuiw)',                                           &
     &       'Write out ice-water momentum transfer coefficient.'
          IF (Sout(idT0mk,ng)) WRITE (out,170) Sout(idT0mk,ng),         &
     &       'Sout(idT0mk)',                                            &
     &       'Write out temperature of molecular sublayer under ice.'
          IF (Sout(idS0mk,ng)) WRITE (out,170) Sout(idS0mk,ng),         &
     &       'Sout(idS0mk)',                                            &
     &       'Write out salinity of molecular sublayer under ice.'
          IF (Sout(idWfr,ng)) WRITE (out,170) Sout(idWfr,ng),           &
     &       'Sout(idWfr)',                                             &
     &       'Write out frazil ice growth rate.'
          IF (Sout(idWai,ng)) WRITE (out,170) Sout(idWai,ng),           &
     &       'Sout(idWai)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Sout(idWao,ng)) WRITE (out,170) Sout(idWao,ng),           &
     &       'Sout(idWao)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Sout(idWio,ng)) WRITE (out,170) Sout(idWio,ng),           &
     &       'Sout(idWio)',                                             &
     &       'Write out ice growth/melt rate.'
          IF (Sout(idWro,ng)) WRITE (out,170) Sout(idWro,ng),           &
     &       'Sout(idWro)',                                             &
     &       'Write out ice melt runoff rate.'
#endif
        END DO
      END IF

  30  FORMAT (/,' READ_IcePar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Ice Parameters, Grid: ',i2.2,                   &
     &        /,' ============================',/)
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t28,a)
  80  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t28,a)
 110  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)
 120  FORMAT (10x,l1,2x,a,t28,a)
 170  FORMAT (10x,l1,2x,a,t30,a)
 280  FORMAT (/,' READ_IcePar - variable info not yet loaded, ', a)
      RETURN
      END SUBROUTINE read_IcePar
#endif
