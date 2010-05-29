      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!=======================================================================
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, itrc, ng, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in UMaine CoSiNE biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'reg1') THEN
            Npts=load_r(Nval, Rval, Ngrids, reg1)
          ELSE IF (TRIM(KeyWord).eq.'reg2') THEN
            Npts=load_r(Nval, Rval, Ngrids, reg2)
          ELSE IF (TRIM(KeyWord).eq.'gmaxs1') THEN
            Npts=load_r(Nval, Rval, Ngrids, gmaxs1)
          ELSE IF (TRIM(KeyWord).eq.'gmaxs2') THEN
            Npts=load_r(Nval, Rval, Ngrids, gmaxs2)
          ELSE IF (TRIM(KeyWord).eq.'beta1') THEN
            Npts=load_r(Nval, Rval, Ngrids, beta1)
          ELSE IF (TRIM(KeyWord).eq.'beta2') THEN
            Npts=load_r(Nval, Rval, Ngrids, beta2)
          ELSE IF (TRIM(KeyWord).eq.'akz1') THEN
            Npts=load_r(Nval, Rval, Ngrids, akz1)
          ELSE IF (TRIM(KeyWord).eq.'akz2') THEN
            Npts=load_r(Nval, Rval, Ngrids, akz2)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'amaxs2') THEN
            Npts=load_r(Nval, Rval, Ngrids, amaxs2)
          ELSE IF (TRIM(KeyWord).eq.'parsats1') THEN
            Npts=load_r(Nval, Rval, Ngrids, parsats1)
          ELSE IF (TRIM(KeyWord).eq.'parsats2') THEN
            Npts=load_r(Nval, Rval, Ngrids, parsats2)
          ELSE IF (TRIM(KeyWord).eq.'pis1') THEN
            Npts=load_r(Nval, Rval, Ngrids, pis1)
          ELSE IF (TRIM(KeyWord).eq.'pis2') THEN
            Npts=load_r(Nval, Rval, Ngrids, pis2)
          ELSE IF (TRIM(KeyWord).eq.'akno3s1') THEN
            Npts=load_r(Nval, Rval, Ngrids, akno3s1)
          ELSE IF (TRIM(KeyWord).eq.'akno3s2') THEN
            Npts=load_r(Nval, Rval, Ngrids, akno3s2)
          ELSE IF (TRIM(KeyWord).eq.'aknh4s1') THEN
            Npts=load_r(Nval, Rval, Ngrids, aknh4s1)
          ELSE IF (TRIM(KeyWord).eq.'aknh4s2') THEN
            Npts=load_r(Nval, Rval, Ngrids, aknh4s2)
          ELSE IF (TRIM(KeyWord).eq.'akpo4s1') THEN
            Npts=load_r(Nval, Rval, Ngrids, akpo4s1)
          ELSE IF (TRIM(KeyWord).eq.'akpo4s2') THEN
            Npts=load_r(Nval, Rval, Ngrids, akpo4s2)
          ELSE IF (TRIM(KeyWord).eq.'akco2s1') THEN
            Npts=load_r(Nval, Rval, Ngrids, akco2s1)
          ELSE IF (TRIM(KeyWord).eq.'akco2s2') THEN
            Npts=load_r(Nval, Rval, Ngrids, akco2s2)
          ELSE IF (TRIM(KeyWord).eq.'aksio4s2') THEN
            Npts=load_r(Nval, Rval, Ngrids, aksio4s2)
          ELSE IF (TRIM(KeyWord).eq.'ak1') THEN
            Npts=load_r(Nval, Rval, Ngrids, ak1)
          ELSE IF (TRIM(KeyWord).eq.'ak2') THEN
            Npts=load_r(Nval, Rval, Ngrids, ak2)
          ELSE IF (TRIM(KeyWord).eq.'bgamma0') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma0)
          ELSE IF (TRIM(KeyWord).eq.'bgamma1') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma1)
          ELSE IF (TRIM(KeyWord).eq.'bgamma2') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma2)
          ELSE IF (TRIM(KeyWord).eq.'bgamma3') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma3)
          ELSE IF (TRIM(KeyWord).eq.'bgamma4') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma4)
          ELSE IF (TRIM(KeyWord).eq.'bgamma5') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma5)
          ELSE IF (TRIM(KeyWord).eq.'bgamma6') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma6)
          ELSE IF (TRIM(KeyWord).eq.'bgamma7') THEN
            Npts=load_r(Nval, Rval, Ngrids, bgamma7)
          ELSE IF (TRIM(KeyWord).eq.'wsd') THEN
            Npts=load_r(Nval, Rval, Ngrids, wsd)
          ELSE IF (TRIM(KeyWord).eq.'wsdsi') THEN
            Npts=load_r(Nval, Rval, Ngrids, wsdsi)
          ELSE IF (TRIM(KeyWord).eq.'wsp') THEN
            Npts=load_r(Nval, Rval, Ngrids, wsp)
          ELSE IF (TRIM(KeyWord).eq.'pco2a') THEN
            Npts=load_r(Nval, Rval, Ngrids, pco2a)
          ELSE IF (TRIM(KeyWord).eq.'si2n') THEN
            Npts=load_r(Nval, Rval, Ngrids, si2n)
          ELSE IF (TRIM(KeyWord).eq.'p2n') THEN
            Npts=load_r(Nval, Rval, Ngrids, p2n)
          ELSE IF (TRIM(KeyWord).eq.'o2no') THEN
            Npts=load_r(Nval, Rval, Ngrids, o2no)
          ELSE IF (TRIM(KeyWord).eq.'o2nh') THEN
            Npts=load_r(Nval, Rval, Ngrids, o2nh)
          ELSE IF (TRIM(KeyWord).eq.'c2n') THEN
            Npts=load_r(Nval, Rval, Ngrids, c2n)
          ELSE IF (TRIM(KeyWord).eq.'ro5') THEN
            Npts=load_r(Nval, Rval, Ngrids, ro5)
          ELSE IF (TRIM(KeyWord).eq.'ro6') THEN
            Npts=load_r(Nval, Rval, Ngrids, ro6)
          ELSE IF (TRIM(KeyWord).eq.'ro7') THEN
            Npts=load_r(Nval, Rval, Ngrids, ro7)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,120)                           &
     &                      'idTvar(idbio(', itrc, '))'
                  exit_flag=5
                  RETURN
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTsur(idbio(itrc))
                IF (i.eq.0) THEN
                  IF (Master) WRITE (out,120)                           &
     &                      'idTsur(idbio(', itrc, '))'
                  exit_flag=5
                  RETURN
                END IF                  
                Hout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#ifdef AVERAGES
          ELSE IF (TRIM(KeyWord).eq.'Aout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Ltrc)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                Aout(i,ng)=Ltrc(itrc,ng)
              END DO
            END DO
#endif
          END IF
        END IF
      END DO
  10  IF (Master) WRITE (out,30) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,110) reg1(ng), 'reg1',                           &
     &            'Microzooplankton excretion rate to ammonium',        &
     &            '[1/day].'
            WRITE (out,100) reg2(ng), 'reg2',                           &
     &            'Mesozooplankton excretion rate to ammonium [1/day].'
            WRITE (out,110) gmaxs1(ng), 'gmaxs1',                       &
     &            'Maximum specific growth rate of small phytoplankton',&
     &            '[1/day].'
            WRITE (out,100) gmaxs2(ng), 'gmaxs2',                       &
     &            'Maximum specific growth rate of diatom [1/day].'
            WRITE (out,100) beta1(ng), 'beta1',                         &
     &            'Microzooplankton maximum grazing rate [1/day].'
            WRITE (out,100) beta2(ng), 'beta2',                         &
     &            'Mesozooplankton maximum grazing rate [1/day].'
            WRITE (out,110) akz1(ng), 'akz1',                           &
     &            'Half saturation constant for microzooplankton',      &
     &            'grazing [mmol_N/m3].'
            WRITE (out,110) akz2(ng), 'akz2',                           &
     &            'Half saturation constant for mesozooplankton',       &
     &            'grazing [mmol_N/m3].'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) amaxs2(ng), 'amaxs2',                       &
     &            'Initial slope of P-I curve of small',                &
     &            'phytoplankton [1/(Watts/m2)/day].'
            WRITE (out,110) parsats1(ng), 'parsats1',                   &
     &            'PAR saturation onset parameter of',                  &
     &            'small phytoplankton [Watts/m2].'
            WRITE (out,110) parsats2(ng), 'parsats2',                   &
     &            'PAR saturation onset parameter of diatom',           &
     &            '[Watts/m2].'
            WRITE (out,110) pis1(ng), 'pis1',                           &
     &            'Ammonium inhibition parameter for small',            &
     &            'phytoplankton [mmol_N/m3].'
            WRITE (out,110) pis2(ng), 'pis2',                           &
     &            'Ammonium inhibition parameter for diatom',           &
     &            '[mmol_N/m3].'
            WRITE (out,110) akno3s1(ng), 'akno3s1',                     &
     &            'Half saturation concentration for nitrate',          &
     &            'uptake by small phytoplankton [mmol_N/m3].'
            WRITE (out,110) akno3s2(ng), 'akno3s2',                     &
     &            'Half saturation concentration for nitrate',          &
     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) aknh4s1(ng), 'aknh4s1',                     &
     &            'Half saturation concentration for ammonium',         &
     &            'uptake by small phytoplankton [mmol_N/m3].'
            WRITE (out,110) aknh4s2(ng), 'aknh4s2',                     &
     &            'Half saturation concentration for ammonium',         &
     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,110) akpo4s1(ng), 'akpo4s1',                     &
     &            'Half saturation concentration for phosphate',        &
     &            'uptake by small phytoplankton [mmol_P/m3].'
            WRITE (out,110) akpo4s2(ng), 'akpo4s2',                     &
     &            'Half saturation concentration for phosphate',        &
     &            'uptake by diatom [mmol_P/m3].'
            WRITE (out,110) akco2s1(ng), 'akco2s1',                     &
     &            'Half saturation concentration for co2',              &
     &            'uptake by small phytoplankton [mmol_C/m3].'
            WRITE (out,110) akco2s2(ng), 'akco2s2',                     &
     &            'Half saturation concentration for co2',              &
     &            'uptake by diatom [mmol_C/m3].'
            WRITE (out,110) aksio4s2(ng), 'aksio4s2',                   &
     &            'Half saturation constant for silicate',              &
     &            'uptake by diatom [mmol_N/m3].'
            WRITE (out,100) ak1(ng), 'ak1',                             &
     &            'Light attenuation coefficient of water [1/m].'
            WRITE (out,110) ak2(ng), 'ak2',                             &
     &            'Specific light attenuation coefficient for',         &
     &            'phytoplankton [1/m/(mmol_N/m3)].'
            WRITE (out,100) bgamma0(ng), 'bgamma0',                     &
     &            'Mesozooplankton specific mortality rate [1/day].'
            WRITE (out,110) bgamma1(ng), 'bgamma1',                     &
     &            'Grazing efficiency of microzooplankton',             &
     &            '[nondimensional].'
            WRITE (out,110) bgamma2(ng), 'bgamma2',                     &
     &            ' Grazing efficiency of mesozooplankton',             &
     &            '[nondimensional].'
            WRITE (out,100) bgamma3(ng), 'bgamma3',                     &
     &            'Death rate of small phytoplankton [1/day].'
            WRITE (out,100) bgamma4(ng), 'bgamma4',                     &
     &            'Death rate of large phytoplankton [1/day].'
            WRITE (out,100) bgamma5(ng), 'bgamma5',                     &
     &            'Decay rate of detritus [1/day].'
            WRITE (out,100) bgamma6(ng), 'bgamma6',                     &
     &            ' '
            WRITE (out,100) bgamma7(ng), 'bgamma7',                     &
     &            'Nitrafication rate [1/day].'
            WRITE (out,100) wsd(ng), 'wsd',                             &
     &            'Sinking velocity of detritus [m/day].'
            WRITE (out,100) wsdsi(ng), 'wsdsi',                         &
     &            'Sinking velocity of detritus silicate [m/day].'
            WRITE (out,100) wsp(ng), 'wsp',                             &
     &            'Sinking velocity of large phytoplankton [m/day].'
            WRITE (out,100) pco2a(ng), 'pco2a',                         &
     &            'Air pCO2 [ppmv].'
            WRITE (out,100) si2n(ng), 'si2n',                           &
     &            'Silicate to nitrogen ratio [mol_Si/mol_N].'
            WRITE (out,100) p2n(ng), 'p2n',                             &
     &            'Phosphorus to nitrogen ratio [mol_P/mol_N].'
            WRITE (out,100) o2no(ng), 'o2no',                           &
     &            'Oxygen to nitrate ratio [mol_O2/mol_NO3].'
            WRITE (out,100) o2nh(ng), 'o2nh',                           &
     &            'Oxygen to ammonium ratio [mol_O2/mol_NH4].'
            WRITE (out,100) c2n(ng), 'c2n',                             &
     &            'Carbon to nitrogen ratio [mol_C/mol_N].'
            WRITE (out,100) ro5(ng), 'ro5',                             &
     &            'Grazing preference for diatom [nondimensional].'
            WRITE (out,110) ro6(ng), 'ro6',                             &
     &            'Grazing preference for mesozooplankton',             &
     &            '[nondimensional].'
            WRITE (out,100) ro7(ng), 'ro7',                             &
     &            'Grazing preference for detritus [nondimensional].'
#ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#endif
#ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTsur(i),ng)) WRITE (out,60)                    &
     &            Hout(idTsur(i),ng), 'Hout(idTsur)',                   &
     &            'Write out tracer flux ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#ifdef AVERAGES
            WRITE (out,'(1x)')
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Aout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Aout(idTvar(i),ng), 'Aout(idTvar)',                   &
     &            'Write out averaged tracer ', i,                      &
     &            TRIM(Vname(1,idTvar(i)))
            END DO
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Rescale biological tracer parameters
!-----------------------------------------------------------------------
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
      DO ng=1,Ngrids
        DO itrc=1,NBT
          i=idbio(itrc)
          tnu4(i,ng)=SQRT(ABS(tnu4(i,ng)))
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
          IF (Tnudg(i,ng).gt.0.0_r8) THEN
            Tnudg(i,ng)=1.0_r8/(Tnudg(i,ng)*86400.0_r8)
          ELSE
            Tnudg(i,ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' read_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' UMaine CoSiNE Model Parameters, Grid: ',i2.2,       &
     &        /,  ' =================================',/)  
  50  FORMAT (1x,i10,2x,a,t30,a)
  60  FORMAT (10x,l1,2x,a,t30,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t30,a)
  80  FORMAT (f11.3,2x,a,t30,a,/,t32,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t30,a,/,t32,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t30,a)
 110  FORMAT (1p,e11.4,2x,a,t30,a,/,t32,a)
 120  FORMAT (/,' read_BioPar - variable info not yet loaded, ',        &
     &        a,i2.2,a)
 130  FORMAT (/,' read_BioPar - variable info not yet loaded, ',a)
 140  FORMAT (10x,l1,2x,a,t30,a,1x,a)

      RETURN
      END SUBROUTINE read_BioPar
