#include "cppdefs.h"
#ifdef STOPERTURB
      SUBROUTINE set_stoperturb(ng, tile)
!
!=======================================================================
!  Copyright (c) 2002-2014 ROMS/TOMS Group                             !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beginning of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
      USE mod_stepping
      USE mod_stoperturb

      END SUBROUTINE set_stoperturb
#endif
