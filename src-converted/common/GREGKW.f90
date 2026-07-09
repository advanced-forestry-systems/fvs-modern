!  GREGKW -- keyword-selected Greg codes (DGDRIVER / MORTDRVR keywords).
!    IDGDRV   -- DG site driver code, -1 unset; 0 none,1 elev,2 cspi,3 bgi,4 esi,5 emt
!    IMORTDRV -- mort coefficient TIER code (NOT a driver family; additive
!                gompit tiers), -1 unset; 0 crown-only,1 size,2 size+BGI
INTEGER IDGDRV, IMORTDRV
COMMON /GREGKW/ IDGDRV, IMORTDRV
