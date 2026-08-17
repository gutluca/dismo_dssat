C=======================================================================
C  Disease Impact and Severity Module
C  DISMO, Subroutine, Gustavo de Angelo Luca, Izael Martins Fattori Jr, Fábio Ricardo Marin
C  Luiz de Queiroz College of Agriculture (ESALQ), University of São Paulo, Piracicaba, Brazil
C  
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  07/17/2023 Written.
C  11/15/2024 Revised.
C  05/20/2025 Fungicide logic. 
C  08/10/2025 Logic/robustness fixes (cohorts, IR, LAF, LWD, cum. LAI)
C  11/05/2025 Write func for "DISMO.OUT" improved
C  11/17/2025 Virtual lesions factor added.
C  12/09/2025 Severity calculation in output file
C  03/05/2026 Removed hardcoded paths and added logic to find parameter file in DSSAT folders.
C  03/10/2026 Moved DISMO.for from Plant\CROPGRO to Plant\Generic-Pest
C  06/25/2026 Improved output file formatting
C  07/17/2026 Added defoliation/senescence logic (disease-induced senescence)
C  07/28/2026 Added monocyclic disease support (NCYCLE parameter: M/P)
C  08/16/2026 Added pre-plant environmental inoculum reconstruction.
C-----------------------------------------------------------------------
      SUBROUTINE DISEASE_LEAF (DYNAMIC,
     &    CONTROL, ISWITCH, Tmin, Tmax, RH, LAI_TOTAL,    ! Input
     &    WTLF, SLDOT,
     &    ESP_LAT_HIST, SUP_INF_LIST, LAI_INF_LIST,       ! Input/State
     &    YRDOY, YREMRG, NVEG0, YREND,                    ! Input
     &    DISEASE_LAI, VIRTUAL_PHOTO_FACTOR,              ! Output
     &    DISEASE_SEN_RATE)                               ! Output
C-----------------------------------------------------------------------
          USE ModuleDefs     
          IMPLICIT NONE
          EXTERNAL F_IR, F_DS, F_CANSPO, F_LR, F_LS, F_PS, F_LAF,
     &             F_LAR, F_PPSR,
     &             F_PPSR_POP, APPLY_FUNGICIDE, CALC_DVIP,
     &             READ_DISEASE_PARAMETERS, F_VIRTUAL_LESIONS,
     &             F_SEVERITY, GETLUN, F_DEFOLIATION,
     &             DISMO_PRESEASON, DISMO_UPDATE_ENVIRONMENT
          SAVE
C-----------------------------------------------------------------------
C  Switches / constants
          LOGICAL, PARAMETER :: USE_FUNGICIDE  = .FALSE.
          LOGICAL, PARAMETER :: USE_WTH_RH     = .TRUE.
          REAL,    PARAMETER :: LAI_MIN_START  = 0.5
          INTEGER, PARAMETER :: MAXDAYS        = 200
          REAL,    PARAMETER :: EPS            = 1.0E-6
C-----------------------------------------------------------------------
          
          INTEGER DAS, DYNAMIC
          
          REAL RH, LWD
          
          REAL FT_D, FT, T, Tmin, Tmax
          REAL FT_G
          REAL TMIN_G, TOT_G, TMAX_G
          REAL TMIN_D, TOT_D, TMAX_D
      
          REAL IR
          REAL FSS, LAI
          REAL LR
          REAL IS
          REAL LA, LAF, LESIONAGEOPT
          REAL PREV_IS
          REAL LAI_TOTAL, DISEASE_LAI
          
          INTEGER k, DAE, DAE_START
          INTEGER KMAX, DAE_IDX, PREV_IDX
          
          INTEGER YRDOY, YREMRG, PLANT_LIVE, NVEG0, YREND
          INTEGER DISEASE_LIVE
          INTEGER IYEAR, IDOY, IDAP
      
          REAL NDS, LESION_S, KVERHULST, RVERHULST
          REAL YMAX, COF_A, COF_B
          REAL LDMIN, LESLIFEMAX
          REAL Lesion_Rate
          REAL ESP_INOC_SEC
          REAL HEALTH_LAI, HEALTH_LAI_AVAIL, NEW_LOSS_TODAY
          REAL HEALTH_LAI_EPI
          REAL INF_AREA_K
          REAL s, beta, fvl, VIRTUAL_PHOTO_FACTOR
          REAL WTLF, SLDOT, DISEASE_SEN_RATE
          REAL rrds
          CHARACTER(LEN=1), SAVE :: NCYCLE        ! 'M' or 'P' from parameter file
          LOGICAL,          SAVE :: IS_MONOCYCLIC ! .T
          
C--------- Primary inoculum build-up -----------          
          REAL DAILY_IP
          REAL, DIMENSION(60), SAVE :: REGIONAL_HISTORY
          INTEGER, SAVE :: HIST_IDX
          REAL, SAVE :: SOURCE_PRESSURE
          REAL, SAVE :: SPOR_CLOUD
          REAL, PARAMETER :: SPOR_DECAY = 0.7937
          LOGICAL, SAVE :: PRESEASON_DONE
          
          REAL, SAVE :: SEC_SPORE_CLOUD
          REAL, SAVE :: SEC_SPORES_PENDING
          REAL, PARAMETER :: SEC_DECAY = 0.7937
          REAL, PARAMETER :: SEC_RELEASE = 1.0
          
          REAL PRI_CLOUD, SEC_CLOUD, CLOUD_TOTAL
          REAL DS_TOTAL, DS_PRI, DS_SEC
          REAL LS_TODAY
          
C--------- Population (individuals) & potential rate per area -----------
          REAL INF_COUNT_PREV, NPREV_POP
          REAL POT_SPO_PER_AREA, PS_K
          
          REAL, DIMENSION(MAXDAYS,5) :: ESP_LAT_HIST 
          REAL, DIMENSION(MAXDAYS)   :: SUP_INF_LIST
          REAL, DIMENSION(MAXDAYS)   :: LAI_INF_LIST
          
          REAL,    DIMENSION(MAXDAYS) :: ADMITTED_AREA
          
          INTEGER DVIP_pts(7)  
          INTEGER idx, SUM7, BufferDays, ResidualDays, NSprays
          LOGICAL FungActive
          INTEGER DVIP_today
          REAL FUNG_EFFICIENCY
          REAL, SAVE :: LAI_PEAK_SEASON
          REAL, SAVE :: SEVERITY_PCT
          INTEGER, SAVE :: LUN_OUT
          LOGICAL, SAVE :: HDR_DONE
          
!-----------------------------------------------------------------------
!         Constructed types
          
          TYPE (ControlType) CONTROL
          TYPE (SwitchType)  ISWITCH
          
          !DYNAMIC = CONTROL % DYNAMIC
          DAS     = CONTROL % DAS

!***********************************************************************
!  RUNINIT — called once (season start)
!***********************************************************************
          
      IF (DYNAMIC .EQ. RUNINIT) THEN

          IS = 0.0
          PLANT_LIVE = 0
          DISEASE_LIVE = 0
          DAE = 0

          ESP_LAT_HIST = 0.0
          SUP_INF_LIST = 0.0
          LAI_INF_LIST = 0.0
          ADMITTED_AREA= 0.0
          
          DISEASE_LAI  = 0.0
          DISEASE_SEN_RATE = 0.0
          
          idx = 1
          DVIP_pts(:) = 0
          SUM7 = 0
          BufferDays = 0
          ResidualDays = 0
          FungActive = .FALSE.
          NSprays = 0
          FUNG_EFFICIENCY = 0.723
          LAI_PEAK_SEASON = 0.0
          SEVERITY_PCT    = 0.0
          REGIONAL_HISTORY = 0.0
          HIST_IDX = 1
          SOURCE_PRESSURE  = 0.0
          SPOR_CLOUD = 0.0
          PRESEASON_DONE = .FALSE.
          
          SEC_SPORE_CLOUD   = 0.0
          SEC_SPORES_PENDING = 0.0

          NCYCLE      = 'P'
          IS_MONOCYCLIC = .FALSE.
          VIRTUAL_PHOTO_FACTOR = 1.0
          HDR_DONE = .FALSE.
          
          CALL READ_DISEASE_PARAMETERS(CONTROL,NDS, LESION_S, 
     &              KVERHULST,YMAX, COF_A, COF_B, RVERHULST, TMIN_G,
     &              TOT_G,TMAX_G,TMIN_D, TOT_D, TMAX_D, LDmin,
     &              LESIONAGEOPT,LESLIFEMAX,DAE_START, beta,rrds,NCYCLE)
               IS_MONOCYCLIC = (NCYCLE .EQ. 'M')
               
               
          IF (CONTROL%RUN .EQ. 1) THEN
              CALL GETLUN('DISOUT',LUN_OUT)
              OPEN(LUN_OUT, FILE='DISMO.OUT',
     &             STATUS='REPLACE')
              WRITE(LUN_OUT,'(A)')
     &             '*DISEASE IMPACT AND SEVERITY MODULE OUTPUT FILE'
          ELSE
              OPEN(LUN_OUT, FILE='DISMO.OUT',
     &             STATUS='OLD', ACCESS='APPEND')
           ENDIF  
!***********************************************************************
!  SEASINIT -- reconstruct environmental inoculum before planting
!***********************************************************************

      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN

          IF (.NOT. PRESEASON_DONE) THEN
              CALL DISMO_PRESEASON(CONTROL, NDS, TMIN_G, TOT_G,
     &             TMAX_G, TMIN_D, TOT_D, TMAX_D, USE_WTH_RH,
     &             SPOR_DECAY, SOURCE_PRESSURE, SPOR_CLOUD,
     &             REGIONAL_HISTORY, HIST_IDX)
              PRESEASON_DONE = .TRUE.
          END IF
          
!***********************************************************************
!  RATE — called every day
!***********************************************************************
          
      ELSEIF (DYNAMIC .EQ. RATE) THEN
          
          !----- Cloud spore dynamic and primary inoculum build-up --------

          CALL DISMO_UPDATE_ENVIRONMENT(Tmin, Tmax, RH, USE_WTH_RH,
     &         NDS, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &         SPOR_DECAY, DAILY_IP, SOURCE_PRESSURE, SPOR_CLOUD,
     &         T, LWD, FT, FT_D, FT_G, REGIONAL_HISTORY, HIST_IDX)
          
          
          SEC_SPORE_CLOUD = SEC_SPORE_CLOUD * SEC_DECAY +
     &                      SEC_SPORES_PENDING * SEC_RELEASE

          SEC_SPORES_PENDING = 0.0

!----- Emergence gate ---------------------------------------------------
          
          IF ((CONTROL%DAS .GT. NVEG0) .AND. (PLANT_LIVE .EQ. 0)) THEN
               DAE = 1
               PLANT_LIVE   = 1
               DISEASE_LIVE = 0
               
               IS = 0.0
               ESP_LAT_HIST = 0.0
               SUP_INF_LIST = 0.0
               LAI_INF_LIST = 0.0
               ADMITTED_AREA= 0.0
           END IF
              
!----- Disease activation (starts after DAE_START) ----------------------
          
          IF ((PLANT_LIVE .EQ. 1) .AND. (DAE .GE. DAE_START)) THEN 
              DISEASE_LIVE = 1
          END IF
          
          IF (DISEASE_LIVE .EQ. 1) THEN
              
!----- Healthy LAI from yesterday’s cumulative loss ---------------------
              
              IF (DAE .GT. 1) THEN
                 PREV_IDX   = MIN(DAE-1,MAXDAYS)
                 HEALTH_LAI = LAI_TOTAL - LAI_INF_LIST(PREV_IDX)
              ELSE
                 HEALTH_LAI = LAI_TOTAL
                 PREV_IDX   = 1
              END IF
              HEALTH_LAI = MAX(0.0, HEALTH_LAI)
              
              IF (LAI_PEAK_SEASON .GT. EPS .AND. DAE .GT. 1) THEN
                  HEALTH_LAI_EPI = MAX(0.0, LAI_PEAK_SEASON -
     &                                 LAI_INF_LIST(PREV_IDX))
              ELSE
                  HEALTH_LAI_EPI = LAI_TOTAL
              END IF

!----- Fungicide decision module (optional) -----------------------------
              
              CALL CALC_DVIP(LWD, T, DVIP_today)
              CALL APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7,
     &            BufferDays, FungActive, ResidualDays, NSprays,
     &            USE_FUNGICIDE, HEALTH_LAI)
              
!----- Deposition of primary and secondary inoculum ---------------------
              ! Primary deposition requires the preseason environmental gate.
              ! Secondary inoculum is local and therefore does not use that gate.
              PRI_CLOUD = 0.0
              IF (SOURCE_PRESSURE .GE. (NDS / MAX(YMAX, EPS))) THEN
                  PRI_CLOUD = SPOR_CLOUD
              END IF
              SEC_CLOUD   = SEC_SPORE_CLOUD
              CLOUD_TOTAL = PRI_CLOUD + SEC_CLOUD
              DS_TOTAL = 0.0
              DS_PRI   = 0.0
              DS_SEC   = 0.0
              FSS      = 0.0

              ! Both sources compete for the healthy leaf area present today.
              IF (CLOUD_TOTAL .GT. EPS .AND. HEALTH_LAI .GT. EPS) THEN
                  CALL F_CANSPO(HEALTH_LAI, CLOUD_TOTAL, LESION_S,FSS)
                  CALL F_DS(FSS, CLOUD_TOTAL, DS_TOTAL)
                  DS_PRI = DS_TOTAL * PRI_CLOUD / CLOUD_TOTAL
                  DS_SEC = DS_TOTAL * SEC_CLOUD / CLOUD_TOTAL
                  SPOR_CLOUD = MAX(SPOR_CLOUD - DS_PRI, 0.0)
                  SEC_SPORE_CLOUD = MAX(SEC_SPORE_CLOUD - DS_SEC,0.0)
              END IF
              
              CALL F_IR(FT, LWD, YMAX, COF_A, COF_B, IR)
              
              IF (FungActive) IR = IR * (1.0 - FUNG_EFFICIENCY)

!----- No healthy leaf -> no successful infections today ----------------
              
              IF (HEALTH_LAI .LE. 0.0) THEN
                  IR = 0.0
                  DS_TOTAL = 0.0
              END IF
              
              CALL F_LS(IR, DS_TOTAL, LS_TODAY)

!----- Current day accumulators ------------------------------------------
              
              DAE_IDX = MIN(DAE,MAXDAYS)
              IS             = 0.0      ! infectious surface (m2 m-2) today
              ESP_INOC_SEC   = 0.0      ! total secondary spores today
              NEW_LOSS_TODAY = 0.0      ! newly activated surface today

!----- Yesterday carry-over for population/logistic ---------------------
              
              PREV_IS  = 0.0
              IF (DAE .GT. 1) PREV_IS = SUP_INF_LIST(PREV_IDX)

!----- Lesion aging rate (constant within the day) ----------------------
              
              CALL F_LAR (FT_D, LESLIFEMAX, Lesion_Rate)

!----- Build yesterday's population (individuals) -----------------------
              
          POT_SPO_PER_AREA = 0.0
          IF (.NOT. IS_MONOCYCLIC) THEN
                  
              INF_COUNT_PREV = 0.0
              IF (DAE .GT. 1) THEN
                 INF_COUNT_PREV = PREV_IS / MAX(LESION_S, EPS)
              END IF

              ! Only infectious lesions are a sporulating population.
              NPREV_POP = MAX(INF_COUNT_PREV,0.0)

!----- Potential spore rate per unit sporulating area (yesterday) -------
              
              CALL F_PPSR_POP(KVERHULST, RVERHULST, NPREV_POP,
     &                        PREV_IS, HEALTH_LAI_EPI, POT_SPO_PER_AREA)
              
          END IF

!----- Cohort loop -------------------------------------------------------
              
              KMAX = MIN(DAE,MAXDAYS)
              DO k = 1, KMAX

                  CALL F_LR (FT_D, LDmin, LR)
                  ESP_LAT_HIST(k, 2) = ESP_LAT_HIST(k, 2) + LR
                  
!--------- Latent -> infectious (happens once) --------------------------
                  
                  IF ((ESP_LAT_HIST(k, 2) .GE. 1.0) .AND.
     &                (ESP_LAT_HIST(k, 3) .EQ. 0.0)) THEN
                       ESP_LAT_HIST(k, 3) = 1.0
                  ! allocate effective area respecting remaining healthy leaf
                       INF_AREA_K = ESP_LAT_HIST(k,1) * LESION_S
                       INF_AREA_K = MAX(INF_AREA_K, 0.0)
                       HEALTH_LAI_AVAIL = MAX(HEALTH_LAI_EPI - 
     &                                         NEW_LOSS_TODAY,0.0)
                       ADMITTED_AREA(k)=MIN(INF_AREA_K,HEALTH_LAI_AVAIL)
                       NEW_LOSS_TODAY=NEW_LOSS_TODAY + ADMITTED_AREA(k)
                  END IF
                  
!--------- After infectious ------------------------------------------------
     
                  IF (ESP_LAT_HIST(k, 3) .EQ. 1.0) THEN 
                      ESP_LAT_HIST(k,4) = ESP_LAT_HIST(k,4) +Lesion_Rate
                      LA = ESP_LAT_HIST(k,4)

                      CALL F_LAF(LA, LAF, LESIONAGEOPT)

!--------- Infectious area of the cohort today (only while LA<=1) -------
                      
                      INF_AREA_K = 0.0
                      IF (LA .LE. 1.0) INF_AREA_K = ADMITTED_AREA(k)
                      IF (LAI_TOTAL .LT. LAI_MIN_START) INF_AREA_K = 0.0

                      IS = IS + INF_AREA_K
                      
!--------- Secondary production while LA<=1 and substrate remains -------
                      
                      IF (LA .LE. 1.0 .AND. .NOT. IS_MONOCYCLIC) THEN
                          LAI  = MAX(HEALTH_LAI_EPI-NEW_LOSS_TODAY, 0.0)
                          IF(LAI.GT.0.0.AND.POT_SPO_PER_AREA.GT.0.0)THEN
                              PS_K = POT_SPO_PER_AREA * INF_AREA_K * 
     &                        FT_D * LAF
                              ESP_INOC_SEC =ESP_INOC_SEC + MAX(PS_K,0.0)
                          END IF
                      END IF
                  END IF
              END DO

!----- Store infectious surface and route secondary inoculum to cloud ---
              
              SUP_INF_LIST(DAE_IDX) = IS
              
          IF (.NOT. IS_MONOCYCLIC) THEN
              ! Emitted spores become available for deposition tomorrow.
              SEC_SPORES_PENDING = MAX(ESP_INOC_SEC, 0.0)
              ESP_LAT_HIST(DAE_IDX, 5) = SEC_SPORES_PENDING
          ELSE
              SEC_SPORES_PENDING = 0.0
          END IF

!----- Register infections only after all existing cohorts advance ------
              ESP_LAT_HIST(DAE_IDX, 1) = ESP_LAT_HIST(DAE_IDX, 1) +
     &                                    LS_TODAY
              
!----- Cumulative removed LAI (m2 m-2) ----------------------------------
              
              IF (DAE .GT. 1) THEN
                 LAI_INF_LIST(DAE_IDX) = LAI_INF_LIST(PREV_IDX) + 
     &                                   NEW_LOSS_TODAY
              ELSE
                 LAI_INF_LIST(DAE_IDX) = NEW_LOSS_TODAY
              END IF
              
              IF (LAI_PEAK_SEASON .GT. EPS) THEN
                  LAI_INF_LIST(DAE_IDX) = MIN(LAI_INF_LIST(DAE_IDX),
     &                                        LAI_PEAK_SEASON)
              END IF

!----- Internal var uses cm2 m-2, output will divide by 10000 --
              
              !----- Effective DISEASE_LAI: scale by severity fraction applied to
              !      current LAI. Natural senescence removes diseased and healthy
              !      tissue proportionally.
              IF (LAI_PEAK_SEASON .GT. EPS) THEN
                  DISEASE_LAI = (LAI_INF_LIST(DAE_IDX) / 
     &                           LAI_PEAK_SEASON) *
     &                          LAI_TOTAL * 10000.0
                  DISEASE_LAI = MIN(DISEASE_LAI,
     &                              LAI_TOTAL * 10000.0)
              ELSE
                  DISEASE_LAI = 0.0
              END IF
              DISEASE_LAI = MAX(DISEASE_LAI, 0.0)
              
!----- Calculate disease severity --
              
              CALL F_SEVERITY(LAI_TOTAL, LAI_INF_LIST(DAE_IDX), 
     &                        LAI_PEAK_SEASON, SEVERITY_PCT)
              
!----- Virtual Lesion logic
          ! Daily fractional severity (only today's new necrosis)
          IF (HEALTH_LAI .GT. EPS) THEN
              s = NEW_LOSS_TODAY / HEALTH_LAI
          ELSE
              s = 0.0
          END IF
          s = MAX(0.0, MIN(s,1.0))

          ! Call subroutine to compute virtual lesion reduction factor (0..1)
          CALL F_VIRTUAL_LESIONS(s, beta, fvl)

          ! Store result in a global variable to be exported to CROPGRO
          VIRTUAL_PHOTO_FACTOR = fvl
          
!----- Defoliation / Senescence Logic
          IF (DISEASE_LIVE .EQ. 1) THEN
              CALL F_DEFOLIATION (WTLF, SLDOT, SEVERITY_PCT, rrds,
     &                           DISEASE_SEN_RATE)
          ELSE
              DISEASE_SEN_RATE = 0.0
          END IF
          
!----- ------------------------------------------------------          
          ELSE
              HEALTH_LAI = LAI_TOTAL
          END IF
          
!----- Advance time ------------------------------------------------------
          
          IF (PLANT_LIVE .EQ. 1) DAE = DAE + 1
          IF (YREND .EQ. YRDOY)  PLANT_LIVE = 0
          
!***********************************************************************
!  OUTPUT — called on output events
!***********************************************************************
          
      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN

!----- Decompose YRDOY into YEAR and DOY --------------------------------
            IYEAR = YRDOY / 1000
            IDOY  = YRDOY - IYEAR * 1000
            IF (IYEAR .LT. 100) THEN
                IF (IYEAR .GT. 50) THEN
                    IYEAR = IYEAR + 1900
                ELSE
                    IYEAR = IYEAR + 2000
                END IF
            END IF

!----- Days after emergence -----------------------------------------------
            IDAP = MAX(DAE - 1, 0)

!----- Write run header once per run -------------------------------------
            IF (.NOT. HDR_DONE) THEN
              HDR_DONE = .TRUE.
              WRITE(LUN_OUT,'(A)') ' '
              WRITE(LUN_OUT,'(A,I4,A,A,A,A,1X,A,I5)')
     &         '*RUN',CONTROL%RUN,'        : ',
     &         CONTROL%ENAME,
     &         '                      ',
     &         CONTROL%MODEL,
     &         CONTROL%TRTNUM
              WRITE(LUN_OUT,'(A,A)')
     &         ' MODEL          : ', CONTROL%MODEL
              WRITE(LUN_OUT,'(A,A)')
     &         ' EXPERIMENT     : ', CONTROL%FILEX
              WRITE(LUN_OUT,'(A)')
     &         ' DATA PATH      :'
              WRITE(LUN_OUT,'(A,I2,A,A,A,A)')
     &         ' TREATMENT',CONTROL%TRTNUM,
     &         '    : ',CONTROL%ENAME,
     &         '                      ',CONTROL%MODEL
              WRITE(LUN_OUT,'(A)') ' '
              WRITE(LUN_OUT,'(A)') '!'
              WRITE(LUN_OUT,'(A)') '!'
              WRITE(LUN_OUT, 25)
   25         FORMAT('@YEAR  DOY   DAS   DAE',
     &         '      LAIH      LWDh      RHU%      FTMP      LAIT',
     &         '      SUM7 NSPRAYS     FACT     SEV%')
            ENDIF

            WRITE(LUN_OUT,
     &       '(I5,I5,I6,I6,5F10.3,F10.1,I8,L9,F9.2)')
     &        IYEAR, IDOY, CONTROL%DAS, IDAP,
     &        HEALTH_LAI,
     &        LWD, RH, FT, LAI_TOTAL, REAL(SUM7),
     &        NSprays, FungActive, SEVERITY_PCT

!***********************************************************************
!  SEASEND — called once (season end)
!***********************************************************************
            
      ELSEIF (DYNAMIC .EQ. SEASEND) THEN
            PLANT_LIVE   = 0
            DISEASE_LIVE = 0
            ESP_LAT_HIST = 0.0
            SUP_INF_LIST = 0.0
            LAI_INF_LIST = 0.0
            ADMITTED_AREA= 0.0
            DISEASE_LAI  = 0.0
            DAE = 1
            idx = 1
            DVIP_pts(:) = 0
            SUM7 = 0
            BufferDays = 0
            ResidualDays = 0
            !FungActive = .FALSE.
            NSprays = 0
            LAI_PEAK_SEASON = 0.0
            SEVERITY_PCT    = 0.0
            DISEASE_SEN_RATE = 0.0
            
            SEC_SPORE_CLOUD   = 0.0
            SEC_SPORES_PENDING = 0.0
            CLOSE(LUN_OUT)
            	 
      ENDIF

      END SUBROUTINE DISEASE_LEAF

!-----------------------------------------------------------------------
!  SUBROUTINES
!-----------------------------------------------------------------------

! ------ READ PARAMETERS
!  keep parameter file consistent with columns order below.
     
      SUBROUTINE READ_DISEASE_PARAMETERS(CONTROL,NDS,LESION_S,KVERHULST,
     &                                  YMAX, COF_A, COF_B, 
     &                                  RVERHULST, TMIN_G, TOT_G,TMAX_G,
     &                                  TMIN_D, TOT_D, TMAX_D,    
     &                                  LDmin, LESIONAGEOPT, LESLIFEMAX,
     &                                  DAE_START, beta, rrds, NCYCLE)

          USE ModuleDefs
          IMPLICIT NONE
          EXTERNAL GETLUN
          
          TYPE (ControlType) CONTROL
          
          REAL    NDS, LESION_S, KVERHULST, RVERHULST
          REAL    YMAX, COF_A, COF_B
          REAL    TMIN_G, TOT_G, TMAX_G
          REAL    TMIN_D, TOT_D, TMAX_D
          REAL    LDmin, LESIONAGEOPT, LESLIFEMAX, beta, rrds
          INTEGER DAE_START

          CHARACTER(LEN=8)   ID        
          CHARACTER(LEN=20)  VRNAME    
          CHARACTER(LEN=120) DISFIL
          CHARACTER(LEN=400) LINE
          CHARACTER(LEN=20)  TARGET_DISEASE
           CHARACTER(LEN=1) NCYCLE
          INTEGER LUN_DIS, IOS, STATE
          LOGICAL FEXIST, FOUND_DISEASE

!         Only look in the local (current working) directory
          DISFIL = 'disease_parameters.txt'
          INQUIRE(FILE=TRIM(DISFIL), EXIST=FEXIST)
          IF (.NOT. FEXIST) RETURN

          CALL GETLUN('DISINP', LUN_DIS)
          
          OPEN(LUN_DIS, FILE=TRIM(DISFIL), 
     &         STATUS='OLD', ACTION='READ', IOSTAT=IOS)

          IF (IOS /= 0) RETURN

!         --- State-machine driven parsing ---
!         STATE 0 : searching for section tags
!         STATE 1 : found *DISEASE CONTROL; waiting for @TARGET_DISEASE header
!         STATE 2 : next valid line is the target disease name value
!         STATE 3 : found *DISEASE DATABASE; waiting for @VAR# column header
!         STATE 4 : reading database records

          TARGET_DISEASE = ' '
          STATE          = 0
          FOUND_DISEASE  = .FALSE.
          NCYCLE         = 'P' !default to policyclic 

          DO WHILE (.TRUE.)
              READ(LUN_DIS, '(A)', IOSTAT=IOS) LINE
              IF (IOS /= 0) EXIT          ! EOF or read error

              LINE = ADJUSTL(LINE)
              IF (LEN_TRIM(LINE) .EQ. 0)  CYCLE   ! skip blank lines
              IF (LINE(1:1) .EQ. '!')     CYCLE   ! skip comment lines

              SELECT CASE (STATE)

              CASE (0)  ! searching for section tags
                  IF (INDEX(LINE,'*DISEASE CONTROL') .GT. 0) THEN
                      STATE = 1
                  ELSEIF (INDEX(LINE,'*DISEASE DATABASE') .GT. 0) THEN
                      STATE = 3
                  END IF

              CASE (1)  ! found *DISEASE CONTROL; wait for @TARGET_DISEASE header
                  IF (LINE(1:1) .EQ. '@') STATE = 2

              CASE (2)  ! next valid line is the target disease name
                  TARGET_DISEASE = TRIM(LINE)
                  STATE = 0             ! continue searching for *DISEASE DATABASE

              CASE (3)  ! found *DISEASE DATABASE; wait for @VAR# column header
                  IF (LINE(1:1) .EQ. '@') STATE = 4

              CASE (4)  ! reading database records
                  READ(LINE, *, IOSTAT=IOS) ID, VRNAME,
     &                NDS, LESION_S, KVERHULST, RVERHULST,
     &                YMAX, COF_A, COF_B,
     &                TMIN_G, TOT_G, TMAX_G,
     &                TMIN_D, TOT_D, TMAX_D,
     &                LDmin, LESIONAGEOPT, LESLIFEMAX,
     &                DAE_START, beta, rrds, NCYCLE
                  IF (IOS .EQ. 0 .AND.
     &                TRIM(VRNAME) .EQ. TRIM(TARGET_DISEASE)) THEN
                      FOUND_DISEASE = .TRUE.
                      EXIT
                  END IF

              END SELECT
          END DO

          CLOSE(LUN_DIS)

          IF (.NOT. FOUND_DISEASE) THEN
              WRITE(*,'(A,A)')
     &            'WARNING (DISMO): disease not found in database: ',
     &            TRIM(TARGET_DISEASE)
          END IF

       END SUBROUTINE READ_DISEASE_PARAMETERS

!-----------------------------------------------------------------------
!  PRE-PLANT ENVIRONMENTAL INOCULUM RECONSTRUCTION
!
!  CROPGRO does not call PEST during RATE before planting.  To keep the
!  DISMO integration self-contained, this routine replays weather from
!  SDATE through the day before PDATE at SEASINIT.  The resolved weather
!  file and planting date are read from DSSAT48.INP, which is generated
!  by the DSSAT input module for the active treatment.
!-----------------------------------------------------------------------

      SUBROUTINE DISMO_PRESEASON(CONTROL, NDS, TMIN_G, TOT_G, TMAX_G,
     &    TMIN_D, TOT_D, TMAX_D, USE_WTH_RH, SPOR_DECAY,SOURCE_PRESSURE,
     &    SPOR_CLOUD, REGIONAL_HISTORY, HIST_IDX)

      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL GETLUN, DISMO_UPDATE_ENVIRONMENT,
     &         DISMO_NORMALIZE_DATE, DISMO_READ_WTH_FIELD

      TYPE (ControlType) CONTROL

      REAL NDS, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D
      REAL SPOR_DECAY, SOURCE_PRESSURE, SPOR_CLOUD
      REAL WTMAX, WTMIN, WRH, DAILY_IP, T, LWD, FT, FT_D, FT_G
      REAL, DIMENSION(60) :: REGIONAL_HISTORY
      INTEGER HIST_IDX
      LOGICAL USE_WTH_RH, FOUND_HEADER, IN_PLANTING
      LOGICAL OK_TMAX, OK_TMIN, OK_RH
      INTEGER LUN_IO, LUN_WTH, IOS, WTH_DATE, FULL_DATE
      INTEGER POS_TMAX, POS_TMIN, POS_RAIN, POS_RHUM
      INTEGER PRESEASON_PDATE
      CHARACTER(LEN=400) LINE, WTH_HEADER
      CHARACTER(LEN=30) WTH_NAME
      CHARACTER(LEN=120) WTH_PATH
      CHARACTER(LEN=240) WTH_FILE

      PRESEASON_PDATE = -99
      WTH_NAME = ' '
      WTH_PATH = ' '
      IN_PLANTING = .FALSE.

!----- Read resolved planting date and weather location -----------------

      CALL GETLUN('DSMINP', LUN_IO)
      OPEN(LUN_IO, FILE=TRIM(CONTROL%FILEIO), STATUS='OLD',
     &     ACTION='READ', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
          WRITE(*,'(A)')
     &      'WARNING (DISMO): Cannot open DSSAT input file; '
     &      //'preseason reconstruction skipped.'
          RETURN
      END IF

      DO WHILE (.TRUE.)
          READ(LUN_IO, '(A)', IOSTAT=IOS) LINE
          IF (IOS .NE. 0) EXIT
          LINE = ADJUSTL(LINE)
          IF (LEN_TRIM(LINE) .EQ. 0) CYCLE

          IF (INDEX(LINE, 'WEATHERW') .EQ. 1) THEN
              READ(LINE(9:), *, IOSTAT=IOS) WTH_NAME, WTH_PATH
              CYCLE
          END IF

          IF (INDEX(LINE, '*PLANTING DETAILS') .EQ. 1) THEN
              IN_PLANTING = .TRUE.
              CYCLE
          END IF

          IF (IN_PLANTING) THEN
              IF (LINE(1:1) .EQ. '*') THEN
                  IN_PLANTING = .FALSE.
              ELSE
                  READ(LINE, *, IOSTAT=IOS) PRESEASON_PDATE
                  IF (IOS .EQ. 0 .AND. PRESEASON_PDATE .GT. 0) THEN
                      IN_PLANTING = .FALSE.
                  END IF
              END IF
          END IF
      END DO
      CLOSE(LUN_IO)

      IF (LEN_TRIM(WTH_NAME) .EQ. 0 .OR.
     &    PRESEASON_PDATE .LE. CONTROL%YRSIM) THEN
          WRITE(*,'(A)')
     &      'WARNING (DISMO): Preseason weather or planting date '
     &      //'was not found; reconstruction skipped.'
          RETURN
      END IF

      IF (LEN_TRIM(WTH_PATH) .EQ. 0) THEN
          WTH_FILE = TRIM(WTH_NAME)
      ELSEIF (WTH_PATH(LEN_TRIM(WTH_PATH):LEN_TRIM(WTH_PATH))
     &         .EQ. CHAR(92) .OR.
     &         WTH_PATH(LEN_TRIM(WTH_PATH):LEN_TRIM(WTH_PATH))
     &         .EQ. '/') THEN
          WTH_FILE = TRIM(WTH_PATH)//TRIM(WTH_NAME)
      ELSE
          WTH_FILE = TRIM(WTH_PATH)//CHAR(92)//TRIM(WTH_NAME)
      END IF

!----- Read daily weather with positions defined by the WTH header ------

      CALL GETLUN('DSMWTH', LUN_WTH)
      OPEN(LUN_WTH, FILE=TRIM(WTH_FILE), STATUS='OLD',
     &     ACTION='READ', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
          WRITE(*,'(A,A)')
     &      'WARNING (DISMO): Cannot open weather file: ',
     &      TRIM(WTH_FILE)
          RETURN
      END IF

      FOUND_HEADER = .FALSE.
      POS_TMAX = 0
      POS_TMIN = 0
      POS_RAIN = 0
      POS_RHUM = 0

      DO WHILE (.TRUE.)
          READ(LUN_WTH, '(A)', IOSTAT=IOS) LINE
          IF (IOS .NE. 0) EXIT
          LINE = ADJUSTL(LINE)

          IF (.NOT. FOUND_HEADER) THEN
              IF (INDEX(LINE, '@DATE') .EQ. 1) THEN
                  WTH_HEADER = LINE
                  POS_TMAX = INDEX(WTH_HEADER, 'TMAX')
                  POS_TMIN = INDEX(WTH_HEADER, 'TMIN')
                  POS_RAIN = INDEX(WTH_HEADER, 'RAIN')
                  POS_RHUM = INDEX(WTH_HEADER, 'RHUM')
                  IF (POS_TMAX .GT. 0 .AND. POS_TMIN .GT. 0
     &                .AND. POS_RAIN .GT. 0) THEN
                      IF (USE_WTH_RH .AND. POS_RHUM .EQ. 0) THEN
                          WRITE(*,'(A)')
     &                    'WARNING (DISMO): RHUM is missing from WTH; '
     &                    //'preseason reconstruction skipped.'
                          CLOSE(LUN_WTH)
                          RETURN
                      END IF
                      FOUND_HEADER = .TRUE.
                  END IF
              END IF
              CYCLE
          END IF

          IF (LEN_TRIM(LINE) .EQ. 0 .OR. LINE(1:1) .EQ. '!') CYCLE
          READ(LINE(1:5), '(I5)', IOSTAT=IOS) WTH_DATE
          IF (IOS .NE. 0) CYCLE

          CALL DISMO_NORMALIZE_DATE(WTH_DATE, FULL_DATE)
          IF (FULL_DATE .LT. CONTROL%YRSIM .OR.
     &        FULL_DATE .GE. PRESEASON_PDATE) CYCLE

          CALL DISMO_READ_WTH_FIELD(LINE, POS_TMAX, POS_TMIN-1,
     &         WTMAX, OK_TMAX)
          CALL DISMO_READ_WTH_FIELD(LINE, POS_TMIN, POS_RAIN-1,
     &         WTMIN, OK_TMIN)
          IF (USE_WTH_RH) THEN
              CALL DISMO_READ_WTH_FIELD(LINE, POS_RHUM, LEN(LINE),
     &             WRH, OK_RH)
          ELSE
              WRH = 0.0
              OK_RH = .TRUE.
          END IF

          IF (OK_TMAX .AND. OK_TMIN .AND. OK_RH) THEN
              CALL DISMO_UPDATE_ENVIRONMENT(WTMIN, WTMAX, WRH,
     &             USE_WTH_RH, NDS, TMIN_G, TOT_G, TMAX_G,
     &             TMIN_D, TOT_D, TMAX_D, SPOR_DECAY, DAILY_IP,
     &             SOURCE_PRESSURE, SPOR_CLOUD, T, LWD, FT, FT_D, FT_G,
     &             REGIONAL_HISTORY, HIST_IDX)
          END IF
      END DO
      CLOSE(LUN_WTH)

      IF (.NOT. FOUND_HEADER) THEN
          WRITE(*,'(A,A)')
     &      'WARNING (DISMO): Cannot read WTH header from: ',
     &      TRIM(WTH_FILE)
      END IF

      END SUBROUTINE DISMO_PRESEASON

!-----------------------------------------------------------------------
!  DAILY ENVIRONMENTAL INOCULUM UPDATE
!-----------------------------------------------------------------------

      SUBROUTINE DISMO_UPDATE_ENVIRONMENT(TMIN, TMAX, RH,
     &    USE_WTH_RH, NDS, TMIN_G, TOT_G, TMAX_G,
     &    TMIN_D, TOT_D, TMAX_D, SPOR_DECAY, DAILY_IP,
     &    SOURCE_PRESSURE, SPOR_CLOUD, T, LWD, FT, FT_D, FT_G,
     &    REGIONAL_HISTORY, HIST_IDX)

      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL F_TAVG, F_DEW, F_RH, F_LWD, T_DEV

      REAL TMIN, TMAX, RH, NDS, TMIN_G, TOT_G, TMAX_G
      REAL TMIN_D, TOT_D, TMAX_D, SPOR_DECAY
      REAL DAILY_IP, SOURCE_PRESSURE, SPOR_CLOUD, T, LWD, FT, FT_D, FT_G
      REAL TDEW, ES, E, RH_LOCAL
      LOGICAL USE_WTH_RH
      
      REAL, DIMENSION(60) :: REGIONAL_HISTORY
      INTEGER :: HIST_IDX

      CALL F_TAVG(TMAX, TMIN, T)
      RH_LOCAL = RH
      IF (.NOT. USE_WTH_RH) THEN
          CALL F_DEW(T, TMIN, TMAX, TDEW)
          CALL F_RH(RH_LOCAL, TDEW, ES, E, T)
      END IF
      CALL F_LWD(RH_LOCAL, LWD)
      CALL T_DEV(T, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D,
     &    TMAX_D, FT, FT_D, FT_G)

      DAILY_IP = NDS * FT_G * (LWD / 24.0)
      SOURCE_PRESSURE = MAX(SOURCE_PRESSURE - 
     &                  REGIONAL_HISTORY(HIST_IDX), 0.0)
      
      REGIONAL_HISTORY(HIST_IDX) = DAILY_IP
      SOURCE_PRESSURE = SOURCE_PRESSURE + DAILY_IP
      
      HIST_IDX = HIST_IDX + 1
      IF (HIST_IDX .GT. 60) HIST_IDX = 1
      
      SPOR_CLOUD = (SPOR_CLOUD * SPOR_DECAY) + DAILY_IP

      END SUBROUTINE DISMO_UPDATE_ENVIRONMENT

!-----------------------------------------------------------------------
!  NORMALIZE DSSAT weather dates (YYDDD or YYYYDDD) to YYYYDDD.
!-----------------------------------------------------------------------

      SUBROUTINE DISMO_NORMALIZE_DATE(INPUT_DATE, FULL_DATE)

      IMPLICIT NONE
      INTEGER INPUT_DATE, FULL_DATE, IYEAR, IDOY

      IYEAR = INPUT_DATE / 1000
      IDOY = INPUT_DATE - IYEAR * 1000
      IF (IYEAR .LT. 100) THEN
          IF (IYEAR .GT. 50) THEN
              IYEAR = IYEAR + 1900
          ELSE
              IYEAR = IYEAR + 2000
          END IF
      END IF
      FULL_DATE = IYEAR * 1000 + IDOY

      END SUBROUTINE DISMO_NORMALIZE_DATE

!-----------------------------------------------------------------------
!  Read one fixed-width value from a DSSAT WTH record.
!-----------------------------------------------------------------------

      SUBROUTINE DISMO_READ_WTH_FIELD(LINE, FIRST, LAST, VALUE, OK)

      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER FIRST, LAST, LAST_POS, IOS
      REAL VALUE
      LOGICAL OK

      OK = .FALSE.
      VALUE = 0.0
      IF (FIRST .LT. 1 .OR. LAST .LT. FIRST) RETURN
      IF (FIRST .GT. LEN(LINE)) RETURN
      LAST_POS = MIN(LAST, LEN(LINE))
      READ(LINE(FIRST:LAST_POS), *, IOSTAT=IOS) VALUE
      IF (IOS .EQ. 0) OK = .TRUE.

      END SUBROUTINE DISMO_READ_WTH_FIELD

! ------ AVERAGE TEMPERATURE

      SUBROUTINE F_TAVG(Tmax, Tmin, T)
          USE ModuleDefs
          IMPLICIT NONE
          REAL Tmax, Tmin, T
          T = (Tmax + Tmin) / 2.0
      END SUBROUTINE F_TAVG
      
! ------ DEW POINT TEMPERATURE (empirical)
      SUBROUTINE F_DEW(T, Tmin, Tmax, Tdew)
          USE ModuleDefs
          IMPLICIT NONE
          REAL T, Tmin, Tmax, Tdew
          Tdew = (-0.036*T) + (0.9679*Tmin)+(0.0072*(Tmax-Tmin)) +1.0111
      END SUBROUTINE F_DEW
      
! ------ RELATIVE HUMIDITY (Tetens form)
      SUBROUTINE F_RH(RH, Tdew, Es, E, T)
          USE ModuleDefs
          IMPLICIT NONE
          REAL RH, Tdew, Es, E, T
          Es = EXP(17.625*T    /(243.04+T))
          E  = EXP(17.625*Tdew /(243.04+Tdew))
          RH = MAX((E / Es) * 100.0, 0.0)
      END SUBROUTINE F_RH
      
! ------ LEAF WETNESS DURATION 
      SUBROUTINE F_LWD(RH, LWD)
          USE ModuleDefs
          IMPLICIT NONE
          REAL RH, LWD
          REAL RH_LOC
          
          RH_LOC = RH
          IF (RH_LOC .LE. 1.0) RH_LOC = RH_LOC * 100.0
          RH_LOC = MIN(MAX(RH_LOC, 0.0), 100.0)
          
          LWD = MAX(31.31 / (1.0 + EXP(-((RH_LOC - 85.17) / 9.13))),0.0)
          LWD = MIN(LWD, 24.0)
      END SUBROUTINE F_LWD

! ------ TEMPERATURE DEVELOPMENT FUNCTION (beta response)
      SUBROUTINE T_DEV (T, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &                  FT, FT_D, FT_G)
      USE ModuleDefs
      IMPLICIT NONE
      REAL FT_D, FT_G, FT, T
      REAL TMIN_D, TOT_D, TMAX_D, TMIN_G, TOT_G, TMAX_G
 
      ! ---- FT_G ---- (spores germination)
      FT_G = 0.0
      IF ((TMAX_G-T) > 0.0 .AND. (TMAX_G-TOT_G) /= 0.0 .AND.
     &    (T-TMIN_G) > 0.0 .AND. (TOT_G-TMIN_G) /= 0.0) THEN
          FT_G = ((TMAX_G-T)/(TMAX_G-TOT_G)) *
     &    ((T-TMIN_G)/(TOT_G-TMIN_G))**((TOT_G-TMIN_G)/(TMAX_G-TOT_G))
          
      END IF

      ! ---- FT_D ---- (general fungal development)
      IF (T .GT. TMIN_D .AND. (TMAX_D-TOT_D) /= 0.0 .AND. 
     &      (TOT_D-TMIN_D) /= 0.0) THEN
          
         IF ((TMAX_D-T) > 0.0 .AND. (T-TMIN_D) > 0.0) THEN
             FT_D = ((TMAX_D-T)/(TMAX_D-TOT_D)) *
     &      ((T-TMIN_D)/(TOT_D-TMIN_D))**((TOT_D-TMIN_D)/(TMAX_D-TOT_D))
             
         ELSE
             FT_D = 0.0
         END IF
         
      ELSE
         FT_D = 0.0
      END IF

      ! --------
      FT = MAX(MIN(FT_G * FT_D, 1.0), 0.0)
      
      END SUBROUTINE T_DEV

! ------ INFECTION RATE
      SUBROUTINE F_IR(FT, LWD, Ymax, A, B, IR)
          USE ModuleDefs
          IMPLICIT NONE
          REAL Ymax, A, B, IR, LWD, FT
          IR = Ymax * FT * (1.0 - EXP(-(A * LWD)**B))
      END SUBROUTINE F_IR
    
! ------ CANOPY SUPPORTED SPORES (capacity fraction)
      SUBROUTINE F_CANSPO(LAI, NDS, LESION_S, FSS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FSS, LAI, NDS, LESION_S
          IF (NDS .LE. 0.0) THEN
              FSS = 0.0
          ELSE
              FSS = LAI / (NDS * LESION_S)
              FSS = MIN(MAX(FSS, 0.0), 1.0)
          END IF
      END SUBROUTINE F_CANSPO

! ------ DEPOSITED SPORES
      SUBROUTINE F_DS(FSS, NDS, DS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FSS, NDS, DS
          DS = FSS * NDS
      END SUBROUTINE F_DS

! ------ LATENCY RATE 
      SUBROUTINE F_LR(FT_D, LDmin, LR)
          USE ModuleDefs
          IMPLICIT NONE
          REAL LR, FT_D, LDmin
          LR = MAX(FT_D, 0.0) / MAX(LDmin, 1.0E-6)
      END SUBROUTINE F_LR
    
! ------ LATENT SPORES
      SUBROUTINE F_LS(IR, DS, LS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL IR, DS, LS
          LS = MAX(IR * DS, 0.0)
      END SUBROUTINE F_LS

! ------ SPORES PRODUCTION (secondary inoculum) 
      SUBROUTINE F_PS(PPSR, LESION_S, FT_D, LAF, PS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL PS, PPSR, LESION_S, FT_D, LAF
          PS = PPSR * LESION_S * FT_D * LAF
      END SUBROUTINE F_PS 
    
! ------ LESION AGE RATE
      SUBROUTINE F_LAR (FT_D, LESLIFEMAX, Lesion_Rate)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FT_D, LESLIFEMAX, Lesion_Rate
          IF (FT_D .LE. 0.0) THEN
              Lesion_Rate = 0.0
          ELSE
              Lesion_Rate = FT_D / MAX(LESLIFEMAX, 1.0E-6)
          END IF
      END SUBROUTINE F_LAR
          
! ------ LESION AGE FACTOR (0..1)
      SUBROUTINE F_LAF(LA, LAF, LESIONAGEOPT)
          USE ModuleDefs
          IMPLICIT NONE
          REAL LAF, LA, LESIONAGEOPT, TMP
          IF (LA .LT. LESIONAGEOPT) THEN
              TMP = LA / MAX(LESIONAGEOPT, 1.0E-6)
          ELSEIF (LA .EQ. LESIONAGEOPT) THEN
              TMP = 1.0
          ELSE 
              TMP = 1.0 - (LA-LESIONAGEOPT)/MAX(1.0-LESIONAGEOPT,1.0E-6)
          END IF
          LAF = MAX(MIN(TMP, 1.0), 0.0)
      END SUBROUTINE F_LAF

! ------ POTENTIAL SPORE RATE 
      SUBROUTINE F_PPSR(KVERHULST, RVERHULST, INFECTIOUS_S, PREV_IS,LAI,
     &                  PPSR)
          USE ModuleDefs
          IMPLICIT NONE
          REAL PPSR, KVERHULST, RVERHULST, INFECTIOUS_S, PREV_IS, LAI
          REAL CARRY
          IF (INFECTIOUS_S .EQ. 0.0) THEN
              PPSR = 0.0   
          ELSE
              CARRY = KVERHULST * LAI
              PPSR = (PREV_IS + (RVERHULST * PREV_IS *(CARRY-PREV_IS)) /
     &               MAX(KVERHULST, 1.0E-6)) / INFECTIOUS_S
              PPSR = MAX(PPSR, 0.0)
          END IF
      END SUBROUTINE F_PPSR

! ------ POTENTIAL SPORE RATE PER UNIT SPORULATING AREA 

      SUBROUTINE F_PPSR_POP(KVERHULST, RVERHULST, NPREV, INF_SURF_PREV,
     &                       LAI_SUSC, POT_SPO_PER_AREA)
          USE ModuleDefs
          IMPLICIT NONE
          REAL KVERHULST, RVERHULST, NPREV, INF_SURF_PREV, LAI_SUSC
          REAL POT_SPO_PER_AREA
          REAL KTOTAL, DN, NEW_SPORES
          REAL EPSL
          EPSL = 1.0E-6

          POT_SPO_PER_AREA = 0.0

          IF (LAI_SUSC .LE. 0.0) RETURN
          KTOTAL = MAX(KVERHULST * LAI_SUSC, EPSL)

          IF (NPREV .LE. 0.0) RETURN

          DN         = RVERHULST * NPREV * (KTOTAL - NPREV) / KTOTAL
          NEW_SPORES = MAX(DN, 0.0)

          IF (INF_SURF_PREV .GT. 0.0) THEN
             POT_SPO_PER_AREA = NEW_SPORES / INF_SURF_PREV
          ELSE
             POT_SPO_PER_AREA = 0.0
          END IF
      END SUBROUTINE F_PPSR_POP
      
! ------ DVIP (Daily infection-probability index) 
! Implementation of Beruski et al. (2020) Plant Disease.

      SUBROUTINE CALC_DVIP(LWD, T, DVIP)
          USE ModuleDefs
          IMPLICIT NONE
          
          REAL, INTENT(IN)     :: LWD, T
          INTEGER, INTENT(OUT) :: DVIP
          
          ! Local variables
          REAL NLes, TERM_LWD, TERM_TEMP, EXP_VAL 
          
          REAL, PARAMETER :: A_COEFF = 12.611
          
          REAL, PARAMETER :: LWD_OPT = 20.0
          REAL, PARAMETER :: LWD_SIG = 9.0
          
          REAL, PARAMETER :: T_OPT   = 23.0
          REAL, PARAMETER :: T_SIG   = 5.0
          
          DVIP = 0
          NLes = 0.0
          
          IF (LWD .LE. 1.0) RETURN

          ! --- Gaussian Model Calculation ---
          TERM_LWD  = ((LWD - LWD_OPT) / LWD_SIG)**2.0
          TERM_TEMP = ((T   - T_OPT)   / T_SIG)**2.0
          
          EXP_VAL = -2.5 * (TERM_LWD + TERM_TEMP)
          
          IF (EXP_VAL .LT. -20.0) THEN
             NLes = 0.0
          ELSE
             NLes = A_COEFF * EXP(EXP_VAL)
          END IF

          ! --- Classification ---
          IF (NLes .LE. 0.5) THEN
             DVIP = 0
          ELSEIF (NLes .LE. 3.0) THEN
             DVIP = 1
          ELSEIF (NLes .LE. 6.0) THEN
             DVIP = 2
          ELSE
             DVIP = 3
          END IF

      END SUBROUTINE CALC_DVIP

! ------ FUNGICIDE DECISION/APPLICATIONS
      SUBROUTINE APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7, 
     &             BufferDays, FungActive, ResidualDays, NSprays,
     &             USE_FUNGICIDE, HEALTH_LAI)
          USE ModuleDefs
          IMPLICIT NONE
          INTEGER DVIP_today, DVIP_pts(7), idx, SUM7
          INTEGER BufferDays, ResidualDays, NSprays
          LOGICAL FungActive, USE_FUNGICIDE
          REAL HEALTH_LAI
          LOGICAL CAN_SPRAY
          
          ! 1. Update Rolling Sum
          SUM7 = SUM7 - DVIP_pts(idx) + DVIP_today
          DVIP_pts(idx) = DVIP_today
          idx = MOD(idx,7) + 1
          
          ! 2. Decrement Buffer
          IF (BufferDays .GT. 0) BufferDays = BufferDays - 1
          
          ! 3. Decrement Residual (Protection)
          IF (FungActive) THEN
              ResidualDays = ResidualDays - 1
              IF (ResidualDays .LE. 0) THEN 
                  FungActive = .FALSE.
                  ! SUM7 = 0  
              END IF
          END IF
          
          ! 4. Spray Decision
          CAN_SPRAY = (HEALTH_LAI .GE. 0.1)
          
          IF (SUM7 .GE. 6 .AND. BufferDays .EQ. 0 .AND. CAN_SPRAY) THEN
              NSprays = NSprays + 1
              
              BufferDays = 16
              
              IF (USE_FUNGICIDE) THEN
                  FungActive   = .TRUE.
                  ResidualDays = 14
              ELSE
                  FungActive   = .FALSE.
                  ResidualDays = 0 
              END IF
          END IF
      END SUBROUTINE APPLY_FUNGICIDE
      
! ------ VIRTUAL LESIONS (Primiano & Amorim, 2020)
      SUBROUTINE F_VIRTUAL_LESIONS(s, beta, fvl)
          USE ModuleDefs
          IMPLICIT NONE
          
          REAL s      ! fractionary severity (0-1)
          REAL beta   ! β parameter for virtual lesions 
          REAL fvl    ! multiplier factor 
          
          fvl = (1.0 - MAX(0.0, MIN(s,1.0)))**beta
          fvl = MAX(0.0, MIN(fvl, 1.0))
          
      END SUBROUTINE F_VIRTUAL_LESIONS

! ------ SEVERITY
      SUBROUTINE F_SEVERITY(LAI_CUR, LAI_DIS_ACCUM, LAI_PEAK, SEV_MONO)
          USE ModuleDefs
          IMPLICIT NONE
 
          REAL, INTENT(IN)    :: LAI_CUR        
          REAL, INTENT(IN)    :: LAI_DIS_ACCUM  
          
          REAL, INTENT(INOUT) :: LAI_PEAK       
          REAL, INTENT(INOUT) :: SEV_MONO
          
          REAL :: RAW_SEV, PEAK_SAFE

          LAI_PEAK = MAX(LAI_PEAK, LAI_CUR)
          PEAK_SAFE = MAX(LAI_PEAK, 1.0E-6)
          RAW_SEV = (LAI_DIS_ACCUM / PEAK_SAFE) * 100.0
          SEV_MONO = MAX(SEV_MONO, RAW_SEV)
    
          SEV_MONO = MIN(SEV_MONO, 100.0)

      END SUBROUTINE F_SEVERITY
      
! ------ DEFOLIATION
      SUBROUTINE F_DEFOLIATION(WTLF, SLDOT, SEVERITY_PCT, rrds,
     &                         DISEASE_SEN_RATE)
          USE ModuleDefs
          IMPLICIT NONE
          
          REAL, INTENT(IN)    :: WTLF, SLDOT, SEVERITY_PCT, rrds
          REAL, INTENT(OUT)   :: DISEASE_SEN_RATE
          
          REAL rrsen, rrsenD
          REAL, PARAMETER :: EPS = 1.0E-6
          
          IF (WTLF .GT. EPS) THEN
              ! relative rate of senescence calculated by dssat
              rrsen = SLDOT / WTLF
              
              ! relative rate of senescence due to disease
              rrsenD = rrds * (SEVERITY_PCT / 100.0)
              
              !physical mass to be removed (avoiding double counting in the same area)]
              DISEASE_SEN_RATE = (rrsenD - (rrsen * rrsenD)) * WTLF
              DISEASE_SEN_RATE = MAX(0.0, DISEASE_SEN_RATE)
          ELSE
              DISEASE_SEN_RATE = 0.0
          END IF
          
          END SUBROUTINE F_DEFOLIATION
              
!=======================================================================

!***********************************************************************
!  Variable listing 
!***********************************************************************
! --------------------------- Arguments --------------------------------
! DYNAMIC           : DSSAT phase flag (RUNINIT/RATE/OUTPUT/SEASEND)
! CONTROL           : ControlType with %DAS (days after sowing), %RUN, etc.
! ISWITCH           : SwitchType (not used here; kept for interface parity)
! Tmin, Tmax (°C)   : Daily min/max air temperature
! RH (%)            : Relative humidity (either from weather or computed)
! LAI_TOTAL (m2 m-2): Canopy leaf area index (total)
! ESP_LAT_HIST(:,:) : State array (MAXDAYS x 5) for cohorts; columns:
!                     (1)=latent amount (spores m-2),
!                     (2)=latent progress (0..1),
!                     (3)=infectious flag (0/1),
!                     (4)=lesion age (d),
!                     (5)=secondary inoculum credited today (spores m-2)
! SUP_INF_LIST(:)   : Infectious surface by day (m2 m-2) → LA_INFECT
! LAI_INF_LIST(:)   : Cumulative removed LAI by day (m2 m-2) → LA_DISEASE
! YRDOY             : Current date (YYDOY)
! YREMRG            : Date of emergence (YYDOY) [not used internally]
! NVEG0             : DAS threshold for emergence gate
! YREND             : Harvest/end date (YYDOY)
! DISEASE_LAI       : OUTPUT — cumulative diseased area (cm2 m-2; printed as m2 m-2)

! --------------------------- Parameters -------------------------------
! USE_FUNGICIDE     : Toggle for applying fungicide (logic retained)
! USE_WTH_RH        : Use RH from weather (.TRUE.) or compute via dew point
! LAI_MIN_START     : Minimal LAI to allow activity (substrate safeguard)
! MAXDAYS           : Max internal history length (days)
! EPS               : Small epsilon for safe divisions

! --------------------------- Locals (scalars) -------------------------
! DAS               : Days after sowing (from CONTROL)
! T (°C)            : Daily mean air temperature
! FT                : Combined temperature response (0..1)
! FT_D, FT_G        : Post-infection and germination temperature responses
! IR (0..1)         : Infection rate for the day
! FSS (0..1)        : Canopy fraction supporting spores (capacity fraction)
! DS_TOTAL (spores m-2): Total primary + secondary spores deposited today
! DS_PRI, DS_SEC (spores m-2): Deposited primary and secondary portions
! LR (d-1)          : Latency progress rate
! LS_TODAY (spores m-2): New infections added after cohort progression
! IS (m2 m-2)       : Infectious surface today
! LA (d)            : Lesion age of a cohort
! LAF (0..1)        : Lesion-age factor
! PREV_IS (m2 m-2)  : Infectious surface yesterday
! LWD (h)           : Leaf wetness duration (capped at 24 h)
! Tdew (°C)         : Dew point temperature (for optional RH calc)
! Es, E             : Saturation and actual vapor pressure (Tetens; RH calc)
! HEALTH_LAI (m2 m-2)      : Healthy/susceptible LAI at day start
! HEALTH_LAI_AVAIL (m2 m-2): Remaining healthy LAI available for new lesions
! NEW_LOSS_TODAY (m2 m-2)  : Newly activated (removed) LAI today
! INF_AREA_K (m2 m-2)      : Infectious area allocated to cohort k today
! Lesion_Rate (d-1)        : Lesion aging rate
! ESP_INOC_SEC (spores m-2): Total secondary inoculum produced today
! SEC_SPORES_PENDING (spores m-2): Today's emission, available tomorrow
! SEC_SPORE_CLOUD (spores m-2): Airborne secondary inoculum after decay
! INF_COUNT_PREV (#)       : Count of infectious lesions yesterday
! NPREV_POP (#)            : Infectious lesion population used for sporulation
! POT_SPO_PER_AREA (sp m-2 d-1 per m2 m-2):
!                     Potential secondary spores per unit sporulating area
! PS_K (spores m-2)  : Secondary spores contributed by cohort k today
! NDS (spores m-2)   : Primary dispersed spores (from parameter file)
! LESION_S (m2)      : Average lesion surface area
! KVERHULST (#/LAI)  : Carrying capacity per unit LAI (population space)
! RVERHULST (d-1)    : Intrinsic logistic rate (population space)
! YMAX, COF_A, COF_B : Infection response parameters
! TMIN_G/TOT_G/TMAX_G (°C): Temperature response for germination
! TMIN_D/TOT_D/TMAX_D (°C): Temperature response for development
! LDMin (d)          : Minimum latency (used to scale LR)
! LESLIFEMAX (d)     : Max lesion lifespan (caps aging)
! DAE                : Days after emergence (internal counter)
! DAE_START          : DAE at which disease becomes active
! KMAX               : Upper bound for cohort loop (≤ MIN(DAE,MAXDAYS))
! DAE_IDX, PREV_IDX  : Indices for today and yesterday in history arrays
! PLANT_LIVE         : 1 while crop is alive (after emergence gate)
! DISEASE_LIVE       : 1 after disease activation gate (DAE ≥ DAE_START)
! DVIP_today         : Daily infection-probability class (0..3)
! idx (1..7)         : Circular index for 7-day DVIP buffer
! SUM7               : Sum of last seven DVIP classes
! BufferDays (d)     : Spray buffer to avoid back-to-back applications
! ResidualDays (d)   : Fungicide residual protection counter
! NSprays (#)        : Number of sprays applied
! FungActive (L)     : Whether fungicide residual is currently active
! FUNG_EFFICIENCY    : Proportional reduction in IR while active
! s    : Daily fractional severity (0-1). Computed as NEW_LOSS_TODAY / HEALTH_LAI. Represents today's proportion of newly necrosed leaf area.
! beta : Empirical parameter read from DISEASE_PARAMETERS.TXT. Controls the intensity of physiological reduction in green tissue.
! fvl  : Virtual lesion reduction factor (0-1). Computed by F_VIRTUAL_LESIONS as (1 - s)**beta. Reduces the effective daily C-assimilation potential.
! VIRTUAL_PHOTO_FACTOR  : Exported variable to CROPGRO (0-1). Receives fvl at RATE stage and is used inside CROPGRO to reduce PGAVL (PGAVL = PGAVL * VIRTUAL_PHOTO_FACTOR).
! LAI_PEAK_SEASON (m2 m-2): State variable (SAVE). Tracks the maximum value of 
!                           LAI_TOTAL observed throughout the current season. 
!                           Used as the denominator to normalize severity, preventing
!                           false 100% values during senescence.
! rrds (d-1) : Relative rate of senescence due to disease, read from input. Used in F_DEFOLIATION to compute the daily mass of leaf area to be removed due to disease.
!
! SEVERITY_PCT (%): OUTPUT variable (SAVE). Monotonous disease severity.
!                   Calculated as (LA_DISEASE / LAI_PEAK_SEASON) * 100.
!                   Contains the "blinded" value that never decreases,
!                   ensuring consistent AUDPC calculation.
! NCYCLE (CHARACTER*1) : Disease cycle type read from disease_parameters.txt.
!                       'P' = polycyclic (secondary inoculum produced; default).
!                       'M' = monocyclic (no secondary inoculum; single cycle).
! IS_MONOCYCLIC (L)   : Logical derived from NCYCLE. When .TRUE., the three
!                       secondary-inoculum blocks (population build-up,
!                       cohort sporulation, ESP_LAT_HIST re-injection) are
!                       skipped entirely. All other model logic (latency,
!                       lesion aging, DISEASE_LAI, VIRTUAL_PHOTO_FACTOR,
!                       DISEASE_SEN_RATE) runs unchanged for both modes.
      
! --------------------------- Locals (arrays) --------------------------
! ESP_LAT_HIST(MAXDAYS,5):
!   (1) latent amount (spores m-2) created on day k
!   (2) latent progress (0..1) accumulated for cohort k
!   (3) infectious flag for cohort k (0 or 1)
!   (4) lesion age (d) for cohort k (advances only after infectious)
!   (5) secondary inoculum tallied for cohort k on its creation day
!
! SUP_INF_LIST(MAXDAYS):
!   Infectious surface per day (m2 m-2); yesterday’s value used to derive
!   infectious counts and as denominator for population-based logistic.
!
! LAI_INF_LIST(MAXDAYS):
!   Cumulative removed/sick LAI per day (m2 m-2).
!
! ADMITTED_AREA(MAXDAYS):
!   Infectious area actually admitted to cohort k (m2 m-2), capped by
!   remaining healthy LAI at activation time.

! --------------------------- Subroutines ------------------------------
! F_TAVG      : Mean temperature
! F_DEW       : Empirical dew point (for optional RH path)
! F_RH        : RH from Tetens (uses T and Tdew)
! F_LWD       : Leaf wetness duration from RH (0..24 h)
! T_DEV       : Beta-type temperature responses (FT_G, FT_D, FT)
! F_IR        : Infection rate vs FT and LWD (Ymax, A, B)
! F_CANSPO    : Canopy spore support fraction (capacity)
! F_DS        : Deposited spores = FSS * NDS
! F_LR        : Latency progress rate (FT_D / LDMin)
! F_LS        : Latent spores = IR * DS (≥ 0)
! F_LAF       : Lesion-age factor (triangular shape peaking at LESIONAGEOPT)
! F_LAR       : Lesion aging rate (FT_D / LESLIFEMAX; 0 if FT_D≤0)
! F_PPSR      : Legacy area-based logistic (kept for compatibility)
! F_PS        : Legacy secondary inoculum per lesion (compatibility)
! F_PPSR_POP  : Population-space logistic → spores per unit sporulating area
! CALC_DVIP   : Daily infection-probability index class (0..3)
! APPLY_FUNGICIDE:
!   7-day risk sum trigger; optional spray; residual IR reduction window.
! F_VIRTUAL_LESIONS: Simulates the green leaf area around the necrotic area that does less photosynthesis
! F_DEFOLIATION: Computes the daily mass of leaf area to be removed due to disease, avoiding double counting with natural senescence.
