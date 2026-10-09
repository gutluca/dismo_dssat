# Disease Impact and Severity Module - DISMO

**Authors**: Gustavo de A. Luca, Izael M. Fattori Jr., Willingthon Pavan, Fábio R. Marin

**Affiliation**: University of São Paulo – ESALQ/USP

---

## Abstract

This repository contains a modified version of the [DSSAT]([https://dssat.net/](https://github.com/DSSAT/dssat-csm-os.git)) (Decision Support System for Agrotechnology Transfer) source code, integrating a generic epidemiological model for foliar fungal diseases into the Generic-Pest routine, via damage-coupling points.

The new module, **DISMO** - Disease Impact and Severity Module (`DISEASE_LEAF` subroutine), simulates the progression of leaf diseases such as Asian Soybean Rust (*Phakopsora pachyrhizi*) using a daily cohort-based state-transition model. It calculates diseased leaf area from environmental drivers and pathogen dynamics, and reduces canopy photosynthesis via the existing Generic-Pest (`PEST`) module.

The module is activated by a single switch in the DSSAT experiment file (FILEX): **`ISWDIS = 'Y'`**.
DISMO will only activate if  ISWDIS = 'Y' AND the input file (disease_parameters.txt) is in the current simulation directory
(for example: if you are simulating soybean, the input parameters should be in 'C:\DSSAT48\Soybean').

> **Architecture note (March 2026):** `DISMO.for` was moved from `Plant\CROPGRO` to `Plant\Generic-Pest` and is now invoked exclusively through the **Generic-Pest (`PEST`) module**. No direct calls to `DISEASE_LEAF` remain in `CROPGRO.for`.

---

## Scientific Background

### 1. The Epidemiological Model (`DISEASE_LEAF` — `DISMO.for`)

The model follows a **cohort-based polycyclic approach**. Each day a new cohort of latent spores is deposited on the canopy; cohorts progress through latency, become infectious, sporulate, and die according to temperature-driven rates. Secondary inoculum produced by infectious cohorts creates new latent cohorts, resulting in the exponential epidemic growth typical of rust pathogens.
If desired, DISMO also supports monocyclic disease simulations. The cycle mode is defined by the NCYCLE column in the disease parameter file: use 'P' for Polycyclic (default) or 'M' for Monocyclic.

#### Disease Cycle — Processes and Functions

| Step | Function | Formula / Logic |
|---|---|---|
| Canopy spore capacity | `F_CANSPO` | `FSS = min(LAI / (NDS × LESION_S), 1)` |
| Spore deposition | `F_DS` | `DS = FSS × NDS` |
| Infection rate | `F_IR` | `IR = Ymax × FT × (1 − exp(−(A × LWD)^B))` |
| Latent spores per cohort | `F_LS` | `LS = IR × DS` |
| Latency progress rate | `F_LR` | `LR = FT_D / LDmin` |
| Latent → infectious | cohort loop | when cumulative `LR >= 1`, cohort becomes infectious |
| Lesion aging rate | `F_LAR` | `Lesion_Rate = FT_D / LESLIFEMAX` |
| Sporulation shape | `F_LAF` | Triangular: rises to 1 at `LESIONAGEOPT`, falls to 0 at age = 1 |
| Secondary inoculum | `F_PPSR_POP` | Verhulst logistic in population space; spores per unit infectious area |
| Virtual lesion factor | `F_VIRTUAL_LESIONS` | `fvl = (1 − s)^β`, where `s = NEW_LOSS_TODAY / HEALTH_LAI` |
| Disease induces senescence | `F_DEFOLIATION` | `DISEASE_SEN_RATE = (rrsenD − (rrsen × rrsenD)) × WTLF`, where `rrsen = SLDOT / WTLF` and `rrsenD = rrds × (SEVERITY_PCT / 100)` |

#### Temperature Response

Two independent **beta-shaped** responses are computed from daily mean temperature `T = (Tmax + Tmin) / 2`:

- **`FT_G`** — spore germination, parameterised by `TMIN_G`, `TOT_G`, `TMAX_G`
- **`FT_D`** — general fungal development (latency, aging, sporulation), parameterised by `TMIN_D`, `TOT_D`, `TMAX_D`
- **`FT = FT_G × FT_D`** — combined response used for the infection rate

#### Leaf Wetness Duration

`LWD` (hours, 0–24) is derived from relative humidity via a logistic function:

```
LWD = 31.31 / (1 + exp(-((RH − 85.17) / 9.13)))
```

When `USE_WTH_RH = .TRUE.` (default), `RH` comes directly from the DSSAT weather file. When `.FALSE.`, `RH` is computed internally from an empirical dew-point formula using `Tmin` and `Tmax`.

#### Healthy LAI Tracking

At each `RATE` call, healthy LAI is computed as:

```
HEALTH_LAI = max(LAI_TOTAL − LAI_INF_LIST(DAE−1), 0)
```

where `LAI_INF_LIST` is the cumulative removed LAI array updated daily. New infectious cohorts can only occupy remaining healthy area (`HEALTH_LAI_AVAIL`), preventing the model from exceeding the physical leaf area.

#### Virtual Lesion Effect on Photosynthesis

The virtual lesion concept (Primiano & Amorim, 2020) accounts for the physiological impairment of still-green tissue surrounding necrotic lesions:

1. `s = NEW_LOSS_TODAY / HEALTH_LAI` — daily fraction of newly necrosed leaf area
2. `fvl = (1 − s)^β` — `β` is read from the parameter file; higher values = stronger reduction
3. `VPHOTF_DISMO = fvl` is stored via `PUT('PLANT','VPHOTF',...)` and applied to `PGAVL` inside `ASMDM` in `PEST.for`

#### Disease Induces Senescence (Defoliation)

The accelerated leaf senescence logic (Willocquet et al., 2025) accounts for early defoliation by integrating natural crop senescence with disease severity, preventing overlapping leaf death:
---

### 2. Fungicide Decision Module

An optional rule-based fungicide model, toggled by `USE_FUNGICIDE` in `DISMO.for` (default `.FALSE.`):

| Component | Value | Description |
|---|---|---|
| Risk index (DVIP) | 0–3 per day | Gaussian function of `LWD` and `T` (Beruski et al., 2020) |
| Spray trigger | `SUM7 >= 6` AND `HEALTH_LAI >= 0.1` | 7-day rolling DVIP sum exceeds threshold with sufficient canopy |
| Buffer period | 16 days | Minimum interval between consecutive sprays |
| Residual protection | 14 days | Efficacy window after each application |
| Efficacy | `FUNG_EFFICIENCY = 0.723` | Proportional reduction in `IR` during residual period |

When `USE_FUNGICIDE = .FALSE.`, spray decisions are still tracked (for diagnostics in the output file), but the infection rate is not reduced.

---

### 3. Integration Architecture

#### Previous architecture (before March 2026)

`DISEASE_LEAF` was called **directly inside `CROPGRO.for`** at every simulation phase (`RUNINIT`, `RATE`, `INTEGR`, `OUTPUT`, `SEASEND`). State arrays (`ESP_LAT_HIST`, `SUP_INF_LIST`, `LAI_INF_LIST`) and weather variables were declared locally in `CROPGRO.for`. The `VIRTUAL_PHOTO_FACTOR` was also applied there directly as `PGAVL = PGAVL * VIRTUAL_PHOTO_FACTOR`.

#### Current architecture

`DISEASE_LEAF` is called **exclusively inside `PEST.for`**. `CROPGRO.for` calls `PEST` unchanged — only minimal additions were made to CROPGRO (`PUT` calls for `XLAID`, `NVEG0D`, `YREMGD`, `YRENDD`, and the `ISWDIS='Y'` condition in the `PEST` call guard).

```
CROPGRO.for
  └── CALL PEST(...)                           ← ISWDIS = 'Y'
        └── PEST.for
              ├── [RUNINIT]   IF RUN_DISMO → CALL DISEASE_LEAF(RUNINIT)
              ├── [SEASINIT]  IF ISWDIS='Y' → CALL IPPROG (classic time-series)
              │               CALL PESTCP / ASMDM / SEEDDM / VEGDM / ROOTDM / OPPEST
              │               IF RUN_DISMO → CALL DISEASE_LEAF(SEASINIT)
              │               PUT('PLANT','VPHOTF', VPHOTF_DISMO)
              ├── [RATE]      IF RUN_DISMO → GET weather & LAI via GET/PUT
              │                           → CALL DISEASE_LEAF(RATE)
              │                           → DISLA += DISEASE_LAI_DISMO
              │                           → PUT('PLANT','VPHOTF', VPHOTF_DISMO)
              │               IF PCN > 0  → classic LINDM/PESTCP/SEEDDM/VEGDM/ROOTDM
              ├── [INTEGR]    CALL ASMDM(PGAVL, PPSR, TPSR, VPHOTF_DISMO → ASMDOT)
              │               IF PCN > 0  → SEEDDM / VEGDM / ROOTDM
              └── [OUTPUT/SEASEND]  CALL OPPEST
                              IF RUN_DISMO → CALL DISEASE_LEAF(OUTPUT/SEASEND)
```

**Key design points:**

- **Classic pest routines** (`PESTCP`, `SEEDDM`, `VEGDM`, `ROOTDM`) run at `RATE` only when `PCN > 0`. With `ISWDIS = 'D'` and no FILET loaded, `PCN = 0` and those routines are skipped.
- **`ASMDM`** always runs at `INTEGR`, regardless of `ISWDIS`. It now accepts `VPHOTF_DISMO` as a fourth input and applies the virtual lesion reduction to `PGAVL` before subtracting classic pest assimilate damage.
- **`DISLA`** (cm² m⁻²) is zeroed

#### Data Exchange Between Modules

| Variable | Producer | Consumer | Mechanism | Description |
|---|---|---|---|---|
| `TMIN`, `TMAX`, `RHUM` | CROPGRO weather | PEST → DISMO | `GET(WEATHER_PEST)` inside PEST | Daily meteorological inputs |
| `XLAI` | CROPGRO | PEST → DISMO | `PUT('PLANT','XLAID',XLAI)` at end of CROPGRO | Total canopy LAI |
| `NVEG0` | CROPGRO/PHENOL | PEST → DISMO | `PUT('PLANT','NVEG0D',NVEG0)` | Emergence DAS gate in DISMO |
| `YREMRG` | CROPGRO/PHENOL | PEST → DISMO | `PUT('PLANT','YREMGD',YREMRG)` | Emergence date |
| `YREND` | CROPGRO | PEST → DISMO | `PUT('PLANT','YRENDD',YREND)` | Harvest date; resets `PLANT_LIVE` |
| `DISEASE_LAI_DISMO` | DISMO | PEST | local variable in PEST | Added to `DISLA`; passed to `VEGDM` and `GROW` |
| `VPHOTF_DISMO` | DISMO | PEST → ASMDM | `PUT('PLANT','VPHOTF',...)` | Reduces `PGAVL` via `ASMDM` at `INTEGR` |
| `DISLA` | PEST | CROPGRO | subroutine output argument | Diseased leaf area used by `GROW` |
| `ASMDOT` | PEST/ASMDM | CROPGRO | subroutine output argument | Assimilate damage subtracted from `PGAVL` |

---

## How to Use

### 1. Compilation

1. Clone this repository.
2. In your DSSAT-CSM Visual Studio project, add `DISMO.for` as a source file under the `Plant\Generic-Pest` filter.
   - Register it in `ALL_BUILD.vcxproj` and `ALL_BUILD.vcxproj.filters` if not already present.
3. Replace `PEST.for` with the modified version from this repository (in `Plant\Generic-Pest\`).
4. The modified `CROPGRO.for` is already in place — it requires no additional action beyond a clean rebuild.
5. Recompile the full DSSAT-CSM solution.

### 2. Activating DISMO: `ISWDIS = 'Y'` in the FILEX AND input file present in current dirrectory.

The module is activated by the **`DISES`** field in the `*SIMULATION CONTROLS` section of the DSSAT experiment file (`.SBX`, `.PNX`, or equivalent FILEX):

```
*SIMULATION CONTROLS
...
@N WATER NITRO SYMBI PHOS POTAS DISES TILL  ...
 1 Y     Y     Y     N    N     Y     N     ...
```

| `DISES` | Behavior |
|---|---|
| `N` | No pest or disease simulation |
| `Y` without the 'disease_parameters.txt'| Classic DSSAT pest module — requires a FILET pest time-series file |

When `ISWDIS = 'Y'` and 'disease_parameters.txt' inside the current simulation directory:
- `CROPGRO.for` calls `PEST` (same entry point as for `'Y'`)
- Inside `PEST.for`, `RUN_DISMO = .TRUE.` activates `DISEASE_LEAF` at every simulation phase
- `IPPROG` is **not** called, so `PCN = 0` and classic coupling-point routines (`PESTCP`, `SEEDDM`, `VEGDM`, `ROOTDM`) are skipped at `RATE` and `INTEGR`
- Weather and LAI data are automatically retrieved inside `PEST` via the DSSAT `GET`/`PUT` exchange mechanism

### 3. Input Files

#### `disease_parameters.txt`

This file contains all pathogen-specific parameters. It is searched in this order:

1. **Current working directory** — the folder from which DSSAT is executed (same location as the `.SBX` experiment file and DSSAT output files)
2. **DSSAT species data path** — resolved at runtime via `PATH('PSD', CONTROL%DSSATP, ...)`, typically `C:\DSSAT48\Soybean\`

If the file is not found in either location, `READ_DISEASE_PARAMETERS` exits silently and parameters retain their default (zero) values, which will suppress all disease dynamics.

**File format:** 5 header/comment lines followed by exactly **one data line** with the parameters in the order below, space-separated:

| # | Parameter | Type | Units | Description |
|---|---|---|---|---|
| 1 | `ID` | CHAR(6) | — | Pathogen identifier code |
| 2 | `VRNAME` | CHAR(20) | — | Pathogen or variant name |
| 3 | `NDS` | REAL | spores m⁻² | Primary dispersed spores per infection cycle |
| 4 | `LESION_S` | REAL | m² | Average area of a single lesion |
| 5 | `KVERHULST` | REAL | — | Carrying capacity per unit LAI (population space) |
| 6 | `RVERHULST` | REAL | d⁻¹ | Logistic intrinsic growth rate |
| 7 | `YMAX` | REAL | 0–1 | Maximum infection rate under optimal conditions |
| 8 | `COF_A` | REAL | — | Infection response coefficient A (LWD sensitivity) |
| 9 | `COF_B` | REAL | — | Infection response coefficient B (LWD shape) |
| 10 | `TMIN_G` | REAL | °C | Minimum temperature for spore germination |
| 11 | `TOT_G` | REAL | °C | Optimum temperature for spore germination |
| 12 | `TMAX_G` | REAL | °C | Maximum temperature for spore germination |
| 13 | `TMIN_D` | REAL | °C | Minimum temperature for fungal development |
| 14 | `TOT_D` | REAL | °C | Optimum temperature for fungal development |
| 15 | `TMAX_D` | REAL | °C | Maximum temperature for fungal development |
| 16 | `LDmin` | REAL | d | Minimum latency period; scales `LR = FT_D / LDmin` |
| 17 | `LESIONAGEOPT` | REAL | 0–1 | Normalised lesion age at peak sporulation |
| 18 | `LESLIFEMAX` | REAL | d | Maximum lesion lifespan; scales aging rate |
| 19 | `DAE_START` | INT | d | Days after emergence when disease becomes active |
| 20 | `beta` | REAL | — | Virtual lesion exponent β; higher = stronger photosynthesis reduction |
| 21 | rrds | REAL | d⁻¹ | Relative rate of disease-induced senescence per disease severity unit |

#### Internal Switches in `DISMO.for`

Two compile-time `PARAMETER` constants control optional behaviour:

| Constant | Default | Effect when `.TRUE.` |
|---|---|---|
| `USE_FUNGICIDE` | `.FALSE.` | Applies fungicide when `SUM7 >= 6`; reduces `IR` by `FUNG_EFFICIENCY` for 14 days |
| `USE_WTH_RH` | `.TRUE.` | Uses `RH` from the DSSAT weather file directly; when `.FALSE.` computes `RH` from empirical dew-point |

### 4. Output Files

#### `DISMO.OUT`

- **Location:** written to the **current working directory** (same folder as DSSAT output files such as `PLANTGRO.OUT`)
- **When:** file is opened (with `STATUS='replace'`) at `CONTROL%RUN = 1` during `RUNINIT`; one row is written per day at every `OUTPUT` call; the file is **closed** at `SEASEND`
- **Multi-season runs:** a new file is created at the start of each run sequence (file is replaced at `RUN = 1`); subsequent seasons within the same run append rows

| Column | Units | Description | Notes |
|---|---|---|---|
| `FILEX` | — | Experiment file name (8 chars) | From `CONTROL%FILEX` |
| `RUN` | — | Simulation run number | From `CONTROL%RUN` |
| `YYDOY` | YYDDD | Current date | e.g. `25031` = year 2025, day 31 |
| `DAS` | d | Days after sowing | From `CONTROL%DAS` |
| `LAI_HEALTH` | m² m⁻² | Healthy LAI at start of day | `LAI_TOTAL − LAI_INF_LIST(prev)`, bounded ≥ 0 |
| `LA_DISEASE` | m² m⁻² | Cumulative diseased/removed LAI | Internally stored as cm² m⁻²; divided by 10 000 for output |
| `LA_INFECT` | m² m⁻² | Instantaneous infectious (sporulating) surface today | Daily total; **not** cumulative; sum across all active cohorts |
| `NEW_LOSS` | m² m⁻² d⁻¹ | Newly activated diseased area today | Summing over all days reconstructs `LA_DISEASE` |
| `LWD` | h | Leaf wetness duration | 0–24 h; derived from `RH` via logistic function |
| `RH` | % | Relative humidity | From weather file (`USE_WTH_RH=.TRUE.`) or computed |
| `FAT_TEMP` | 0–1 | Combined temperature response `FT = FT_G × FT_D` | Used for infection rate |
| `LAI_TOTAL` | m² m⁻² | Total canopy LAI from CROPGRO | Retrieved via `GET('PLANT','XLAID',...)` |
| `SUM7` | 0–21 | 7-day rolling sum of DVIP risk classes | Each day contributes 0–3 |
| `NSPRAYS` | count | Cumulative fungicide applications | Zero if `USE_FUNGICIDE = .FALSE.` |
| `FUNG_ACT` | T / F | Fungicide residual protection active | When `T`, `IR` is reduced by `FUNG_EFFICIENCY` |
| `Severity%` | % | Monotone disease severity | `(LA_DISEASE / LAI_PEAK_SEASON) × 100`; never decreases within a season |

> **Unit note on `DISEASE_LAI`:** the variable is passed from DISMO to PEST as `DISEASE_LAI_DISMO` in **cm² m⁻²**, matching the units of `DISLA` used by `VEGDM` and `GROW`. The output file column `LA_DISEASE` divides by 10 000 to report in m² m⁻².

---

## Citation

If you use this code or data in your research, please cite this repository as follows:

Luca, G. de A., Fattori Junior, I. M., Del Ponte, E. M., & Marin, F. R. (2025). *Code and data for: Process-Based Simulation of Soybean Rust in Brazil: A DSSAT-Coupled Approach* (Version v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.17266278

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17266278.svg)](https://doi.org/10.5281/zenodo.17266278)

---
- **DISMO calibration and validation**:
  - Luca, G. de A., Fattori Junior, I. M., Del Ponte, E. M., & Marin, F. R. (2026). Process-based simulation of soybean rust in Brazil: A DSSAT-coupled approach. European Journal of Agronomy, 178, 128139.https://doi.org/10.1016/j.eja.2026.128139*

- **Virtual Lesions**:
  - Primiano, I. V., & Amorim, L. (2020). *Tropical Plant Pathology*
    
- **Defoliatoin logic**
  - Willocquet L, Bregaglio S, Ferrise R, Kim K, Savary S (2025) DYNAMO-A: a generic simulation model coupling crop growth and disease epidemic. PLoS One 20(4): e0321261. 
  https://doi.org/10.1371/journal.pone.0321261

- **Fungicide Decision Index (DVIP)**:
  - Beruski, N. D., et al. (2020). *Plant Disease*

- **Fungicide Recommendations**:
  - Godoy, C. V., et al. (2024). *Embrapa Technical Circular*

- **Epidemiological Model Basis**:
  - Caubel, J., et al. (2017). *European Journal of Agronomy*, 90, 53–66.

- **DSSAT Framework**:
  - Jones, J. W., et al. (2003). *The DSSAT Cropping System Model*. *European Journal of Agronomy*, 18(3–4), 235–265.
