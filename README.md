# CosmicTB Analysis Code

Analysis code for the CosmicTB testbeam experiment at Irfu/CEA, developed during the 2025 internship. The code covers data decoding, ray tracking, detector alignment, efficiency analysis and Fe55 gain measurements for a P2 detector equipped with DREAM ASICs.

For the full operational pipeline (data taking, directory structure, file naming, what to keep), see the companion document `README_pipeline.md`.

---

## Repository structure

```
avril_code_ctb/
│
├── CodeDamienNew/          # P2 data decoder
│   ├── DreamDataReader.cxx # Main decoder — reads .fdf files, outputs .root
│   ├── FeuReadOut.h        # Feu readout structures
│   ├── dreamdataline.h     # Data line definitions
│   ├── Pedestal.C          # Pedestal computation
│   ├── RMSPedestaux.C      # RMS pedestal calculation
│   ├── CommonNoisePedSubstr.C  # Common noise + pedestal subtraction
│   ├── Display.C           # Event display
│   ├── Makefile            # Build file
│   └── bash_files/         # Run scripts (RunAllP2.sh, RunAll.sh, ...)
│
├── m3_tracking/            # Ray tracking code (Feu 1)
│   ├── src/                # C++ source files
│   ├── include/            # Header files
│   ├── m3_bash_files/      # Run scripts (run_tracking_P2.sh, ...)
│   ├── config*.json        # Configuration files
│   └── Makefile            # Build file
│
├── TBanalysisP2_2025_mod.C # Main cosmic analysis macro (recommended)
├── TBanalysisP2_2025.C     # Main cosmic analysis macro (previous version)
├── TBanalysisP2_2025_v2.C  # Main cosmic analysis macro (v2)
├── TBanalysisP2_Fe55.C     # Fe55 analysis macro
├── TBanalysisP2_2024.C     # P2 analysis macro from 2024
├── TBanalysis.C            # Generic analysis macro
├── TBanalysisRD3.C         # RD3 detector analysis
├── TBanalysisTPOT_2022.C   # TPOT analysis (2022)
│
├── rotation.C              # Alignment macro — residuals, outputs outtree_aligned.root
├── efficacite.C            # Efficiency map (requires mapping + aligned data)
├── efficacite_no_mapping.C # Efficiency map without mapping
├── dead_ch.C               # Dead channel map (requires mapping + aligned data)
├── eff_vs_volt.C           # Efficiency vs HV scan
├── gain_2D.C               # 2D gain map for Fe55
├── connector_dead_ch.C     # Pedestal / dead channel check per connector
├── readRMS.C               # Read and display RMS pedestal file
├── producerkp2pdf.C        # PDF producer utility
├── SetStyle.C              # ROOT style settings — include in all macros
│
├── outTreeAnalysis/        # Additional output tree analysis macros
├── macro_odl/              # Older alignment and mapping macros
├── bash_files/             # Top-level run scripts (analysisP2.sh, RunAna.sh)
├── cfgfile/                # Configuration files for data taking
└── config_cos.json         # CosmicBench configuration
```

---

## Dependencies

- **ROOT** (CERN) — all `.C` macros and the decoder output rely on ROOT
- **C++ compiler** — g++ with C++11 support minimum (see `cpp11.patch` in m3_tracking)
- **make** — to build `CodeDamienNew` and `m3_tracking`

The code was developed and tested on CentOS 7 (`sedipcaa28.extra.cea.fr`). Compatibility with other systems has not been verified.

---

## Compilation

### CodeDamienNew (P2 decoder)

```bash
cd CodeDamienNew/
# edit FeuIDs.txt and DreamDataReader.cxx to set the number of Feus (Nfeu)
# and the Feu ID mapping:
#   index 0 → Feu 3 (ID 40)
#   index 1 → Feu 4 (ID 22)
#   index 2 → Feu 7 (ID 21)
make
```

This produces the `DreamDataReader` executable.

### m3_tracking (ray tracker)

```bash
cd m3_tracking/
make
```

This produces the `tracking` executable (and `AutoAlign`, `DataReader`, `wrapper`, `MultiCluster`, `carac_all`).

---

## Running the code

### P2 decoding

Edit `CodeDamienNew/bash_files/RunAllP2.sh` to set the input path, filenames, number of files and Feus, then:

```bash
# run inside a screen session — can take a long time
. CodeDamienNew/bash_files/RunAllP2.sh
```

### Ray tracking

Edit `m3_tracking/m3_bash_files/run_tracking_P2.sh` to set the input path and filenames, then:

```bash
. m3_tracking/m3_bash_files/run_tracking_P2.sh
```

### ROOT analysis macros

All macros are run from the ROOT interpreter. Edit the input file path inside the macro before running.

```bash
# Full cosmic analysis
root -l 'TBanalysisP2_2025_mod.C++("/path/to/work/folder/")'

# Fe55 analysis
root -l 'TBanalysisP2_Fe55.C++("/path/to/work/folder/")'

# Alignment (run before efficiency macros)
root -l rotation.C

# Efficiency map (requires aligned data + mapping)
root -l efficacite.C

# Dead channel map (requires aligned data + mapping)
root -l dead_ch.C

# Efficiency vs HV
root -l eff_vs_volt.C

# Fe55 gain map
root -l gain_2D.C

# Pedestal / dead channel check (input must be a .thr.aux file)
root -l
root [0] .L connector_dead_ch.C+
root [1] connector_dead_ch()
```

Alternatively, use the shell script for file-by-file analysis:

```bash
# edit bash_files/analysisP2.sh first
. bash_files/analysisP2.sh
```

---

## Notes for future developers

- `SetStyle.C` sets the ROOT style — load it at the start of any new macro with `.L SetStyle.C` for consistent plots.
- `macro_odl/` contains older versions of alignment and mapping macros, kept for reference.
- `outTreeAnalysis/` contains macros for further analysis of the output trees produced by the main analysis macros.
- When adding a new Feu, update both `FeuIDs.txt` and the `Nfeu` variable in `DreamDataReader.cxx`, then recompile with `make`.
- The `config*.json` files in `m3_tracking/` control the tracking configuration — `config_cos.json` is the reference configuration for CosmicTB.

---

## Author

Avril Gaurot — Irfu/CEA internship 2025-2026
Supervisor contact: [Maxence Vandebroucke - Alexandra Kallitspoulou]
```
