# romsDiag
romsDiag is a plotting routine to call romsObj and create like-for-like ROMS and validation figures.
    
## Table of Contents

- [Updates](#updates)
- [Getting started](#getting-started)
- [Using romsDiag](#using-romsDiag)
- [Generating PDFs](#generating-pdfs)
- [Support](#support)

Requires MATLAB 2013 or above.

## Updates
* 07/17/2024 -- First commit of romsObj, romsDiag, romsComp, and romsOpt 
* 07/29/2024 -- Updated routines to make PDF generation automatic 

## Getting started
#### Update diag_options.m to add new BGC variables (if any)
    - As of 07/29/2024, there are 7 choices to pick from:
        - plots(1) = Shallow depth-slice comparisons (z-slices)
        - plots(2) = Transect comparisons
        - plots(3) = Surface field comparisons
        - plots(4) = Oxygen-minimum-zone (OMZ) thickness comparisons
        - plots(5) = Sinking particulate organic carbon (POC) flux comparisons
        - plots(6) = Gridded biogeochemical tracers (O2, NO2, NH4, N2O, NO3) vs ROMS
        - plots(7) = Deep depth-slice comparisons (similar to plots(1)) 

    NOTE: I've also provided the file 'diag_overrides.m' under the example simName directories (e.g. ./peru_chile_0p1/diag_overrides.m)
    You can update any of the settings in 'diag_options' here, and these will be applied if you're calling that 'simName' in 'run_diag.m'
    Useful for setting permanent overrides for specific grids (i.e., which lons/lats you wish to create transect comparisons)
    Useful for setting temporary overrides (i.e., turning off comparisons of certain variables via plots(#).opt)

## Using romsDiag
#### Update run_diag.m to specify simName, runName, and which ROMS files to use 
    - Routine only works (or, was designed) for:  
        - 1 annually-averaged ROMS output file
        - 12 monthly-averaged ROMS output files
    - If daily-output (365 or 366 netcdf files), use NCO-tools to create an annually-averaged file first

#### Call run_diag.m to generate diagnostic figures (may take a while) 
    - in MATLAB:
    >> run_diag;

## Generating PDFs 
    - Typically, I create the diagnostic PDFs on my local machine, using 'scp' to copy figures from remote-->local
    - Recommended instructions:
        - (1) Download contents of 'PDF' folder from GitHub to local machine 
        - (2) Navigate to folder where contents are held 
        - (3) Update 'get_figs_make_pdf.sh' to edit:
            - simName (same as in 'run_diag.m')
            - runName (same as in 'run_diag.m')
            - simTitle (LaTeX friendly name of simulation, e.g. Pacific-wide 25km)
            - runTitle (LaTeX friendly name of run, e.g. Test)
            - diagPath (location of romsObj figures on remote server)
            - scp command for copying files from remote --> local
        - (4) Update 'string_replacements.tex' if any variables were added to 'diag_options.m'
            - This sets the captions under each figure in the PDF
        - (5) Execute the shell script:
            - ./get_figs_make_pdfs.sh
            - NOTE: May require loosened permissions via 'chmod' (e.g., chmod 777 get_figs_make_pdfs.sh)
            - NOTE: Assumes local terminal has access to 'pdflatex' command
                - If not, download .Tex software and compile the PDF seperately after running ./get_figs_make_pdfs.sh 
        - (6) Rename output PDF ('diagnostics_template_auto.pdf') since it is always overwritten

## Support
Contact Daniel McCoy at the Carnegie Institution for Science (dmccoy@carnegiescience.edu) 

