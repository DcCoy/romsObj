#!/bin/sh

# Simulation and run information
simName="peru_chile_0p1"    
runName="microbes_eth_obligate_tune0" # runName

# Alternative titles (for LaTeX document)
simTitle="Peru Chile 10km"
runTitle="Tune0" 

# Get diagPath for plots on remote server
# NOTE: UPDATE THIS PATH FOR YOUR OWN SERVER
diagPath=/data/project1/demccoy/romsObj/diag/${simName}/plots/${runName}/

# Get output directories
mkdir ${simName}
mkdir ${simName}/${runName}/

# Directory containing the files
DIR="${simName}/${runName}"

# Remove old figures
rm ${DIR}/*png
sleep 5

# Download new figures from remote server
# NOTE: UPDATE THIS COMMAND FOR YOUR OWN SERVER 
scp demccoy@poseidon:${diagPath}*png ${DIR}/.

######################################################################
# Start the LaTeX file and save number of variables and variable names
######################################################################
# Generate LaTeX file with variables
OUTPUT_FILE="variables.tex"
rm $OUTPUT_FILE
cd ${DIR}

####################################
# Get simulation and run information
####################################

# Get simName, runName, simTitle, runTitle 
echo "\\\\newcommand*{\\\\simName}{${simName}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"
echo "\\\\newcommand*{\\\\runName}{${runName}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"
echo "\\\\newcommand*{\\\\simTitle}{${simTitle}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"
echo "\\\\newcommand*{\\\\runTitle}{${runTitle}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

#######################################
# Get surface (2D) vars, headers, diags
#######################################
count=0
files=$(ls *_2D_roms.png)
prefixes=$(echo "$files" | sed 's/...//' | sed 's/_2D_roms.png//')
echo "\\\\newcommand*{\sfcVars}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    echo "{$prefix}" >> "$OUTPUT_FILE"
    count=$((count + 1))
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get headers
headers=$(echo "$files" | cut -d'_' -f1)
echo "\\\\newcommand*{\sfcHdr}{" >> "$OUTPUT_FILE"
for this_header in $headers; do
    echo "{$this_header}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numsfcVars
echo "\\\\newcommand*{\\\\numsfcVars}{${count}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get number of sfcdiags
echo "\\\\newcommand*{\sfcDiags}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    file_count=$(ls *_${prefix}_2D_diag*.png 2>/dev/null | wc -l)
    echo "{$file_count}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

############################################
# Get zshallow (zslice) vars, headers, diags
############################################
count=0
files=$(ls *_zshallow_1_roms.png)
prefixes=$(echo "$files" | sed 's/...//' | sed 's/_zshallow_1_roms.png//')
echo "\\\\newcommand*{\shallowVars}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    echo "{$prefix}" >> "$OUTPUT_FILE"
    count=$((count + 1))
    numshallowZ=$(ls *_${prefix}*_zshallow_*_roms.png 2>/dev/null | wc -l)
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numshallowZ
echo "\\\\newcommand*{\\\\numshallowZ}{${numshallowZ}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get headers
headers=$(echo "$files" | cut -d'_' -f1)
echo "\\\\newcommand*{\shallowHdr}{" >> "$OUTPUT_FILE"
for this_header in $headers; do
    echo "{$this_header}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numshallowVars
echo "\\\\newcommand*{\\\\numshallowVars}{${count}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get number of diags
echo "\\\\newcommand*{\shallowDiags}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    file_count=$(ls *_${prefix}_zshallow_1_diag*.png 2>/dev/null | wc -l)
    echo "{$file_count}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

###############################
# Get deep vars, headers, diags
###############################
count=0
files=$(ls *_zdeep_1_roms.png)
prefixes=$(echo "$files" | sed 's/...//' | sed 's/_zdeep_1_roms.png//')
echo "\\\\newcommand*{\deepVars}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    echo "{$prefix}" >> "$OUTPUT_FILE"
    count=$((count + 1))
    numdeepZ=$(ls *_${prefix}*_zdeep_*_roms.png 2>/dev/null | wc -l)
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numdeepZ
echo "\\\\newcommand*{\\\\numdeepZ}{${numdeepZ}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get headers
headers=$(echo "$files" | cut -d'_' -f1)
echo "\\\\newcommand*{\deepHdr}{" >> "$OUTPUT_FILE"
for this_header in $headers; do
    echo "{$this_header}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numdeepVars
echo "\\\\newcommand*{\\\\numdeepVars}{${count}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get number of diags
echo "\\\\newcommand*{\deepDiags}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    file_count=$(ls *_${prefix}_zdeep_1_diag*.png 2>/dev/null | wc -l)
    echo "{$file_count}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

##############################################
# Get transect vars, locations, headers, diags
##############################################
count=0
files=$(ls *_transect_1_roms.png)
prefixes=$(echo "$files" | sed 's/...//' | sed 's/_transect_1_roms.png//')
echo "\\\\newcommand*{\sectVars}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    echo "{$prefix}" >> "$OUTPUT_FILE"
    count=$((count + 1))
    numSects=$(ls *_${prefix}*_transect_*_roms.png 2>/dev/null | wc -l)
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numSects
echo "\\\\newcommand*{\\\\numSects}{${numSects}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get headers
headers=$(echo "$files" | cut -d'_' -f1)
echo "\\\\newcommand*{\sectHdr}{" >> "$OUTPUT_FILE"
for this_header in $headers; do
    echo "{$this_header}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numsectVars
echo "\\\\newcommand*{\\\\numsectVars}{${count}}" >> "$OUTPUT_FILE"
echo " " >> $OUTPUT_FILE

# Get number of diags
echo "\\\\newcommand*{\sectDiags}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    file_count=$(ls *_${prefix}_transect_1_diag*.png 2>/dev/null | wc -l)
    echo "{$file_count}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

#############################################
# Get gridded vars, headers, number of depths
#############################################
count=0
files=$(ls *_gridded_obs*_1.png)
prefixes=$(echo "$files" | sed 's/.*_obs_\([a-zA-Z]*\)_.*/\1/')
echo "\\\\newcommand*{\griddedVars}{" >> "$OUTPUT_FILE"
for prefix in $prefixes; do
    echo "{$prefix}" >> "$OUTPUT_FILE"
    count=$((count + 1))
    numgriddedZ=$(ls *gridded_obs_${prefix}*.png 2>/dev/null | wc -l)
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numgriddedZ
echo "\\\\newcommand*{\\\\numgriddedZ}{${numgriddedZ}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get headers
headers=$(echo "$files" | cut -d'_' -f1)
echo "\\\\newcommand*{\griddedHdr}{" >> "$OUTPUT_FILE"
for this_header in $headers; do
    echo "{$this_header}" >> "$OUTPUT_FILE"
done
echo "}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

# Get numgriddedVars
echo "\\\\newcommand*{\\\\numgriddedVars}{${count}}" >> "$OUTPUT_FILE"
echo " " >> "$OUTPUT_FILE"

#############################################
# Finalize PDF 
#############################################
# Move variables to working directory
mv ${OUTPUT_FILE} ../../.
cd ../../.
pdflatex diagnostics_template_auto.tex
sleep 5
# Remove other output files (if they exist)
rm texput.log
rm diagnostics_template_auto.lof
rm diagnostics_template_auto.log
rm diagnostics_template_auto.aux
rm diagnostics_template_auto.synctex*
echo "DONE! Check diagnostics_template_auto.pdf"
