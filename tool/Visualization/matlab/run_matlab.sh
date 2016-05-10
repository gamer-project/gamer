export MATLABPATH=YOUR_SCRIPT_PATH
SCRIPT=YOUR_SCRIPT_NAME

matlab -nodisplay -nodesktop -nojvm -nosplash -r "$SCRIPT; exit;"
