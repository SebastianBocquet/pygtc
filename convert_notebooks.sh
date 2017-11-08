#!/bin/bash
# Call this with a list of jupyter notebooks to convert to rst. It will automatically put them in the docs folder with links to _static
for var in "$@"
do
  jupyter nbconvert --to rst "$var" --output-dir docs --NbConvertApp.output_files_dir="_static/{notebook_name}_files"
done
