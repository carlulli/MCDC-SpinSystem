# shell script that reads in parameter file and runs main

grep -v '^#' < ./params.txt | while IFS=$'\t' read -r -a line
do
  echo $line
done
