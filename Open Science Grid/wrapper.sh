# Written by Jose Alejandro Medina. 
#!/bin/sh
# $1 = pdb file
# $2 = run number
pdb_file=$1
run_number=$2

# First, make the dnaseqs file
linesperrun=10
# Find what line number to start at
start_line=$((linesperrun*run_number + linesperrun))

# Create the expected scwrl directory
mkdir -p scwrl4_lin
cp Scwrl4 Scwrl4.ini bbDepRotLib.bin scwrl4_lin/ 

echo "Before"
ls -l

cat dnaseqs | head -n $start_line | tail -n $linesperrun > dnaseqs.new

mv dnaseqs.new dnaseqs

./calculaterobustnessLOOP $pdb_file

rm dnaseqs
