# Written by Jose Alejandro Medina. 
#!/bin/sh
pdb_dir=$1

# test the input to be correct
if [ "x" = "x$pdb_dir" ]; then
   echo "Please specify a pdb directory"
   exit 1
fi

# Make the condor output directory
mkdir -p condor_out

# First, create base submit file
cat > submit.condor << EOF
universe = vanilla
arguments = \$(PDB_FILE_BASENAME) \$(PROCESS)
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
output = condor_out/output.\$(PDB_FILE_BASENAME).\$(PROCESS)
error = condor_out/error.\$(PDB_FILE_BASENAME).\$(PROCESS)
executable = wrapper.sh
transfer_input_files = \$(PDB_FILE), calculaterobustnessLOOP, dnaseqs, scwrl4_lin/Scwrl4, scwrl4_lin/Scwrl4.ini, scwrl4_lin/bbDepRotLib.bin, vendruscolomatrix
transfer_output_remaps = "logfile=logfile.\$(PDB_FILE_BASENAME).\$(PROCESS); robustout=robustout.\$(PDB_FILE_BASENAME).\$(PROCESS)"
request_memory = 1024
request_cpus = 1
log = condor.log
EOF

# Count the number of lines in aastable and divide by 10
lines=`wc -l dnaseqs  | awk '{print $1}'`
num_submits=$(($lines / 10))

# Now, for each pdb file in the directory given as the argument, create a submit line
for pdb_file in `find $pdb_dir -type f`; do
cat >> submit.condor << EOF
PDB_FILE=$pdb_file
PDB_FILE_BASENAME=`basename $pdb_file`
queue $num_submits
EOF
done

# Print a helpful message
echo "Created a condor submit file submit.condor"
echo "To submit the job, issue the command:"
echo "condor_submit submit.condor"
