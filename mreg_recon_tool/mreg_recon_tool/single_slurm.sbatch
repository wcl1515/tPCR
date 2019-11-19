#!/bin/bash

#SBATCH --job-name=single
# #SBATCH --workdir=<directory> 	# set the working dir of the batch script to directory before executed
#SBATCH --output=single.out	
#SBATCH --error=single.err
#SBATCH --nodes=1
#SBATCH --partition=engine 		# Options: [testmode, engine]
#SBATCH --nodelist= 	# Options: [yeti,nessie,ram], empty: no node specification
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00  		# Format: [hh:mm:ss]
#SBATCH --mem=100 			# suffix [K,M,G,T], default [M]
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user=<your-mail>

# Create a local work directory
# mkdir -p ./$SLURM_JOB_ID

# any linux terminal command is possible here
echo test

# start matlab and run your matlab-script (*.m)
# /usr/local/MATLAB/R2017b/bin/matlab -singleCompThread -nodisplay -nosplash -r "try, reconstruction_off_mo(101:1:200), catch ME, disp(ME.message), exit, end"
# matlab -singleCompThread -nodisplay -nosplash -r "try, simu(141:1:150,5e-5,'l1',48,1,1), catch ME, disp(ME.message), exit, end"
 matlab -singleCompThread -nodisplay -nosplash -r "try, PC_recon_slurm, catch ME, disp(ME.message), exit, end"
# as soon as matlab is finished the next lines will be executed. For example:
echo -e 'This is not a matlab command/n Matlab is closed after finishing your script.'
# have a look in your outputfile.
