##fineSTRUCTURE##

#######################################################################
# This script is largely derived from Rui's script however it has additional comments and starts with
# phasing with a reference panel and installing fineSTRUCTURE (assumes plink format)
# If a recombination map is not available please generate a uniform recombination map
#######################################################################


###### Notes before starting #######
# 1.)It is essential to phase all files for any chromopainter run simultaneously to avoid bias
#    shapeit works well for these purposes and is faster than Beagle etc
#    as fineSTRUCTURE is quite robust to minor phasing error shapeit is sufficient. 
#
# 2.)chromopainter is not well optimised for dealing with missing data so I would reccommend removing
#    missingness in your data before running.
#####################################


################################################
# Step 1
# Phasing with Shapeit:
################################################

#1) Splitting into chromosomes
for chr in $(seq 1 22); do
plink --bfile [file] --chr $chr --make-bed --out chr.$chr
done


#2.)Shapeit basic allignment to reference panel (Here I use the 1000GP phase 3; available )
#   This will give you a list of snps to exclude due to strand error etc (see chr1.alignments.snp.strand.exclude)
for chr in $(seq 1 22); do
shapeit -check --input-bed chr.$chr.bed chr.$chr.bim chr.$chr.fam  --input-map 1000GP_Phase3/genetic_map_chr"$chr"_combined_b37.txt --input-ref 1000GP_Phase3/1000GP_Phase3_chr"$chr".hap.gz 1000GP_Phase3/1000GP_Phase3_chr"$chr".legend.gz 1000GP_Phase3/1000GP_Phase3.sample --output-log chr$chr.alignments
done &
  
#3.) Shapeit phasing once alligned
#    NB: Make sure all of the files have a chr$chr.alignments.snp.strand.exclude file, 
#    if not make an decoy empty file with the title chr$chr.alignments.snp.strand.exclude
#    this line will throw an error if the file does not exist

for chr in $(seq 1 22); do
shapeit --input-bed chr.$chr.bed chr.$chr.bim chr.$chr.fam --input-map 1000GP_Phase3/genetic_map_chr"$chr"_combined_b37.txt  --input-ref 1000GP_Phase3/1000GP_Phase3_chr"$chr".hap.gz 1000GP_Phase3/1000GP_Phase3_chr"$chr".legend.gz 1000GP_Phase3/1000GP_Phase3.sample --exclude-snp chr$chr.alignments.snp.strand.exclude --output-max phasedref.$chr.haps phasedref.$chr.sample --thread 10
done 

#################################################
# Step 2:
# Other preprocessing: .phase and recomb files
#################################################

# 1.) Convert to chromopainter.phase file format 
for chr in $(seq 1 22); do 
perl impute2chromopainter.pl -r 1000GP_Phase3/genetic_map_chr"$chr"_combined_b37.txt  phasedref.$chr.haps chromopainter.$chr.phased -f -J & done


#If the first line is a zero in your phase file, remove it to make it chromopainterV2 compatible 
#Not sure if latest version of impute2chromopainter.pl does this for you, but the version I use leaves it in 
for line in $(seq 1 22); do
sed '1d' chromopainter.$line.phased.phase > new_chromopainter.$line.phased.phase
done

#2.) Making recombination map files for chromopainter
##Make genetic map have 4 columns. You can do this in many ways. I use:
for i in $(seq 1 22); do  
cat 1000GP_Phase3/genetic_map_chr"$i"_combined_b37.txt|awk 'BEGIN {}{print "$i" " " $0}END{print "$i" " " $0}' > recomb.$i.txt; done

##make recombination maps specific to your phase files
for chr in $(seq 1 22); do 
perl convertrecfile.pl -M hap chromopainter.$chr.phased.phase recomb.$chr.txt recomb_map.$chr.txt 
done

#3.) Make an idfile
#making ids file from a fam file is easy and sufficient for basic fineSTRUCTURE usage
#Just need a list of identifiers in order matching your phase file for now
cat chr.1.fam  | cut -f1,2 -d' '  | tr ' ' '_' > idfile

#once/if you have population lables and want to use chromopainterV2 use the format:
#  col1: id  col2: population lable  col3: 1 for inclusion, 0 for exclusion

#######
#Running finestructure
######


#transfer to ICHEC
# (at fionn)
#Make an analysis folder ie 
mkdir Analysis
mkdir Analysis/chromopainter_out

echo "Analyisis" > list
##Note we use list file in case you want to run many analyses at once (just add their folder name to the list.)
##Populate Analysis/chromopainter_out with the "new_chromopainter*.phased.phase and recomb_map.*.txt files and idfile(use SCP or SFTP)


###############
Install finestructure v2 
###############
mkdir fineSTRUCTUREv2; cd fineSTRUCTUREv2
wget http://www.maths.bris.ac.uk/~madjl/finestructure/fs-2.0.8.tar.gz
tar xvfz fs-2.0.8.tar.gz; cd fs-2.0.8/
module load dev gcc/5.1.0 libs gsl/gcc/1.16
./configure CC=gcc CXX=g++ CFLAGS="-O3 -mtune=native" CXXFLAGS="-O3 -mtune=native" --prefix=/ichec/work/tclif[your extension]/fineSTRUCTUREv2/fs-2.0.8
make
make install

##Nb 


#load modules
##these are the modules for finestructure. (fs-2.0.8)
module load dev gcc/5.1.0 libs gsl/gcc/1.16
module load dev intel/2013-sp1
module load dev java/1.7.0_45
module load apps taskfarm/latest
module load libs boost/intel/1.55.0


##############################
#                            #
#         Stage 1            #
#                            #
##############################

# # # # estimating the parameters mu and Ne # # # # 
# (mutation rate, effective population size)
##finestructure will write a commandfile for the EM step which can be taskfarmed on HPC
##This will include a line for each individual in your phasefile chr"$chr"_cp/commandfiles/commandfile1.txt


###Important: edit the pbs file to include your project name ie: -A tclif025c, email (-M email) and an appropriate wall time and number of nodes.

cat list | while read i
do 
echo "Processing Sample $i"
cd $i
for chr in $(seq 22)
do
#STAGE1
#pre-stage1
fs  chr"$chr"_cp.cp -n -phasefiles chromopainter_out/new_chromopainter.$chr.phased.phase -recombfiles chromopainter_out/recomb_map.$chr.txt -hpc 1 -idfile chromopainter_out/idfile -writes1 
done

#make pbs for stage 1
## Remember to set your pbs wall time appropriately for the size of your data; 
##If unsure try running your largest chromosome for 24 hours on say 4 nodes. (Needed for 3000 plus inds; 2 nodes are sufficient for ~1000 inds)

for chr in $(seq 1 22)
do 
echo "
#!/bin/bash
#PBS -l nodes=2:ppn=24
#PBS -l walltime=24:00:00
#PBS -A [tclif-your-project]
#PBS -j oe
#PBS -m bea
# send it to this address
#PBS -M [your email]
# inherit the current environment
#(necessary if pre-loading modules for example)
# This job's working directory
#echo Working directory is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR
# For a MPI run on the Stokes system
module load dev gcc/5.1.0 libs gsl/gcc/1.16
module load dev intel/2013-sp1
module load dev java/1.7.0_45
module load apps taskfarm/latest
module load libs boost/intel/1.55.0
export TASKFARM_PPN=24
taskfarm ./chr"$chr"_cp/commandfiles/commandfile1.txt" > run_S1_chr"$chr".pbs
done

echo "Finished writes1 Sample $i"

cd ..
done

###
# running stage 1
###
cat list | while read i
do 
echo "Processing Sample $i"
cd $i

#run stage 1

for chr in $(seq 1 22)
do
qsub run_S1_chr"$chr".pbs & done
cd ..
done

####################
#When stage 1 is done, combine files:
###################

cat list |while read i
do 
echo "Processing Sample $i"
cd $i

for chr in $(seq 22)
do

fs  chr"$chr"_cp.cp -n -phasefiles chromopainter_out/new_chromopainter.$chr.phased.phase -recombfiles chromopainter_out/recomb_map.$chr.txt -hpc 1 -idfile chromopainter_out/idfile -combines1 
done
cd ..
done

##############################
#                            #
#      Stage 2               #
#                            #
##############################
# # # # Estimating 'c' and creating the genome-wide chromopainter output for all individuals.  # # # # 

#This generally takes up less time than stage 1 so you can reduce your walltime
#run the same command as above, if -go is selected fs will find out what needs to be done
#Note that for large datasets with many inds or snps however this step is memory intensive
#In this case change the number of cores running at a time by editing the line export TASKFARM_PPN=n
#Also change the number of nodes needed as with stage 1


cat list | while read i
do 
echo "Processing Sample $i"
cd $i

for chr in $(seq 1 22)
do
#STAGE2
fs  chr"$chr"_cp.cp -n -phasefiles chromopainter_out/new_chromopainter.$chr.phased.phase -recombfiles chromopainter_out/recomb_map.$chr.txt -hpc 1 -idfile chromopainter_out/idfile -go
done


#something like this should show up
#Inferred Ne=43.0635 and mu=0.00164003
#HPC mode: commands for stage 2 written to file chr22_cp_commandfile2.txt. Rerun when those commands have been completed and the results copied back to this directory.
#For local parallel execution try "cat chr22_cp_commandfile2.txt | parallel"
#IMPORTANT: The run ended with a requirement to run commands externally from file "chr22_cp_commandfile2.txt".

cd ..
done

#make pbs for stage 2


cat list |while read i
do 
echo "Processing Sample $i"
cd $i


for chr in $(seq 1 22)
do 
echo "
#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=24:00:00
#PBS -A [project-name]
#PBS -j oe
#PBS -m bea
# send it to this address
#PBS -M [email]
# inherit the current environment
#(necessary if pre-loading modules for example)
# This job's working directory
#echo Working directory is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR
# For a MPI run on the Stokes system
module load dev gcc/5.1.0 libs gsl/gcc/1.16
module load dev intel/2015-u2
module load dev java/1.7.0_45
module load apps taskfarm/latest
module load libs boost/intel/1.55.0
export TASKFARM_PPN=24
taskfarm ./chr"$chr"_cp/commandfiles/commandfile2.txt" > run_S2_chr"$chr".pbs
done

#run stage 2
for chr in $(seq 1 22)
do
qsub run_S2_chr"$chr".pbs & done
cd ..
done


##once run we can repeat -go command 
cat list | while read i
do 
echo "Processing Sample $i"
cd $i
#Combine stage 2
for chr in $(seq 1 22)
do
fs  chr"$chr"_cp.cp -n -phasefiles chromopainter_out/new_chromopainter.$chr.phased.phase -recombfiles chromopainter_out/recomb_map.$chr.txt -hpc 1 -idfile chromopainter_out/idfile -go
done
echo "Combined stage 2 for Sample $i"


# first, you need to combine the chromosomes to do whole genome analysis
fs chromocombine -o pooled chr1_cp_linked chr2_cp_linked chr3_cp_linked chr4_cp_linked chr5_cp_linked chr6_cp_linked chr7_cp_linked chr8_cp_linked chr9_cp_linked chr10_cp_linked chr11_cp_linked chr12_cp_linked chr13_cp_linked chr14_cp_linked chr15_cp_linked chr16_cp_linked chr17_cp_linked chr18_cp_linked chr19_cp_linked chr20_cp_linked chr21_cp_linked chr22_cp_linked
echo "Finished Chromocombine for Sample $i"
cd ..
done

# Note: Inspect value of rescaling factor "c". If c<0.05 read this
# c is a summary of the data, something like effective number of chunks

# For interpretation of the empirical evaluation presented, we note that when c
# is ‘too large’, the effective number of chunks is reduced and therefore any mistakes
# in population assignment will tend to be under-split, i.e. we will not distinguish
# efficiently between similar populations. When c is ‘too small’ our model believes
# it has more independent chunks than is true and therefore will tend to over-split
# populations. The smallest c that does not over-split is called efficient, and larger
# c are called conservative.
# see this, for more info
# http://www.maths.bris.ac.uk/~madjl/finestructure/manualse12.html


####################
#                  #
#      Stage 3     #
#                  #
####################

#this stage is finestructure mcmc and making the tree 
#Note our server has an old version of fineSTRUCTURE on it so if you want to follow the online manual you will need to install the 
#Latest version locally.
#In my experience they are very consistent (99.8%) for basic functionality

# -x <num>        Number of burn in iterations for MCMC method.
# -y <num>        Number of sample iterations for MCMC method.
# -z <num>        Thin interval in the output file, for MCMC method.


##If your individual names start with a number they will not run. 
#In this case attempt the following hack to add an A to all the ind names

cat pooled.chunkcounts.out|sed 's/^/A/g'|sed '2s/[^[:blank:]]\{1,\}/A&/g'|sed '2s/ARecipient/Recipient/g'|sed '2s/ARecipient/Recipient/g'|sed '1s/A#Cfactor/#Cfactor/g' > pooled.chunkcounts_fixed.out


##Need to run this locally. Not feasible to run on ICHEC
#run fs small iter to make sure data is sensible
fs fs -x 100000 -y 100000 -z 100 pooled.chunkcounts.out smallmcmc
fs fs -x 10000 -m T pooled.chunkcounts_fixed.out smallmcmc smalltree

##For visualisation I would reccomend using lawson's script and not finegui, but if you must
finegui -c pooled.chunkcounts.out -m smallmcmc -t smalltree ### The tree often won't load, use Lawson's script


#run fs big iter; This may take weeks so do in screen or tmux in case you disconnect
#Best to run twice to test for convergence
##note choose the -t switch to match the number of combinations you want the merge tree step to try (often the default of 1500 is too small for larger datasets)
fs fs -x 3000000 -y 1000000 -z 100  pooled.chunkcounts.out bigmcmc1
fs fs -x 10000 -m T -t 10000 pooled.chunkcounts.out big_mcmc bigtree1

fs fs -x 3000000 -y 1000000 -z 100  pooled.chunkcounts.out bigmcmc2
fs fs -x 10000 -m T -t 10000 pooled.chunkcounts.out big_mcmc bigtree2

##For a visual comparison of MCMC convergence see my hacks script.

####Note if the mcmc doesn't converge you can pick it up using 
fs fs -x 0 -y [num] -z [num] [chunkcounts] [old_out.mcmc.xml] [out.mcmc.longer.xml]

##Rerunning tree building with maximum concordance state
##This is the improvement for cluster assignment that Leslie et al. use in POBI paper
##locally install latest version
./fs fs -x 10000 -m T -T 1 -v -k 2 -t 10000 pooled.chunkcounts.out big_mcmc bigtree_maxconcordance


# Note: there is a lot of experimentation to be done before trusting a result
# and documentation is really poor.
# If you have doubts about something, best thing is to email Dan Lawson












