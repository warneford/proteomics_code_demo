# Bash script strict mode (helps with debugging)
#!/bin/zsh
#set -eo pipefail #Refer to http://redsymbol.net/articles/unofficial-bash-strict-mode/ for details
IFS=$'\n\t'
#record start time
#Extract 5' umi scstart=$(date +%s)
echo "\n\n\n\n\ndia-NN LC-MS/MS database search"
echo "updated 16May2022 / Robert Warneford-Thomson"

# #Initialize
zstart=$(date +%s)

HOME=/home/rwt
wd=`pwd`
raw_dir=$HOME/shared/RBRID_DIA
threads=120

mkdir -p $wd/output
mkdir -p $wd/logs
mkdir -p $wd/temp
mkdir -p $wd/mzML; cd $wd/mzML

# Convert .raw files to open source format .mzML using ThermoRawFileParser

pstart=$(date +%s) #record start
conda activate py3
cd $raw_dir; ls *raw | parallel -j$threads "mono \
/home/rwt/packages/anaconda3/envs/py3/bin/ThermoRawFileParser.exe \
-i $raw_dir/{} \
-o $wd/mzML/ \
--format=2 \
--gzip"

pend=$(date +%s) #record end
runtime=$(($pend-$pstart))

echo $(printf '%dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))) \
"runtime to convert .raw files to indexed .mzML\n"

#unzip mzML files
cd $wd/mzML; mkdir -p $wd/unzipped_mzML  
ls *gz | parallel "unpigz -k -p $threads {}"
mv *mzML $wd/unzipped_mzML 

#set up dia-NN database search
pstart=$(date +%s) #record start

cd $wd

# set up DiaNN command
cat << EOF > $wd/dia_nn.cfg
--dir $wd/unzipped_mzML \
--lib   \
--threads $threads \
--verbose 1  \
--qvalue 0.01  \
--matrices  \
--out $wd/output/report.tsv  \
--out-lib $wd/logs/dia_nn.report-lib.tsv \
--temp $wd/temp \
--gen-spec-lib  \
--predictor  \
--fasta $HOME/proteomes/uniprot/human/uniprot_human_sp_tr_plus_isoforms_rel_2022_01_up000005640.fasta
--fasta-search  \
--min-fr-mz 200  \
--max-fr-mz 3000  \
--met-excision  \
--cut K*,R*  \
--missed-cleavages 1  \
--min-pep-len 6  \
--max-pep-len 40  \
--min-pr-mz 300  \
--max-pr-mz 1800  \
--min-pr-charge 2  \
--max-pr-charge 4  \
--unimod4  \
--var-mods 1  \
--var-mod UniMod:35,15.994915,M  \
--var-mod UniMod:1,42.010565,*n  \
--reanalyse  \
--smart-profiling  \
--pg-level 1  \
--peak-center  
EOF

#run DiaNN
diann-1.8.1 `cat $wd/dia_nn.cfg` 

pend=$(date +%s) #record end
runtime=$(($pend-$pstart))

echo $(printf '%dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))) \
"runtime to perform dia-NN database search\n"
