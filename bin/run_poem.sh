#!/bin/bash

# get the scripts path
SCRIPT=`realpath -s $0`
SCRIPTPATH=`dirname $SCRIPT`

# set the path for the python
#lpython=/home/xiaoh/Downloads/compiler/intel/intelpython27/bin/python
#python=pypy
lpython=python
python=python

#######################################
# parse the args
######################################
while [ $# -gt 1 ]
do
    key="$1"
    case $key in
        -f|--fasta)
        fas="$2"
        shift # past argument
        ;;
        -a|--assembly)
        asm="$2"
        shift # past argument
        ;;
        -p|--predict)
        gpd="$2"
        shift # past argument
        ;;
        -l|--fasta)
        lr="$2"
        shift # past argument
        ;;
        *)
            # unknown option
        ;;
    esac
    shift # past argument or value
done

if [ ! $fas ]; then

	echo '#######################################'
	echo '#'
	echo '# usage:'
    echo 'for genome|assembly|contig'
	echo '$ bash this_script.sh -f genome.fsa -a n -p prokka'
    echo ''
    echo 'for short reads (length <= 600)'
	echo '$ bash this_script.sh -f reads.fsa -a y -p prokka -l n'

    echo 'for short reads (length > 600)'
	echo '$ bash this_script.sh -f reads.fsa -a y -p prokka -l y'

	echo '#'
	echo '#######################################'
	echo ''
	exit 1

fi


# make an temp dir for the intermediate results
temp=$fas\_output
mkdir -p $temp


########################
# assembly with IDBA_ud
########################
if [[ $asm == "Y" ]] || [[ $asm == "y" ]]
then
    echo "assembly mode"

    if [[ $lr == "n" ]] || [[ $lr == "N" ]]; then
        idba_ud -r $fas -o $temp/assembly --pre_correction > $temp/asm.log
    else
        idba_ud -l $fas -o $temp/assembly --pre_correction > $temp/asm.log
    fi
    # check the header of the fasta
    awk -F" " '{if($1~/^>/){print $1}else{print $0}}' $temp/assembly/contig.fa > $temp/input.fsa
else
    echo "no assembly mode"
    # check the header of the fasta
    awk -F" " '{if($1~/^>/){print $1}else{print $0}}' $fas > $temp/input.fsa
fi

fasta=$temp/input.fsa


#exit 1

############################################
# gene prediction by Metagenemark or prokka
############################################
echo '##########################################################################'
echo 'step 1: Gene prediction'
echo '##########################################################################'
echo ''


if [[ $gpd == "gmk" ]] || [[ $gpd == "genemark" ]]; then

    #gmhmmp=/home/xiaoh/Downloads/genome/evaluator/quast-3.1/libs/genemark/linux_64/gmhmmp 
    gmhmmp=gmhmmp 
    $gmhmmp -A $fasta\_gmk_aa.fsa -p 0 -f G -m $SCRIPTPATH/../config/MetaGenemark/MetaGeneMark_v1.mod $fasta

elif [[ $gpd == "prokka" ]] || [[ $gpd == "pka" ]]; then

    #/usr/bin/perl /home/xiaoh/Downloads/compiler/intel/intelpython27/bin/prokka --quiet --fast --prefix prokka_out --metagenome --force --outdir $fasta\_prokka $fasta
    #prokka --quiet --fast --prefix prokka_out --metagenome --centre X --compliant --force --outdir $fasta\_prokka $fasta
    # prokka needs clean contig names:
    $python $SCRIPTPATH/../lib/clean_header.py $fasta > $fasta\_new.fsa
    prokka --quiet --fast --prefix prokka_out --metagenome --force --outdir $fasta\_prokka $fasta\_new.fsa
    #$python $SCRIPTPATH/../lib/pka2gmk.py $fasta $fasta\_prokka/prokka_out.faa $fasta\_prokka/prokka_out.gff > $fasta\_gmk_aa.fsa
    $python $SCRIPTPATH/../lib/pka2gmk.py $fasta\_new.fsa $fasta\_prokka/prokka_out.faa $fasta\_prokka/prokka_out.gff > $fasta\_gmk_aa.fsa

else
    exit 1
fi


#########################################
# re-id the aa predicted by metagenemark
# the features is sep by ||(not |)
#########################################
echo '##########################################################################'
echo 'step 2: Modify the header of predicted gene...'
echo '##########################################################################'
echo ''

$python $SCRIPTPATH/../lib/reid.py $fasta\_gmk_aa.fsa > $fasta\_aa.fsa

##########################
# cd-hit to find ortholog
##########################
echo '##########################################################################'
echo 'step 3: Cdhit to do clustering and filter the close and remote species...'
echo '##########################################################################'
echo ''

cd-hit -c .70 -aL .85 -M 12000 -T 0 -d 268435456 -i $fasta\_aa.fsa -o $fasta\.cdhit > $temp/run.log

##############################
# only keep the .7 < id < .98
##############################
$python $SCRIPTPATH/../lib/cdhit_parse_filter.py $fasta\.cdhit.clstr $fasta\_aa.fsa > $fasta\.flt

#################
# cog annotation
#################
echo '##########################################################################'
echo 'step 4: Use blast or Diamond to find homologous in cog database...'
echo '##########################################################################'
echo ''
$python $SCRIPTPATH/../lib/cog_annot.py $fasta\.flt $SCRIPTPATH/../database/cog/cog2014 diamond | uniq | sort -rk1 | uniq  > $fasta\.cog


########################################################
# the nr annotation and binning?, use diamond very fast
# diamond
#
# use the networkx to store and find operon
########################################################
echo '##########################################################################'
echo 'step 5: Operonic adjacency prediction...'
echo '##########################################################################'
echo ''
# convert the gene predict to the format the operon predictor need.
$python $SCRIPTPATH/../lib/to_list.py $fasta\_aa.fsa > $fasta\.locus

echo "$lpython $SCRIPTPATH/../lib/deep_operon.py predict $fasta $fasta\.locus $SCRIPTPATH/../config/Operon_Predictor/model.hdf5 > $fasta\.adjacency"
$lpython $SCRIPTPATH/../lib/deep_operon.py predict $fasta $fasta\.locus $SCRIPTPATH/../config/Operon_Predictor/model.hdf5 > $fasta\.adjacency

$lpython $SCRIPTPATH/../lib/adj2operon.py $fasta\.adjacency $fasta\.cog > $fasta\.operon



########################################################
# the nr annotation and binning?, use diamond very fast
# diamond
#
# use the networkx to store and find operon
########################################################
echo '##########################################################################'
echo 'step 6: find the core genes of operon...'
echo '##########################################################################'
echo ''

#$python $SCRIPTPATH/../lib/find_core.py $fasta\.cog $fasta\.adjacency 2 $SCRIPTPATH/../database > $fasta\.cog_adjacency.txt
#$python $SCRIPTPATH/../lib/find_core.py $fasta\.cog $fasta\.adjacency 5 > $fasta\.cog_adjacency.txt
$python $SCRIPTPATH/../lib/find_core.py $fasta\.cog $fasta\.adjacency 5 > $fasta\.cog_adjacency



##########################################################
# parse the results to cytoscape's node, edge and network
##########################################################
echo '##########################################################################'
echo 'step 7: Output the result for cytoscape...'
echo '##########################################################################'
echo ''

#$python $SCRIPTPATH/../lib/cytoscape.py $fasta\.core_cog_adjacency $SCRIPTPATH/database/cog/cog2014/cognames2003-2014.tab
$python $SCRIPTPATH/../lib/cytoscape.py $fasta\.core_cog_adjacency



