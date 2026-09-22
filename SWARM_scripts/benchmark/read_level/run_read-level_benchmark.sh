module load singularity

SIF=tensorflow_24.01-tf2-py3-pysam.sif
FASTA=IVT-m_data/I-M_combined.fasta
SCRIPT=../../SWARM_read_level.py ;

mkdir outputs

for target_mod in "m6A" "m5C" "pU"; do
        for sample_mod in "m6A" "m5C" 'pU' "NM"; do
                BLOW5=IVT-m_data/blow5/IVT-m1-m2.$sample_mod.blow5 ;
                SAM=IVT-m_data/events/IVT-m1-m2.$sample_mod.events.sam ;
                OUT=outputs/IVT-m1-m2.sample-${sample_mod}.target-${target_mod} ;
                singularity exec --nv $SIF python3 $SCRIPT -m $target_mod --sam $SAM --fasta $FASTA --raw $BLOW5 -o $OUT ;
        done ;

done
