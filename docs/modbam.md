# Modification visualization

## Generate mod.sam

Run `SWARM_read_level.py` with --modsam tag.

*Current implementation only supports single-threaded preprocessing which may result in slower prediction compared to default for models with Mini architecture (Ψ and m5C). 

## Process modbam 

### Filter coordinates
Run postprocess/filter_modbam_sites.py to filter mod.bam files for significant sites of interest.

```bash
# first convert mod.sam to mod.bam
samtools view -b $OUT.mod.sam > $OUT.mod.bam

# run filtering
python3 postprocess/filter_modbam_sites.py input.mod.bam out.sig.mod.bam sig.sites.bed
```

### Merge different modifications

Run postprocess/merge_modbams.py to merge mod.bam files with different modifications predicted in same reads.

--merge takes a csv file with each row being a bam for a given modification, and two columns: 
1) path to mod.bam
2) modification <pU, m6A, m5C>

```bash
python3 postprocess/merge_modbams.py --merge merge.input.csv -o merged.sig.mod.bam
```

### Filter modification tags

Genome browsers visualise raw modification calls using colour scales that exaggerate the confidence in modification presence, having highly visible colours even for modification calls with probability under 0.5.  

Run postprocess/merge_modbams.py with -t tag for setting minimun read-level probability to be visualised
```bash
python3 postprocess/merge_modbams.py --merge merge.input.csv -o merged.sig.mod.bam -t 0.9
```


