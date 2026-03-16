# Modification visualization

## Generate mod.bam

Run `SWARM_read_level.py` with --modbam tag.

*Current implementation only supports single-threaded preprocessing which may result in slower prediction compared to default for models with Mini architecture (Ψ and m5C). 

## Process modbam 

Use `process_modbam.py` to process modbam files for more informative visualisation of modifications. 

### Merge different modifications

```bash
process_modbam.py --inputs <INPUTS.tsv> -o <OUT>
```


### Filter modification tags

Genome browsers visualise raw modification calls using colour scales that exaggerate the confidence in modification presence, having highly visible colours even for modification calls with probability under 0.5.  

```bash
process_modbam.py -i <mod.bam> -o <OUT> -r <0.5> 
```


### Filter coordinates

It is often useful to highlight modification trends in long reads only for a particular set of positions. 

```bash
process_modbam.py -i <mod.bam> -o <OUT> --positions <BED>
```

