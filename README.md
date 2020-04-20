# OLC

## Installation


### Getting the latest source code with git

```
# Get a local copy of MTG-Link source code
git clone --recursive https://github.com/anne-gcd/OLC.git
```


## User Manual

### Usage

The OLC command line interface is composed of multiple parameters. You can get a summary of all available parameters by running:
```
./olc.py --help

usage: olc.py -in <input_sequence> -reads <reads_file> -s <seed_size> -o <minimum_overlap_size> -l <maximum_assembly_length> [options]
                                
Gapfilling, using an Overlap-Layout-Consensus (OLC) method

optional arguments:
  -h, --help            show this help message and exit
  -reads READS          File of reads
  -s SEED_SIZE          Seed size used for indexing the reads (bp)
  -o MIN_OVERLAP        mMinimum overlapping size (bp)
  -a [ABUNDANCE_MIN [ABUNDANCE_MIN ...]]
                        Minimal abundance(s) of reads used for gapfilling ; 
                        extension's groups having less than this number of reads are discarded from the graph
  -l MAX_LENGTH         Maximum assembly length (bp) (it could correspond to the length of the gap to fill (+length input sequences)
                        OR it could be a very high length ; to prevent for searching indefinitely
  -out OUTDIR           Output directory for the results' files

```



