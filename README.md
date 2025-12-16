# NAGNAG 3 Splice Sites
Code developed and used for analyses presented in _The Role of NAGNAG 3' Splice Sites in the Flow of Genetic Information_

### Code

Scripts for data parsing and analysis are in `nagnags`
To run all analyses, call `./run.sh`

Because HGMD data is not publically available, it is not provided here. Thus, for related analyses to run, you must provide them yourself. The program assumes that there is a folder named `HGMD_pro` in `src` containing tab-delimited files `splice.txt` and `hg38_coords.txt`. Alternatively, you may specify their locations with command-line options `-s` and `-c`, respectively. To skip HGMD analyses, use the `-H` flag.

## Support
If you have questions, email hpearce@unc.edu or alain@unc.edu
