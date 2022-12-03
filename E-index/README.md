Codes for generating data/figures used in Fig. 6, and S13.

`build.end.model` is designed to build the end model from control samples;
`calc.E-index.multi-thread` is designed to calculate E-index for testing samples;
`hg38.info` contains the size information for each chromosome in human genome (hg38);
`ENCODE.blacklist.hg38.bed` contains problematic regions collected from
[Amemiya et al. Scientific Reports 2019, 9:9354](https://www.nature.com/articles/s41598-019-45839-z)
(used in `calc.E-index.multi-thread` as `mask.bed`).

For `build.end.model` and `calc.E-index.multi-thread`, run the program without parameters will show its usage.