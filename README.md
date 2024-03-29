### Codes and programs used in [An et al. Nature Communications 2023](https://www.nature.com/articles/s41467-023-35959-6 "An Nature Communications 2023").

---

Distributed under the [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/ "CC BY-NC-ND")
license for **personal and academic usage only**.

Please refer to each directory (where `wk.sh` is the main script) for detailed information.

Notes:
1. Codes to process the raw fastq files (e.g., preprocessing, alignment, and duplicate-removal) are NOT covered here;
you may need to install [Ktrim](https://github.com/hellosunking/Ktrim/),
[bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
[Msuite2](https://github.com/hellosunking/Msuite2/),
and [samtools](https://www.htslib.org/) to perform these procedures;
2. [bedtools](https://github.com/arq5x/bedtools2) and [R](https://www.r-project.org/) are required to run these codes;
3. Nucleosome track is obtained from [NucMap](https://ngdc.cncb.ac.cn/nucmap/),
here is the [link to the file used in this work](https://download.cncb.ac.cn/nucmap/organisms/v1/Homo_sapiens/byDataType/Nucleosome_peaks_DANPOS/Homo_sapiens.hsNuc0390101.nucleosome.DANPOSPeak.bed.gz).