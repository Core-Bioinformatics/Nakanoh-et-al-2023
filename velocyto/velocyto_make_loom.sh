#!/bin/bash

refs="/path/to/genomes/CellRanger_references/"
gtf="${refs}/refdata-cellranger-GRCh38-5.0.0/genes/genes.gtf"

refpath="/path/to/3.cellranger"
path0="${refpath}/SLX-21814_SITTA3_H7T7CDRX2/"
path1="${refpath}/SLX-21814_SITTB3_H7T7CDRX2/"
path2="${refpath}/SLX-21814_SITTC3_H7T7CDRX2/"
path3="${refpath}/SLX-21814_SITTD3_H7T7CDRX2/"

velocyto run10x  ${path0} ${gtf} &
velocyto run10x  ${path1} ${gtf} &
velocyto run10x  ${path2} ${gtf} &
velocyto run10x  ${path3} ${gtf} &
