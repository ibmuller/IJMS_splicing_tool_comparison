miso --run SE_indexmap/ data/mySTARdir/CEMWT_r30dm/CEMWTAligned.sortedByCoord.out.bam --read-len 101 --output-dir CEMWT_SE

index_gff --index data/gff/hg38_miso/SE.hg38.gff3 SE_indexmap
