cat data/myFLUXdir/stringtie_CEMWT_sorted.gtf | awk '{if($3 == "transcript"){print substr($18,2,(length($18)-3))*30}}' > CEMWT_expressions_reads.txt
