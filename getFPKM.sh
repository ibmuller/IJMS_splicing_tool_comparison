cat stringtie_CEMWT_sorted.gtf | awk '{if($3=="transcript"){print substr($12,2,15) "\t" substr($18,2,(length($18)-3))}}' >  sorted_FPKM.txt

