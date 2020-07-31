bamloc=$1
for file in ${bamloc}*.bam; do filename=`echo $file | cut -d "." -f 1`; samtools view -H $file | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > ${filename}_chr.bam; done
