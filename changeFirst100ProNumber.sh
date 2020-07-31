pro=$1
enst=$2
out=$3
cat ${pro} | awk 'BEGIN{count = 0} {if(count <= 100 && $5!=0 && $6!=0){$5=($5*100);$6=($6*100);count++}print $0}'  OFS='\t' > ${out} 
