#! /bin/bash
outlog=$1
start=`date +%s`
printf "Memory\t\tCPU\tTime\n"
printf "Memory\t\tCPU\tTime\n" > ${outlog}.log
while [ True  ]; do
end=`date +%s`
MEMORY=$(free -m | awk 'NR==2{printf "%.2f%%\t\t", $3*100/$2 }')
CPU=$(top -bn1 | grep load | awk '{printf "%.2f%%\t\t\n", $(NF-2)}')
RUNTIME=$((end-start))
echo "$MEMORY$CPU$RUNTIME" >> ${outlog}.log
echo "$MEMORY$CPU$RUNTIME"
sleep 1
done


