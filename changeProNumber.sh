#!/bin/bash
pro=$1
expression=$2
echo "" > test.pro
while read -u 3 -r file1 && read -u 4 -r file2; do
	echo "${file1}" | awk -v file2=${file2} '{$6=(file2);printf "%s\t%s\t%s\t%i\t%f\t%.0f\n", $1, $2, $3, $4, $5, $6}' >> test.pro
done 3<$1 4<$2
