#!/bin/sh
touch stochTimeTests.txt
rm stochTimeTests.txt
touch stochTimeTests.txt
echo "Filename\tLength\tNormal\tminus10\tminus8\tminus6" >> stochTimeTests.txt
for((i = 100; i <= 4000;i=i+25))
do
    for((j = 1; j <= 8; j=j+1))
    do
	echo ${i}-${j}
	python random_seq.py $i > ./randseqs/${i}-${j}.seq;
	../src/time_test ./randseqs/${i}-${j}.seq >> stochTimeTests.txt;
    done
done
