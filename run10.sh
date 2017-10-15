#!/bin/sh

for i in `seq -w 0 19`
do
	./build/Scheduler STG/50/rand00${i}.stg
done
