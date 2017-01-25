#!/bin/bash

for n in $1/tdc*.root
do
	./Proto $n
done
