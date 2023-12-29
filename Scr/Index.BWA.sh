#!/bin/bash

RefDir=$1"/Ref"

for Ind in $RefDir/*fa
do
  bwa index $Ind &
done 

wait
