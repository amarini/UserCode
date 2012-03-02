#!/bin/bash

[ -f "$1" ] || exit 1;

cat "$1" | grep "JetPt\|Q$" | sed 's/Q$//g'|sed 's/}/Quark&/g' > Quark.txt
cat "$1" | grep "JetPt\|G$" | sed 's/G$//g'|sed 's/}/Gluon&/g' > Gluon.txt


NUM=$( (cat Quark.txt | wc -l | tr -d '\n'; echo "/3-1") | bc )
echo "[quark]"> nCharged.txt
cat Quark.txt | grep -A $NUM nCharged  >> nCharged.txt
echo "[quark]"> nNeutral.txt
cat Quark.txt | grep -A $NUM nNeutral  >> nNeutral.txt
echo "[quark]"> ptD.txt
cat Quark.txt | grep -A $NUM ptD >> ptD.txt
echo "[gluon]">> nCharged.txt
cat Gluon.txt | grep -A $NUM nCharged  >> nCharged.txt
echo "[gluon]">> nNeutral.txt
cat Gluon.txt | grep -A $NUM nNeutral  >> nNeutral.txt
echo "[gluon]">> ptD.txt
cat Gluon.txt | grep -A $NUM ptD  >> ptD.txt
