#!/bin/bash

[ -f "$1" ] || exit 1;
[ -z "$2" ] && export DIR="./"
[ -n "$2" ] && export DIR="$2"

cat "$1" | grep "JetPt\|Q$" | sed 's/Q$//g'|sed 's/}/_quark&/g' > $DIR/Quark.txt
cat "$1" | grep "JetPt\|G$" | sed 's/G$//g'|sed 's/}/_gluon&/g' > $DIR/Gluon.txt


NUM=$( (cat $DIR/Quark.txt | wc -l | tr -d '\n'; echo "/3-1") | bc )
echo "[quark]"> $DIR/nCharged.txt
cat $DIR/Quark.txt | grep -A $NUM nCharged  >> $DIR/nCharged.txt
echo "[quark]"> $DIR/nNeutral.txt
cat $DIR/Quark.txt | grep -A $NUM nNeutral  >> $DIR/nNeutral.txt
echo "[quark]"> $DIR/ptD.txt
cat $DIR/Quark.txt | grep -A $NUM ptD >> $DIR/ptD.txt
echo "[gluon]">> $DIR/nCharged.txt
cat $DIR/Gluon.txt | grep -A $NUM nCharged  >> $DIR/nCharged.txt
echo "[gluon]">> $DIR/nNeutral.txt
cat $DIR/Gluon.txt | grep -A $NUM nNeutral  >> $DIR/nNeutral.txt
echo "[gluon]">> $DIR/ptD.txt
cat $DIR/Gluon.txt | grep -A $NUM ptD  >> $DIR/ptD.txt
