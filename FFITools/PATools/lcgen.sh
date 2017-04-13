#!/bin/bash
INLIST=fits.ls
LC=LC/
for frame in $(cat $INLIST);
do 
    cat FITS/$frame.fiphot;
done | grcollect - --col-base 1 --prefix $LC --extension rlc
