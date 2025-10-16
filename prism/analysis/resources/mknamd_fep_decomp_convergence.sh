#!/bin/bash

if [[ "$#" -eq 0 ]]; then
    echo ">Usage: mknamd_fep_decomp_convergence alchemy.fepout [-multiple (Optional)]"
    echo ">If you don't have fep.namd in this folder"
    echo ">Usage: mknamd_fep_decomp_convergence alchemy.fepout [alchEquilSteps] [numSteps] [alchOutFreq] [numPoints] [-multiple (Optional)]"
    exit 0
fi
# bash mknamd_fep_decomp_convergence.sh *.fepout 10000 510000 500 10
# bash $mfdc *.fepout 20000 510000 5 10
# bash $mfdc *.fepout 10000 510000 5 10

numPoints=$5

if [[ $2 == "-multiple" ]]; then
    tag_multiple=true
else
    tag_multiple=false
fi

if $tag_multiple; then

    triallist=$(find . -maxdepth 1 -name 'trial[0-9]*' | sort -V)
    for trial in $triallist; do
        [ ! -f "$trial/$1" ] && {
            echo -e "$trial/$1 does not exist!"
            exit 1
        }
        echo "Computing for $trial"
        outdir=$trial/outdecomp
        if [ -d "$outdir" ]; then
            rm -r $outdir
        fi
        mkdir $outdir
        grep "^#Free energy change" $trial/$1 | awk '{print NR, $9, $12, $19}' >$outdir/raw
        grep -n "^#STARTING COLLECTION" $trial/$1 | awk -F: '{print $1+1}' >$outdir/l1.txt
        grep -n "^#Free energy change" $trial/$1 | awk -F: '{print $1}' >$outdir/l2.txt
        paste $outdir/l1.txt $outdir/l2.txt >$outdir/l3.txt
        awk '{print NR,$1,($2-$1)}' $outdir/l3.txt >$outdir/nlist.txt
        while read -r nn ntail nhead; do
            tail -n +$ntail $trial/$1 | head -$((nhead + 1)) >> $outdir/raw$nn.txt
            tail -n +$ntail $trial/$1 | head -$((nhead)) | awk 'NR%1==0' | awk '{print ($4 - $3), ($6 - $5), $7, $9}' >> $outdir/data$nn.txt
        done <$outdir/nlist.txt
    done

else

    triallist="1"
    [ ! -f "$1" ] && {
        echo -e "$1 does not exist!"
        exit 1
    } # $1: filename
    echo 'Start processing'
    outdir=outdecomp
    outfile=decompose_summary
    if [ -d "$outdir" ]; then
        rm -r $outdir
    fi
    mkdir $outdir
    grep "^#Free energy change" $1 | awk '{print NR, $9, $12, $19}' >$outdir/raw.txt # λ2, dG, net change
    grep -n "^#STARTING COLLECTION" $1 | awk -F: '{print $1+1}' >$outdir/l1.txt      # line number plus 1. one number per line. a window per line
    grep -n "^#Free energy change" $1 | awk -F: '{print $1}' >$outdir/l2.txt         # line number. the range of collection
    paste $outdir/l1.txt $outdir/l2.txt >$outdir/l3.txt                              # it's like zip(list1, list2) in Python
    awk '{print NR,$1,($2-$1)}' $outdir/l3.txt >$outdir/nlist.txt                    # new line number, start, # of lines, respectively
    while read -r nn ntail nhead; do                                                 # from nlist.
        tail -n +$ntail $1 | head -$((nhead + 1)) >> $outdir/raw$nn.txt               # each window, all "start collection" data and the "#Free energy change" line
        tail -n +$ntail $1 | head -$((nhead)) | awk 'NR%1==0' | awk '{print ($4 - $3), ($6 - $5), $7, $9}' >> $outdir/data$nn.txt
        # from raw$nn.txt. ΔElec, ΔVdW, dE, temperature. awk 'NR%1==0' means excluding line 0
        echo finished extracting data for window $nn
    done <$outdir/nlist.txt
fi

nhead=$(head -n 1 $outdir/nlist.txt | awk '{print $3}')
if [ -f "$outfile.dat" ]; then
    rm $outfile.dat
fi

printf "%12s,%11s,%11s,%11s,%11s\n" fraction sum_dg sum_delec sum_dvdw sum_couple >$outfile.dat # table head

for trial in $triallist; do

    if $tag_multiple; then
        outdir=$trial/outdecomp
    else
        outdir=outdecomp
    fi

    for i in $(seq 1 $numPoints); do
        rm $outdir/decomp.txt $outdir/error.txt                              # the next timestep does not need that for the last timestep
        lastline=$((i * (nhead - 1) / numPoints + (i == numPoints) * ((nhead - 1) % numPoints))) # last line to extract, the cumulative time
        while read -r nn v1 v2; do
            # each window. nn, v1, v2, which are from nlist
            # $1~4: from data$nn.txt.
            # k_b=-0.01987 K^{-1}=R_kJmol*kJ2kcal=8.314 J/mol/K / 1000 kJ/mol * 1/4.182
            # extract time series data. get all data at the last point.
            head -n $lastline $outdir/data$nn.txt >$outdir/data_temp.txt
            awk '{sum_temp += $4; \
            esum_elec += (exp($1 / -0.001987 / $4)); \
            esum_vdw += (exp($2  / -0.001987 / $4)); \
            esum_total += (exp($3 / -0.001987 / $4)); \
            sum_elec += $1; sum2_elec += $1*$1; \
            sum_vdw += $2; sum2_vdw += $2*$2; \
            sum_couple += $1*$2; n++;} END \
            {printf "%3d %5.2f % 6.2f % 10.2f % 6.2f % 6.2f % 6.2f % 6.2f % 10.2f % 6.2f % 6.2f % 6.2f % 6.2f % 10.2f\n", \
                n, \
                sum_temp/n, \
                (-0.001987*sum_temp/n*log(esum_total/n)), \
                (-0.001987*sum_temp/n*log(esum_elec/n)), \
                (sum_elec/n) - (sum2_elec/n - sum_elec*sum_elec/n/n)/1.24 - (sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                (sum_elec/n), \
                -(sum2_elec/n - sum_elec*sum_elec/n/n)/1.24, \
                -(sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                (-0.001987*sum_temp/n*log(esum_vdw/n)), \
                (sum_vdw/n) - (sum2_vdw/n - sum_vdw*sum_vdw/n/n)/1.24 - (sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                (sum_vdw/n), \
                -(sum2_vdw/n - sum_vdw*sum_vdw/n/n)/1.24, \
                -(sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
                -0.001987*sum_temp/n*log(esum_total/n) + 0.001987*sum_temp/n*log(esum_elec/n) + 0.001987*sum_temp/n*log(esum_vdw/n) \
                }' <$outdir/data_temp.txt >>$outdir/decomp.txt # one line per window. I don't know what that is but it's an important intermediate file
            awk '{sum_temp += $4; \
            sum1 += (exp($3 / -0.001987 / $4)); \
            sum2 += (exp(2 * $3 / -0.001987 / $4)); \
            sum_total += $3; sum2_total += $3*$3; n++;} END \
            {printf "%3d %5.2f % 6.2f %6.2f % 6.2f %6.2f %6.2f %7.0f %3.2e %3.2e %3.2e\n", \
                n, \
                sum_temp/n, \
                (-0.001987*sum_temp/n*log(sum1/n)), \
                sqrt(((n*sum2/sum1/sum1-1)/0.3844)), \
                (sum_total/n), \
                sqrt(sum2_total/n - sum_total*sum_total/n/n), \
                (sum_total/n + 0.001987*sum_temp/n*log(sum1/n)), \
                (exp((sum_total/n - (-0.001987*sum_temp/n*log(sum1/n)))*n/sum_temp/0.001987)), \
                (sum1/n), \
                (sum1*sum1/n/n), \
                (sum2/n)}' <$outdir/data_temp.txt >>$outdir/error.txt
            # echo finished processing window $nn
        done <$outdir/nlist.txt
        echo finished processing timestep $i/$numPoints

        # print to screen
        awk '{sum_dg += $3; \
            sum_delec += $4; \
            sum_dvdw += $9; \
            sum_couple += $14; n++;} END \
            {printf "mknamd> Decomp gives you: %.2f\n\telec   % 7.2f\n\tvdw    % 7.2f\n\tcoupl  % 7.2f\n",
            sum_dg, sum_delec, sum_dvdw, sum_couple}' <$outdir/decomp.txt
        # output each window
        # awk '{printf "mknamd> Decomp gives you: %.4f\n\telec   % 7.4f\n\tvdw    % 7.4f\n\tcoupl  % 7.4f\n",
        #    $3,  $4,  $9,  $14}' < $outdir/decomp.txt

        # record data
        if $tag_multiple; then
            awk -v "prefix=${trial#./trial}" -v "fraction=$((10 ** 3 * lastline / (nhead - 1)))e-3" \
                '{sum_dg += $3; \
                sum_delec += $4; \
                sum_dvdw += $9; \
                sum_couple += $14; n++;} END \
                {printf "%12d,%11.2f,%11.2f,%11.2f,%11.3f\n",
                prefix,sum_dg,sum_delec,sum_dvdw,sum_couple,fraction}' <$outdir/decomp.txt >>${outfile}.dat
            #TODO: not sure how to put the result of multiple
        else
            awk -v "fraction=$((10 ** 3 * lastline / (nhead - 1)))e-3" \
                '{sum_dg += $3; \
                sum_delec += $4; \
                sum_dvdw += $9; \
                sum_couple += $14; n++;} END \
                {printf "%12.3f,%11.2f,%11.2f,%11.2f,%11.2f\n",
                fraction,sum_dg,sum_delec,sum_dvdw,sum_couple}' <$outdir/decomp.txt >>${outfile}.dat
        fi

    done

    if $tag_multiple; then # total dG
        echo "mknamd> <$trial>"
        echo -n "mknamd>    FEP gives you: "
        tail -1 $trial/$1 | awk '{print $19}'
    else
        echo -n "mknamd>    FEP gives you: "
        tail -1 $1 | awk '{print $19}'
    fi

done

# Clean up intermediate outdecomp folder to save disk space
# The summary data has been saved to decompose_summary.dat
if $tag_multiple; then
    for trial in $triallist; do
        outdir=$trial/outdecomp
        if [ -d "$outdir" ]; then
            echo "Cleaning up $outdir..."
            rm -r $outdir
        fi
    done
else
    outdir=outdecomp
    if [ -d "$outdir" ]; then
        echo "Cleaning up $outdir..."
        rm -r $outdir
    fi
fi

echo "✓ Processing complete. Summary saved to $outfile.dat"

#http://membrane.urmc.rochester.edu/?page_id=126
#https://en.wikipedia.org/wiki/Bootstrapping_(statistics)
