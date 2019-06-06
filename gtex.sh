qR.rb  -t "06:00:00" -a CWALLACE-SL2-CPU -p skylake-himem -j gtex -y 1-22 -n chr -r gtex.R

## catchup
for i in `seq 1 22`; do
    ./gtex.R --args chr=$i
done

## this will generate files to run GUESSFM, which are called now
ls /home/cew54/scratch/gtex-gfm/chr*/runme.sh

for f in /home/cew54/scratch/gtex-gfm/chr*/runme.sh; do
    qlines_asarray.rb -a CWALLACE-SL2-CPU -c 1 -t "04:00:00" -r $f
done

## collate output - run on queue to make intermediate files
qR.rb  -t "02:00:00" -a CWALLACE-SL2-CPU -c 6 -p skylake-himem -j gtex -y 1-22 -n chr -r gtex-collate.R
qR.rb  -t "02:00:00" -a CWALLACE-SL2-CPU -c 6 -p skylake-himem -j gtex -y 1-22 -n chr -r gtex-collate.R

## run interactively to generate plots
./gtex-collate.R


