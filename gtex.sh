qR.rb  -t "06:00:00" -a CWALLACE-SL2-CPU -p skylake-himem -j gtex -y 1-22 -n chr -r gtex.R

## catchup
for i in `seq 1 22`; do
    ./gtex.R --args chr=$i
done


for i in `seq 1 20`; do
    # qlines.rb -c 1 -t "02:00:00" -r  /home/cew54/scratch/gtex-gfm/chr${i}/runme.sh
    qlines_asarray.rb -c 1 -t "04:00:00" -r  /home/cew54/scratch/gtex-gfm/chr${i}/runme.sh
done
