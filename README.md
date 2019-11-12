# coloc-cond-mask

Code to run simulations for coloc-prior-mask project.

## Simulations for p12 variability

SPECIAL corresponds to hypothesis, H0...H4

```{sh}
for ld in lowld highld; do
	for NSNP in 200 400 700; do
		for n in 100 200 500 1000 2000 5000 10000; do
			qR.rb -r -y 0-4 -n SPECIAL ./coloc-vary-p12-sim.R --args NSIM=200 N=$n NSNP=$NSNP -j p12-$n ld=$ld;
		done
	done
done

## wait...

./coloc-vary-p12-collate-v2.R
```

## Simulations to compare conditioning, masking, varying r2thr

Files: $COLOCINDEP/cvaryr2-v2*

SPECIAL depends on NCV.

**NCV=4: simulate two causal variants per trait**
5. share 1, weakest effect for one, strongest for other
4. share 2, opposite effects
3. share 1, equal effect (weakest) + indep each
2. share 2, equal effects
1. share 1, equal effect (strongest), + indep each
0. share 0

**NCV=3: simulate one causal variant for one trait, two for other**
3. share 1, stronger effect
2. share 1, weaker effect
1. share 1, equal effect (strongest), + indep each
0. share 0

for ld in lowld highld; do
```{sh}
for N in 1000 2000; do
for ld in highld lowld; do
for NSNP in 1000; do
    for ncv in 3 4; do
        for s in 0 1 2 3 4 5; do 
        qR.rb -c 2 -a CWALLACE-SL2-CPU -p skylake-himem -y 1-20 \
		-j cm$s-$ncv-$N-$ld -r \
		./coloc-mask-vary-r2-sim.R --args \
		NSIM=50 SPECIAL=$s NCV=$ncv NSNP=$NSNP N=$N ld=$ld; 
        done
    done
done
done
done

## wait...

./coloc-mask-vary-r2-collate-v2.R
```
done

## GTeX analysis - estimating the proportion of SNPs which are eQTL causal variants

qR.rb -r -t "04:00:00" ./gtex-finemap.R

