# coloc-cond-mask

simulate two causal variants per trait
 simulations options
 5 = share 1, weakest effect for one, strongest for other
 4 = share 2, opposite effects
 3 = share 1, equal effect (weakest) + indep each
 2 = share 2, equal effects
 1 = share 1, equal effect (strongest), + indep each
 0 = share 0
./coloc-cond-vs-indep-sim.R --args NSIM=100 SPECIAL=0


quick check
```{sh}
qR.rb -r -y 0-5 -n SPECIAL ./coloc-cond-vs-indep-sim.R --args NSIM=100
```

run
```{sh}
for s in 0 1 2 3 4 5; do qR.rb -r -y 1-50 ./coloc-cond-vs-indep-sim.R --args NSIM=200 SPECIAL=$s NCV=4; done
for s in 0 1 2 3; do qR.rb -r -y 1-50 ./coloc-cond-vs-indep-sim.R --args NSIM=200 SPECIAL=$s NCV=3; done
```
```{sh}
for s in 0 1 2 3 4 5; do qR.rb -r -y 1-10 ./coloc-cond-vs-indep-sim-v2.R --args NSIM=100 SPECIAL=$s NCV=4; done
for s in 0 1 2; do qR.rb -r -y 1-10 ./coloc-cond-vs-indep-sim-v2.R --args NSIM=100 SPECIAL=$s NCV=3; done
```


plot results
./coloc-cond-vs-indep-collate.R
./coloc-cond-vs-indep-collate-v2.R



simulate one causal variant for one trait, two for other 
simulations options
 3 = share 1, stronger effect
 2 = share 1, weaker effect
 1 = share 1, equal effect (strongest), + indep each
 0 = share 0

```{sh}
 for s in 0 1 2 3; do qR.rb -r -y 1-50 ./coloc-cond-vs-indep-onetwo-sim.R --args SPECIAL=$s; done
```

plot results
./coloc-cond-vs-indep-onetwo-collate.R


## Simulations for p12 variability

```{sh}
for ld in lowld highld; do
	for NSNP in 200 400 700; do
		for n in 100 200 500 1000 2000 5000 10000; do
			qR.rb -r -y 0-5 -n SPECIAL ./coloc-vary-p12-sim.R --args NSIM=200 N=$n NSNP=$NSNP -j p12-$n ld=$ld;
		done
	done
done
```

## Simulations for varying r2thr
```{sh}
for ld in lowld highld; do
	for ncv in 3 4; do
		for s in 0 1 2 3 4 5; do 
			qR.rb -r -y 1-100 ./coloc-mask-vary-r2-sim.R --args NSIM=20 SPECIAL=$s NCV=$ncv ld=$ld; 
		done
	done
done
```

## Model averaging - abandonned
```{sh}
 for s in 0 1 2 3 4; do
 for n in 500 1000 2000; do qR.rb -r -y 1-10 ./coloc-avg-prior-sim.R --args SPECIAL=$s NSIM=1000 N=$n -j avg$s-$n;
 done
 done
```
