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
for s in 0 1 2 3 4 5; do qR.rb -r -y 1-50 ./coloc-cond-vs-indep-sim.R --args NSIM=200 SPECIAL=$s; done
```

plot results
./coloc-cond-vs-indep-collate.R



simulate one causal variant for one trait, two for other 
simulations options
 3 = share 1, stronger effect
 2 = share 1, weaker effect
 1 = share 1, equal effect (strongest), + indep each
 0 = share 0

: for s in 0 1 2 3; do qR.rb -r -y 1-50 ./coloc-cond-vs-indep-onetwo-sim.R --args SPECIAL=$s; done

plot results
./coloc-cond-vs-indep-onetwo-collate.R


