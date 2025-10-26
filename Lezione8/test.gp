set terminal qt font "CMU Serif, 12" enhanced #size 650,412 background rgb "#151515"

exp_sin(x) = exp(sin(x-1))

plot "data/test_data.dat", exp_sin(x)

pause mouse close
