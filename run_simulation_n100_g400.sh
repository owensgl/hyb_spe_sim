number=$1
perl ./make_f1s.pl 100 | perl generations.pl 400 >> simulated_genotypes_n100_g400_r100.${number}.txt
