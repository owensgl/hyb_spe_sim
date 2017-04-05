#Run using cat simulation_stats.txt | parallel -j 10 --colsep ' ' bash ./full_sim_pipeline.sh {1} {2} {3} {4}
n_reps=10
pop_size=$1
generations=$2
n_loci=$3
selection=$4
bin=/home/owens/bin/hyb_spe_sim
for rep in `seq 1 $n_reps`;
do
perl $bin/make_f1s.pl $pop_size | perl $bin/generations_withselection.pl $generations $n_loci $selection | \
grep -v "#" | head -n 35 | tee >(perl $bin/measure_true_ancestry.pl > $bin/simulations/EST.XRQ.GATK.simtrueancestry_${pop_size}_${generations}_${n_loci}_${selection}_${rep}.txt) | perl $bin/simulate_genotypes.pl /home/owens/RNAseq/EST.XRQ.GATK.noref.ctab /home/owens/RNAseq/Helianthus.hybrid.poplist.txt | perl /home/owens/bin/reformat/SNPtable2hmp.pl | perl /home/owens/bin/pop_gen/hmp2hyblik_v2.pl /home/owens/RNAseq/Helianthus.hybrid.poplist.forsims.txt | perl /home/owens/bin/pop_gen/hybcurve_summarizer_bywindow_nocm.pl > $bin/simulations/EST.XRQ.GATK.simancestry_${pop_size}_${generations}_${n_loci}_${selection}_${rep}.txt
#grep -v "#" | head -n 35  | perl $bin/measure_true_ancestry.pl > tmp.trueancestry.txt

done
