file=$1
zcat simulations/${file}.txt.gz | perl format_simulations_for_correlations.pl  > simulations/${file}.summarized.txt
