make; ./main -N ../../graphs/astro-ph.graph -R 1 -T 1 -K 16707 > astro_output.txt; sed -n 33416,50122p astro_output.txt > astro_output_static.txt ; sed -n 16710,33415p astro_output.txt > astro_output_dynamic.txt ; diff astro_output_static.txt astro_output_dynamic.txt ; cat astro_output.txt | grep tally

