#!/bin/bash

#Plot real and imaginary amplitudes nomralized by the analytic result
rm sample_plot.eps;  ./cuhre.gnu; mv sample_plot.eps plots/ddAAA_cuhre.eps
rm sample_plot.eps;  ./vegas.gnu; mv sample_plot.eps plots/ddAAA_vegas.eps

#Plot the collinears, soft and uv limit
TITLE_NAMES=( "\$p_1\$ Collinear Limit" "\$p_2\$ Collinear Limit" "Soft Limit" "UV Limit" )
LIMIT_FILES=( "explore_Collinear_p1.csv" "explore_Collinear_p2.csv" "explore_Soft1.csv" "explore_UV.csv" )
OUTPUT_FILES=( "ddAAA_collp1.eps" "ddAAA_collp2.eps" "ddAAA_soft.eps" "ddAAA_uv.eps" )
echo ${TITLE_NAMES[@]}
for i in {0..3}; do
	sed -e "s/TITLE_NAME/${TITLE_NAME[$i]}/g" -e "s/LIMIT_FILE/${LIMIT_FILES[$i]}/g" limit.gnu > tmp.gnu
	gnuplot tmp.gnu; mv sample_plot.eps plots/${OUTPUT_FILES[$i]}
	
done


