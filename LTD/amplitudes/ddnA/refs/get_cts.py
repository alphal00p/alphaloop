#!/bin/bash

echo "name,ct_re,ct_im"

for ID in {4..8}; do
    AMP_NAME=dd${ID}A
    RUN_LTD="~/BitBucket/alphaloop/rust_backend/target/release/ltd -f ./hyperparameters.yaml --amplitudes ./amplitudes_v2.yaml -l ./topologies_v2.yaml --amplitude ${AMP_NAME} integrated_ct"
    
    eval "$RUN_LTD" | grep "Finite" | sed 's/^.*://g' | tr -d ' i' | sed -r "s/([0-9.+\-]+e[0-9\-][^+-]*)([0-9.+\-]+e[0-9\-][^+-]*)/${AMP_NAME},\1,\2/g" \

done
