


AMP_NAME=$1
RUN_LTD="./../../../../rust_backend/target/release/ltd -f ./hyperparameters.yaml --amplitudes ./amplitudes_v2.yaml -l ./topologies_v2.yaml --amplitude ${AMP_NAME} integrated_ct"

eval "$RUN_LTD" \
    | grep -E "(UPlus|Uminus|UBarPlus|UBarMinus|VPlus|VMinus|VBarPlus|VBarMinus|APlus|AMinus)" \
    | sed -r "s/([0-9.+\-]+e[0-9\-][^+-]*)([0-9.+\-]+e[0-9\-][^+-]*)/[\1,\2],/g" \
    | sed -r "s/0.[0]+e0/0/g" \
    | tr -d '|i ' | sed 's/,$//g'\
    | sed -e "s/(p0.*=/    spinor_u:\n      - /g" \
          -e "s/(p1.*=/    spinor_vbar:\n      - /g" \
          -e "s/(p2.*=/    cpol:\n      - /g" \
          -e "s/(p.*=/      - /g" \

