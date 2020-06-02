#!/usr/bin/env bash
for a in epem_a_hhttx_no_h_*.yaml; do echo $a; ./select_numerator.py $1 $a; done;
