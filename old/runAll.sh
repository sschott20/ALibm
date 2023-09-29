#!/bin/bash

./polygen ${BF16INTERVALPATH}/bf16LogIntervals ${TF32INTERVALPATH}/tf32LogIntervals ${FLOAT34INTERVALPATH}/Float34ROLogIntervals configs/log.config 4
./polygen ${BF16INTERVALPATH}/bf16Log2Intervals ${TF32INTERVALPATH}/tf32Log2Intervals ${FLOAT34INTERVALPATH}/Float34ROLog2Intervals configs/log2.config 1
./polygen ${BF16INTERVALPATH}/bf16Log10Intervals ${TF32INTERVALPATH}/tf32Log10Intervals ${FLOAT34INTERVALPATH}/Float34ROLog10Intervals configs/log10.config 4
./polygen ${BF16INTERVALPATH}/bf16ExpIntervals ${TF32INTERVALPATH}/tf32ExpIntervals ${FLOAT34INTERVALPATH}/Float34ROExpIntervals configs/exp.config 4
./polygen ${BF16INTERVALPATH}/bf16Exp2Intervals ${TF32INTERVALPATH}/tf32Exp2Intervals ${FLOAT34INTERVALPATH}/Float34ROExp2Intervals configs/exp2.config 1
./polygen ${BF16INTERVALPATH}/bf16Exp10Intervals ${TF32INTERVALPATH}/tf32Exp10Intervals ${FLOAT34INTERVALPATH}/Float34ROExp10Intervals configs/exp10.config 4
./polygen ${BF16INTERVALPATH}/bf16SinhForSinhIntervals ${TF32INTERVALPATH}/tf32SinhForSinhIntervals ${FLOAT34INTERVALPATH}/Float34ROSinhForSinhIntervals configs/sinhforsinh.config 1
./polygen ${BF16INTERVALPATH}/bf16CoshForSinhIntervals ${TF32INTERVALPATH}/tf32CoshForSinhIntervals ${FLOAT34INTERVALPATH}/Float34ROCoshForSinhIntervals configs/coshforsinh.config 1
./polygen ${BF16INTERVALPATH}/bf16SinhForCoshIntervals ${TF32INTERVALPATH}/tf32SinhForCoshIntervals ${FLOAT34INTERVALPATH}/Float34ROSinhForCoshIntervals configs/sinhforcosh.config 1
./polygen ${BF16INTERVALPATH}/bf16CoshForCoshIntervals ${TF32INTERVALPATH}/tf32CoshForCoshIntervals ${FLOAT34INTERVALPATH}/Float34ROCoshForCoshIntervals configs/coshforcosh.config 1
./polygen ${BF16INTERVALPATH}/bf16SinpiForSinpiIntervals ${TF32INTERVALPATH}/tf32SinpiForSinpiIntervals ${FLOAT34INTERVALPATH}/Float34ROSinpiForSinpiIntervals configs/sinpiforsinpi.config 1
./polygen ${BF16INTERVALPATH}/bf16CospiForSinpiIntervals ${TF32INTERVALPATH}/tf32CospiForSinpiIntervals ${FLOAT34INTERVALPATH}/Float34ROCospiForSinpiIntervals configs/cospiforsinpi.config 1
./polygen ${BF16INTERVALPATH}/bf16SinpiForCospiIntervals ${TF32INTERVALPATH}/tf32SinpiForCospiIntervals ${FLOAT34INTERVALPATH}/Float34ROSinpiForCospiIntervals configs/sinpiforcospi.config 1
./polygen ${BF16INTERVALPATH}/bf16CospiForCospiIntervals ${TF32INTERVALPATH}/tf32CospiForCospiIntervals ${FLOAT34INTERVALPATH}/Float34ROCospiForCospiIntervals configs/cospiforcospi.config 1
