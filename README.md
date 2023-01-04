# CDFLSystem
## Download

Model weights and datasets:

https://pan.baidu.com/s/1i9mB4iUxoaf-B-I5bmlpyQ?pwd=unm6

password: unm6

## CMGNet_code

python test_infer_molecular_formula_mul_pr_rerank_copy.py \
--weights_path ./CMGNet_code/weight_smiles_decoder/20220404_bart_3stage_long_1/epoch_199_loss_0.067575.pth \
--use_dict

## SampRNN_code

python transfer_learning.py \
--weight ./SampRNN_code/weight/pretain_5epoch/pretrain_5_epoch.ckpt \
--smiles_file ./SampRNN_code/case.smi

python sample.py \
--filename ./SampRNN_code/weight/case.ckpt

## CReSS_code

python example_search_library_list.py \
--pretrain_model_path ./CReSS_code/models/2_5_w_model/8.pth \
--json_list ./CReSS_code/RNN_G+C/G+C_lib/case_lib.json

python example_calculate_cosdistance.py
