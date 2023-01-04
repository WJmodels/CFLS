# CDFLSystem
## Download

Model weights and datasets:

https://pan.baidu.com/s/1i9mB4iUxoaf-B-I5bmlpyQ?pwd=unm6

password: unm6

## CMGNet_code

python test_infer_molecular_formula_mul_pr_rerank_copy.py \
--weights_path ./CMGNet/weight_smiles_decoder/20220404_bart_3stage_long_1/epoch_199_loss_0.067575.pth \
--use_dict

## SampRNN

python transfer_learning.py \
--weight ./SampRNN/weight/pretain_5epoch/pretrain_5_epoch.ckpt
--smiles_file case.smi

## CReSS

