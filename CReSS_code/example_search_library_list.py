import torch
import json
from tqdm import tqdm
from infer import ModelInference


def get_feature_from_json(json_list,
                          save_name,
                          model_inference,
                          n=256,
                          type="smiles",
                          flag_get_value=False):
    context = []
    print("start parse json")
# import training set from json_list
    for file in json_list:
        with open(file, "r") as f:
            context_tmp = json.loads(f.read())
            context_tmp = [
                i[type][0] for i in tqdm(context_tmp) if len(i[type]) > 0
            ]
        context += context_tmp
#    print("#######", context)
    print("Size of the library: ", len(context))
    # if flag_get_value == "only":
    #     return context
    # if type == "nmr":
    #     fn = model_inference.nmr_encode
    if type == "smiles":
        fn = model_inference.smiles_encode
    contexts = []
    print("start load batch")
    for i in range(0, len(context), n):
        contexts.append(context[i:i + n])
    print("start encode batch")
    result = [fn(i).to("cuda:0") for i in tqdm(contexts)]
    result = torch.cat(result, 0)
    if flag_get_value is True:
        if save_name is not None:
            torch.save((result, context), save_name)
        return result, context

    if save_name is not None:
        torch.save(result, save_name)
    return result


def get_topK_result(nmr_feature, smiles_feature, topK):
    indices = []
    scores = []
    with torch.no_grad():
        for i in tqdm(nmr_feature):
            nmr_smiles_distances_tmp = (
                i.unsqueeze(0) @ smiles_feature.t()).to("cuda:0")
        scores_, indices_ = nmr_smiles_distances_tmp.topk(topK,
                                                          dim=1,
                                                          largest=True,
                                                          sorted=True)
        indices.append(indices_)
        scores.append(scores_)
    indices = torch.cat(indices, 0)
    scores = torch.cat(scores, 0)
    return indices, scores


if __name__ == "__main__":
    # Load the model
    config_path = "models/2_5_w_model/8.json"
    pretrain_model_path = "models/2_5_w_model/8.pth"
    model_inference = ModelInference(config_path=config_path,
                                     pretrain_model_path=pretrain_model_path,
                                     device="cuda:0")

    res_out = []


    # json_list=["data/40W_lib/test_set_notex_lib_6471.json", "data/40W_lib/all_training_set_lib_370898.json", "data/40W_lib/test_set_ex_lib_41494.json"]
    json_list=["RNN_G+C/G+C_lib/5493_RNN_1000mol_lib.json"]
    # json_list=["data/190W_lib/Chembl_190W_lib_1829125.json"]
    
    smiles_feature, smiles_list = get_feature_from_json(
        json_list=json_list,
        model_inference=model_inference,
        n=64,
        save_name=None,
        type="smiles",
        flag_get_value=True)
    print("json_list_feature_over")

    with open("RNN_G+C/case_input/case.json", "r") as f:
        n = 0
        for line in f.readlines():
            n += 1
            print(n)
            content = json.loads(line)
            smiles = content['smiles']       
            nmr_list = content['nmr'] 

            # Extract NMR spectral feature vector 
            nmr_feature = model_inference.nmr_encode(nmr_list)

            # Construct a reference library by extracting structural feature vectors from SMILES strings
            # This might take a long time
            


            # Get top100 candidates by searching library 
            indices, scores = get_topK_result(nmr_feature, smiles_feature, 100)
        
            
            res_all = []
            # Print the result
            for (sco, idx) in zip(scores, indices):
                for ii, i in enumerate(idx): 
                    si_item = {}
                    si_item["top"] = ii
                    si_item["scores"] = sco[ii].item()
                    si_item["smiles"] = smiles_list[i]
                    res_all.append(si_item)
            res = {}
            res["smiles"] = smiles
            res["nmr"] = nmr_list
            res["result"] = res_all
            res_out.append(res)

    with open('RNN_G+C/case_output/case_out.json', 'a') as fo:
        json.dump(res_out, fo, indent = 2)