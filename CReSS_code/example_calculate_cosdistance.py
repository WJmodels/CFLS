from infer import ModelInference
import json

config_path = "models/2_5_w_model/8.json"
pretrain_model_path = "models/2_5_w_model/8.pth"

model_inference = ModelInference(config_path=config_path,
                                 pretrain_model_path=pretrain_model_path,
                                 device="cuda:0")

# smiles_string = [
#     "C(=O)C1(O)C(C=O)=CCC2C(C)(C)CCCC21C", #smiles1
#      "C1(C)(C)C=[N+]([O-])C(C)(C)N1C" #smiles2
# ]

# nmr_list = [
#         [17.6, 18.3, 22.6, 26.5, 31.7, 33.5, 41.8, 42.0, 
#         42.2, 78.34, 140.99, 158.3, 193.4, 203.0], #nmr1
#         [23.3, 23.5, 26.1, 60.5, 90.0, 132.1] #nmr2
#         ]

# smiles_feature = model_inference.smiles_encode(smiles_string)
# nmr_feature = model_inference.nmr_encode(nmr_list)
# print(model_inference.get_cos_distance(smiles_feature, nmr_feature))

lis = []
with open('case.json', 'r') as fi:
    contents = json.load(fi)
    
    i = 0
   
    for content in contents:
        
        i += 1
        print(i,"/",len(contents))
    
        smiles_string = content['result']       
        nmr_list = content['nmr'] 

        smiles_feature = model_inference.smiles_encode(smiles_string)
        nmr_feature = model_inference.nmr_encode(nmr_list)
        cos_distance_gpu = model_inference.get_cos_distance(smiles_feature, nmr_feature)
        cos_distance_cpu = cos_distance_gpu.cpu()
        cos_distance = cos_distance_cpu.detach().numpy().tolist()

        dic = {}
        dic['index'] = content['index']
        dic['cos_distance'] = cos_distance
        lis.append(dic)
        

with open('case_cos.json', 'a') as fo:
    json.dump(lis, fo, indent = 2)