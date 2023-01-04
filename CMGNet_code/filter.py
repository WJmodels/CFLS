from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import rdMolDraw2D
# from rdkit import Chem as chem_tmp
from rdkit.Chem import Descriptors
from collections import Counter
import json
import signal

def set_timeout(num, callback):
    def wrap(func):
        def handle(signum, frame):
            raise RuntimeError

        def to_do(*args, **kwargs):
            try:
                signal.signal(signal.SIGALRM, handle)
                signal.alarm(num)
                print('start alarm signal.')
                r = func(*args, **kwargs)
                print('close alarm signal.')
                signal.alarm(0)
                return r
            except RuntimeError as e:
                callback()
        return to_do
    return wrap


def after_timeout():
    return 1

periodic_table_of_elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
                              "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

def get_molecular_formula(mol_):
    try:
        if mol_ is None:
            return ""
        mol = Chem.AddHs(mol_)
        dict_ = dict(Counter(atom.GetSymbol() for atom in mol.GetAtoms()))
        str_ = ""
        for i in periodic_table_of_elements:
            value = dict_.pop(i) if i in dict_.keys() else None
            if value is not None:
                str_ = str_ + i + str(value)
        if dict_:
            for k, v in dict_.items():
                str_ = str_ + k + str(v)
            print(dict_)
            # raise
        return str_
    except:
        return ""

@set_timeout(10, after_timeout)
def judge_unqualified_a(smiles, molecular_formula, fragment):
    a = 1    

    mol = Chem.MolFromSmiles(smiles)

    if use_dict["molecular_formula"] is True:
        if molecular_formula == "":
            molecular_formula = None
        if molecular_formula is not None and get_molecular_formula(mol) != molecular_formula:
            a = False
        print ('a0')
        print (a)

#     if use_dict["molecular_weight"] is True:
#         if molecular_weight == -1:
#             molecular_weight = None
#         if molecular_weight is not None and abs(Descriptors.ExactMolWt(mol) - molecular_weight) > molecular_weight_range:
#             return False
    if a == 1:
        if use_dict["fragment"] is True:
            if isinstance(fragment, str):
                fragment_ = Chem.MolFromSmiles(fragment)
            else:
                fragment_ = fragment
            print ('a1')
            if fragment_ is not None and not mol.HasSubstructMatch(fragment_):
                a = False
            print ('a2')        
    else:
        None
#     if mol_weight is None:
#         return False
    return a

def judge_unqualified_b(smiles):
    b = 1    
    try:
        mol = Chem.MolFromSmiles(smiles)
        print ('b0')
        mol_weight = Descriptors.MolWt(mol)
        print ('b1')
    except:
        b = False
        print ('b1')

#     if mol_weight is None:
#         return False

            

    return b
    
def judge_unqualified(smiles, molecular_formula, fragment):
    if judge_unqualified_b(smiles) == False:
        return False
    else:
        return judge_unqualified_a(smiles, molecular_formula, fragment)

def removeUnqualified(pi, po):

    
    # piQ = 'yz_data/realc_v.json'
    # fiQ = open(piQ, 'r')
    # liQ = fiQ.readlines()
    piQ = 'ACD40w_v/ACD40w_v.json'
    fiQ = open(piQ, 'r')
    liQ = fiQ.readlines()



    fo = open(po, 'a')
    fi = open(pi, 'r', encoding='utf8')


    context = json.load(fi)
    for m in range(10058, len(context)):
    # for m in range(10058, 10059):
        print (m)
        contex = context[m]
    #     contex = context[10057]
        result = contex['result']
        smiQuery = contex['smiles'][0]
    #     print (smiQuery)
        i = contex['index']
        Query = json.loads(liQ[i])
        resultb = []
        for k in range(0, len(result)):
            print (k)
            smiles = result[k] 
            print (smiles)
            if len(Query['fragments']) > 0 and contex['fragment_idx'] != None:
    #             fragment = Query['fragments'][contex['fragment_idx']]
                fragment = Query['fragments'][contex['fragment_idx']].replace('*', '')
                print (fragment)
            else:
                fragment = None
            molecular_formula = Query['molecular_formula']
    #         for smiles in result:
    #             print (smiles)
    #             if smiles == 'C(C(=O)OC)=C(c1c[nH]c2ccccc12)c1c[nH]c1':
    #                 print ('Find it')
    #             else:
    #                 None

            if judge_unqualified(smiles, molecular_formula, fragment) == False:
                None
    #                 print ('False')
    #                 j = result.index(smiles)
    #                 print (j)
    #                 print (result[j])
    #                 del result[j]
            else:
                resultb.append(smiles)

        dic = {}
        dic['smiles'] = contex['smiles']
        dic['nmr'] = contex['nmr']
        dic['index'] = contex['index']
        dic['result'] = resultb
#         print (dic['result'])
        json.dump(dic, fo)
        fo.write('\n')

    fo.close()
    fi.close()
    # fiQ.close()
    fiQ.close()

if __name__ == "__main__":
    # from yz_getRank_byTanimo import getRank_byTanimo
    # from yz_getRank_byTanimo import getRank_byTanimob
    # from yz_getCMC import getCMC

    use_dict = {"C13-NMR": True, "molecular_formula": True,
                "fragment": True, "SMILES": True}
    po = '20220404_bart_3stage_long_1_result/epoch_199_loss_0.067575_1_frg_max_ex_2.json'    
    pi = '20220404_bart_3stage_long_1_result/epoch_199_loss_0.067575_1_frg_max_ex.json'
    removeUnqualified(pi, po)

