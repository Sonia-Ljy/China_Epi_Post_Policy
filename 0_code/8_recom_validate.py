
import pandas as pd
from Bio import SeqIO


def get_pp_loc(mut):
    if "del" in mut or "ins" in mut:
        loc = int(mut.split("_")[1][1:-3])
    elif "stop" in mut:
        loc = int(mut.split("_")[1][1:-4])
    else:
        loc = int(mut.split("_")[1][1:-1])
    return loc


def get_prim_pp(first_m):
    prim_pp = first_m.split("_")[1][0]
    return prim_pp


def trans_del_ins(first_m,pp,prim_pp,loc):
    import re
    if pp == "S":
        new_pp = "Spike_"
    elif pp in ["ORF3","ORF6","ORF7","ORF8","ORF10","ORF3a","ORF6a","ORF8a","ORF10a","ORF3b","ORF6b","ORF8b","ORF10b"]:
        new_pp = "NS"+re.findall(r"\d+\.?\d*",pp)[0]+"_"
    elif pp == "ORF7a":
        new_pp = "NS7a_"
    elif pp == "ORF7b":
        new_pp = "NS7b_"
    else:
        new_pp = pp+"_"
        
    if first_m[-1] == "-":
        second_m = new_pp+prim_pp+loc+"del"
    elif first_m[-1] == "*":
        second_m = new_pp+prim_pp+loc+"stop"
    else:
        second_m = new_pp+prim_pp+loc+first_m.split("_")[1][-1]
    return second_m

def regular_mut(mut_list):
    if "NS6_D61L" in mut_list:
        mut_list.remove("NS6_D61L")
    if "Spike_R158G" in mut_list:
        mut_list.remove("Spike_R158G")
        mut_list.append("Spike_R158del")
    if "Spike_E156del" in mut_list:
        mut_list.remove("Spike_E156del")
        mut_list.append("Spike_E156G")
    if "Spike_L242del" in mut_list:
        mut_list.remove("Spike_L242del")
        mut_list.append("Spike_L244del")
    if "NS7b_F13del" in mut_list:
        mut_list.remove("NS7b_F13del")
    if "NS8_D119del" in mut_list:
        mut_list.remove("NS8_D119del")
    if "NS8_F120del" in mut_list:
        mut_list.remove("NS8_F120del")
    if "NS8_T87I" in mut_list:
        mut_list.remove("NS8_T87I")
    if "NS8_D107E" in mut_list:
        mut_list.remove("NS8_D107E")
    if "NS8_I121L" in mut_list:
        mut_list.remove("NS8_I121L")
    if "NS8_F120V" in mut_list:
        mut_list.remove("NS8_F120V")
    if "NS8_L118V" in mut_list:
        mut_list.remove("NS8_L118V")
    if "NS8_I121del" in mut_list:
        mut_list.remove("NS8_I121del")
    if "NS8_F120L" in mut_list:
        mut_list.remove("NS8_F120L")
        
    return mut_list



def sort_epi_mutation(epiV):
    genome_proteim = ["NSP1","NSP2","NSP3","NSP4","NSP5","NSP6","NSP7","NSP8","NSP9","NSP10","NSP12","NSP13","NSP14","NSP15","NSP16","Spike","NS3","E","M","NS6","NS7","NS8","N","NS10"]
    sort_epiV = []
    for pp in genome_proteim:
        temp_mut = []
        for mut in epiV:
            pp_epi = mut.split("_")[0]
            if pp_epi == pp:
                temp_mut.append(mut)
            else:
                continue
        if len(temp_mut) >= 1:
            mut_loc = {}
            for n in temp_mut:
                if "del" in n or "ins" in n:
                    loc = int(n.split("_")[1][1:-3])
                elif "stop" in n:
                    loc = int(n.split("_")[1][1:-4])
                else:
                    loc = int(n.split("_")[1][1:-1])
                mut_loc[n] = loc
            mut_loc_sorted = sorted(mut_loc.items(), key=lambda x: x[1], reverse=False)
            mut_sort = [k[0] for k in mut_loc_sorted]
            sort_epiV.extend(mut_sort)
        else:
            continue
    return sort_epiV



# 查看有多少是被nextclade判定为重组的序列
dirpath = "/home/soniali/Desktop/02_china_recom_github/0_raw_data/GISAID_nextclade/"
next_recom = []
with open(dirpath+"nextclade.tsv") as f:
    next(f)
    for row in f.readlines():
        if row.split("\t")[2] == "recombinant":
            next_recom.append(row.split("\t")[1])

print(len(next_recom))

# with open("/home/soniali/Desktop/02_china_recom_github/0_raw_data/GISAID_nextclade/nextclade.aligned.fasta","r") as f:
#     for record in SeqIO.parse(f,"fasta"):
#         if record.id in next_recom:
#             with open("/home/soniali/Desktop/02_china_recom_github/3_recom/aligned_china_next_recom_105.fasta","a+") as h:
#                 h.write(">"+str(record.id)+"\n")
#                 h.write(str(record.seq)+"\n")

meta_file = "/home/soniali/Desktop/02_china_recom_github/0_raw_data/Qualified_china_meta_merged.txt"
df_meta = pd.read_csv(meta_file) 
df_meta2 = df_meta[df_meta["Accession_ID"].isin(next_recom)]
df_meta2.to_csv("/home/soniali/Desktop/02_china_recom_github/3_recom/"+"meta_105.csv",index=None,sep = ",")

seq_id_fas = {}
with open("/home/soniali/Desktop/02_china_recom_github/3_recom/aligned_china_next_recom_105.fasta","r") as f:
    for record in SeqIO.parse(f,"fasta"):
        seq_id_fas[str(record.id)] = str(record.seq)
        

# 根据GISAID获取的全球meta文件中所有序列的谱系及氨基酸变异信息，计算各谱系的特征变异，以10%或75%为阈值
import pandas as pd
df = pd.DataFrame(columns=["strain","date","country","province","lineage","mutation"])
with open("/Data4/GISAID_Data/20231026/metadata.tsv") as file:
    next(file)
    num = 0
    for i in file:
        date = i.strip().split("\t")[5]
        if ("2022-12-01" <= date <= "2023-10-30"):
            country = i.strip().split("\t")[6].split("/")[1].strip()
            if country == "China":
                num+=1
                df.loc[num,"strain"] = i.strip().split("\t")[4]
                df.loc[num,"date"] = date
                df.loc[num,"country"] = country
                df.loc[num,"lineage"] = i.strip().split("\t")[13]
                df.loc[num,"mutation"] = i.strip().split("\t")[16][1:-1]
                try:
                    df.loc[num,"province"] = i.strip().split("\t")[6].split("/")[2].strip()
                except:
                    df.loc[num,"province"] = ""

                if num % 10000 == 0 :
                    print(num)
                    df.to_csv("/Data4/GISAID_Data/20231026/metadata_china_"+str(num)+".tsv",index = None)
                    print("df_shape_",str(df.shape[0]))
                    df = pd.DataFrame(columns=["strain","date","country","province","lineage","mutation"])
        else:
            continue

df.to_csv("/Data4/GISAID_Data/20231026/metadata_china_"+str(num)+".tsv",index = None)

# 将上述文件合并为"metadata_china.csv"

# 根据GISAID获取的全球meta文件中所有序列的谱系及氨基酸变异信息，计算各谱系的特征变异，以10%或75%为阈值
DIRPATH = '/home/soniali/Desktop/02_china_recom_github/3_recom/'
variant_surveillance_path = DIRPATH + 'metadata_china.csv'
lineage_10_path = DIRPATH + 'lineagesFM/lineage_10_china.txt'
lineage_75_path = DIRPATH + 'lineagesFM/lineage_75_china.txt'
mutation_num_path = DIRPATH + 'lineagesFM/mutation_num_china.txt'
output_file = DIRPATH + 'putative_recombination.csv'

col_names = ['sample_id','lineage_X', 'lineage_Y', \
        'mutation_pattern', "more_mut","raw_p_value","adjusted_p_value","X_mutations", "Y_mutations", "shared_mutations", "denovo_mutations"]
with open(output_file, "w") as file_epi:
    for c in col_names[0:-1]:
        file_epi.write(c + ",")
    file_epi.write(col_names[-1] + "\n")
    

mutations = []
lineages = {}
with open(variant_surveillance_path,"r") as f:
    next(f)
    num = 0
    for row in f.readlines():
        num+=1
        l = row.strip().split(",")[4]
        if l not in lineages:
            lineages[l] = {'count': 1, 'mutations': {}}
        else:
            lineages[l]['count'] += 1
        
        substitutions = row.strip().split(",")[5:][1:-1]
        substitutions.append(row.strip().split(",")[5:][0].split('"')[-1])
        substitutions.append(row.strip().split(",")[5:][-1].split('"')[0])
        for s in substitutions:
            if 'ins' not in s and s != '':
                mutations.append(s)
                lineages[l]['mutations'][s] = lineages[l]['mutations'].get(s, 0) + 1

mutations = list(set(mutations))
with open(mutation_num_path, 'w') as f:
    f.write('Num'+'\n')
    f.write(str(len(mutations))+'\n')
    
lineages = dict(sorted(lineages.items(), key=lambda x: x[0]))
defining_SNPs_75 = {}
# defining_SNPs_10 = {}
for l in lineages:
    if l != 'None' and l != 'XA' :
        for m in lineages[l]['mutations']:
            if lineages[l]['mutations'][m] / lineages[l]['count'] > 0.1:
                # defining_SNPs_10.setdefault(l, []).append(m)
                if l not in defining_SNPs_75:
                    defining_SNPs_75[l] = []
            if lineages[l]['mutations'][m] / lineages[l]['count'] > 0.75:
                defining_SNPs_75[l].append(m)

# with open(lineage_10_path, 'w') as f:
#     for l in defining_SNPs_10:
#         f.write(l + ',' + str(lineages[l]['count']) + ',' + ','.join(defining_SNPs_10[l]) + '\n')

with open(lineage_75_path, 'w') as f:
    for l in defining_SNPs_75:
        f.write(l + ',' + str(lineages[l]['count']) + ',' + ','.join(defining_SNPs_75[l]) + '\n')


# 提取被nextclade判定为重组的序列
dirpath = "/home/soniali/Desktop/02_china_recom_github/0_raw_data/GISAID_nextclade/"
next_recom_mutation = {}
next_recom = []
with open(dirpath+"nextclade.tsv") as f:
    next(f)
    for row in f.readlines():
        if row.split("\t")[2] == "recombinant":
            info = (row.strip().split("\t"))
            epi = info[1]
            next_recom.append(epi)
            mut_info = info[29].split(",")
            del_info = info[30].split(",")
            mut_info.extend(del_info)
            next_recom_mutation[epi] = mut_info

print(len(next_recom)) #105

query_seq = {}
with open(variant_surveillance_path,"r") as f:
    next(f)
    for row in f.readlines():
        epi = row.strip().split(",")[0].split('"')[1]
        if {epi} - set(next_recom) == set():
            substitutions = row.strip().split(",")[5:][1:-1]
            substitutions.append(row.strip().split(",")[5:][0].split('"')[-1])
            substitutions.append(row.strip().split(",")[5:][-1].split('"')[0])
            query_seq[epi] = substitutions

len(query_seq) #77

num = 0
for epi in next_recom:
    if epi not in query_seq.keys():
        num+=1
        mut_info  = next_recom_mutation[epi]
        mut_list = []
        for mut in mut_info:
            first_m = mut.replace(":","_")
            pp = first_m.split("_")[0]
            prim_pp = get_prim_pp(first_m)
            if pp == "S":
                loc = str(get_pp_loc(first_m))
                second_m = trans_del_ins(first_m,pp,prim_pp,loc)
                mut_list.append(second_m)
            elif pp in ["E","M","ORF3","ORF6","ORF7","ORF8","N","ORF10","ORF3a","ORF6a","ORF7a","ORF8a","ORF10a","ORF3b","ORF6b","ORF7b","ORF8b","ORF10b"]:
                loc = str(get_pp_loc(first_m))
                second_m = trans_del_ins(first_m,pp,prim_pp,loc)
                mut_list.append(second_m)
            elif pp == "ORF1a":
                df1 = df[df["gene_full"] == pp]
                df2 = df1[df1["peptidePos"] == get_pp_loc(first_m)]
                prim_pp = get_prim_pp(first_m)
                df3 = df2[df2["aa"] == prim_pp]
                if len(set(df3["product"].tolist())) == 1:
                    new_pp = df3["product"].tolist()[0]
                    new_loc = str(df3["aaPos"].tolist()[0])
                    second_m = trans_del_ins(first_m,new_pp,prim_pp,new_loc)
                    # print(mut,"     ",second_m)
                    mut_list.append(second_m)
                    # epi_mut_list[epi].extend(mut_list)
                else:
                    print(epi,first_m)
            elif pp == "ORF1b":
                df1 = df[df["gene_full"] == pp]
                df2 = df1[df1["nextclade_aaPos"] == get_pp_loc(first_m)]
                prim_pp = get_prim_pp(first_m)
                df3 = df2[df2["aa"] == prim_pp]
                if len(set(df3["product"].tolist())) == 1:
                    new_pp = df3["product"].tolist()[0]
                    new_loc = str(df3["aaPos"].tolist()[0])
                    second_m = trans_del_ins(first_m,new_pp,prim_pp,new_loc)
                    mut_list.append(second_m)
                else:
                    print(epi,first_m)

        query_seq[epi] = regular_mut(mut_list)
        
print(num)

from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import hypergeom
from collections import Counter

df_china_meta = pd.read_csv("/home/soniali/Desktop/02_china_recom_github/0_raw_data/Qualified_china_meta_merged.txt")
df_china_meta = df_china_meta[~df_china_meta["Accession_ID"].isin(next_recom)]
df_china_meta['date'] = pd.to_datetime(df_china_meta['Sample_Collection_Date']) # date转为时间格式
element_counts = Counter(df_china_meta["Lineage"].tolist())

element_counts_dict = dict(element_counts)
element_counts_dict_sort = dict(sorted(element_counts_dict.items(), key=lambda item: item[1], reverse=True))

Lineage_v = defining_SNPs_75
feature_mutations = mutations
mutaions_num = int(len(mutations))
len_UXY = 3
calculate_num = 0
for epi in next_recom:
    calculate_num+=1
    print(calculate_num,"   ",epi)
    epiV = query_seq[epi]
    epi_feat = len(set(epiV) & set(feature_mutations))
    # P-value for Non-recombination
    epi_record = {}
    aftertime_lin = []
    for lin_A in Lineage_v:
        aftertime_lin.append(lin_A)
        all_AA = len(Lineage_v[lin_A]) 
        all_AA_epi = len(set(Lineage_v[lin_A]) & set(epiV))
        pVal = hypergeom.sf(all_AA_epi - 1, mutaions_num, all_AA, epi_feat)
        epi_record[str(lin_A) + "_" + str(lin_A)] = pVal

    # the least p-value for the Non-recombinant
    min_AA = min(epi_record, key = epi_record.get)
    A_already = []
    for A in aftertime_lin:
        A_already.append(A)
        A_epi = set(Lineage_v[A]) & set(epiV)
        if len(A_epi) < len_UXY:
            continue
        else:
            afterA_linB = set(aftertime_lin) - set(A_already)
            for B in afterA_linB:
                B_epi = set(Lineage_v[B]) & set(epiV)
                if len(B_epi) < len_UXY:
                    continue
                else:
                    unique_A = A_epi - B_epi
                    unique_B = B_epi - A_epi
                    if len(unique_A) < len_UXY or len(unique_B) < len_UXY:
                        continue
                    else:
                        all_AB = len(set(Lineage_v[A]) | set(Lineage_v[B]))  
                        all_AB_epi = len(set(set(Lineage_v[A]) | set(Lineage_v[B])) & set(epiV)) 
                        pVal = hypergeom.sf(all_AB_epi - 1, mutaions_num, all_AB, epi_feat)
                        epi_record[str(A) + "_" + str(B)] = pVal

    raw_pvals = list(epi_record.values())
    rejected, p_adjusted, _, alpha_corrected = multipletests(raw_pvals, alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)

    lin_adjP = {}
    for p in range(len(p_adjusted)):
        lin_pair = list(epi_record.keys())[p]
        lin_adjP[lin_pair] = p_adjusted[p]
    
    # 判断亲本谱系存在
    sorted_dict_by_value = dict(sorted(lin_adjP.items(), key=lambda item: item[1]))
    lin_count = {}
    for pair in sorted_dict_by_value:
        a_lin, b_lin = pair.split("_")[0], pair.split("_")[1]
        if a_lin != b_lin:
            a,b = pair.split("_")[0], pair.split("_")[1]
            if (a in element_counts_dict_sort) == True and (b in element_counts_dict_sort) == True:
                lin_count[pair] = [element_counts_dict_sort[a],element_counts_dict_sort[b]]

    if len(lin_count) >= 1:
        raw_ab = 0
        for ab in lin_count:
            if sum(lin_count[ab]) >= raw_ab:
                temp_ab = ab
                raw_ab = sum(lin_count[ab])

        min_adjp_pair = temp_ab
        
        if min_adjp_pair == min_AA or lin_adjP[min_adjp_pair] >= 0.05: 
            continue
        else:
            print(epi,"   ",min_adjp_pair,"    ",min_AA)

            lin_A_draw, lin_B_draw = min_adjp_pair.split("_")[0],min_adjp_pair.split("_")[1]

            feature_SNPA = Lineage_v[lin_A_draw]
            feature_SNPB = Lineage_v[lin_B_draw]
            A_B_shared = set(feature_SNPA) & set(feature_SNPB)
            UA_mutate = (set(feature_SNPA) & set(epiV)) - set(A_B_shared)
            UB_mutate = (set(feature_SNPB) & set(epiV)) - set(A_B_shared)
            sample_special = set(epiV) - (set(feature_SNPA) | set(feature_SNPB))

            UA_mutate_unique = []
            UB_mutate_unique = []
            shared_mut = []
            denovo_mut = []

            lin_record = ""
            epiV = sort_epi_mutation(epiV)
            for j in epiV:
                if j in A_B_shared:
                    shared_mut.append(j)
                elif j in UA_mutate:
                    UA_mutate_unique.append(j)
                    lin_record = lin_record + "X"
                elif j in UB_mutate:
                    UB_mutate_unique.append(j)
                    lin_record = lin_record + "Y"
                elif j in sample_special:
                    denovo_mut.append(j)
            
            moremut = set(UA_mutate_unique + UB_mutate_unique) - set(Lineage_v[min_AA.split("_")[0]])
            # print("------------------------------------")
            
            if moremut == set():
                continue
            else:
                print("SUCCESS ONE",epi)
                with open(output_file, "a+") as file_epi:
                    file_epi.write(epi + "," +lin_A_draw + "," + lin_B_draw + "," + lin_record + "," +"/".join(moremut)+","+\
                        str(epi_record[min_adjp_pair])+","+str(lin_adjP[min_adjp_pair])+","+"/".join(UA_mutate_unique) + "," + "/".join(UB_mutate_unique) +\
                            "," + "/".join(shared_mut) + "," + "/".join(denovo_mut) + "\n")

# 观察"putative_recombination.csv"结果
putative_recom = ["EPI_ISL_18289734","EPI_ISL_18105656","EPI_ISL_18289815","EPI_ISL_18284946","EPI_ISL_18289774","EPI_ISL_18401744"]
df_china_meta['date'] = pd.to_datetime(df_china_meta['Sample_Collection_Date'])
df_china_meta1 = df_china_meta[df_china_meta["date"]<="2023-08-01"]
df_china_meta2 = df_china_meta1[df_china_meta1["date"]>="2023-07-01"]
df_china_meta3 = df_china_meta2[df_china_meta2['province'].isin(['Shanghai', 'Anhui',"Gansu"])]
element_counts = Counter(df_china_meta3["Lineage"].tolist())
element_counts_dict = dict(element_counts)
element_counts_dict_sort = dict(sorted(element_counts_dict.items(), key=lambda item: item[1], reverse=True))

df_china_105 = pd.read_csv("/home/soniali/Desktop/02_china_recom_github/3_recom/meta_105.csv")
epi_meta = []
for epi in putative_recom:
    df_epi = df_china_105[df_china_105["Accession_ID"] == epi]
    print(list(df_epi["Sample_Collection_Date"])[0], ",",list(df_epi["province"])[0], ",",list(df_epi["Lineage"])[0])
    
df =  pd.read_csv(output_file)
for i in df.index:
    if "XXXX" in df.loc[i,"mutation_pattern"] and "YYYY" in df.loc[i,"mutation_pattern"]:
        print(df.loc[i,"sample_id"])
        
# "FR.1.1"，"EG.5.1.1"

df_china_meta['date'] = pd.to_datetime(df_china_meta['Sample_Collection_Date']) 

df_china_meta1 = df_china_meta[df_china_meta["date"]<="2023-08-01"]
df_china_meta2 = df_china_meta1[df_china_meta1["date"]>="2023-07-01"]
df_china_meta3 = df_china_meta2[df_china_meta2['province'].isin(['Shanghai', 'Anhui',"Gansu"])]

df_EG511 = df_china_meta3[df_china_meta3["Lineage"] == "EG.5.1.1"]
df_EG511_anhui = df_EG511[df_EG511["province"] == "Anhui"] #EPI_ISL_18289730
df_EG511_sh = df_EG511[df_EG511["province"] == "Shanghai"] #EPI_ISL_18105651 *
df_EG511_gansu = df_EG511[df_EG511["province"] == "Gansu"] #EPI_ISL_18114672

df_FR11 = df_china_meta3[df_china_meta3["Lineage"] == "FR.1.1"]
df_FR11_anhui = df_FR11[df_FR11["province"] == "Anhui"] #
df_FR11_sh = df_FR11[df_FR11["province"] == "Shanghai"] #EPI_ISL_18108456 *
df_FR11_gansu = df_FR11[df_FR11["province"] == "Gansu"] #

# set(Lineage_v["FR.1.1"]) - set(Lineage_v["EG.5.1.1"])

with open("/home/soniali/Desktop/02_china_recom_github/0_raw_data/GISAID_nextclade/nextclade.aligned.fasta","r") as f:
    for record in SeqIO.parse(f,"fasta"):
        if record.id in ["EPI_ISL_18101315","EPI_ISL_18289734","EPI_ISL_18105656","EPI_ISL_18289815","EPI_ISL_18284946","EPI_ISL_18289774","EPI_ISL_18401744","EPI_ISL_18108456"]:
            with open("/home/soniali/Desktop/01_china_recom/0_raw_data/next_recom/snipit/FR.1.1_XCN_EG.5.1.1.fasta","a+") as h:
                h.write(">"+str(record.id)+"\n")
                h.write(str(record.seq)+"\n")


infile = "/home/soniali/Desktop/02_china_recom_github/3_recom/snipit"
temp_file = "FR.1.1_EG.5.1.1"
#### left
with open(temp_file+".fasta", "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        ID_name = str(record.id)
        seq = str(record.seq)
        seq_masked = ""
        for i in range(0,len(seq)):
            if i+1 > 22629:
                seq_masked = seq_masked+"n"
            elif i+1 <= 22629:
                seq_masked = seq_masked+seq[i:i+1]
        
        with open(infile+temp_file+"_masked_left.fasta", "a+") as h:
            h.write(">"+ID_name+"\n")
            h.write(seq_masked+"\n")


#### right
with open(infile+temp_file+".fasta", "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        ID_name = str(record.id)
        seq = str(record.seq)
        seq_masked = ""
        for i in range(0,len(seq)):
            if i+1 < 22664:
                seq_masked = seq_masked+"n"
            elif i+1 >= 22664:
                seq_masked = seq_masked+seq[i:i+1]
        
        with open(infile+temp_file+"_masked_right.fasta", "a+") as h:
            h.write(">"+ID_name+"\n")
            h.write(seq_masked+"\n")