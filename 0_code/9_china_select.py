
# 根据GISAID获取的全球meta文件中所有序列的谱系及氨基酸变异信息
import pandas as pd
df = pd.DataFrame(columns=["strain","date","country","province","lineage","mutation"])
with open("/Data4/GISAID_Data/20231118/metadata_20231118.tsv") as file:
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

df.to_csv("/Data4/GISAID_Data/20231118/metadata_china_"+str(num)+".tsv",index = None)
# 将上述文件合并为"metadata_china.csv"