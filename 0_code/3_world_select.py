import pandas as pd

df = pd.DataFrame(columns=["strain","date","country","lineage"])
with open("/Data4/GISAID_Data/20231118/metadata_20231118.tsv") as file:
    next(file)
    num = 0
    for i in file:
        date = i.strip().split("\t")[5]
        if ("2022-12-01" <= date <= "2023-10-30"):
            country = i.strip().split("\t")[6].split("/")[1].strip()
            if country in ["USA","China","Japan","United Kingdom"]:
                num+=1
                df.loc[num,"strain"] = i.strip().split("\t")[4]
                df.loc[num,"date"] = date
                df.loc[num,"country"] = country
                df.loc[num,"lineage"] = i.strip().split("\t")[13]
                if num % 10000 == 0 :
                    print(num)
                    df.to_csv("/Data4/GISAID_Data/20231118/metadata_enrolled_"+str(num)+".tsv",index = None)
                    print("df_shape_",str(df.shape[0]))
                    df = pd.DataFrame(columns=["strain","date","country","lineage"])
        else:
            continue

df.to_csv("/Data4/GISAID_Data/20231026/metadata_enrolled_."+str(num)+"tsv",index = None)


