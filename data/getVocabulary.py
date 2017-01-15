#coding: utf-8
# 檔案內容讀取範例
檔案 = open("WebsterUnabridgedDictionary.txt",encoding="UTF-8")
count = 0
所有資料 = 檔案.readlines()
print("這個檔案共有:"+str(len(所有資料))+"行")

for 行數 in range(len(所有資料)):
    # 若該行的所有字母都是大寫, 表示為單字開頭
    行資料 = 所有資料[行數].rstrip()
    if(行資料.isupper()):
        # 列出單字
        #print(行資料.encode("UTF-8"))
        count += 1
    else:
        # 列出單字的所屬解釋
        #print("解釋", end="")
        #print(行資料.encode("UTF-8"))
        continue

print("單字數共有:"+str(count)+"個")
檔案.close() # 關閉檔案
