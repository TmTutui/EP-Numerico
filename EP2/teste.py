
f = open("teste.txt", "r")
for i in f:
    print(i, type(i))
    for a in i:
        print(a)
        continue
    break
    
f.close()