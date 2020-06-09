
lista = []
for k1 in range(10):
    lista_aux = []
    for k2 in range(10):
        lista_aux.append(k1*k2)
    lista.append(lista_aux)

for i in range(len(lista)):
    print(lista[i])
