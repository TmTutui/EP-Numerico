#!/usr/bin/env python3
def main():
    
    tarefa = int(input("Digite a tarefa (1 ou 2): "))
    
    if(tarefa == 1):
        alternative = input("Type an alternative (a,b ou c):")
        if(alternative.lower() == "a"):      
            from Primeira_tarefa import letra_a
            letra_a

        elif(alternative.lower() == "b"):
            from Primeira_tarefa import letra_b
            letra_b
            
        elif(alternative.lower() == "c"):
            from Primeira_tarefa import letra_c
            letra_c
            
        else:
            print("Não existe esse item, tente novamente!")


    elif(tarefa == 2):
        alternative = input("Type an alternative (b ou c):")

        if(alternative.lower() == "b"):
            from Segunda_tarefa import item_b
            item_b

        elif(alternative.lower() == "c"):
            from Segunda_tarefa import item_c
            item_c
        

main()
