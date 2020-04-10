def main():
    escolhido = False
    
    while (escolhido == False):
        alternative = input("Type an alternative (a,b ou c):")
        if(alternative.lower() == "a"):      
            escolhido = True
            letra_a()

        # elif(alternative.lower() == "b"):
        #     escolhido = True
        #     letra_b()

        # elif(alternative.lower() == "c"):
        #     escolhido = True
        #     letra_c()

        else:
            print(" You did not type an existing alternative! ")

main()