# Um problema inverso para obtenção de distribuição de Temperatura 
# EP - MAP3121

Em um problema inverso, a partir do efeito tentamos determinar a causa. Neste EP (Exercício Programa) resolvemos um problema ligado à equação do calor. Ele está dividido em duas tarefas.

## Primeira Tarefa

Nesta fase, vamos tratar de um problema direto, sobre como determinar a evolução da distribuição de temperatura em uma barra sujeita a fontes de calor, a partir de uma dada distribuição inicial.

Nela, implementamos o método (11) disponível no enunciado, de tal maneira que as variáveis N é definida em tempo de execução. 
Para os três itens, usamos lambda = 0.25 , lambda = 0.50 e lambda = 0.51 .

## Segunda Tarefa
Nesta segunda fase vamos utilizar um método implícito (Euler implícito (item b) ou método de Crank-Nicolson)

__________________

## Getting Started

### Pré-requisitos
```
Python 3.5+
Numpy == 1.18.2
Matplotlib == 3.2.1
tqdm == 4.45.0
```

#### Instalação
```
$ pip install numpy
$ pip install matplotlib
$ pip install tqdm
```
## Running
Estando na pasta raiz do projeto:
```$ python3 main.py```

Em seguida, serão solicitadas as seguintes informações:
    1- Qual a tarefa deseja rodar (1 ou 2);
    2- Qual item deseja rodar (**a**, **b** ou **c** para **tarefa 1**; **b** ou **c** para **tarefa 2**)
    3- Digite o valor almejado de N

Os gráficos que correspondem à saída do programa serão salvos em ```path/to/root/NÚMERO_ESCOLHIDO_PARA_tarefa/figuras_ALTERNATIVA_ESCOLHIDA```.


# Autores:
    - Geraldo Marques de Sousa Junior - 10771480
    - Túlio Navarro Tutui - 10290208



