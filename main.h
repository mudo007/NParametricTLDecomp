
/*
 Decomp_nparametrica: "Reorganizes a multivariable equation so it is
 suitable to be realized onto a singluar translinear loop analogue 
 current mode circuit."
 
 Copyright (C) 2014  Diogo Andrade
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 email:diogo007@gmail.com
 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

//definição de debug
//#define DEBUG_LEXICO
//#define DEBUG_EXPR
//#define DEBUG_EXPAND
//#define DEBUG_SIMPLIFY
//#define DEBUG_POLYDIV
//#define DEBUG_PARFRAC
//#define DEBUG_VECTOR_GEN
//#define DEBUG_MEMORIA

//Mensagens de erro
#define ERRO_001    "ERRO: Memoria Insuficiente"
#define ERRO_002    "ERRO: Equacao muito longa"
#define ERRO_003    "ERRO: Equacao errada"
#define ERRO_004    "ERRO: Limite de variaveis violado (50)"

//definicoes boleanas
#define TRUE    1
#define FALSE   0

//tipos
//estrutura de tokens para parse inicial
typedef struct _token
{
    int tipo;           //dependendo do tipo do token, uma das vari·veis abaixo È preenchida
    char *literal;      //string do literal
    char codigo;        //caracter do código correspondente à string do literal, ou no caso de operadores e abre_fecha pareneteses, o caracter correspondente
    int parametro;      //valor numerico do parametro
    int operador;       //operador matem·tico
    int abre_fecha;     //abre ou fecha parenteses
    struct _token *proximo_token;  //ponteiro para o proximo token da fila.
}token;

//estrutura para a lista de literais e seus códigos, limitados a 50 no total
typedef struct _tabela_literais
{
    char *literal;
    char codigo;     //o codigo é só uma letra de A-Z maiuscula ou minuscula
    struct _tabela_literais *proximo_codigo;
} tabela_literais;

//Estrutura da lista que representará a expressão expandida
typedef struct _lista_expr
{
	char	*codigos_numerador;	//As variáveis contidas no numerador do termo
	struct _lista_expr 	*codigos_denominador;	//As variáveis contidas no numerador do termo
	//char 	sinal;		//se é positivo ou negativo
	double 	parametro;	//Inteiro que multiplica o termo
	struct _lista_expr	*proximo;
	struct _lista_expr	*anterior;
} lista_expr;

//Estrutura da árvore de expressões:
typedef struct _arvore_expr
{
    token  *elemento;               //ponteiro para o token que será inserido na árvore
    struct _arvore_expr *esquerda;  //ponteiro para o elemento a direita
    struct _arvore_expr *direita;   //ponteiro para o elemento a esquerda
} arvore_expr;

//Estrutura de pilha para conversão da expressão em RPN
typedef struct _pilha_expr
{
    token  *elemento;               //ponteiro para o token que será inserido na pilha   
    struct _pilha_expr *proximo;
} pilha_expr;

//Estrutura das pilhas de construção da árvore
typedef struct _pilha_arvore
{
    arvore_expr *node;
    struct _pilha_arvore *proximo;
} pilha_arvore;

typedef struct _polinomio
{
    lista_expr *P;
    int         id;
    
}polinomio;

//estrutura de vetor de combinação linear das variáveis de entrada
typedef struct _vetor_polinomios
{
    polinomio *polinomio;
    struct _vetor_polinomios *proximo_polinomio;  
	struct _vetor_polinomios *polinomio_anterior;  
    
}vetor_polinomios;

//estrutura que guarda combinações de polinomios que podem ser semente de uma decomposição
typedef struct _vetor_sementes
{
    polinomio P1;
    polinomio P2;
    lista_expr *quociente;
    lista_expr *R1;
    lista_expr *R2;
    struct _vetor_sementes *conjunto_prox;  
	struct _vetor_sementes *conjunto_ant;  
    
}vetor_sementes;



//estrutura que guarda decomposição paramétrica encontrada
typedef struct _vetor_decomp
{
    vetor_polinomios *poly_pares;
    vetor_polinomios *poly_impares;
    lista_expr *resto_impar;
    lista_expr *resto_par;
    
    struct _vetor_decomp *prox_decomp;
    struct _vetor_decomp *ant_decomp;
    
    
} vetor_decomp;

//variavel global para debug
int global = 0;

//definições para a struct

#define tipo_literal     1
#define tipo_parametro  2
#define tipo_operador   3
#define tipo_abre_fecha 4

//tipos de operadores
#define operador_mais           1
#define operador_menos          2
#define operador_multiplicacao  3
#define operador_potenciacao    4
#define operador_negacao        5
//TODO #define operador_divis„o 4

//tipos de abre_fecha
#define abre_parenteses     0
#define fecha_parenteses    1

//protótipos de funções de análise léxica
void erro (char* );
token *le_tokens(char *);
/*int isdigit(char );
int isalpha (char );
int isspace(char );*/
int isoperator (char );
int isabre_fecha (char );

token *cria_token(token **);
void destroi_lista( token *);
char insere_tabela_literais(tabela_literais **, char*);
void destroi_tabela_literais(tabela_literais *);
void constroi_tabela_literais(tabela_literais **, token *);


//protótipo das funções de construtor de árvore de expressões
token *retira_pilha (pilha_expr **);
int insere_pilha(pilha_expr **, token *);
void destroi_pilha(pilha_expr *);
int insere_lista_expr(token **, token *);
token *constroi_lista_expr(token *);
int prioridade_operador(token *);
arvore_expr *constroi_arvore_expr(token *);
arvore_expr *cria_no_arvore(token *);
int insere_pilha_arvore(pilha_arvore **, arvore_expr *);
arvore_expr *retira_pilha_arvore (pilha_arvore **);
void destroi_arvore_expr(arvore_expr *);
void destroi_lista_expr( token *);



//protótipo das funções de expansor de expressões
void string_sort(char **);
lista_expr *constroi_lista_expressoes_exp(arvore_expr *);
lista_expr *multiplica_expr( lista_expr *, lista_expr *);
void destroi_lista_expr_expandida(lista_expr *);

//protótipo das funções de simplificador de expressıes
lista_expr *simplifica_expr_expandida(lista_expr *);
lista_expr *remove_lista_expr(lista_expr *, lista_expr *);
lista_expr *copia_lista_expr(lista_expr *);


//protótipo de funções de divisão polinomial
int lexdeg(char *, char *);
lista_expr *lexdegbubblesort(lista_expr *);
int polydiv(lista_expr *, lista_expr *, lista_expr **, lista_expr **);
lista_expr *divide_monomio(lista_expr *,lista_expr *);
void subtrai_expr(lista_expr **, lista_expr *);
void soma_expr(lista_expr *, lista_expr *);
lista_expr *constroi_elemento_zerado(void);

//protótipo das funcoes de expansao em fracoes parciais
int partial_fraction_expansion(lista_expr *, lista_expr *, lista_expr *,lista_expr **, lista_expr **, lista_expr **);
lista_expr *substitui_var(lista_expr *, lista_expr *, char );


//protótipo das funcoes de geracao do vetor de polinomios
vetor_polinomios *gera_vetor(vetor_polinomios *, lista_expr *, lista_expr *, int , int );
vetor_polinomios *elimina_zero(vetor_polinomios *);
vetor_polinomios *remove_polinomio(vetor_polinomios *);
vetor_polinomios *remove_polinomio_retorna_anteror(vetor_polinomios *);
lista_expr *gera_polinomio_base(tabela_literais *);
vetor_polinomios *remove_polinomios_negativos (vetor_polinomios *);
vetor_polinomios *remove_polinomios_redundantes (vetor_polinomios *);

//prototipo das funcoes que implementam o algoritmo propriamente dito
int deg(lista_expr *);
vetor_sementes *gera_vetor_semente(vetor_polinomios *, lista_expr *);
vetor_sementes *novo_vetor_semente(void);
void destroi_lista_sementes(vetor_sementes *);
vetor_polinomios *novo_vetor_polinomios(void);
vetor_decomp *novo_vetor_decomp(void);
void insere_polinomio(vetor_polinomios **,polinomio *);
void destroi_decomp(vetor_decomp *);
void destroi_vetor_polinomios(vetor_polinomios *);
void destroi_vetor_decomp(vetor_decomp *);
vetor_decomp *encontra_decomp(vetor_sementes *, int, lista_expr *, tabela_literais *);
int prova_real(vetor_decomp *, lista_expr *);
vetor_decomp *insere_lista_decomp(vetor_decomp *, vetor_decomp *);
vetor_decomp *copia_vetor_semente(vetor_sementes *);
void encontra_decomp_recursiva(vetor_decomp *, vetor_decomp *, vetor_decomp **, int, lista_expr *, tabela_literais *, int *);
vetor_decomp *copia_semente(vetor_sementes *);
vetor_decomp *copia_decomp(vetor_decomp *);
void elimina_decomp_redundantes(vetor_decomp *);
vetor_polinomios *ordena_polinomios(vetor_polinomios *);
int compara_decomp(vetor_decomp *, vetor_decomp *);

//funcoes de tamanho de memoria
int decomp_size(vetor_decomp *);
int poly_size(vetor_polinomios *);
int expr_size(lista_expr *);

//funcoes de contagem
void encontra_decomp_dummie(vetor_sementes *, int);
void encontra_decomp_recursiva_dummie(vetor_decomp *, vetor_decomp *, int , int *, int);
void encontra_decomp_dummie_mulder(vetor_sementes *,vetor_polinomios *, int);
void encontra_decomp_recursiva_dummie_mulder(vetor_decomp *, vetor_polinomios *, int , int *, int);




//protótipo das funcoes de interface
void imprime_lista_expr_expandida(lista_expr *, tabela_literais *);
void print_monomio(char *, tabela_literais *);
void imprime_arvore_expr(arvore_expr *);
void imprime_tokens(token* );
void imprime_decomposicao(vetor_decomp *, tabela_literais *);



//erro - imprime uma mensagem de erro na tela
void erro (char* n_erro)
{
    printf("\n\r\t%s\n",n_erro);
}

//le_tokens - realiza o parse da equaÁ„o de entrada para construir uma lista de tokens
token *le_tokens(char *equacao_entrada)
{
    //o parser é uma maquina de estados que percorre a string da equaÁ„o caracter a caracter
    int estado = 0;                             //variavel de estado
    char *inicio_literal = NULL;                //ponteiro temporario para manipulacao de literais
    token *lista_token = NULL;                  //ponteiro para a lista de tokens a ser retornada
    token *token_atual = NULL;                  //proximo token a ser preenchido;
    
    //variaveis para calculos temporarios
    int n_chars;    //calculo de numero de caracteres de um literal
    int cont = 0;   //contador de caracteres       
    //repetir atÈ encontrar o caracter nulo, e continuar o laço se estiver no meio de um reconhecimento de literal.
    while (*(equacao_entrada + cont) != '\0' || estado == 2)
    {
#if defined DEBUG_LEXICO
    printf("\n%c\n",*(equacao_entrada + cont));
#endif 
        switch(estado)
        {
            case 0: //estado inicial
                if (isspace(*(equacao_entrada + cont)))
                {
                    //nada È feito
                    break;
                }
                //operadores matematicos
                if (isoperator(*(equacao_entrada + cont)))
                {

                    //cria um token na lista
                    token_atual = cria_token(&lista_token);
                    
                    //preenche o tipo
                    token_atual->tipo = tipo_operador;
                    
                    //descobre qual o operador e preenche na estrutura
                    switch(*(equacao_entrada + cont))
                    {
                        case '+':
                            token_atual->operador = operador_mais;
                            token_atual->codigo = '+';
                            break;
                        case '-':
                            //a diferenciação entre subtração e a operação de negação será tratada mais tarde.
                            token_atual->operador = operador_menos;
                            token_atual->codigo = '-';
                            break;
                        case '*':
                            token_atual->operador = operador_multiplicacao;
                            token_atual->codigo = '*';
                            break;
                        case '^':
                            token_atual->operador = operador_potenciacao;
                            token_atual->codigo = '^';
                            break;     
                    }
                    break;
                }
                //parametros
                if (isdigit(*(equacao_entrada + cont)))
                {
                    //cria um token na lista
                    token_atual = cria_token(&lista_token);
                    
                    //preenche o tipo
                    token_atual->tipo = tipo_parametro;
                    
                    //comeÁa a preencher o numero
                    token_atual->parametro = (int)(*(equacao_entrada + cont) - '0');

                    
                    //muda o estado
                    estado = 1;                    
                    break;
                }
                //literais
                if (isalpha(*(equacao_entrada + cont)))
                {

                    //cria um token na lista
                    token_atual = cria_token(&lista_token);
                    
                    //preenche o tipo
                    token_atual->tipo = tipo_literal;
                    
                    //marca o inicio do literal
                    token_atual->literal = (equacao_entrada + cont);
                    
                    //muda o estado
                    estado = 2;
                    
                    //há um erro aqui. quando está lendo o ultimo caractere, e atualiza cont após o break, o laço é interrompido, pois cehga no \0
                    
                    break;
                }
                //abre e fecha parenteses
                if (isabre_fecha(*(equacao_entrada + cont)))
                {

                    //cria um token na lista
                    token_atual = cria_token(&lista_token);
                    
                    //preenche o tipo
                    token_atual->tipo = tipo_abre_fecha;
                     //descobre qual o operador e preenche na estrutura
                    switch(*(equacao_entrada + cont))
                    {
                        case '(':
                            token_atual->abre_fecha = abre_parenteses;
                            token_atual->codigo = '(';
                            break;
                        case ')':
                            token_atual->abre_fecha = fecha_parenteses;
                            token_atual->codigo = ')';
                            break; 
                    }
                    break;
                }
                
                //se n„o encontrar nenhum dos caracteres esprados, a equaÁ„o est· errada
                destroi_lista(lista_token);
                erro(ERRO_003);
                //imprimir tamÈm a equacao a partir do erro
                printf("-> %s\n",(equacao_entrada + cont));
                return(NULL);
                break;
                
            case 1: //numero
                //se o proximo caracter È um digito, multiplicar o parametro atual por 10 e adicionar o digito

                if (isdigit(*(equacao_entrada + cont)))
                {  
                    token_atual->parametro *= 10;
                    token_atual->parametro += (int)(*(equacao_entrada + cont) - '0');  
                }
                else
                {
                     //qualquer outra coisa significa que o numero ja foi lido
                    //volta o caracter, para que seja lido no proximo loop
                    cont--;
                    estado = 0;
                }
                break;
                
            case 2: //literal
                //se o proximo for digito ou letra, ou _, continuar lendo
                if(isdigit(*(equacao_entrada + cont)) || isalpha(*(equacao_entrada + cont)) || (*(equacao_entrada + cont) == '_'))
                    break;
                else
                {
                    //salva a posiÁ„o do inico do literal
                    inicio_literal = token_atual->literal;
                    
                    //calcula o numero de caracteres necessarios
                    n_chars = (int)((equacao_entrada + cont) - inicio_literal);
                    
                    //aloca memoria para o literal encontrado com espaÁo para o terminador nulo
                    if((token_atual->literal = (char *)malloc((n_chars + 1)*sizeof(char))) == NULL)
                    {
                        erro(ERRO_001);
                        destroi_lista(lista_token);
                        return (NULL);
                    }
                    //copia a string do literal;
                    if((strncpy ( token_atual->literal, inicio_literal, n_chars )) == NULL)
                    {
                        erro(ERRO_001);
                        destroi_lista(lista_token);
                        return (NULL);
                    }
                   
                    //coloca o terminador nulo
                    *(token_atual->literal + n_chars) = '\0';

                    //volta o caracter, para que seja lido no proximo loop
                    cont--;
                    estado = 0;
                    break;  
                }
                break;
                
        }
        //proximo caracter a ser lido
        cont++;
    }
    //retorna a lista de tokens
    return (lista_token);
}

//isdigit - verifica se um caracter È um dÌgito
/*int isdigit(char digito)
{
    if( (digito >= '0') && (digito <= '9'))
        return TRUE;
    else
        return FALSE;  
}

//isalpha - verifica se um caracter È uma letra do alfabeto ou underline
int isalpha (char digito)
{
    if (((digito >= 'a') && (digito <= 'z')) || ((digito >= 'A') && (digito <= 'Z')) || (digito == '_'))
            return TRUE;
    else
        return FALSE;
}
 //isspace - verifica se o caracter È um espaÁo em brtanco ou pula-linha
 int isspace(char digito)
 {
 if ((digito == ' ') || (digito == '\t') || (digito == '\n'))
 {
 
 return TRUE;
 }
 else
 return FALSE;
 }*/

//isoperator - verifica se o caracter È um operador matem·tico
int isoperator (char digito)
{
    if ((digito == '+') || (digito == '-') || (digito == '*') || (digito == '^'))
            return TRUE;
    else
        return FALSE;
}

//isoabre_fecha - verifica se o caracter È um abre ou fecha parenteses
int isabre_fecha (char digito)
{
    if ((digito == '(') || (digito == ')'))
            return TRUE;
    else
        return FALSE;
}


//cria_token - aloca memoria para um novo token no fim da fila e retorna o ponteiro para ele
token *cria_token(token **lista_token)
{
    token *token_atual;     //variavel para percorrer a lista
    
    //verifica se a lista está vazia e a cria se necess·rio
    if ((*lista_token) == NULL)
    {
        if (((*lista_token) = (token *)malloc(sizeof(token))) != NULL)
        {
            (*lista_token)->proximo_token = NULL;
            (*lista_token)->tipo = 0;
            (*lista_token)->literal = NULL;
            (*lista_token)->parametro = 0;
            (*lista_token)->operador = 0;
            (*lista_token)->abre_fecha = 0;
            
            return(*lista_token);
        }
    }
    else
    {
        token_atual = (*lista_token);
        while(token_atual->proximo_token != NULL)
        {
            token_atual = token_atual->proximo_token;
        }
        if ((token_atual->proximo_token = (token *)malloc(sizeof(token))) != NULL)
        {
            token_atual->proximo_token->proximo_token = NULL;
            token_atual->proximo_token->tipo = 0;
            token_atual->proximo_token->literal = NULL;
            token_atual->proximo_token->parametro = 0;
            token_atual->proximo_token->operador = 0;
            token_atual->proximo_token->abre_fecha = 0;
            return(token_atual->proximo_token);
        }
    }
    //se ocorreu um erro apenas retorna NULL
    return (NULL);
}

//destroi lista - funÁ„o para destruir todos os ponteiros da lista em caso de erro ou fim do programa
void destroi_lista( token *lista)
{
    //percorre a lista recursivamente atÈ encontrar o ultimo elemento
    if(lista->proximo_token != NULL)
    {
        destroi_lista(lista->proximo_token);
    }
    //destroi os ponteiros dos literais e dos proprios tokens nos retornos das chamadas recursivas.
    if (lista->literal != NULL)
    {
        free(lista->literal); 
    }
    lista->literal = NULL;
    lista->proximo_token = NULL;
    free(lista);
}


//percorre a lista imprimindo os tokens encontrados
void imprime_tokens(token *lista_token)
{
    int contador = 0;
    token   *p_lista = lista_token;
    
    while(p_lista != NULL)
    {
        switch(p_lista->tipo)
        {
            case tipo_literal:
                printf("#%d literal \"%s\"\n\r ",contador, p_lista->literal);
                contador++;
                break;
                
            case tipo_parametro:
                printf("#%d parametro \"%d\"\n\r ",contador, p_lista->parametro);
                contador++;
                break;
            case tipo_operador:
                printf("#%d operador \"%c\"\n\r ",contador, p_lista->codigo);
                contador++;
                break;
            case tipo_abre_fecha:
                printf("#%d abre_fecha \"%c\"\n\r ",contador, p_lista->codigo);
                contador++;
                break;
                
        }
        p_lista = p_lista->proximo_token;
    }
}

//cria_tablea_literais - aloca memoria para um novo literal e retorna o seu codigo
char insere_tabela_literais(tabela_literais **lista_literais, char *literal)
{
    tabela_literais *codigo_atual;     //variavel para percorrer a lista
    
    //verifica se a lista est· vazia e a cria se necess·rio
    if ((*lista_literais) == NULL)
    {
        
        if (((*lista_literais) = (tabela_literais *)malloc(sizeof(tabela_literais))) != NULL)
        {
            //alocar memoria e copiar a string
            (*lista_literais)->literal = (char *)malloc(strlen(literal)+1); //TODO : erro alocacao memoria
            strcpy((*lista_literais)->literal, literal);
            (*lista_literais)->codigo = 'A'; //o primeiro 
            (*lista_literais)->proximo_codigo = NULL;
            
            return((*lista_literais)->codigo);
        }
    }
    else
    {
        codigo_atual = (*lista_literais);
        
        while(1)
        {
            //procura se o cÛdigo j· n„o existe
            if ( !strcmp(codigo_atual->literal, literal))
                return(codigo_atual->codigo);   
            
            if(codigo_atual->proximo_codigo == NULL)
                break;
            else
               codigo_atual = codigo_atual->proximo_codigo; 
        }
        
        
        if ((codigo_atual->proximo_codigo = (tabela_literais *)malloc(sizeof(tabela_literais))) != NULL)
        {
            codigo_atual->proximo_codigo->proximo_codigo = NULL;
            //alocar memoria e fazer uma stringcopy
            codigo_atual->proximo_codigo->literal = (char *)malloc(strlen(literal)+1); //TODO : erro alocacao memoria
            strcpy(codigo_atual->proximo_codigo->literal, literal);
            if (codigo_atual->codigo == 'Z') //passar para as letras minusculas
                codigo_atual->proximo_codigo->codigo = 'a';
            else if (codigo_atual->codigo == 'z') //estouro de codigos
            {
                erro(ERRO_004);
                return(0);
            }
            else
                codigo_atual->proximo_codigo->codigo = codigo_atual->codigo + 1;
            return(codigo_atual->proximo_codigo->codigo);
        }
    }
    //se ocorreu um erro apenas retorna NULL
    return ('\0');
}
//destroi lista - funÁ„o para destruir todos os ponteiros da lista em caso de erro ou fim do programa
void destroi_tabela_literais(tabela_literais *lista_literais) 
{
    //percorre a lista recursivamente atÈ encontrar o ultimo elemento
    if(lista_literais->proximo_codigo != NULL)
    {
        destroi_tabela_literais(lista_literais->proximo_codigo);
    }
    //destroi os elementos de forma recursiva.
    free(lista_literais->literal);
    free(lista_literais);
}

//percorre a lista imprimindo os tokens encontrados
void constroi_tabela_literais(tabela_literais **lista_literais, token *lista_token)
{
    token   *p_lista = lista_token;
    char codigo;

    while(p_lista != NULL)
    {
        if(p_lista->tipo == tipo_literal)
        {
            codigo = insere_tabela_literais(lista_literais, p_lista->literal);
            //tem que trocar a string do token apenas pelo codigo na lista de tokens, senão isso não faz sentido.
            p_lista->codigo = codigo;
#if defined DEBUG_LEXICO
            printf("#codigo:%c literal \"%s\"\n\r ",codigo,  p_lista->literal);
#endif
        }
        p_lista = p_lista->proximo_token;
    }
}


//Funcoes de manipulação da pilha

//destruicao de pilhas
void destroi_pilha(pilha_expr *pilha_destruida)
{
    if (pilha_destruida->proximo != NULL)
        destroi_pilha(pilha_destruida->proximo);
    free(pilha_destruida);
}

//inserir um elemento na pilha
int insere_pilha(pilha_expr **pilha, token *elemento)
{
    pilha_expr *topo_pilha;
    
    //aloca memoria para o proximo elemento
    topo_pilha = (pilha_expr*)malloc(sizeof(pilha_expr));
    if(topo_pilha == NULL)
        return(0); //houve erro de falta de memoria
    
    //aponta o elemento da pilha para o token adequado
    topo_pilha->elemento = elemento;
    
    //posiciona o novo elemento no topo da pilha
    topo_pilha->proximo = *pilha;
    *pilha = topo_pilha;
    
    return (1); //sucesso
    
}

//retirar um elemento da pilha
token *retira_pilha (pilha_expr **pilha)
{
    token *elemento;
    pilha_expr *topo_pilha;
    
    //guarda o token que será retornado
    elemento = (*pilha)->elemento;
    
    //guarda o novo topo da pilha
    topo_pilha = (*pilha)->proximo;
    
    //destroi o ponteiro a ser removido
    free(*pilha);
    
    //atualiza o topo da pilha
    *pilha = topo_pilha;
    
    //retorna o elemento
    return(elemento);
    
}

//construir pilha de operandos e operadores
token *constroi_lista_expr(token *lista_token)
{   
    
    //crio uma variavel para percorrer a lista de tokens
	token *p_lista = lista_token;
    //crio a pilha de operadores
    pilha_expr *pilha_operadores = NULL;
    //crio a lista de nós da arvore de expressoes
    token *lista_expressao = NULL;
    //flag para detecção de oeprador negativo unario no inicio da expressao
    int flag_unario = 1;
    //enquanto não chegar ao final da lista, iterar:

    while(p_lista != NULL)
    {
        switch(p_lista->tipo)
        {
            case tipo_parametro:
            case tipo_literal: //pilha de operandos
                insere_lista_expr(&lista_expressao,p_lista);
                break;

            case tipo_operador:
                //unario: caso especial. Se encontrou um operador qualquer, e o proximo token é o operador '-', significa que é uma operação unaria,
                // e devo modificar o token seguinte
                if (p_lista->proximo_token != NULL)
                    if (p_lista->proximo_token->tipo == tipo_operador) 
                        if (p_lista->proximo_token->operador == operador_menos)
                            p_lista->proximo_token->operador = operador_negacao;
                //unario: devo verificar a flag_unario, caso seja 1, e for um operador '-', significa que é um operador unario no inicio da expressao
                if (flag_unario && (p_lista->operador == operador_menos)) 
                    p_lista->operador = operador_negacao;
                //devo verificar a precedencia do operador (quanto menor o valor, maior a precedencia) e inserir na pilha diretamente se o novo estiver prioridade maior em relação ao topo da pilha.
                if ((pilha_operadores != NULL) && (prioridade_operador(p_lista) < prioridade_operador(pilha_operadores->elemento)))
                {
                    insere_pilha(&pilha_operadores, p_lista);
                }
                else 
                //caso contrario, devo remover todos os operadores da pilha de maior prioridade e inserir o novo operador na pilha
                {
                    while ((pilha_operadores != NULL) && (prioridade_operador(p_lista) >= prioridade_operador(pilha_operadores->elemento)))
                       insere_lista_expr(&lista_expressao,retira_pilha(&pilha_operadores)); 
                    //inserir o novo operador na pilha
                    insere_pilha(&pilha_operadores, p_lista);

                }
                break;
                
            case tipo_abre_fecha: //precedencia de parenteses
                //unario: caso especial. Se o proximo token é o operador '-', após abre parenteses, significa que é uma operação unaria,
                // e devo modificar o token seguinte
                //abre parenteses é inserido na pilha, mas jamais será colocado na lista de expressao, serve apenas como controle
                if (p_lista->abre_fecha == abre_parenteses)
                {
                    if (p_lista->proximo_token != NULL)
                        if (p_lista->proximo_token->tipo == tipo_operador) 
                            if (p_lista->proximo_token->operador == operador_menos)
                                p_lista->proximo_token->operador = operador_negacao;
                    insere_pilha(&pilha_operadores, p_lista);
                }
                else if (p_lista->abre_fecha == fecha_parenteses)
                //se for fecha parenteses, tenho que retirar todos os operadores até encontrar o abre-parenteses na pilha, e colocá-los na lista
                {
                    while(pilha_operadores->elemento->tipo != tipo_abre_fecha)
                    {
                        //retiro um operador da pilha e coloco no fim da lista de expressao
                        insere_lista_expr(&lista_expressao,retira_pilha(&pilha_operadores));
                        
                        //se retirou o ultimo elemento da pilha sem encontrar abre-fecha, erro de parentees
                        if (pilha_operadores == NULL)
                        {
                            //erro parenteses mal colocados
                            return NULL;
                        }
                    }
                    //retiro o abre parenteses que sobrou na pilha
                    if(pilha_operadores->elemento->abre_fecha == abre_parenteses)
                        retira_pilha(&pilha_operadores);
                }
                break;
                    
        }
        //limpa o flag_unario, pois já não estamos mais no primeiro token.
        flag_unario = 0;
        p_lista = p_lista->proximo_token;
        
    }
    //desempilhar todos os operadores restantes
    while(pilha_operadores != NULL)
    {
        insere_lista_expr(&lista_expressao,retira_pilha(&pilha_operadores));
    }
    return(lista_expressao);
}


void imprime_lista_expr(token *lista)
{
    token   *p_lista = lista;
    
    //titulo
    printf("\n A expressao em formato RNP e: ");
    while(p_lista != NULL)
    {
        switch(p_lista->tipo)
        {
            case tipo_literal:
                printf("%s ", p_lista->literal);
                break;
                
            case tipo_parametro:
                printf("%d ", p_lista->parametro);
                break;
                
            case tipo_operador:
            case tipo_abre_fecha:

                printf("%c ", p_lista->codigo);
                break;
                
        }
        p_lista = p_lista->proximo_token;
    }
}


//inserir um elemento na lista encadeada dupla d expressões
int insere_lista_expr(token **lista, token *elemento)
{
    token *p_lista = *lista;
    token *fim_lista = NULL;
    
    //aloca memoria para o novo elemento
    fim_lista = (token *)malloc(sizeof(token));
    if(fim_lista == NULL)
        return(0); //houve erro de falta de memoria
    
    //preenche os campos do elemento
    memcpy(fim_lista,elemento,sizeof(token));
    fim_lista->proximo_token = NULL;
    
    //se a lista estiver vazia, retornar o elemento criado
    if (*lista == NULL)
    {
        *lista = fim_lista;
        return(1); //sucesso
    }
    else 
    {
        //percorre a lista até o final
        while(p_lista->proximo_token != NULL)
        {
            p_lista = p_lista->proximo_token;
        }
        
        //aponta o ponteiro anterior do novo elemento para o final da lista, e vice-versa
        p_lista->proximo_token = fim_lista;
        
        return (1); //sucesso
    }
    
}

//prioridade_operador -> funcao que retorna a precedencia do oeprador, sendo que um numero menor signigfica maior prioridade
int prioridade_operador(token *elemento)
{
    //se tentar comparar com um token inexistente
    if(elemento == NULL)
        return(99);
    //parenteses sempre tem prioridade maior
    if (elemento->tipo == tipo_abre_fecha)
        return (99);
    else
        //caso sejam operadores
        switch(elemento->operador)
        {
            case operador_mais:           
            case operador_menos:  
                return (4);
            case operador_multiplicacao: 
                return (3);
            case operador_negacao:
                return(2);
            case operador_potenciacao:
                return(1);
        }
    return(0);
}

//destroi lista - função para destruir todos os ponteiros da lista em caso de erro ou fim do programa
void destroi_lista_expr( token *lista)
{
    //percorre a lista recursivamente atÈ encontrar o ultimo elemento
    if(lista->proximo_token != NULL)
    {
        destroi_lista_expr(lista->proximo_token);
    }
    lista->literal = NULL;
    lista->proximo_token = NULL;
    free(lista);
}

  
//funcao que constroi a arvore binaria de expressoes a partir da expressao RPN gerada
arvore_expr *constroi_arvore_expr(token *lista)
{
    token *p_lista = lista;
    pilha_arvore *pilha = NULL;
    arvore_expr *temp_node = NULL;
    //percorre a lista da expressao RPN construindo a pilha de árvore
    while(p_lista != NULL)
    {
        switch(p_lista->tipo)
        {
            case tipo_literal:
            case tipo_parametro:
                //para literais ou parametros, os tokens sao inseridos diretamente na pilha
                insere_pilha_arvore(&pilha, cria_no_arvore(p_lista));
                break;
            case tipo_operador:
                //para operadores, retirar os ultimos 2 nos de arvore, criar uma arvore a partir deles e re-inserir na pilha
                
                //construir o node do operador
                temp_node = cria_no_arvore(p_lista);
                
                //unario: Para operacoes negacao unaria, preencho apenas o nó da direita
                if (p_lista->operador == operador_negacao)
                {
                    temp_node->direita = retira_pilha_arvore(&pilha);
                    temp_node->esquerda = NULL;
                }
                else
                {
                    //retira os ultimos 2 elementos da pilha para construir a expressao do no
                    temp_node->direita = retira_pilha_arvore(&pilha);
                    temp_node->esquerda = retira_pilha_arvore(&pilha);
                }
                
                //re-insere o no na pilha
                insere_pilha_arvore(&pilha, temp_node);
                
                //reseta o ponteiro de temp_node
                temp_node = NULL;
                
                break;
        }
        p_lista = p_lista->proximo_token;
    }
    
    //ao final do processo, devera sobrar apenas um elemento na pilha, contnedo a arvore completa da expressao
    return(retira_pilha_arvore(&pilha));
}

arvore_expr *cria_no_arvore(token *elemento)
{
    arvore_expr *novo_no;
    //aloca memoria para o proximo elemento
    novo_no = (arvore_expr*)malloc(sizeof(arvore_expr));
    if(novo_no == NULL)
        return(0); //houve erro de falta de memoria
    //aponta o no da pilha para o elemento e inicializa o resto da estrutura node
    novo_no->elemento = elemento;
    novo_no->direita = NULL;
    novo_no->esquerda = NULL;
    
    return(novo_no);
}

//inserir um elemento na pilha
int insere_pilha_arvore(pilha_arvore **pilha, arvore_expr *elemento)
{
    pilha_arvore *topo_pilha = NULL;
    
    topo_pilha = (pilha_arvore *)malloc(sizeof(pilha_arvore));
    if (topo_pilha == NULL)
        return (0); //erro
    //cria um novo elemento na pilha
    
    //posiciona o novo elemento no topo da pilha
    topo_pilha->node = elemento;
    topo_pilha->proximo = *pilha;
    *pilha = topo_pilha;
    return (1); //sucesso
    
}

//retirar um elemento da pilha
arvore_expr *retira_pilha_arvore (pilha_arvore **pilha)
{
    arvore_expr *elemento;
    pilha_arvore *topo_pilha;
    
    if (*pilha != NULL) 
    {
        //guarda o token que será retornado
        elemento = (*pilha)->node;
        
        //guarda o novo topo da pilha
        topo_pilha = (*pilha)->proximo;
        
        //desaloca memoria da estrutura da pilha
        free(*pilha);
        
        //atualiza o topo da pilha
        *pilha = topo_pilha;
        
        //retorna o elemento
        return(elemento);
    }
    else
        return (NULL);
}

//funcao que desaloca a memoria da arvore 
void destroi_arvore_expr(arvore_expr *arvore)
{
    
    if (arvore!=NULL)
    {
        //percorre sub-arvore esquerda
        destroi_arvore_expr(arvore->esquerda);
        //percorre sub-arvore direita
        destroi_arvore_expr(arvore->direita);
        
        //apos ter desalocado todas as sub-arvores, desaloca o proprio nó
        free((void*)arvore);
    }
}
        
        
        

//funcao que imprime a expressao da arvore
void imprime_arvore_expr(arvore_expr *arvore)
{
    if (arvore!=NULL)
    {
        if(arvore->elemento->tipo == tipo_operador)
        {
            //cada sub-expressao estará contida em parenteses
            printf("( ");
        }
        //percorre a sub-arvore da esquerda
        imprime_arvore_expr(arvore->esquerda);
        //imprime o caracter do token
        switch(arvore->elemento->tipo)
        {
            case tipo_literal:
                printf("%s", arvore->elemento->literal);
                break;
                
            case tipo_parametro:
                printf("%d", arvore->elemento->parametro);
                break;
                
            case tipo_operador:
                printf(" %c ", arvore->elemento->codigo);
                break;
        }
        //percorre a sub-arvore da direita
        imprime_arvore_expr(arvore->direita);
         if(arvore->elemento->tipo == tipo_operador)
        {
            //cada sub-expressao estará contida em parenteses
            printf(" )");
        }
    }
}


//função que reordena um astring em ordem alfabética
void string_sort(char **input)
{
	char *p_input;
	char swap;
	int flag = 0;	
	//o loop só para quando não houve reordenação alguma
	do
	{	
		//aponta para o inicio da string
		p_input = *input;
        
		//reset na flag de controle do loop
		flag = 0;
		
		//percorre a string até chegar ao final
		while(*(p_input + 1) != '\0')
		{

			//compara termos adjacentes
			if ((*p_input) > *(p_input+1))
			{
				//realiza a troca entre os termos adjacentes, caso o anterior seja maior
				swap = *p_input;
				*p_input = *(p_input+1);
				*(p_input + 1) = swap;
				
				//sinaliza que houve troca
				flag = 1;
			}
            //incrementa o ponteiro para proxima comparacao
            p_input++;
		}
	} while(flag);
    
}

//funcao que imprime a expressao da arvore
lista_expr *constroi_lista_expressoes_exp(arvore_expr *arvore)
{
    lista_expr *no_esq; //toda a expressão representada pelo no esquerdo
    lista_expr *no_dir; //toda a expressao representada pelo no direito
    lista_expr *elemento_lista = NULL; //elemento novo a ser alocado
    lista_expr *elemento_ant = NULL; //referencia para o ultimo elemento criado
    lista_expr *p_elemento_dir; //ponteiro para percorrer a lista da sub-arvore direita
    int i;  //indice de proposito geral
    
    if (arvore!=NULL)
    {
        //percorre a sub-arvore da esquerda
        no_esq = constroi_lista_expressoes_exp(arvore->esquerda);
        //percorre a sub-arvore da direita
        no_dir = constroi_lista_expressoes_exp(arvore->direita);
        //Realiza operacoes dependendo do tipo do nó
        switch(arvore->elemento->tipo)
        {
            //se chegou até um literal ou um parametro, é uma folha da arvore, portanto devemos alocar o elemento da lista de expressoes expandidas
            case tipo_literal:
                if((elemento_lista = (lista_expr *)malloc(sizeof(lista_expr))) == NULL)
                {
                    //erro memoria insuficiente
                    return(NULL);
                }
                //aloca apenas 2 posicoes, uma para a variavel e outra para o \0
                if((elemento_lista->codigos_numerador = (char *)malloc(2*sizeof(char))) == NULL)
                {
                    //erro memoria insuficiente
                    return(NULL);
                }
                //copia o codigo e o \0 na string
                elemento_lista->codigos_numerador[0] = arvore->elemento->codigo;
                elemento_lista->codigos_numerador[1] = '\0';
                elemento_lista->codigos_denominador = NULL;
                //por padrao, todos os termos sao positivos
                //elemento_lista->sinal = operador_mais; 
                //a variavel inicia-se multiplicada por 1
                elemento_lista->parametro = 1.0;	
                elemento_lista->proximo = NULL;
                elemento_lista->anterior = NULL;
                
                break;
 
            case tipo_parametro:
                if((elemento_lista = (lista_expr *)malloc(sizeof(lista_expr))) == NULL)
                {
                    //erro memoria insuficiente
                    return(NULL);
                }
                //escalares nao tem nenhuma lista de codigos associado
                elemento_lista->codigos_numerador = NULL;	
                elemento_lista->codigos_denominador = NULL;
                //por padrao, todos os termos sao positivos
               // elemento_lista->sinal = operador_mais; 
                //insere-se o parametro já transformado para inteiro
                elemento_lista->parametro = (double)(arvore->elemento->parametro);	
                elemento_lista->proximo = NULL;
                elemento_lista->anterior = NULL;
                
                break;
            
            //aqui estará a maior complexidade do codigo, onde as expansoes serao de fato realizadas
            case tipo_operador:
                switch(arvore->elemento->operador)
                {   
                    case operador_mais:  
                        //o operador positivo apenas concatena a lista do no esquerdo com a lista do no direito
                        //se o nó esquerdo for nulo, significa que há apenas uma operação unaria, então não deve haver concatenação
                        //aponta para o final da lista do no esquerdo
                        elemento_lista = no_esq;
                        while (elemento_lista->proximo != NULL) 
                            elemento_lista = elemento_lista->proximo;
                        
                        //concatenação das expressões
                        elemento_lista->proximo = no_dir;
                        no_dir->anterior = elemento_lista;
                        
                        //recupera o inicio da lista
                        while (elemento_lista->anterior != NULL) 
                            elemento_lista = elemento_lista->anterior;

                        break;
                        
                    case operador_menos:  
                        //o operador negativo troca todos os sinais dos elementos da sub-arvore direita e depois concatena como na adição
                        p_elemento_dir = no_dir;
                        while(p_elemento_dir != NULL)
                        {
                            /*if(p_elemento_dir->sinal == operador_mais)
                                p_elemento_dir->sinal = operador_menos;
                            else 
                                p_elemento_dir->sinal = operador_mais;*/
                            p_elemento_dir->parametro *= -1.0;
                            
                            p_elemento_dir = p_elemento_dir->proximo;
                        }
                        
                        //aponta para o final da lista do no esquerdo
                        elemento_lista = no_esq;
                        while (elemento_lista->proximo != NULL) 
                            elemento_lista = elemento_lista->proximo;
                        
                        //concatenação das expressões
                        elemento_lista->proximo = no_dir;
                        no_dir->anterior = elemento_lista;
                        
                        //recupera o inicio da lista
                        while (elemento_lista->anterior != NULL) 
                            elemento_lista = elemento_lista->anterior;
                        
                        break;
                        
                    case operador_negacao:
                        //o operador negação apenas troca todos os sinais dos elementos da sub-arvore direita
                        p_elemento_dir = no_dir;
                        while(p_elemento_dir != NULL)
                        {
                           /* if(p_elemento_dir->sinal == operador_mais)
                                p_elemento_dir->sinal = operador_menos;
                            else 
                                p_elemento_dir->sinal = operador_mais;*/
                            p_elemento_dir->parametro *= -1.0;
                            
                            p_elemento_dir = p_elemento_dir->proximo;
                        }
                        
                        //como nao há operando na nó esquerdo, não há o que concatenar
                        elemento_lista = no_dir;
                        
                        break;
                        
                    case operador_multiplicacao: 
                        
                        // o procedimento foi encapsulado para ser utilizado na potenciação
                        elemento_lista = multiplica_expr(no_esq,no_dir);
                        
                        //desalocar as listas originais
                        destroi_lista_expr_expandida(no_dir);
                        destroi_lista_expr_expandida(no_esq);
                        
                        //erro
                        if (elemento_lista == NULL) 
                        {
                            //erro
                            return(NULL);
                        }
                        break;
                        
                    case operador_potenciacao:
                        //primeiramente verificar se o expoente é um numero inteiro
                        if (no_dir->codigos_numerador != NULL) 
                        {
                            //erro - expressão no expoente
                            return(NULL);
                        }
                        
                        //TODO divisao: tambem retornar erro se o expoente for negativo
                       // if (no_dir->sinal == operador_menos) 
                        if (no_dir->parametro < 0.0) 
                        {
                            //erro - expoente negativo
                            return(NULL);
                        }
                        
                        // se o expoente for 0, retornar 1
                        if (no_dir->parametro == 0.0) 
                        {
                            if((elemento_lista = (lista_expr *)malloc(sizeof(lista_expr))) == NULL)
                            {
                                //erro memoria insuficiente
                                return(NULL);
                            }
                            //escalares nao tem nenhuma lista de codigos associado
                            elemento_lista->codigos_numerador = NULL;	
                            elemento_lista->codigos_denominador = NULL;
                            //por padrao, todos os termos sao positivos
                            //elemento_lista->sinal = operador_mais; 
                            //insere-se o numero 1 no campo do parametro
                            elemento_lista->parametro = 1.0;	
                            elemento_lista->proximo = NULL;
                            elemento_lista->anterior = NULL;
                            
                            //desalocar as listas do no direito e esquerdo
                            destroi_lista_expr_expandida(no_dir);
                            destroi_lista_expr_expandida(no_esq);

                        } //se o expoente for 1, retornar o proprio no esquerdo
                        else if (no_dir->parametro == 1.0) 
                        {
                            elemento_lista = no_esq;
                            
                            //desalocar o nó direito
                            destroi_lista_expr_expandida(no_dir);
                            
                        }
                        else //multiplicar o no esquerdo por si mesmo quantas vezes for o parametro do no direito
                        {
                            elemento_ant = no_esq;
                            for (i=0;i < (int)(no_dir->parametro - 1.0);i++)
                            {
                                elemento_lista = multiplica_expr(elemento_ant, no_esq);
                                
                                //desaloca o elemento anterior, mas apenas se não for o no_Esq original
                                if (elemento_ant != no_esq) 
                                    destroi_lista_expr_expandida(elemento_ant);
                                
                                //aponta o elemento anterior para o atual
                                elemento_ant = elemento_lista;
                            }
                        }
                    
                        //desaloca as listas originais
                        destroi_lista_expr_expandida(no_dir);
                        destroi_lista_expr_expandida(no_esq);
                        break;
                }   
                break;
        }


    }
    return(elemento_lista);
}

//realiza multiplicacao entre duas listas de expressao
lista_expr *multiplica_expr( lista_expr *no_esq, lista_expr *no_dir)
{

    lista_expr *elemento = NULL; //elemento novo a ser alocado
    lista_expr *elemento_ant = NULL; //referencia para o ultimo elemento criado
    lista_expr *p_elemento_esq; //ponteiro para percorrer a lista da sub-arvore esquerda
    lista_expr *p_elemento_dir; //ponteiro para percorrer a lista da sub-arvore direita
    int temp;
    
    
    
    //aqui deve-se combinar todos os elementos da lista esquerda com os da lista direita, gerando uma nova lista
    //no processo que deverá ser desalocada ao final da operação.
    p_elemento_esq = no_esq;
    p_elemento_dir = no_dir;
    
    //caso algum dos termos seja um ponteiro nulo
    if (no_esq == NULL || no_dir == NULL)
        return NULL;

    while (p_elemento_esq != NULL) 
    {
        while (p_elemento_dir != NULL) 
        {
            //reseta o ponteiro elemento para a etapa
            elemento = NULL;
            //aloco cada elemento novo
            if((elemento = (lista_expr *)malloc(sizeof(lista_expr))) == NULL)
            {
                //erro memoria insuficiente
                return(NULL);
            } 
            
            //inicializar ponteiros
            elemento->proximo = NULL;
            elemento->anterior = NULL;
            
            //TODO: DIVISAO implementar as operacoes cabiveis com os denominadores também
            //vejo a quantidade de variaveis em cada termo da multiplicacao e aloco uma string adequada
            if ((p_elemento_esq->codigos_numerador != NULL) || (p_elemento_dir->codigos_numerador != NULL)) 
            { 
                temp = 0; //acumulador de tamanho de string
                if (p_elemento_esq->codigos_numerador != NULL) 
                {
                    temp += strlen(p_elemento_esq->codigos_numerador);
                }
                
                if (p_elemento_dir->codigos_numerador != NULL) 
                {
                    temp += strlen(p_elemento_dir->codigos_numerador);
                }
                
                elemento->codigos_numerador = (char *)malloc(sizeof(char)*(1 + temp));
                *(elemento->codigos_numerador) = '\0';
                elemento->codigos_denominador = NULL;
                //copio a string das variaveis do primeiro elemento
                if (p_elemento_esq->codigos_numerador != NULL) 
                    strcpy(elemento->codigos_numerador,p_elemento_esq->codigos_numerador);
                //e concateno com os do segundo
                if (p_elemento_dir->codigos_numerador != NULL) 
                    strcat(elemento->codigos_numerador,p_elemento_dir->codigos_numerador);
                //finalmente, reordeno a nova string
                string_sort(&(elemento->codigos_numerador));
            }
            else //multiplicacao entre 2 parametros
            {
                elemento->codigos_numerador = NULL;
            }   
        
            //inicializar o ponteiro do denominador
            elemento->codigos_denominador = NULL;
        
            //definir o sinal do elemento (equivale a um xor)
           /* if(p_elemento_esq->sinal == p_elemento_dir->sinal)
                elemento->sinal = operador_mais;
            else 
                elemento->sinal = operador_menos;*/
                
            //multiplicar os parametros
            elemento->parametro = p_elemento_esq->parametro * p_elemento_dir->parametro;
                
            //atualizar os ponteiros
            if(elemento_ant == NULL)
                elemento_ant = elemento;
            else 
            {
                elemento_ant->proximo = elemento;
                elemento->anterior = elemento_ant;
                elemento_ant = elemento;
            }
        
            //incrementa o ponteiro da lista do nó da direita
            p_elemento_dir = p_elemento_dir->proximo;
        }
    
        //avanço o ponteiro da lista do nó esquerdo
        p_elemento_esq = p_elemento_esq->proximo;
    
        //reseto o ponteiro da sub-arvore esquerda
        p_elemento_dir = no_dir;
    }

    //aproveitando que a lista é duplamente encadeada, eu posso recuperar o inicio atraves dos ponteiros para anterior
    while (elemento->anterior != NULL)
        elemento = elemento->anterior;
    
    return (elemento);
}

//desalocar um lista de expressoes
void destroi_lista_expr_expandida(lista_expr *lista)
{
    //teste de consistencia
    if (lista == NULL) 
        return;
    
    //chamadas recursivas
    if (lista->proximo != NULL)
       destroi_lista_expr_expandida(lista->proximo);
    
    //desalocar cada estrutura da lista
    if (lista->codigos_denominador != NULL)
        destroi_lista_expr_expandida(lista->codigos_denominador);
    if (lista->codigos_numerador != NULL)
        free(lista->codigos_numerador);
    //desalocar o nó
    free(lista);
    
}

void imprime_lista_expr_expandida(lista_expr *lista, tabela_literais *tabela)
{
    lista_expr   *p_lista = lista;
    double teste;
    
    //imprime todos os campos do elemento
    while(p_lista != NULL)
    {
        //imprime o sinal
        if (p_lista->parametro > 0.0)
            printf("+");
        else if (p_lista->parametro == -1.0)
            printf("-");
        
        //imprime o parametro
        if (p_lista->parametro != 1.0 && p_lista->parametro != -1.0)
            //checar se o coeficiente é inteiro
            if ((teste = (int)p_lista->parametro - p_lista->parametro) == 0.0)
                printf("%d",(int)(p_lista->parametro));
            else
                printf("%2.2f",p_lista->parametro);

        
        //imprime o numerador
        if (p_lista->codigos_numerador != NULL)
            //printf("%s", p_lista->codigos_numerador);
            print_monomio(p_lista->codigos_numerador, tabela);
        else  if (p_lista->parametro == 1.0 || p_lista->parametro == -1.0) 
            printf("1.00");
        
        //imprime o denominador
        if (p_lista->codigos_denominador != NULL)
            //printf("/(%s)",p_lista->codigos_denominador);
            imprime_lista_expr_expandida(p_lista->codigos_denominador, tabela);
        
        //imprime um espaço para o proximo elemento
        printf(" ");
        
        //imprime o proximo elemento
        p_lista = p_lista->proximo;
    }
}

//funcão que imprime o monomio segundo suas variáveis a partir dos codigos internos
void print_monomio(char *monomio, tabela_literais *tabela)
{
    int i;
    tabela_literais *p_tabela;
    
    //percorrer os codigos dos monomios, contando quantos tem igual, e imprimindo com expoente quando for o caso
    while (*monomio != '\0') 
    {
        i = 1; //primeira variavel
        //testar se os proximos codigos são iguais
        while (*monomio == *(monomio + 1))
        {
            i++; //incrementa o numero de variaveis iguais
            monomio++;
        }
        //inicializa a tabela de literias
        p_tabela = tabela;
        //Pesquisa o código na tabela de literais
        while (p_tabela->codigo != *monomio) 
            p_tabela = p_tabela->proximo_codigo;
        
        //imprime a variável correspondente ao código
        printf("%s",p_tabela->literal);
        
        //imprime o expoente
        if (i > 1)
            printf("^%d",i);
        //proximo codigo
        monomio++;
        
        //adicionar um sinal de multiplicação caso o proximo nao seja o \0
        if (*monomio != '\0') 
            printf("*");
    }
    
}


// função que agrupa elementos iguais da expressão expandida
lista_expr *simplifica_expr_expandida(lista_expr *p_lista)
{
    lista_expr *p_1, *p_2, *p_remove, *p_inicio;  //ponteiros para percorrer a lista
    double resultado;          //avaliação da soma ou subtração dos termos
    
    //realiza uma copia de p_lista
    //TODO: controle de erro
    if (p_lista == NULL)
        return NULL;
    
    p_inicio = copia_lista_expr(p_lista);
    
    //inicializa os ponteiros
    p_1 = p_inicio;
    
    if (p_1->proximo == NULL)
    {
        //a expressão só tem um elemento
        return(p_1);
    }
    
    //loop de comparação
    while (p_1 != NULL) 
    {
        //inicializar p_2 como o segundo elemento
        p_2 = p_1->proximo;
        while (p_2 != NULL) 
        {
            //verificar se o conjunto de variaveis dos dois elementos são iguais
            // TODO divisão: levar em conta o denominador também
            if (p_1->codigos_numerador != NULL && p_2->codigos_numerador != NULL) 
            {
                if (!strcmp(p_1->codigos_numerador, p_2->codigos_numerador)) 
                {
                    //inicializar a variavel de resultado
                    resultado = 0.0;
                    //somar ou subtrair o primeiro termo
                    /* if (p_1->sinal == operador_mais)
                     resultado+= p_1->parametro;
                     else
                     resultado-= p_1->parametro;*/
                    
                    resultado+= p_1->parametro;
                    //somar ou subtrair o segundo termo 
                    /*if (p_2->sinal == operador_mais)
                     resultado+= p_2->parametro;
                     else
                     resultado-= p_2->parametro;*/
                    
                    resultado+= p_2->parametro;
                    
                    //BUG: arredondar quando o resultado é residual.Por exemplo, 0.9 - 0.9 = 11.1E-15
                    
                    if ((resultado <= 0.0000000000001 && resultado >= 0.0) || (resultado >= -0.0000000000001 && resultado <= 0.0)) 
                        resultado = 0.0;
                    //guardar o resultado em p_1 sem sinal
                    p_1->parametro = resultado;
                    
                    //avaliar o sinal
                    /* if (resultado >= 0)
                     p_1->sinal = operador_mais;
                     else
                     p_1->sinal = operador_menos;*/
                    
                    //ja atualiza p_2
                    p_remove = p_2;
                    p_2 = p_2->proximo;
                    //remover p_2 da lista
                    //TODO: implementar a função
                    p_inicio = remove_lista_expr(p_inicio, p_remove);
                }
                else //atualiza p_2
                    p_2 = p_2->proximo;
            }
            else //atualiza p_2
                p_2 = p_2->proximo;
            
        }
        //atualiza p_1
        p_1 = p_1->proximo;
    }
    //percorre a lista novamente elimienando todos os resultados 0, deixando somente um caso não sobre elementos na lista
    p_1 = p_inicio;
    while (p_1 != NULL) 
    {
        //se sobrou apenas um elemento (p_1 == p_inicio) e  (p_1->proximo == NULL) então não descartá-lo
        if ((p_1->parametro == 0.0) && (p_1 == p_inicio) && (p_1->proximo == NULL))
        {
            //desalocar os codigos se houver algum
            if(p_1->codigos_numerador != NULL)
            {
                free(p_1->codigos_numerador);
                p_1->codigos_numerador = NULL;
            }
            if(p_1->codigos_denominador != NULL)
            {
                //free(p_1->codigos_denominador);
                destroi_lista_expr_expandida(p_1->codigos_denominador);
                p_1->codigos_denominador = NULL;
            }
            
            p_1 = p_1->proximo;
 
        }
        //caso o parametro seja 0 em outras circunstancias
        else if (p_1->parametro == 0.0)
        {
            //atualiza p_1
            p_remove = p_1;
            p_1 = p_1->proximo;
            //remove p_1 da lista
            p_inicio = remove_lista_expr(p_inicio, p_remove);
        }
        else
            p_1 = p_1->proximo; //nao é zero
    }
    
    return p_inicio;
}
  
//funcao que remove um ponteiro da lista de expressoes expandida
lista_expr *remove_lista_expr(lista_expr *p_inicio, lista_expr *p_remove)
{
    lista_expr *p_lista; //ponteiro para percorrer a lista
    
    //inicializa p_lista
    p_lista = p_inicio;
    
    //percorre a lista a procura de p_remove
    while (p_lista != p_remove && p_lista != NULL) 
        p_lista = p_lista->proximo;
    
    //se não encontrou o elemento, retorna simplesmente o inicio da lista denovo
    if (p_lista == NULL) 
    {
        return p_inicio;
    }
    
    //se o ponteiro a ser removido esta no inicio da lista atualizar p_inicio antes da remocao
    if (p_lista == p_inicio) 
    {
        p_inicio = p_inicio->proximo;
    }
    
    //atualizar o pontiero proximo do elemento anterior a p_remove
    if (p_remove->anterior != NULL)
        (p_remove->anterior)->proximo = p_remove->proximo;
    
    //atualizar o ponteiro anterior do elemento proximo a p_remove
    if (p_remove->proximo != NULL)
        (p_remove->proximo)->anterior = p_remove->anterior;
    
    //liberar memoria dos elementos internos de p_remove
    if (p_remove->codigos_numerador != NULL) 
        free(p_remove->codigos_numerador);
    
    if (p_remove->codigos_denominador != NULL) 
       // free(p_remove->codigos_denominador);
        destroi_lista_expr_expandida(p_remove->codigos_denominador);
    
    //desalocar p_remove
    free(p_remove);
    
    //retornar o inicio da lista
    return p_inicio;
}

lista_expr *copia_lista_expr(lista_expr *lista)
{
    lista_expr *p_lista_prox = NULL, *p_lista;
    
    //chamadas recursivas
    if (lista->proximo != NULL)
        p_lista_prox = copia_lista_expr(lista->proximo);
    
    //alocar memoria para o nó:
    //TODO controle de erro
    p_lista = (lista_expr *)malloc(sizeof(lista_expr));
    
    //inicializa os elementos
    p_lista->anterior = NULL;
    p_lista->proximo = NULL;
    p_lista->codigos_denominador = NULL;
    
    //copia os elementos
    p_lista->parametro = lista->parametro;
    //p_lista->sinal = lista->sinal;
    if (lista->codigos_numerador != NULL) 
    {
        //TODO controle de erro
        p_lista->codigos_numerador = (char *)malloc((strlen(lista->codigos_numerador)+1)*sizeof(char));
        strcpy(p_lista->codigos_numerador, lista->codigos_numerador);
    }
    else
        p_lista->codigos_numerador = NULL;
    
    if (lista->codigos_denominador != NULL) 
    {
        //TODO controle de erro
        //p_lista->codigos_denominador = (char *)malloc((strlen(lista->codigos_denominador)+1)*sizeof(char));
        //strcpy(p_lista->codigos_denominador, lista->codigos_denominador);
        p_lista->codigos_denominador = copia_lista_expr(lista->codigos_denominador);
    }
    else
        p_lista->codigos_denominador = NULL;
    
    //atualiza os ponteiros de anterior e proximo
    p_lista->proximo = p_lista_prox;
    if (p_lista_prox != NULL) 
        p_lista_prox->anterior = p_lista;
    
    return p_lista;
    
}

//reordena um polinomio na forma lexdeg
lista_expr *lexdegbubblesort(lista_expr *poly)
{
    lista_expr *p_atual;    //ponteiro para elemento atual do polinomio
    lista_expr *p_troca;    //ponteiro auxiliar para troca 
    int flag_troca = 0;     //flag que indica se houve alguma troca de elementos
    
    if (poly == NULL)
        return NULL;
    //percorre a lista de polinomios e troca os elementos 1 a 1 quando for necessario. Só para quando a percorrer o polinomio inteiro e não fizer troca alguma
    do
    {
        //reseta o flag de flag_troca
        flag_troca = 0;
        //aponta para o inicio do polinomio
        p_atual = poly;
        //rebobina para o inicio do polinomio
        while (p_atual->anterior != NULL) 
            p_atual = p_atual->anterior;
        //percorre a lista
        while (p_atual->proximo != NULL) 
        {
            p_troca = p_atual->proximo;
            
            //compara dois monomios adjacentes. Se o proximo tem precedencia lexdeg, é feita uma troca entre os monomios na lista do polinomio
            if (lexdeg(p_atual->codigos_numerador,p_troca->codigos_numerador))
            {
                //sinaliza que houve troca
                flag_troca = 1;
                //realiza a troca
                //ajusta os ponteiros de borda
                p_atual->proximo = p_troca->proximo;
                if (p_atual->proximo != NULL)                     
                    p_atual->proximo->anterior = p_atual;
                
                p_troca->anterior = p_atual->anterior;
                if(p_troca->anterior != NULL)
                    p_troca->anterior->proximo = p_troca;
                
                //ajusta os ponteiros entre p_atual e p_troca
                p_troca->proximo = p_atual;
                p_atual->anterior = p_troca;
            }
            else
                p_atual = p_atual->proximo;
        }
    }while (flag_troca == 1); 
    
    //recuperar o inicio do polinomio
    while (poly->anterior != NULL) 
    {
        poly = poly->anterior;
    }
    
    return poly;
}

//lexdeg função que compara dois monomios e retorna 1 caso o segundo argumento tenha precedencia lexdeg em relação ao primeiro
int lexdeg(char *primeiro, char *segundo)
{
    char parametro = 'A'; //parametro de comparacao lexica
    char *p_deg;    //ponteiros para comparação de grau monomial
    int deg_1, deg_2;           //grau dos monomios
    
    //se o segundo monomio for uma constante, retornar 0
    if (segundo == NULL) 
        return 0;
    //se o primeiro for uma constante, retornar 1
    else if (primeiro == NULL)
        return 1;
    
    //loop infinito de procura, o retorno será efetuado dentro das condicionais
    while (1) 
    {
        //procura o parametro no segundo monomio
        if ((p_deg = strchr(segundo, parametro)) != NULL) 
        {
            //procuro o parametro no primeiro
            if ((p_deg = strchr(primeiro, parametro)) != NULL) 
            {
                //se o argumento prioritário foi encontrado em ambos monomios, comparar o grau dos monomios para o mesmo parametro
                //encontrar o grau dos monomios
                p_deg = primeiro;   //inicializa o ponteiro
                //inicializa o grau
                deg_1 = 0;
                while (p_deg!= NULL) 
                {
                    //a cada vez que encontra o parametro, refazer a busca a partir do proximo caracter
                    p_deg = strchr(p_deg, parametro);
                    if (p_deg != NULL)
                    {
                        deg_1 += 1;
                        //atualiza o ponteiro da string
                        p_deg+=1;
                    }
                }
                
                //repetir o procedimento para encontrar o grau do segundo monomio para o presente parametro
                p_deg = segundo;   //inicializa o ponteiro
                //inicializa o grau
                deg_2 = 0;
                while (p_deg!= NULL) 
                {
                    //a cada vez que encontra o parametro, refazer a busca a partir do proximo caracter
                    p_deg = strchr(p_deg, parametro);
                    if (p_deg != NULL)
                    {
                        deg_2 += 1;
                        //atualiza o ponteiro da string
                        p_deg+=1;
                    }
                }
                
                //de posse dos graus dos monomios, comparar
                if (deg_2 > deg_1) 
                {
                    //no caso do grau do segundo monomio ser maior, retornar 1
                    return 1;
                }
                else if (deg_2 < deg_1)
                {
                    //no caso do grau do segundo monomio ser menor que o primeiro, retornar 0
                    return 0;
                }
                
                //caso contrario, se o grau for identico para a mesma variavel em ambos monomios, não faço mais nada e atualizo o parametro 
                //para o proximo loop
                
            }
            //se o argumento prioritario foi encontrado apenas no segundo monomio, retornar 1
            else
                return 1;
            
        }
        else if ((p_deg = strchr(primeiro, parametro)) != NULL) 
        {
            //se encontrou o parametro prioritario apenas no primeiro monomio, retornar 0.
            return 0;
        }
        
        //caso nao tenha encontrado o parametro prioriatario em nenhum monomio, ou caso o grau da variavel relatica ao parametro atual for igual
        //em ambos monomios, refazer a comparação com o próximo parametro.
        
        //se já está na letra Z, passar a busca para "a" minusculo
        if (parametro == 'Z') 
        {
            parametro = 'a';
        }
        else
            parametro += 1;
        
        //caso nao tenha encontrado variavel alguma valida, retornar 0
        if (parametro > 'z') 
        {
            return(0);
            //erro 
        }
        
    }
    
}
//funcao que realiza a divisao polinomial.
int polydiv(lista_expr *dividendo, lista_expr *divisor, lista_expr **quociente, lista_expr **resto)
{
    lista_expr *dividendo_cpy = NULL;
    lista_expr *resultado_monomio = NULL;
    lista_expr *poly_ptr;
    lista_expr *resto_temp;
    lista_expr *quociente_temp;
    
    //inicializa variaveis de retorno
    *quociente = NULL;
    *resto = NULL;
    
    //se passar um dividendo 0 ou inexistente, devo retornar imediatamente
    if (dividendo == NULL) 
    {
        *quociente = constroi_elemento_zerado();
        *resto = constroi_elemento_zerado();
        
        return 1;
    } //se o dividendo possui apenas um monomio de valor 0
    else if (dividendo->proximo == NULL && dividendo->parametro == 0)
    {
        *quociente = constroi_elemento_zerado();
        *resto = constroi_elemento_zerado();
        
        return 1;
    }

    //inicializar resto e quociente
    quociente_temp = NULL;
    resto_temp = NULL;
    //devo copiar o dividendo a cada operação, visto que ela pode falhar e eu ter que recomeçar
    dividendo_cpy = copia_lista_expr(dividendo); //este passo pode deixar o algoritmo bem pesado. Posso otimizar para que a lista vire um bloco de memoria, e
                                                 //realizar a copia via memcpy.
    
    //1) Compara-se o Lt(dividendo) com o Lt(divisor) -> Leading term.
    //2) Se a string do LT(divisor) está contida na Lt(dividendo), -> cuidado: AAABBBCCC = A^3*B^3*C^3, e ABC. A segunda divide a primeira, mas não há string ABC na primeira! encontro a string que falta e crio um termo novo.
    //o resultado_monomio é a string que falta no divisor para chegar ao dividendo. Se o divisor nao estiver totalmente contido no dividendo, ela retorna null, ou seja, o LT do dividendo já começa a integrar o resto.
    
    //repetir a procedimento até o dividendo ficar vazio
    while (dividendo_cpy != NULL) 
    {
        //divide os LT's do dividendo e divisor, guardando o resultado em resultado_monomio
        resultado_monomio = divide_monomio(dividendo_cpy,divisor);
        
        //condicao de erro
        if (resultado_monomio == NULL) 
        {
            //erro, divisao por 0
            //limpar variaveis
            destroi_lista_expr_expandida(dividendo_cpy);
            //se der merda aqui, tenho que zerar as variaveis de retorno
            *quociente = constroi_elemento_zerado();
            *resto = constroi_elemento_zerado();
            return 0;
        }
        else if (resultado_monomio->parametro == 0.0) // não houve divisão
        {
            //O LT do dividendo passa a integrar o resto
            if (resto_temp == NULL)
            {
                resto_temp = dividendo_cpy;
                poly_ptr = resto_temp; //manter a referencia para por NULL no final
            }
            else
            {
                //inserir o resto no final da lista do resto
                poly_ptr = resto_temp;
                while (poly_ptr->proximo != NULL)  
                    poly_ptr = poly_ptr->proximo;
                poly_ptr->proximo = dividendo_cpy;
                dividendo_cpy->anterior = poly_ptr;
                poly_ptr = poly_ptr->proximo; //manter a referencia para por NULL no final
            }
            //O LT passa a ser o proximo monomio
            dividendo_cpy = dividendo_cpy->proximo;
            if (dividendo_cpy != NULL) 
                dividendo_cpy->anterior = NULL;
            //colocar NULL no final do resto
            poly_ptr->proximo = NULL;
            
            //desalocar o resultado_monomio
            destroi_lista_expr_expandida(resultado_monomio);
            
        }
        else //caso houve divisão 
        {
            
            //eliminar o LT do dividendo, pois sabemos que ele será zerado na operação
            poly_ptr = dividendo_cpy;
            dividendo_cpy = dividendo_cpy->proximo;
            if (dividendo_cpy != NULL)
                dividendo_cpy->anterior = NULL;
            poly_ptr->proximo = NULL;
            destroi_lista_expr_expandida(poly_ptr);
            
            //multiplicar o monomio resultante pelo divisor sem o LT (pois este produto irá cancelar com o LT do dividendo) e subtrair do dividendo
            if (divisor->proximo != NULL) //se o divisor tiver apenas um monomio, pula-se esta etapa.
            {
                poly_ptr = multiplica_expr(resultado_monomio, divisor->proximo);
                
                //subtrai-se poly_ptr do dividendo
                subtrai_expr(&dividendo_cpy, poly_ptr);
            }
            
            //o monomio resultante integra o quociente
            if (quociente_temp == NULL) 
                quociente_temp = resultado_monomio;
            else
            {
                //inserir o monomio resultante ao final do quociente
                poly_ptr = quociente_temp;
                while (poly_ptr->proximo != NULL)
                    poly_ptr = poly_ptr->proximo;
                poly_ptr->proximo = resultado_monomio;
                resultado_monomio->anterior = poly_ptr;
                resultado_monomio->proximo = NULL;
            }
            
        }

    }
    
    //ao final da divisão, reordenar e simplificar o quociente e o resto
    if (quociente_temp != NULL) 
    {
        poly_ptr = simplifica_expr_expandida(quociente_temp);
        destroi_lista_expr_expandida(quociente_temp);
        quociente_temp = lexdegbubblesort(poly_ptr);
        
    }
    else 
    {
        //construir um elemento zerado
        quociente_temp = constroi_elemento_zerado();
        
    }
    if (resto_temp != NULL) 
    {
        poly_ptr = simplifica_expr_expandida(resto_temp);
        destroi_lista_expr_expandida(resto_temp);
        resto_temp = lexdegbubblesort(poly_ptr);
        
    }
    else 
    {
        //construir um elemento zerado
        resto_temp = constroi_elemento_zerado();        
    }
    
    *quociente = quociente_temp;
    *resto = resto_temp;

    return 1;
    
}


//função que subtrai duas expressões
void subtrai_expr(lista_expr **no_esq, lista_expr *no_dir)
{
    lista_expr *poly_ptr;
    
    //troca o sinal dos elementos do nó direito
    poly_ptr = no_dir;
    while(poly_ptr != NULL)
    {
       /* if(poly_ptr->sinal == operador_mais)
            poly_ptr->sinal = operador_menos;
        else 
            poly_ptr->sinal = operador_mais;*/
        
        poly_ptr->parametro *= -1.0;
        
        poly_ptr = poly_ptr->proximo;
    }
    
    //aponta para o final da lista do no esquerdo
    if (*no_esq != NULL) 
    {
        poly_ptr = *no_esq;
        while (poly_ptr->proximo != NULL) 
            poly_ptr = poly_ptr->proximo;
        
        //concatenação das expressões
        poly_ptr->proximo = no_dir;
        if (no_dir != NULL) 
        {
            no_dir->anterior = poly_ptr;
        }
        
    }
    else //o elemento do nó esquerdo é vazio
    {
        *no_esq = no_dir;
    }
    
    
    
}

//função que soma duas expressões
void soma_expr(lista_expr *no_esq, lista_expr *no_dir)
{
    lista_expr *poly_ptr;
    
    //aponta para o final da lista do no esquerdo
    poly_ptr = no_esq;
    while (poly_ptr->proximo != NULL) 
        poly_ptr = poly_ptr->proximo;
    
    //concatenação das expressões
    poly_ptr->proximo = no_dir;
    no_dir->anterior = poly_ptr;
}

//função que realiza a divisão entre LT's -> OTIMIZACAO
lista_expr *divide_monomio(lista_expr *dividendo, lista_expr *divisor)
{
    char *dividendo_ptr;
    char *divisor_ptr;  
    char *div_string = NULL;        //guardar a string do monomio resultante
    lista_expr *div_monomio = NULL; //guardar a estrutur do monomio resultante
    double div_parametros = 1.0;    //guardar o parametro resultante
    int tamanho_string =0 ;             //calcula o tamanho da string do monomio a ser alocada
    
    if (divisor == NULL)
        return NULL;
    //teste de erro
    if (divisor->parametro == 0.0) 
    {
        //erro divisao por 0
        return NULL;
    }
    //teste de apenas parametro
    if (dividendo->codigos_numerador == NULL)
    {
        //erro o monomio é apenas um parametro
        return NULL;
    }
    
    dividendo_ptr = dividendo->codigos_numerador;
    divisor_ptr = divisor->codigos_numerador;
    
    //divide os parametros
    div_parametros = dividendo->parametro / divisor->parametro;
    
    //testar se é uma simples divisão de parâmetros
    if (dividendo->codigos_numerador == NULL && divisor->codigos_numerador == NULL) 
    {
        //nao haverá string de parametros
        div_string = NULL;
    }
    else //há uma string no divisor ou no dividendo ou em ambos TODO: aqui pode dar pau com string nula
    {
        //calcula-se o tamanho da string a ser alocado
        if (dividendo->codigos_numerador != NULL) 
            tamanho_string = (int)strlen(dividendo->codigos_numerador);
        if (divisor->codigos_numerador != NULL) 
            tamanho_string -= (int)strlen(divisor->codigos_numerador);
        
        //se o divisor tiver grau maior que o dividendo, deve-se retornar 0
        if (tamanho_string <0)
        {
            div_parametros = 0.0;
            div_string = NULL;
        }
        //caso o dividendo e divisor tenham grau semelhante
        else if (tamanho_string == 0)
        {
            //se as strings forem diferentes, muda o resultado de div_parametros para 0
            if (strcmp(dividendo->codigos_numerador, divisor->codigos_numerador))
                div_parametros = 0.0;
            //caso contrario, mantém-se o div_parametros calculado anteriormente.
            //em qualquer hipotese, nao há string de variaveis aqui
            //div_string = NULL; -> já está inicializado assim
        }
        //string no dividendo maior que no divisor
        else
        {
            //caso o divisor seja apenas um parametro
            if (divisor->codigos_numerador == NULL) 
            {
                //copia a string do dividendo na string do monomio
                div_string = (char *)malloc((tamanho_string+1)*sizeof(char)); //TODO: tratar erro de alocação de memoria
                *div_string = '\0';
                strcpy(div_string, dividendo->codigos_numerador);
            }
            //há string no divisor, encontrar a string que a completa para que fique igual ao dividendo.
            else
            {
                //aloca a string do resultado
                div_string = (char *)malloc((tamanho_string+1)*sizeof(char)); //TODO: tratar erro de alocação de memoria
                *div_string = '\0';
                
                //busca os caracteres da string do divisor na string do dividendo
                while ((*dividendo_ptr) != '\0') 
                {
                    //caso não encontre o caracter no divisor, adicioná-lo ao resultado
                    if (*dividendo_ptr != *divisor_ptr) 
                    {
                        strncat(div_string,dividendo_ptr,1);
                        //incrementa ponteiro apenas do dividendo
                        dividendo_ptr++;
                    }
                    //caso os caracteres sejam iguais 
                    else
                    {
                        //incrementa ambos ponteiros
                        dividendo_ptr++;
                        divisor_ptr++;
                    }
                    
                }
                
                //ao final do processo, se o divisor não tiver sido esgotado, então os monomios não dividem
                if (*divisor_ptr != '\0')
                {
                    //liberar a memoria do resultado
                    free(div_string);
                    div_string = NULL;
                    //o parametro resultante deverá ser 0
                    div_parametros = 0.0;
                }
                
            }
        }

    }
    
    //com todas as possibilidades cobertas, basta alocar a memoria para o monomio resultante e preencher a estrutura
    div_monomio = (lista_expr *)malloc(sizeof(lista_expr)); // TODO: tratar erro
    
    div_monomio->proximo = NULL;
    div_monomio->anterior = NULL;
    div_monomio->codigos_denominador = NULL;
    div_monomio->codigos_numerador = div_string;
    div_monomio->parametro = div_parametros;
    
    
    //encontra o sinal
   /* if (dividendo->sinal != divisor->sinal)
        div_monomio->sinal = operador_menos;
    else
        div_monomio->sinal = operador_mais;*/
    
    return div_monomio;
        
}

lista_expr *constroi_elemento_zerado(void)
{
    lista_expr *elemento;
    //construir um elemento zerado
    elemento = (lista_expr *)malloc(sizeof(lista_expr)); //TODO: controle de erro
    elemento->codigos_numerador = NULL;
    elemento->codigos_denominador = NULL;
    elemento->proximo = NULL;
    elemento->anterior = NULL;
    //(*quociente)->sinal = operador_mais;
    elemento->parametro = 0.0;  
    
    return elemento;
}


/*

Frações parciais->  P3/(P1 * P2) =  Q + a/P1 + b/P2

- O denominador de cada elemento da lista_expr deverá ser um ponteiro para outra lista_expr.

etapas da expansao em fracao parcial:

*/
//retorna 0 se não encontrar uma expansão certinha, ou 1 se encontrar
int partial_fraction_expansion(lista_expr *P3, lista_expr *P1, lista_expr *P2, lista_expr **quociente,
                                lista_expr **numerador_a, lista_expr **numerador_b)
{
    lista_expr *lista_ptr1, *lista_ptr2, *lista_ptr3, *P3_resto, *resto;
    int achou = 0;  //controle de loop de busca
    char common_var = 0;    //variável em comum encontrada para a expansão
    double coeficiente = 0.0; //coeficiente da variavel comum 
    
    //inicialização de retornos
    *numerador_a = NULL;
    *numerador_b = NULL;
    *quociente = NULL;
    resto = NULL;
    
    //primeiro: dividir P3 por P1*P2 e guardar o Quociente

    //multiplica-se P1 e P2
    lista_ptr1 = multiplica_expr(P1, P2);
    
    //condição de erro
    if (lista_ptr1 == NULL) 
    {
        return 0;
    }
    
    //simplifica, reordena e elimina a lista desordenada
    lista_ptr2 = simplifica_expr_expandida(lista_ptr1);
    destroi_lista_expr_expandida(lista_ptr1);
    lista_ptr1 = lexdegbubblesort(lista_ptr2);
    
    //divide-se P3 por P1 *P2, guardando em lista_ptr3, o resto que vai gerar as frações parciais.
    lista_ptr3 = NULL;
    polydiv(P3, lista_ptr1, quociente, &lista_ptr3); //TODO: controle de erro
    
    //salva o resto para depois
    P3_resto = copia_lista_expr(lista_ptr3);
    
    //elimina-se P1*P2 expandido
    destroi_lista_expr_expandida(lista_ptr1);
    
    //1) escolhe-se uma variavel comum a P1 e P2
    //comparo cada variavel de cada monomio de p1 com todas as variaveis de P2, até encontrar.
    lista_ptr1 = P1;
    
    //compara cada elemento de P1 com todos de P2
    while (lista_ptr1 != NULL && !achou) 
    {
        //reinicializa lista_ptr2
        lista_ptr2 = P2;
        while (lista_ptr2 != NULL && !achou) 
        {   
            //partindo do principio que os polinomios geradores só possuem monomios de primeiro grau, devido
            // a serem formado por combinação linear entre as variáveis, a comparação será simplificada.
            if (*(lista_ptr1->codigos_numerador) == *(lista_ptr2->codigos_numerador)) 
            {
                achou = 1; //fim do loop
                common_var = *(lista_ptr1->codigos_numerador);
            }
            lista_ptr2 = lista_ptr2->proximo;
            
        }
        lista_ptr1 = lista_ptr1->proximo;
    }
    if (!achou) 
    {
        //eleger uma variavel qualquer caso nao haja variaveis em comum
        //neste caso, a primeira variavel de P1
        common_var = *(P1->codigos_numerador);
        //TODO: implementar algum tipo de contador para verificar estes casos e refinar mais tarde
    }
    
    // 2) encontra-se o polinomio que faça P1 ser zero ao substituir na variavel comum
    // -gero um polinomio igual a P1, mas sem a variavel e com sinal trocado. Divido o resultado 
    // pelo coeficiente do monomio que contem a variavel.
    
    lista_ptr1 = copia_lista_expr(P1);
    //salvo o inicio da lista
    lista_ptr2 = lista_ptr1;
    
    //novamente, P1 e P2 são apenas combinações lineares das variaveis, então este passo está optimizado, comparando
    // apenas o caractere do numerador, e procurando apenas uma unica vez.
    while (lista_ptr2 != NULL) 
    {
        if (*(lista_ptr2->codigos_numerador) == common_var) 
        {
            //guardo o coeficiente da variável comumn no monomio a ser removido
            coeficiente = lista_ptr2->parametro;
            //removo o ponteiro
            lista_ptr1 = remove_lista_expr(lista_ptr1, lista_ptr2);
            //interrompo o laço
            break;
        }
        else //do contrario avanço ao proximo.
            lista_ptr2 = lista_ptr2->proximo;
    }
    
    //TODO: caso lista_ptr2 seja NULL, procurar a variavel em P2, condição em que não há variável comum
    //entretatno, se eu sempre eleger uma variável de P1, este problema jamais ocorrerá.
    
    //se lista_ptr1 ficar null, criar um elemento zerado
    if (lista_ptr1 == NULL) 
        lista_ptr1 = constroi_elemento_zerado();
    else
    {
        lista_ptr2 = lista_ptr1;
        while (lista_ptr2 != NULL) 
        {
            //agora que removi common_var dos monomios, devo inverter o sinal dos que sobraram e dividir pelo coeficiente
            //salvo o inicio da lista
            lista_ptr2->parametro *= -1.0;
            lista_ptr2->parametro /= coeficiente;
            lista_ptr2 = lista_ptr2->proximo;
        }

    }
    
    //3) substitui o polinomio encontrado na variavel comum em lista_expr3
    //agora é a parte mais complicada, devo scanear cada variável de cada monomio de lista_expr3, 
    //e se encontrar a variavel comum, eu a excluo, e multiplico o monomio resultante por cada um dos monomios de 
    //lista_ptr1 resultante do passo anterior. continuo escaneando a partir do primeiro monomio modificado,
    //já que a variável pode estar em um grau maior.

    
    lista_ptr2 = substitui_var(lista_ptr3, lista_ptr1, common_var);
    
    //destroi-se opolinomio de substituicao
    destroi_lista_expr_expandida(lista_ptr3);
    
    //atualiza-se o ponteiro
    lista_ptr3 = lista_ptr2;
    
    //4) simplifica e reordena P3_resto
    //BUG: 0,9 - 0,9 = 1.11E-15, ou seja tenho que forçar um zero no braço.
    lista_ptr2 = simplifica_expr_expandida(lista_ptr3);
    
   
    //destroi-se P3_resto sem simplificação
    destroi_lista_expr_expandida(lista_ptr3);
    
    //ordena-se o polinomio
    lista_ptr2 = lexdegbubblesort(lista_ptr2);
    
    //substituo o polinomio de substituição em P2
    lista_ptr3 = substitui_var(P2, lista_ptr1, common_var);
    
    //destruo o polinomio de substituição na variavel comum
    destroi_lista_expr_expandida(lista_ptr1);
    
    //simplifico e reordeno
    lista_ptr1 = simplifica_expr_expandida(lista_ptr3);
    lista_ptr1 = lexdegbubblesort(lista_ptr1);
    destroi_lista_expr_expandida(lista_ptr3);
    lista_ptr3 = NULL;
    
   // 5) divide-se P3_resto por P2_subst e guarda o resultado em numerador_a
    resto = NULL;
    
    if(!polydiv(lista_ptr2, lista_ptr1, numerador_a, &resto))
    {
        //divisao por 0
        destroi_lista_expr_expandida(*numerador_a);
        *numerador_a = NULL;
        destroi_lista_expr_expandida(lista_ptr1);
        destroi_lista_expr_expandida(lista_ptr2);
        destroi_lista_expr_expandida(*quociente);
        destroi_lista_expr_expandida(P3_resto);
        destroi_lista_expr_expandida(resto);
        *quociente = NULL;
        return 0;
        
    
    }
    destroi_lista_expr_expandida(lista_ptr1);
    destroi_lista_expr_expandida(lista_ptr2);
    
    //HIPOTESE  - se o resto for diferente de 0, então obrigatoriamente teremos um numerador em "a" de grau igual a p_3, resto.
    //Eu poderia ter feito toda a conta e encontrado os numeradores para depois testar o grau, mas quero que esta etapa seja eficiente.

    if (resto->parametro != 0.0) 
    {
        //limpar tudo e retornar 0
        destroi_lista_expr_expandida(*numerador_a);
        destroi_lista_expr_expandida(*quociente);
        destroi_lista_expr_expandida(P3_resto);
        destroi_lista_expr_expandida(resto);
        *numerador_a = NULL;
        *quociente = NULL;
        return 0;
    }
    
    //liberar o resto para a proxima etapa
    destroi_lista_expr_expandida(resto);
    resto = NULL;
    
    //foi encontrado um numerador_a decente. Encontrar agora o numerador_b
    // Para achar o segundo numerador da fracao parcial:
    //1) multiplico "a" por P2
    //copio P3
    lista_ptr2 = copia_lista_expr(P3_resto); 
    lista_ptr1 = multiplica_expr(*numerador_a, P2);
    
    //2) subtraio lista_ptr1 de P3_resto, e guardo o resultado em lista_ptr2
    subtrai_expr(&lista_ptr2, lista_ptr1);
    lista_ptr3 = simplifica_expr_expandida(lista_ptr2);
    destroi_lista_expr_expandida(lista_ptr2);
    lista_ptr2 = lexdegbubblesort(lista_ptr3);
    
    //3) divido o resultado por P1
    polydiv(lista_ptr2, P1, numerador_b, &resto);
    
    //HIPOTESE - igual acima. 
    if (resto->parametro != 0.0) 
    {
        //limpar tudo e retornar 0
        destroi_lista_expr_expandida(lista_ptr2);
        destroi_lista_expr_expandida(*numerador_b);
        destroi_lista_expr_expandida(*numerador_a);
        destroi_lista_expr_expandida(*quociente);
        *numerador_a = NULL;
        *numerador_b = NULL;
        *quociente = NULL;
        destroi_lista_expr_expandida(resto);
        destroi_lista_expr_expandida(P3_resto);
        
        return 0;
    }
    
    //se chegou até aqui, a operação foi um sucesso
    destroi_lista_expr_expandida(lista_ptr2);
    destroi_lista_expr_expandida(resto);
    destroi_lista_expr_expandida(P3_resto);
    return 1;
    
}


lista_expr *substitui_var(lista_expr *p_destino, lista_expr *p_fonte, char common_var)
{
    char *str_ptr;
    lista_expr *lista_ptr1, *lista_ptr2, *lista_ptr3, *p_dest;
    p_dest = copia_lista_expr(p_destino);
    while (TRUE) 
    {
        //procuro pela variavel comum no monomio
        if(p_dest->codigos_numerador != NULL)
        {
            if ((str_ptr = strchr(p_dest->codigos_numerador, common_var)) != NULL) 
            {
                //removo o caracter da variavel da string, deslocando o restante dela para a esquerda
                while (*str_ptr != '\0') 
                {
                    *str_ptr = *(str_ptr + 1); 
                    str_ptr++;
                }
                
                
                //inicialização de ponteiros
                lista_ptr1 = NULL;
                lista_ptr2 = NULL;
                
                //agora isolo o monomio, caso P3 tenha mais de 1
                if (p_dest->proximo != NULL) 
                {
                    //salvo o proximo monomio
                    lista_ptr1 = p_dest->proximo;
                    
                    //isolo o monomio
                    p_dest->proximo = NULL;
                }
                
                if (p_dest->anterior != NULL) 
                {
                    //salvo o monomio anterior
                    lista_ptr2 = p_dest->anterior;
                    
                    //pode ser que isso seja desnecessario se eu for optimizar
                    p_dest->anterior = NULL;
                }
                
                //a principio, a string conter apenas um \0 nao impede a multiplicação de funcionar corretamente
                // mas se der bug, devo começar por aqui.
                //realizo a multiplicação
                lista_ptr3 = multiplica_expr(p_fonte, p_dest);
                
                //desaloco o monomio original
                destroi_lista_expr_expandida(p_dest);
                
                //reaproveito o ponteiro para apontar para o proximo monomio a ser avaliado, que é o primeiro 
                //do resultado da multiplicacao, para que esgote-se todos os graus da variável comum
                p_dest = lista_ptr3;
                
                //conecto o polinomio resultante da multiplicacao 
                lista_ptr3->anterior = lista_ptr2;
                if (lista_ptr2 != NULL) 
                    lista_ptr2->proximo = lista_ptr3;
                
                //encontro o ultimo monomio reultante da multiplicacao
                while (lista_ptr3->proximo != NULL)
                    lista_ptr3 = lista_ptr3->proximo;
                
                //e o conecto com o resto de P3
                lista_ptr3->proximo = lista_ptr1;
                if (lista_ptr1!= NULL)
                    lista_ptr1->anterior = lista_ptr3;
                
                //como já determinei qual monomio será avaliado no proximo passo, posso pular para o inicio do laço
                continue;
            }
        }// else: TODO eu simplesmente pulo o monomio se for só um parametro.
        
        //continuo a busca ou termino o laço, "rebobinando" o ponteiro para o inicio da lista
        if (p_dest->proximo == NULL) 
        {
            while (p_dest->anterior != NULL) 
                p_dest = p_dest->anterior;
            
            break;
        }
        else
            p_dest = p_dest->proximo;
        
    }
    return p_dest;
}




//função que gera um conjunto de vetores de entrada
vetor_polinomios *gera_vetor(vetor_polinomios *ultimo_gerado, lista_expr *polinomio_base, lista_expr *variavel_atual, int lim_inferior, int lim_superior)
{
    //variavel que vai dar um identificador aos polinomios
    static int poly_id = 0;
    vetor_polinomios *elemento_atual = NULL;
    int i;
    for (i=lim_superior; i>=lim_inferior; i--) 
    {
        //atribuo a variavel atual o valor de i;
        variavel_atual->parametro = i;
        //percorre o polinomio base até chegar na ultima variável
        if (variavel_atual->proximo != NULL)
            ultimo_gerado = gera_vetor(ultimo_gerado, polinomio_base, variavel_atual->proximo, lim_inferior, lim_superior);
        else //chegou na ultima variavel
        {
            //ao chegar na ultima variavel, gerar vetores variando o parametro da ultima variavel de lim_inferior a lim_superior
            for (i=lim_superior; i>=lim_inferior; i--) 
            {
                variavel_atual->parametro = i;
                elemento_atual = (vetor_polinomios *)malloc(sizeof(vetor_polinomios)); //TODO: controle de erro
                elemento_atual->polinomio = (polinomio *)malloc(sizeof(polinomio)); //TODO: controle de erro
                //copio o polinomio base com os parametros do jeito que estão
                //elimino os elementos zero
                elemento_atual->polinomio->P = simplifica_expr_expandida(polinomio_base);
                
                //atribuo um identificador e incremento
                elemento_atual->polinomio->id = poly_id;
                poly_id++;
                
                //concateno na lista de polinomios
                if (ultimo_gerado != NULL) 
                {
                    ultimo_gerado->proximo_polinomio = elemento_atual;
                    elemento_atual->polinomio_anterior = ultimo_gerado;
                    ultimo_gerado = elemento_atual;
                }
                else 
                {
                    ultimo_gerado = elemento_atual;
                    elemento_atual->proximo_polinomio = NULL;
                    elemento_atual->polinomio_anterior = NULL;
                }
                
                
            }
            
            return ultimo_gerado;
        }
        
        
    }
    
    return ultimo_gerado;
        
}

vetor_polinomios *elimina_zero(vetor_polinomios *lista)
{
    vetor_polinomios *p_lista = lista;
    while (p_lista != NULL) 
    {
        //apenas no polinomio nulo, o primeiro monomio é zero
        if (p_lista->polinomio->P->parametro == 0.0) 
        {
            p_lista = remove_polinomio(p_lista);
            
            //rebobinar a lista
            while (p_lista->polinomio_anterior != NULL) 
                p_lista = p_lista->polinomio_anterior;
            
            return p_lista;
        }
        
        //incrementa a busca
        p_lista = p_lista->proximo_polinomio;
            
    }
    
    //se nada encontrou, retorna o inicio da lista
    return lista;
}
//retorna o elemento anterior ao removido, ou o inico da lista
vetor_polinomios *remove_polinomio(vetor_polinomios *elemento)
{
    vetor_polinomios *esq, *dir;
    
    //gravo as bordas
    esq = elemento->polinomio_anterior;
    dir = elemento->proximo_polinomio;
    
    //ajusto os ponteiros à esquerda
    if (esq != NULL) 
        esq->proximo_polinomio = dir;
    
    //ajusto os ponteiros à direita
    if (dir != NULL) 
        dir->polinomio_anterior = esq;
    
    //eliminar o elemento
    destroi_lista_expr_expandida(elemento->polinomio->P);
    free(elemento->polinomio);
    free(elemento);
    
    //caso tenha esvaziado a lista
    if (dir == NULL && esq == NULL) 
        return NULL;
    
    //se apeans o esquerdo for nulo, ja sabemos que o inicio da lista é o da direita
    if (esq == NULL) 
    {
        return dir;
    }
    else //retornar o esquerdo
    {
        return esq;
    }
        
}

//retorna o anterior, pois pode ser que remova o primeiro elemento
vetor_polinomios *remove_polinomio_retorna_anterior(vetor_polinomios *elemento)
{
    vetor_polinomios *esq, *dir;
    
    //gravo as bordas
    esq = elemento->polinomio_anterior;
    dir = elemento->proximo_polinomio;
    
    //ajusto os ponteiros à esquerda
    if (esq != NULL) 
        esq->proximo_polinomio = dir;
    
    //ajusto os ponteiros à direita
    if (dir != NULL) 
        dir->polinomio_anterior = esq;
    
    //eliminar o elemento
    destroi_lista_expr_expandida(elemento->polinomio->P);
    free(elemento->polinomio);
    free(elemento);
    
    //caso tenha esvaziado a lista
    if (dir == NULL && esq == NULL) 
        return NULL;
    
    //retornar o proximo polinomio
    return esq;
}


//gera o polinomio base
lista_expr *gera_polinomio_base(tabela_literais *lista_literais)
{
    tabela_literais *percorre_tabela_literais = lista_literais;
    lista_expr *polinomio_base, *monomio_atual;
    
    polinomio_base = NULL;
    
    // percorre as variaveis gravadas na tabela de literais, gerando um polinomio que é a combinação linear delas.
    while (percorre_tabela_literais != NULL) 
    {
        monomio_atual = constroi_elemento_zerado();
        monomio_atual->codigos_numerador = (char *)malloc(2*sizeof(char));
        monomio_atual->codigos_numerador[0] = percorre_tabela_literais->codigo;
        monomio_atual->codigos_numerador[1] = '\0';
        monomio_atual->parametro = 1.0;
        monomio_atual->codigos_denominador = NULL;
        if (polinomio_base == NULL) 
        {
            polinomio_base = monomio_atual;
            polinomio_base->anterior = NULL;
            polinomio_base->proximo = NULL;
        }
        else
        {
            polinomio_base->proximo = monomio_atual;
            monomio_atual->anterior = polinomio_base;
            polinomio_base = polinomio_base->proximo;
        }
        percorre_tabela_literais = percorre_tabela_literais->proximo_codigo;
        
    }
    //rebobina o polinomio base
    while (polinomio_base->anterior != NULL) 
    {
        polinomio_base = polinomio_base->anterior;
    }
    
    return polinomio_base;

}

vetor_polinomios *remove_polinomios_negativos (vetor_polinomios *vetor_entrada)
{
    vetor_polinomios *percorre_vetor;
    lista_expr *percorre_expr;
    int teste;
    
    //inicializo o ponteiro de varredura
    percorre_vetor = vetor_entrada;
    
    //condição de erro
    if (percorre_vetor == NULL)
        return NULL;
    //varredura de todas as entradas geradas removendo tododos os polinomios onde todos os coeficientes são negativos
    while (percorre_vetor != NULL) 
    {
        //aponta para o polinomio
        percorre_expr = percorre_vetor->polinomio->P;
        
        //inicializo a variável de teste
        teste = 0;
        
        //percorro o polinomio procurando coeficientes positivos. Se achar algum mudo o estado da variável de teste
        while (percorre_expr != NULL) 
        {
            if (percorre_expr->parametro > 0.0) 
            {
                teste = 1;
                break;
            }
            
            //incremento o ponteiro
            percorre_expr = percorre_expr->proximo;
        }
        
        //se teste for 0, o polinomio é todo negativo, e deve ser removido do vetor
        if (!teste) 
        {
            //a variável de retorno guarda o polinomio anterior, para que o final da lista seja mantido
            percorre_vetor = remove_polinomio_retorna_anterior(percorre_vetor);
            
        }
        //se nao for o ultimo elemento incrementa, caso contrario sai do loop
        else if(percorre_vetor->proximo_polinomio != NULL)
            percorre_vetor = percorre_vetor->proximo_polinomio;
        else
            break;
        
    } 
    
    //rebobina o vetor reduizdo
    while (percorre_vetor->polinomio_anterior != NULL) 
    {
        percorre_vetor = percorre_vetor->polinomio_anterior;
    }
    
    
    return percorre_vetor;
}

vetor_polinomios *remove_polinomios_redundantes (vetor_polinomios *vetor_entrada)
{
    vetor_polinomios *p_vetor1, *p_vetor2;
    lista_expr *expr1, *expr2, *Q, *R;
    
    //inicializo o ponteiro de varredura
    p_vetor1 = vetor_entrada;
    
    //condição de erro
    if (p_vetor1 == NULL)
        return NULL;
    //para cada vetor, eu divido ele por todos os outros polinomoios. se o resto da zero, é porque um é combinacao linear do outro
    while (p_vetor1->proximo_polinomio != NULL) 
    {
        p_vetor2 = p_vetor1->proximo_polinomio;
        while (p_vetor2 != NULL) 
        {
            //copio a expressao de vetor1
            expr1 = copia_lista_expr(p_vetor1->polinomio->P);
            //copio a expressao de vetor2
            expr2 = copia_lista_expr(p_vetor2->polinomio->P);
            
            //divido os dois
            polydiv(expr1, expr2, &Q, &R);
            
            //se o resto for 0, e o quociente for uma constante, devo remover p_vetor_2
            if (R->parametro == 0.0 && Q->codigos_numerador == NULL) 
            {
                //p_vetor2 apontará para o polinomio anterior, logo nao há problema em incrementarmos o ponteiro log a seguir
                p_vetor2 = remove_polinomio(p_vetor2);
            }
            
            //destruo expr
            destroi_lista_expr_expandida(expr1);
            destroi_lista_expr_expandida(expr2);
            destroi_lista_expr_expandida(Q);
            destroi_lista_expr_expandida(R);
            Q = NULL;
            R = NULL;
            
            //incremento o ponteiro
            if (p_vetor2 != NULL) 
                p_vetor2 = p_vetor2->proximo_polinomio;
        }
        
        //incremento a busca principal
        if (p_vetor1->proximo_polinomio!= NULL) 
            p_vetor1 = p_vetor1->proximo_polinomio;
        else
            break;
        
    } 
    
    //rebobina o vetor reduizdo
    while (p_vetor1->polinomio_anterior != NULL) 
    {
        p_vetor1 = p_vetor1->polinomio_anterior;
    }
    
    
    return p_vetor1;
}

//função que retorna o grau de um polinomio
int deg(lista_expr *poly_in)
{
    int grau = 0;
    int grau_max = 0;
    lista_expr *ptr_poly = poly_in;
    
    //percorre monomio a monomio verificando o grau, que nada mais é senão o tamanho da string de codigos do numerador
    while (ptr_poly != NULL) 
    {
        if (ptr_poly->codigos_numerador == NULL) 
            grau = 0;
        else 
            grau = (int)strlen(ptr_poly->codigos_numerador);
        
        if (grau > grau_max) 
            grau_max = grau;
        
        ptr_poly = ptr_poly->proximo;
        
    }
    
    return grau_max;
    
}

vetor_sementes *gera_vetor_semente(vetor_polinomios *vetor_in, lista_expr *eq_entrada)
{
    vetor_polinomios *ptr_vetor_in = vetor_in->proximo_polinomio;
    vetor_polinomios *ptr_vetor_base = vetor_in;
    vetor_sementes *vetor_gerado = NULL;
    
    //resultados
    lista_expr *R1;
    lista_expr *R2;
    lista_expr *Q;
    
    //se por um acaso passar apenas um elemento, retornar nulo
    if (ptr_vetor_in == NULL)
    {
        return NULL;
    }
    //percorre o vetor de entrada, tentando todas as combinações possíveis. Como a ordem não importa, vou sempre incrementando vetor_in também
    while (ptr_vetor_base->proximo_polinomio != NULL) 
    {
        //ptr_vetor_in = vetor_in->proximo_polinomio;
        //não restrição nenhuma em poder combinar consigo mesmo
        ptr_vetor_in = ptr_vetor_base;
        
        while (ptr_vetor_in != NULL) 
        {
            //tentar realizar a expansão em fracoes parciais
            if(partial_fraction_expansion(eq_entrada, ptr_vetor_base->polinomio->P, ptr_vetor_in->polinomio->P, &Q,&R1, &R2))
            {
                //consegui uma expansão em fracoes parciais adequada, preparar vetor de semente, ao final ele será rebobinado
                if (vetor_gerado == NULL) 
                    vetor_gerado = novo_vetor_semente();
                else
                {
                    //crio um novo elemento ja ligado ao anterior
                    vetor_gerado->conjunto_prox = novo_vetor_semente();
                    //realizo o duplo encadeamento
                    vetor_gerado->conjunto_prox->conjunto_ant = vetor_gerado;
                    //atualizo o ponteiro
                    vetor_gerado = vetor_gerado->conjunto_prox;
                }
                
                //preencho os campos da semente, P1 e P2 devem ser copiados, pois ainda serão utilizados
                vetor_gerado->P1.P = copia_lista_expr(ptr_vetor_base->polinomio->P);
                vetor_gerado->P1.id = ptr_vetor_base->polinomio->id;
                vetor_gerado->P2.P = copia_lista_expr(ptr_vetor_in->polinomio->P);
                vetor_gerado->P2.id = ptr_vetor_in->polinomio->id;
                //estes outros 3 não precisam de copia, pois serão utilizados apenas dentro do vetorde sementes
                vetor_gerado->quociente = Q;
                vetor_gerado->R1 = R1;
                vetor_gerado->R2 = R2;
                    
            }
            
            //reseta as variaveis
            R1 = NULL;
            R2 = NULL;
            Q = NULL;
            
            //incremento ponteiro para proximo teste
            ptr_vetor_in = ptr_vetor_in->proximo_polinomio;
        }
        
        //incremento o ponteiro para a proxima bateria de testes
        ptr_vetor_base = ptr_vetor_base->proximo_polinomio; 
    }
    //nao encontrou nenhuma semente
    if (vetor_gerado == NULL)
        return NULL;
    //rebobina o vetor de sementes
    while (vetor_gerado->conjunto_ant != NULL) 
        vetor_gerado = vetor_gerado->conjunto_ant;
    
    //retorna o vetor gerado
    return vetor_gerado;
    
}

//cria elemento de vetor de sementes
vetor_sementes *novo_vetor_semente(void)
{
    vetor_sementes *novo;
    
    novo = (vetor_sementes *)malloc(sizeof(vetor_sementes)); //TODO: conferir erro de alocacao de memoria
    
    novo->P1.P = NULL;
    novo->P1.id = 0;
    novo->P2.P = NULL;
    novo->P2.id = 0;
    novo->quociente = NULL;
    novo->R1 = NULL;
    novo->R2 = NULL;
    novo->conjunto_ant = NULL;
    novo->conjunto_prox = NULL;
    
    return novo;
}

//destroi toda a lista de sementes
void destroi_lista_sementes(vetor_sementes *entrada)
{
    if (entrada != NULL) 
        destroi_lista_sementes(entrada->conjunto_prox);
    else
        return;
     
    destroi_lista_expr_expandida(entrada->P1.P);
    destroi_lista_expr_expandida(entrada->P2.P);
    destroi_lista_expr_expandida(entrada->quociente);
    destroi_lista_expr_expandida(entrada->R1);
    destroi_lista_expr_expandida(entrada->R2);
    free(entrada);
    
}


//cria um novo elemento de vetor de polinomios 
vetor_polinomios *novo_vetor_polinomios(void)
{
    vetor_polinomios *novo;
    
    novo = (vetor_polinomios *)malloc(sizeof(vetor_polinomios));

    novo->polinomio = NULL;
    novo->polinomio_anterior = NULL;
    novo->proximo_polinomio = NULL;
    
    return (novo);
}

//cria elemento de vetor de decomposição
vetor_decomp *novo_vetor_decomp(void)
{
    vetor_decomp *novo;
    novo = (vetor_decomp *)malloc(sizeof(vetor_decomp));
    novo->poly_pares = NULL;
    novo->poly_impares = NULL;
    novo->resto_impar = NULL;
    novo->resto_par = NULL;
    
    novo->prox_decomp = NULL;
    novo->ant_decomp = NULL;
    
    return novo;
}

//insere um polinomio na lista de polinomios
void insere_polinomio(vetor_polinomios **vetor,polinomio *elemento)
{
    vetor_polinomios *p_vetor, *novo_poly;
    
    //cria-se o novo elemento
    novo_poly = novo_vetor_polinomios();
    novo_poly->polinomio = elemento;
    
    //inicializa a lista
    p_vetor = *vetor;
    
    //caso a lista esteja vazia, tornar o novo elemento o primeiro elemento da lista
    if (p_vetor == NULL) 
    {
        *vetor = novo_poly;
    }
    else 
    {
        //procura o fim da fila
        while (p_vetor->proximo_polinomio != NULL) 
            p_vetor = p_vetor->proximo_polinomio;
        
        //insere o elemento
        p_vetor->proximo_polinomio = novo_poly;
        novo_poly->polinomio_anterior = p_vetor;
    }
    
}


void destroi_decomp(vetor_decomp *entrada)
{
    if (entrada->resto_impar != NULL) 
        destroi_lista_expr_expandida(entrada->resto_impar);
    if (entrada->resto_par != NULL) 
        destroi_lista_expr_expandida(entrada->resto_par);
//destroi os vetores de polinomios
    if (entrada->poly_pares != NULL) 
        destroi_vetor_polinomios(entrada->poly_pares);
    if (entrada->poly_impares != NULL) 
        destroi_vetor_polinomios(entrada->poly_impares);
    
    free(entrada);
}

void destroi_vetor_polinomios(vetor_polinomios *entrada)
{
    if (entrada->proximo_polinomio != NULL) 
    {
        destroi_vetor_polinomios(entrada->proximo_polinomio);
    }
    
    //o vetor de polinomios apenas pode referenciar polinomios e nao copiá-los
    free(entrada);
}

void destroi_vetor_decomp(vetor_decomp *entrada)
{
    if (entrada == NULL) 
        return;
    if (entrada->prox_decomp != NULL) 
    {
        destroi_vetor_decomp(entrada->prox_decomp);
    }
    
    destroi_decomp(entrada);
}
//tirar a prova real entre a decomposição encontrada e a equação de entrada
int prova_real(vetor_decomp *decomp, lista_expr *eq_entrada)
{
    lista_expr *acum_par, *acum_impar, *acum, *eq_entrada_copia;
    vetor_polinomios *poly_par_ptr = decomp->poly_pares;
    vetor_polinomios *poly_impar_ptr = decomp->poly_impares;
    
    acum_par = copia_lista_expr(poly_par_ptr->polinomio->P);
    acum_impar = copia_lista_expr(poly_impar_ptr->polinomio->P);
    
    while (poly_par_ptr->proximo_polinomio!= NULL)
    {
        acum = multiplica_expr(acum_par, poly_par_ptr->proximo_polinomio->polinomio->P);
        destroi_lista_expr_expandida(acum_par);
        acum_par = acum;
        
        poly_par_ptr = poly_par_ptr->proximo_polinomio;
    }
    
    //multiplico o resultado pelo lambda, 
    acum = multiplica_expr(acum_par, decomp->resto_impar);
    destroi_lista_expr_expandida(acum_par);
    acum_par = acum;
    
    //repito o mesmo para a parte impar
    while (poly_impar_ptr->proximo_polinomio!= NULL)
    {
        acum = multiplica_expr(acum_impar, poly_impar_ptr->proximo_polinomio->polinomio->P);
        destroi_lista_expr_expandida(acum_impar);
        acum_impar = acum;
        
        poly_impar_ptr = poly_impar_ptr->proximo_polinomio;
    }
    
    //multiplico o resultado pelo lambda
    acum = multiplica_expr(acum_impar, decomp->resto_par);
    destroi_lista_expr_expandida(acum_impar);
    acum_impar = acum;
    
    //somo a parte impar com a parte par
    soma_expr(acum_impar, acum_par);
    
    //subtraio o resultado da equação de entrada
    eq_entrada_copia = copia_lista_expr(eq_entrada);
    
    subtrai_expr(&acum_impar, eq_entrada_copia);
    
    //simplifica 
    acum = simplifica_expr_expandida(acum_impar);
    destroi_lista_expr_expandida(acum_impar);
    acum_impar = acum;
    
    //se o resultado for 0 deu certo
    if(acum_impar->codigos_numerador == NULL && acum_impar->parametro == 0.0 && acum_impar->proximo == NULL)
        return TRUE;
    else 
        return FALSE;
    
}

//insere uma lista de decomposicoes dentro de outra, retornando um ponteiro para o final dela
vetor_decomp *insere_lista_decomp(vetor_decomp *lista_destino, vetor_decomp *lista_origem)
{
    vetor_decomp *ptr_lista_origem, *ptr_lista_destino;
    
    ptr_lista_origem = lista_origem;
    if (lista_destino == NULL) 
    {
        //aponto para o final da lista a ser inserida
        while (ptr_lista_origem->prox_decomp != NULL)     
            ptr_lista_origem = ptr_lista_origem->prox_decomp;
        
        //retorno o fim da lista destino como o fim da lista de origem
        return ptr_lista_origem;
        
    }
    else
    {
        ptr_lista_destino = lista_destino;
        //avançar o ponteiro para fim da lista destino
        while (ptr_lista_destino->prox_decomp != NULL) 
            ptr_lista_destino = ptr_lista_destino->prox_decomp;
        
        //rebobinar o ponteiro da lista de origem
        while (ptr_lista_origem->ant_decomp != NULL) 
            ptr_lista_origem = ptr_lista_origem->ant_decomp;
        
        //ligar as 2
        ptr_lista_destino->prox_decomp = ptr_lista_origem;
        ptr_lista_origem->ant_decomp = ptr_lista_destino;
        
        //avançar o ponteiro até o final
        while (ptr_lista_destino->prox_decomp != NULL) 
            ptr_lista_destino = ptr_lista_destino->prox_decomp;
        
        return ptr_lista_destino;
    }
}

//funcao que pega os dados de um vetor semente e os transfere para um vetor decomp
vetor_decomp *copia_vetor_semente(vetor_sementes *entrada)
{
    vetor_sementes *ptr_sementes;
    vetor_decomp   *ptr_decomp = NULL;
    vetor_decomp   *lista_decomp = NULL;
    
    ptr_sementes = entrada;
    while (ptr_sementes != NULL) 
    {
        ptr_decomp = novo_vetor_decomp();
        ptr_decomp->poly_pares = novo_vetor_polinomios();
        ptr_decomp->poly_pares->polinomio = &(ptr_sementes->P2);
        ptr_decomp->poly_impares = novo_vetor_polinomios();
        ptr_decomp->poly_impares->polinomio = &(ptr_sementes->P1);
        ptr_decomp->resto_par = copia_lista_expr(ptr_sementes->R2) ;
        ptr_decomp->resto_impar = copia_lista_expr(ptr_sementes->R1) ;
        
        //insere o elemento criado na lista
        lista_decomp = insere_lista_decomp(lista_decomp, ptr_decomp);
        
        //incrementa o ponteiro
        ptr_sementes = ptr_sementes->conjunto_prox;
        
    }
    
    //rebobinar a lista_decomp
    while (lista_decomp->ant_decomp != NULL) 
        lista_decomp = lista_decomp->ant_decomp;
    
    return lista_decomp;
}


//função que implementa a decomposição em si - versao recursiva
vetor_decomp *encontra_decomp(vetor_sementes *entrada, int grau, lista_expr *expr_simplificada, tabela_literais *lista_literais)
{
    vetor_sementes *secundario_ptr; //ponteiro para os loops internos
    
    vetor_decomp *decomp_atual; //lista que será gerada dentro do loop atual
    vetor_decomp *decomp_referencia; //lista de refrencia para o loop atual
    vetor_decomp *lista_decomp = NULL; //elemento manipulado dentro do loop
    
    int total_decomp = 0;

    
    //inicialização da lista de referência
    decomp_referencia = copia_vetor_semente(entrada);
    decomp_atual = decomp_referencia;
    decomp_atual = NULL;
    
    //inicialização do ponteiro de busca secundário, pode ser que o par de polinomios inicial possa combinar com si mesmo
    secundario_ptr = entrada;
    
    //loop de construcao das decomposicoes
    while (decomp_referencia != NULL) 
    {
        //para cada vetor incial de decomp referencia, encontro todas as decomposicoes que possam ser geradas a partir dele
        encontra_decomp_recursiva(decomp_referencia, decomp_referencia, &lista_decomp, grau, expr_simplificada, lista_literais, &total_decomp ); 
        
        //apagar o elemento testado e atualizar o ponteiro de referencia
        decomp_atual = decomp_referencia->prox_decomp;
        destroi_decomp(decomp_referencia);
        decomp_referencia = decomp_atual;
        if (decomp_referencia!= NULL)
            decomp_referencia->ant_decomp = NULL;
    }
    
    
    //rebobinar a lista de decomposições
    if (lista_decomp != NULL) 
        while (lista_decomp->ant_decomp != NULL) 
            lista_decomp = lista_decomp->ant_decomp;
    
    //imprimir numero de decomposicoes
    printf("\n\n o numero de decoposicoes validas: %d", total_decomp);
    
    //retornar
    return lista_decomp;
}


//funcao que pega os dados de um vetor semente e os transfere para um vetor decomp
vetor_decomp *copia_semente(vetor_sementes *entrada)
{
    vetor_decomp   *ptr_decomp = NULL;
    
    ptr_decomp = novo_vetor_decomp();
    ptr_decomp->poly_pares = novo_vetor_polinomios();
    ptr_decomp->poly_pares->polinomio = &entrada->P2;
    ptr_decomp->poly_impares = novo_vetor_polinomios();
    ptr_decomp->poly_impares->polinomio = &entrada->P1;
    ptr_decomp->resto_par = copia_lista_expr(entrada->R2);
    ptr_decomp->resto_impar = copia_lista_expr(entrada->R1);
    ptr_decomp->ant_decomp = NULL;
    ptr_decomp->prox_decomp = NULL;
    
    
    
    return ptr_decomp;
}

vetor_decomp *copia_decomp(vetor_decomp *entrada)
{
    vetor_decomp *nova_decomp;
    vetor_polinomios *poly_ptr;
    
    nova_decomp = novo_vetor_decomp();
    //copio os polinomios do primario em nova_decomp, primeiro a parte par
    poly_ptr = entrada->poly_pares;
    while (poly_ptr != NULL) 
    {
        insere_polinomio(&(nova_decomp->poly_pares), poly_ptr->polinomio);
        poly_ptr = poly_ptr->proximo_polinomio;
    }
    
    //copio os polinomios do primario em nova_decomp, depois a parte impar
    poly_ptr =  entrada->poly_impares;
    while (poly_ptr != NULL) 
    {
        insere_polinomio(&nova_decomp->poly_impares, poly_ptr->polinomio);
        poly_ptr = poly_ptr->proximo_polinomio;
    }
    
    //copia os outros parametros
    nova_decomp->resto_impar = NULL; //copia_lista_expr(entrada->resto_impar); estes parametros serao inseridos dentro da funcao recursiva
    nova_decomp->resto_par = NULL; //copia_lista_expr(entrada->resto_par);
    nova_decomp->ant_decomp = NULL;
    nova_decomp->prox_decomp = NULL;
    
    return nova_decomp;

}

int decomp_size(vetor_decomp *entrada)
{
    int count = 0;
    
    //inicia com o esqueleto da estrutura
    count = sizeof(vetor_decomp);
    
    //soma os polinomios pares
    count+= poly_size(entrada->poly_pares);
    
    //soma os polinomios impares
    count+= poly_size(entrada->poly_impares);
    
    //soma os restos
    count+= expr_size(entrada->resto_par);
    count+= expr_size(entrada->resto_impar);
    
    return count;
}

int poly_size(vetor_polinomios *entrada)
{
    int count = 0;
    
    //soma o tamanho em memoria de todos os polinomios afrente
    if (entrada->proximo_polinomio != NULL) 
    {
        count+= poly_size(entrada->proximo_polinomio);
    }
     //soma o proprio tamanho
    count += sizeof(vetor_polinomios);
    
    return count;
    
}

int expr_size(lista_expr *entrada)
{
    int count = 0;
    
    //soma a memoria dos proximos monomios
    if (entrada->proximo != NULL) 
    {
        count+= expr_size(entrada->proximo);
    }
    
    //soma os proprios campos
    count+= sizeof(lista_expr);
    
    if (entrada->codigos_numerador != NULL) 
    {
        count += (strlen(entrada->codigos_numerador) + 1)*sizeof(char);
    }
    
    if (entrada->codigos_denominador != NULL) 
    {
        count+= expr_size(entrada->codigos_denominador);
    }
    
    return count;
    
}







//Funcao que testa se um par de polinomios pode ser combinado com outro para formar uma decomposição, retornando a decomposicao parcial
void encontra_decomp_recursiva(vetor_decomp *primario, vetor_decomp *secundario, vetor_decomp **retorno, int grau, lista_expr *expr_simplificada, tabela_literais *lista_literais, int *total_decomp)
{
    lista_expr *Q = NULL;
    lista_expr *R1= NULL;
    lista_expr *R2= NULL;
    lista_expr *R_dummy = NULL;
    vetor_decomp *decomp_atual = NULL;
    vetor_polinomios *poly_ptr = NULL;
    vetor_decomp *ptr_sementes = NULL;
    vetor_decomp *ptr_decomp = NULL;
    int contador = 0;
    int flag = 0;
    
    //aponto para o inicio do vetor sementes
    ptr_sementes = secundario;
    
    //vou tentar combinar o primario com TODOS os elementos do secundário
    while (ptr_sementes != NULL)
    {
        
        //testar P1 principal com P1 atual
        if(partial_fraction_expansion(primario->resto_impar, primario->poly_impares->polinomio->P, ptr_sementes->poly_impares->polinomio->P, &Q,&R1, &R_dummy))
        {
            //limpar o R_dummy e o &Q
            destroi_lista_expr_expandida(R_dummy);
            destroi_lista_expr_expandida(Q);
            R_dummy = NULL;
            Q = NULL;
            
            //testar P2 principal com P2 atual
            if(partial_fraction_expansion(primario->resto_par, primario->poly_pares->polinomio->P, ptr_sementes->poly_pares->polinomio->P, &Q,&R2, &R_dummy))
            {
                //limpar o R_dummy e o &Q
                destroi_lista_expr_expandida(R_dummy);
                destroi_lista_expr_expandida(Q);
                
                //copio o primario, pois pode ser semente para outras decomposicoes
                decomp_atual = copia_decomp(primario);
                
                //insiro os polinomios que fazem parte da decomposição
                //significa que secundario->P1 é o proximo polinomio Par, e seconudario->P2 o proximo polinomio impar, e é invertido mesmo
                insere_polinomio(&decomp_atual->poly_pares,ptr_sementes->poly_impares->polinomio);
                insere_polinomio(&decomp_atual->poly_impares,ptr_sementes->poly_pares->polinomio);
                
                //atualiza Resto par e Resto impar
                decomp_atual->resto_impar = R1;
                decomp_atual->resto_par = R2;
                
                //quando o numero de polinomios pares e impares for igual ao grau da equacao de entrada pode ser que já tenha terminado
                contador = 0;
                poly_ptr = decomp_atual->poly_pares;
                while (poly_ptr != NULL)
                {
                    ++contador;
                    poly_ptr = poly_ptr->proximo_polinomio;
                }
                if (contador == grau)
                {
                    //testo se a decomposicaO encontrada é valida
                    //se um dos restos for 0, significa que simplesmente fatoramos o polinomio
                    if (decomp_atual->resto_impar->parametro == 0.0 || decomp_atual->resto_par->parametro == 0.0)
                    {
                        destroi_decomp(decomp_atual);
                        decomp_atual = NULL;
                    }
                    else
                    {
                        //tirar a prova real
                        if(prova_real(decomp_atual, expr_simplificada))
                        {
                            //ordenar os polinomios antes de inserir a decomposicão
                            decomp_atual->poly_pares = ordena_polinomios(decomp_atual->poly_pares);
                            decomp_atual->poly_impares = ordena_polinomios(decomp_atual->poly_impares);
                            
                            //procuro por uma decomposicao equivalente no vetor de retorno
                            ptr_decomp = *retorno;
                            
                            //inicializo uma flag de busca
                            flag = 0;
                            while (ptr_decomp != NULL)
                            {
                                //se eu encontro uma decomposicao equivalente, interrompo a busca e nao insiro decomp atual no vetro de retorno
                                if (compara_decomp(decomp_atual, ptr_decomp))
                                {
                                    flag = 1;
                                    break;
                                }
                                ptr_decomp = ptr_decomp->ant_decomp;
                            }
                            //se nao encontrou nenhuma decomposicao equivalente, inserir decomp autal no vetor
                            if (!flag)
                            {
                                *retorno = insere_lista_decomp(*retorno, decomp_atual);
                                //ja imprimo a decomposicao encontrada
                                imprime_decomposicao(decomp_atual, lista_literais);
                                //imcremento as decomp validas
                                (*total_decomp)++;
                                
                            }
                            else
                                destroi_decomp(decomp_atual);
                            
                            //limpa o ponteiro para o proximo teste
                            //TESE eu poderia terminar aqui, pois seria muito dificil encontrar dois ultimos pares que façam parte da mesma decomposicao
                            decomp_atual = NULL;
                            
                        }
                        else
                        {
                            destroi_decomp(decomp_atual);
                            decomp_atual = NULL;
                        }
                    }
                }
                else
                {
                    //caso contrario, deve-se proceder com a decomposicao recursivamente
                    encontra_decomp_recursiva(decomp_atual, ptr_sementes, retorno, grau, expr_simplificada, lista_literais, total_decomp);
                    //como as decomposições validas serao adicionadas no if acima, quando o programa chegar aqui, significa
                    //que nao preciso mais de decomp_atual;
                    destroi_decomp(decomp_atual);
                    decomp_atual = NULL;
                }
                
                //se chegou aqui, signigica que o segundo teste, com os polinomios invertidos, não é necessario
                Q = NULL;
                R1= NULL;
                R2= NULL;
                R_dummy = NULL;
                
                //atualiza o ponteiro
                ptr_sementes = ptr_sementes->prox_decomp;
                
                continue;
            }
            else
            {
                //limpar tudo
                destroi_lista_expr_expandida(R_dummy);
                destroi_lista_expr_expandida(Q);
                destroi_lista_expr_expandida(R1);
                destroi_lista_expr_expandida(R2);
            }
        }
        else
        {
            //limpar tudo
            destroi_lista_expr_expandida(R_dummy);
            destroi_lista_expr_expandida(Q);
            destroi_lista_expr_expandida(R1);
            
        }
        
        //limpar as variaveis de retorno
        Q = NULL;
        R1= NULL;
        R2= NULL;
        R_dummy = NULL;
        //testar também P1 principal com P2 atual
        if(partial_fraction_expansion(primario->resto_impar, primario->poly_impares->polinomio->P, ptr_sementes->poly_pares->polinomio->P, &Q,&R1, &R_dummy))
        {
            //if dummie para por um breakpoint exatamente onde esta dando pau no windows
            //limpar o R_dummy e o &Q
            destroi_lista_expr_expandida(R_dummy);
            destroi_lista_expr_expandida(Q);
            R_dummy = NULL;
            Q = NULL;
            
            //testar P2 principal com P1 atual
            if(partial_fraction_expansion(primario->resto_par, primario->poly_pares->polinomio->P, ptr_sementes->poly_impares->polinomio->P, &Q,&R2, &R_dummy))
            {
                //limpar o R_dummy e o &Q
                destroi_lista_expr_expandida(R_dummy);
                destroi_lista_expr_expandida(Q);
                R_dummy = NULL;
                Q = NULL;
                
                //copio o primario, pois pode ser semente para outras decomposicoes
                decomp_atual = copia_decomp(primario);
                
                //insiro os polinomios que fazem parte da decomposição
                //significa que secundario->P1 é o proximo polinomio Par, e seconudario->P2 o proximo polinomio impar, e é invertido mesmo
                insere_polinomio(&decomp_atual->poly_pares,ptr_sementes->poly_pares->polinomio);
                insere_polinomio(&decomp_atual->poly_impares,ptr_sementes->poly_impares->polinomio);
                
                //atualiza Resto par e Resto impar
                decomp_atual->resto_impar = R1;
                decomp_atual->resto_par = R2;
                
                //quando o numero de polinomios pares e impares for igual ao grau da equacao de entrada pode ser que já tenha terminado
                contador = 0;
                poly_ptr = decomp_atual->poly_pares;
                while (poly_ptr != NULL)
                {
                    ++contador;
                    poly_ptr = poly_ptr->proximo_polinomio;
                }
                if (contador == grau)
                {
                    //testo se a decomposicap encontrada é valida
                    //se um dos restos for 0, significa que simplesmente fatoramos o polinomio
                    if (decomp_atual->resto_impar->parametro == 0.0 || decomp_atual->resto_par->parametro == 0.0)
                    {
                        destroi_decomp(decomp_atual);
                        decomp_atual = NULL;
                    }
                    else
                    {
                        //tirar a prova real
                        if(prova_real(decomp_atual, expr_simplificada))
                        {
                            //ordeno o polinomio antes de inserir
                            decomp_atual->poly_pares = ordena_polinomios(decomp_atual->poly_pares);
                            decomp_atual->poly_impares = ordena_polinomios(decomp_atual->poly_impares);
                            
                            //procuro por uma decomposicao equivalente no vetor de retorno
                            ptr_decomp = *retorno;
                            
                            //inicializo uma flag de busca
                            flag = 0;
                            while (ptr_decomp != NULL)
                            {
                                //se eu encontro uma decomposicao equivalente, interrompo a busca e nao insiro decomp atual no vetro de retorno
                                if (compara_decomp(decomp_atual, ptr_decomp))
                                {
                                    flag = 1;
                                    break;
                                }
                                ptr_decomp = ptr_decomp->ant_decomp;
                            }
                            
                            //se nao encontrou nenhuma decomposicao equivalente, inserir decomp autal no vetor
                            if (!flag)
                            {
                                *retorno = insere_lista_decomp(*retorno, decomp_atual);
                                //ja imprimo a decomposicao encontrada
                                imprime_decomposicao(decomp_atual, lista_literais);
                                //imcremento as decomp validas
                                (*total_decomp)++;
                                
                            }
                            else
                                destroi_decomp(decomp_atual);
                            //limpa o ponteiro para o proximo teste
                            //TESE eu poderia terminar aqui, pois seria muito dificil encontrar dois ultimos pares que façam parte da mesma decomposicao
                            decomp_atual = NULL;
                        }
                        else
                        {
                            destroi_decomp(decomp_atual);
                            decomp_atual = NULL;
                        }
                    }
                }
                else
                {
                    
                    //caso contrario, deve-se proceder com a decomposicao recursivamente
                    encontra_decomp_recursiva(decomp_atual, ptr_sementes, retorno, grau, expr_simplificada,lista_literais, total_decomp);
                    
                    //como as decomposições validas serao adicionadas no if acima, quando o programa chegar aqui, significa
                    //que nao preciso mais de decomp_atual;
                    destroi_decomp(decomp_atual);
                    decomp_atual = NULL;
                }
                
            }
            else
            {
                //limpar tudo
                destroi_lista_expr_expandida(R_dummy);
                destroi_lista_expr_expandida(R2);
                destroi_lista_expr_expandida(R1);
                destroi_lista_expr_expandida(Q);
                
            }
            
        }
        else
        {
            //limpar tudo
            destroi_lista_expr_expandida(R_dummy);
            destroi_lista_expr_expandida(Q);
            destroi_lista_expr_expandida(R1);
        }
        //return decomp_atual;
        Q = NULL;
        R1= NULL;
        R2= NULL;
        R_dummy = NULL;
        
        //atualiza o ponteiro
        ptr_sementes = ptr_sementes->prox_decomp;
    }
    
}

//funcao que elimina as decomp redundantes nao previsiveis
void elimina_decomp_redundantes(vetor_decomp *entrada)
{
    vetor_decomp *ptr_entrada = entrada;
    vetor_decomp *ptr_aux, *ptr_remove;
    int flag;
    
    //primeiro reordeno os polinomios pares e impares de cada decomp
    while (ptr_entrada!= NULL) 
    {
        //reordeno os polinomios pares e impares
        ptr_entrada->poly_pares = ordena_polinomios(ptr_entrada->poly_pares);
        ptr_entrada->poly_impares = ordena_polinomios(ptr_entrada->poly_impares);
        
        ptr_entrada = ptr_entrada->prox_decomp;
    }
    
    ptr_entrada = entrada;
    //percorro toda a lista de decomposicoes, comparando com todas abaixo
    while (ptr_entrada!= NULL) 
    {
        ptr_aux = ptr_entrada->prox_decomp;
        while(ptr_aux != NULL)
        {
            flag = compara_decomp(ptr_entrada, ptr_aux);
            
            //finalmente, se a falg for 1, é porque os polinomios sao redundantes, entao posso excluir o polinomio auxiliar
            if (flag) 
            {
                //marco o ponteiro a ser removido
                ptr_remove = ptr_aux;
                
                //ja incremento o ponteiro auxiliar para a proxima iteração
                ptr_aux = ptr_aux->prox_decomp;
                
                //removo p_remove
                ptr_remove->ant_decomp->prox_decomp = ptr_remove->prox_decomp;
                if (ptr_remove->prox_decomp != NULL) 
                    ptr_remove->prox_decomp->ant_decomp = ptr_remove->ant_decomp;
                destroi_decomp(ptr_remove);
            }
            else
                ptr_aux = ptr_aux->prox_decomp;
        }
        
        ptr_entrada = ptr_entrada->prox_decomp;
        
    }
}
//ordena uma lista de polinomiosem ordem crescente
vetor_polinomios *ordena_polinomios(vetor_polinomios *entrada)
{
    vetor_polinomios *ptr, *aux;
    
    ptr = entrada;
    while (ptr->proximo_polinomio!= NULL) 
    {
        aux = ptr->proximo_polinomio;
        if (aux->polinomio->id < ptr->polinomio->id) 
        {
            //realizo a troca entre ptr e aux
            ptr->proximo_polinomio = aux->proximo_polinomio;
            aux->polinomio_anterior = ptr->polinomio_anterior;
            
            if (ptr->polinomio_anterior != NULL) 
                ptr->polinomio_anterior->proximo_polinomio = aux;
            if (aux->proximo_polinomio != NULL)
                aux->proximo_polinomio->polinomio_anterior = ptr;
            aux->proximo_polinomio = ptr;
            ptr->polinomio_anterior = aux;
            
            //como e um bubblesort, devo voltar ao inicio da lista
            while (ptr->polinomio_anterior != NULL) 
                ptr = ptr->polinomio_anterior;
        }
        else 
            ptr = ptr->proximo_polinomio;
        
    }
    //ao final do processo, basta rebobinar o polinomio
    while (ptr->polinomio_anterior != NULL) 
        ptr = ptr->polinomio_anterior;
    
    return ptr;
}

//funcao que imprime uma decomposicao
void imprime_decomposicao(vetor_decomp *decomposicao, tabela_literais *lista_literais)
{
    vetor_polinomios *percorre_polinomios;
    
    printf("\n %+3f*",decomposicao->resto_par->parametro);
    percorre_polinomios = decomposicao->poly_impares;
    while (percorre_polinomios!= NULL) 
    {
        printf("(");
        imprime_lista_expr_expandida(percorre_polinomios->polinomio->P, lista_literais);
        printf(")");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    
    printf(" %+3f*",decomposicao->resto_impar->parametro);
    percorre_polinomios = decomposicao->poly_pares;
    while (percorre_polinomios!= NULL) 
    {
        printf("(");
        imprime_lista_expr_expandida(percorre_polinomios->polinomio->P, lista_literais);
        printf(")");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
 /*   
    //imprimir os ids dos polinomios
    percorre_polinomios = decomposicao->poly_impares;
    while (percorre_polinomios!= NULL) 
    {
        printf(" %d",percorre_polinomios->polinomio->id);
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    
    printf(" ;");
    percorre_polinomios = decomposicao->poly_pares;
    while (percorre_polinomios!= NULL) 
    {
        printf(" %d",percorre_polinomios->polinomio->id);
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }*/
}

//funcao que retorna 1 caso as decomposicoes sejam equivalentes e 0 caso nao sejam
int compara_decomp(vetor_decomp *decomp1, vetor_decomp *decomp2)
{
    
    vetor_polinomios *ptr_impar1, *ptr_impar2;
    vetor_polinomios *ptr_par1, *ptr_par2;
    
    int flag;
    
    //comparo elemento a elemento de cada conjunto de polinomios para ver se sao iguais
    ptr_impar1 = decomp1->poly_impares;
    ptr_impar2 = decomp2->poly_impares;
    ptr_par1 = decomp1->poly_pares;
    ptr_par2 = decomp2->poly_pares;
        
    //inicializo minha flag 
    flag = 1; 
        
    while (ptr_par1 != NULL) 
    {
        //se algum for diferente, limpa a flag
        if (ptr_par1->polinomio->id != ptr_par2->polinomio->id) 
            flag = 0;
        if (ptr_impar1->polinomio->id != ptr_impar2->polinomio->id) 
            flag = 0;
        
        //atualiza os dois ponteiros
        ptr_par1 = ptr_par1->proximo_polinomio;
        ptr_par2 = ptr_par2->proximo_polinomio;
        ptr_impar1 = ptr_impar1->proximo_polinomio;
        ptr_impar2 = ptr_impar2->proximo_polinomio;
    }
    
    //se a flag for 0, testar com os vetores trocados
    if (!flag)
    {
        //comparo elemento a elemento de cada conjunto de polinomios para ver se sao iguais
        ptr_impar1 = decomp1->poly_impares;
        ptr_impar2 = decomp2->poly_impares;
        ptr_par1 = decomp1->poly_pares;
        ptr_par2 = decomp2->poly_pares;
        
        //inicializo minha flag 
        flag = 1;
        
        while (ptr_par1 != NULL) 
        {
            //se algum for diferente, limpa a flag
            if (ptr_par1->polinomio->id != ptr_impar2->polinomio->id) 
                flag = 0;
            if (ptr_impar1->polinomio->id != ptr_par2->polinomio->id) 
                flag = 0;
            
            //atualiza os dois ponteiros
            ptr_par1 = ptr_par1->proximo_polinomio;
            ptr_par2 = ptr_par2->proximo_polinomio;
            ptr_impar1 = ptr_impar1->proximo_polinomio;
            ptr_impar2 = ptr_impar2->proximo_polinomio;
        }
    }

    return flag;
    
}


//Funcao que testa se um par de polinomios pode ser combinado com outro para formar uma decomposição, retornando a decomposicao parcial
void encontra_decomp_recursiva_dummie(vetor_decomp *primario, vetor_decomp *secundario, int grau, int *total_tries, int profundidade)
{

    vetor_decomp *ptr_sementes = NULL;
    
    //aponto para o inicio do vetor sementes
    ptr_sementes = secundario;
    
    //vou tentar combinar o primario com TODOS os elementos do secundário
    while (ptr_sementes != NULL)
    {
        //testou P1 principal com P1 atual & P2 principal com P2 atual e deu certo, pois estamos vendo o pior caso
        //se ja chamou recursivamente um numero de entradas igual ao grau, é porque ja chegou no final
        if (profundidade == grau)
        {
            *total_tries+=1;
        }
        //caso contrario, fazer uma chamada recursiva
        else
        {
            encontra_decomp_recursiva_dummie(primario, ptr_sementes, grau, total_tries, profundidade + 1);
        }
        
        //atualiza o ponteiro
        ptr_sementes = ptr_sementes->prox_decomp;
    }
    profundidade-=1;
    
}

void encontra_decomp_dummie(vetor_sementes *entrada, int grau)
{
    vetor_sementes *secundario_ptr; //ponteiro para os loops internos
    
    vetor_decomp *decomp_atual; //lista que será gerada dentro do loop atual
    vetor_decomp *decomp_referencia; //lista de refrencia para o loop atual
    
    int total_decomp = 0;
    
    
    //inicialização da lista de referência
    decomp_referencia = copia_vetor_semente(entrada);
    decomp_atual = decomp_referencia;
    decomp_atual = NULL;
    
    //inicialização do ponteiro de busca secundário, pode ser que o par de polinomios inicial possa combinar com si mesmo
    secundario_ptr = entrada;
    
    //loop de construcao das decomposicoes
    while (decomp_referencia != NULL)
    {
        //para cada vetor incial de decomp referencia, encontro todas as decomposicoes que possam ser geradas a partir dele
        encontra_decomp_recursiva_dummie(decomp_referencia, decomp_referencia, grau, &total_decomp, 2 );
        
        //apagar o elemento testado e atualizar o ponteiro de referencia
        decomp_atual = decomp_referencia->prox_decomp;
        destroi_decomp(decomp_referencia);
        decomp_referencia = decomp_atual;
        if (decomp_referencia!= NULL)
            decomp_referencia->ant_decomp = NULL;
    }
    
    //imprimir numero de decomposicoes
    printf("\n\n o numero maximo de combinacoes testadas com meu algoritmo: %d", total_decomp);
    
}

//Funcao que testa se um par de polinomios pode ser combinado com outro para formar uma decomposição, retornando a decomposicao parcial
void encontra_decomp_recursiva_dummie_mulder(vetor_decomp *primario, vetor_polinomios *secundario, int grau, int *total_tries, int profundidade)
{
    
    vetor_polinomios *ptr_sementes = NULL;
    
    //aponto para o inicio do vetor sementes
    ptr_sementes = secundario;
    
    //vou tentar combinar o primario com TODOS os elementos do secundário
    while (ptr_sementes != NULL)
    {
        //testou P1 principal com P1 atual & P2 principal com P2 atual e deu certo, pois estamos vendo o pior caso
        //se ja chamou recursivamente um numero de entradas igual ao grau, é porque ja chegou no final
        if (profundidade == grau)
        {
            *total_tries+=1;
        }
        //caso contrario, fazer uma chamada recursiva
        else
        {
            encontra_decomp_recursiva_dummie_mulder(primario, ptr_sementes, grau, total_tries, profundidade + 1);
        }
        
        //atualiza o ponteiro
        ptr_sementes = ptr_sementes->proximo_polinomio;
    }
    profundidade-=1;
    
}

void encontra_decomp_dummie_mulder(vetor_sementes *entrada, vetor_polinomios *lista_polinomios, int grau)
{
    vetor_sementes *secundario_ptr; //ponteiro para os loops internos
    
    vetor_decomp *decomp_atual; //lista que será gerada dentro do loop atual
    vetor_decomp *decomp_referencia; //lista de refrencia para o loop atual
    
    int total_decomp = 0;
    
    
    //inicialização da lista de referência
    decomp_referencia = copia_vetor_semente(entrada);
    decomp_atual = decomp_referencia;
    decomp_atual = NULL;
    
    //inicialização do ponteiro de busca secundário, pode ser que o par de polinomios inicial possa combinar com si mesmo
    secundario_ptr = entrada;
    
    //loop de construcao das decomposicoes
    while (decomp_referencia != NULL)
    {
        //para cada vetor incial de decomp referencia, encontro todas as decomposicoes que possam ser geradas a partir dele
        encontra_decomp_recursiva_dummie_mulder(decomp_referencia, lista_polinomios, grau, &total_decomp, 2 );
        
        //diminui em 1 o total_decomp, pois ele nao combina com o outro semente.
        total_decomp--;
        
        //realiza denovo a mesma funcao, pois ele separa em dois problemas distintos
        encontra_decomp_recursiva_dummie_mulder(decomp_referencia, lista_polinomios, grau, &total_decomp, 2 );
        
        //diminui em 1 o total_decomp, pois ele nao combina com o outro semente.
        total_decomp--;
        
        //apagar o elemento testado e atualizar o ponteiro de referencia
        decomp_atual = decomp_referencia->prox_decomp;
        destroi_decomp(decomp_referencia);
        decomp_referencia = decomp_atual;
        if (decomp_referencia!= NULL)
            decomp_referencia->ant_decomp = NULL;
    }
    
    //imprimir numero de decomposicoes
    printf("\n\n o numero maximo de combinacoes testadas com metodo mulder: %d", total_decomp);
    
}





