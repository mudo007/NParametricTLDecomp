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
#include "main.h"

int main(int argc, char *argv[])
{
    //Vari‡veis
    char    equacao_entrada[302];
    int i;
    token   *lista_token = NULL;
    tabela_literais *lista_literais = NULL;
    token *expressao_RPN = NULL;
    arvore_expr *arvore = NULL;
    lista_expr  *expr_expandida = NULL;
    lista_expr  *expr_simplificada = NULL;
    lista_expr *polinomio_base;
    vetor_polinomios *lista_polinomios = NULL;
    vetor_polinomios *percorre_polinomios;
    vetor_sementes  *lista_sementes = NULL;
    vetor_sementes  *percorre_sementes = NULL;
    int             contador;
    vetor_decomp    *decomposicoes_encontradas;
    int             lim_inferior;
    int             lim_superior;
    time_t   cronometro;

    
    
//*********
// INICIO *
//*********

    printf("NParametricTLDecomp  Copyright (C) 2014  Diogo Andrade \nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under certain conditions.\n\nThis application performs non-parametric Translinear decomposition onto a homogeneous (all monomials having the same degree) multivariate polynomial.\nThe results are suitable for Translinear analog circuit realization with proper adjustments.\nIf an error occurs, send a brief description with the polynomial inserted to diogo007@gmail.com.\n\r \
           Below, type the polynomial to be decomposed (100 characters maximum). Coefficients should be only integer numbers. \n\r\t Example:\n\r\t \
           x^3 + x^2*y - x*z^2 + y*z^2\n\r  (Hit \"enter\" to use it)");
    
    //leitura da string de entrada
    
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);
    }
    
    //equacao padrao, caso o usuario aperte enter direto. pode ser igual a nova linha ou retorno de carro para tratar como o sistema operacional reconhece o ENTER em varias plataformas diferentes
    if (*equacao_entrada == '\n' || *equacao_entrada == '\r')
        sprintf(equacao_entrada, "x^3 + x^2*y - x*z^2 + y*z^2");
    
//*****************
// ANALISE LEXICA *
//*****************
   
    //leitura dos tokens
    if((lista_token = le_tokens(equacao_entrada)) == NULL)
    {
        //limpeza de ponteiros
        system("PAUSE");	
        return 0;      
    }
    
    //converte os literais em codigos - apenas de exemplo
    constroi_tabela_literais(&lista_literais,lista_token);
           
#if defined DEBUG_LEXICO
    //caso a leitura dos tokens tenha sido correta, imprimir os tokens 
    imprime_tokens(lista_token);
#endif
        
//********************
// ANALISE SINTçTICA *
//********************
    
    //Cria‹o das pilhas de avalia‹o de express‹o
           
    expressao_RPN = constroi_lista_expr(lista_token);    
    arvore = constroi_arvore_expr(expressao_RPN);
#if defined DEBUG_EXPR
    imprime_arvore_expr(arvore);
    imprime_lista_expr(expressao_RPN);
#endif
    
    //Expans‹o da arvore de express›es
    expr_expandida = constroi_lista_expressoes_exp(arvore);
#if defined DEBUG_EXPAND
    printf("\n equacao expandida:");
    imprime_lista_expr_expandida(expr_expandida,lista_literais);
#endif   
 
//********************
// ANALISE SEMåNTICA *
//********************

    //Simplifica‹o da express‹o
    expr_simplificada = simplifica_expr_expandida(expr_expandida);
    destroi_lista_expr_expandida(expr_expandida);
#if defined DEBUG_SIMPLIFY
    printf("\n equacao simplificada:");
    imprime_lista_expr_expandida(expr_simplificada,lista_literais);
#endif 
    
    //ordena‹o lexdeg
   expr_simplificada = lexdegbubblesort(expr_simplificada);
#if defined DEBUG_SIMPLIFY
    printf("\n equacao reordenada lex:");
    imprime_lista_expr_expandida(expr_simplificada,lista_literais);
#endif 
    

/***********************
 Base-Polynomial generation and reduction
 ***************************/

    //imprimir a lista de literais e construir o polinomio base
    
    polinomio_base = gera_polinomio_base(lista_literais);

    //leitura dos limites superiores e inferiores dos coeficientes
    printf("\n Insert the inferior limit of variables coefficients (standard -1)");
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);  
    }
    
    //tratamento do nœmero inserido
    //se apertar enter direto, Ž -1
    if (*equacao_entrada == '\n' || *equacao_entrada == '\r')
        lim_inferior = -1;
    else
    {
        lim_inferior = atoi(equacao_entrada);
    }
    
    //limite superior
    printf("\n Insert the superior limit of variables coefficients (standard +1)");
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);  
    }
    
    //tratamento do nœmero inserido
    //se apertar enter direto, Ž -1
    if (*equacao_entrada == '\n' || *equacao_entrada == '\r')
        lim_superior = 1;
    else
    {
        lim_superior = atoi(equacao_entrada);
    }

    //imprimir os parametros colhidos
    printf("\n Expanded, simplified and re-ordered Polynomial:");
    imprime_lista_expr_expandida(expr_simplificada, lista_literais);
    printf("\n coefficient inferior limit: %d", lim_inferior);
    printf("\n coefficient superior limit: %d\n", lim_superior);
    
    
    //construcao dos vetores, aqui com coeficientes especificados pelo usuario
    lista_polinomios = gera_vetor(lista_polinomios, polinomio_base, polinomio_base,lim_inferior, lim_superior);
    
    //rebobina a lista
    while (lista_polinomios->polinomio_anterior != NULL) 
    {
        lista_polinomios = lista_polinomios->polinomio_anterior;
    }
    
    //contagem de polin™mios T0
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Total number of Base-Polynomials generated (T0 set) is: %d\n",contador);

    //elimina o polinomio nulo
    lista_polinomios = elimina_zero(lista_polinomios);
    
    //elimina os polinomios inteiramente negativos
    lista_polinomios = remove_polinomios_negativos(lista_polinomios);

    //contagem de polin™mios T1
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Number of Base-Polynomials after removing strict negative BP's (T1 set) is: %d\n",contador);

    //elimina os polinomios redundantes
    lista_polinomios = remove_polinomios_redundantes(lista_polinomios);

    //contagem de polin™mios T2
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Number of Base-Polynomials after redundancy removal (T2 set) is: %d\n",contador);

    //gerar vetor T3
    lista_sementes = gera_vetor_semente(lista_polinomios, expr_simplificada);

    //contagem de pares de polin™mios T3
    percorre_sementes = lista_sementes;
    contador = 0;
    while (percorre_sementes != NULL)
    {
        contador++;
        percorre_sementes = percorre_sementes->conjunto_prox;
    }
    printf("\n Number of Base-Polynomial pairs (T3 set) is: %d\n",contador);

    //gerar T4 ap—s T3 ter sido gerado
    lista_polinomios = reconta_polinomios(lista_sementes,lista_polinomios);

    //contagem de T4
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Number of Base-Polynomials after seed generation (reduced set T4) is: %d\n",contador);


/***********************
Recursive Division By BP pairs
***********************/

    //Escolha de qual algoritmo ser‡ utilizado
    //leitura dos limites superiores e inferiores dos coeficientes
    printf("\n Select recursive division Algorithm:[A]ndrade, [M]ulder fast, Mulder Memory [S]afe: (Default: A)");
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);
    }
    
    // Se op‹o estranha for utilizada, rea;lizar algoritmo de diogo
    if (equacao_entrada[0] != 'M' && equacao_entrada[0] != 'm' && equacao_entrada[0] != 'A' && equacao_entrada[0] != 'a' &&
        equacao_entrada[0] != 'S' && equacao_entrada[0] != 's')
    {
        equacao_entrada[0] = 'A';
    }
    
    //preparar para a execu‹o
    global_num_parfrac = 0;
    decomposicoes_encontradas = NULL;
    inicia_cronometro(&cronometro);
    
    //selecionar o algoritmo a ser utilizado
    switch (equacao_entrada[0])
    {
        case 'A':
        case 'a':
            printf("\n Running Andrade's algorithm...");
            decomposicoes_encontradas = encontra_decomp(lista_sementes,deg(expr_simplificada), expr_simplificada, lista_literais);
            break;
            
        case 'M':
        case 'm':
            printf("\n Running Mulder's algorithm (fast version)...");
            decomposicoes_encontradas = encontra_decomp_mulder(lista_polinomios, deg(expr_simplificada), expr_simplificada, lista_literais);
            break;
            
        case 'S':
        case 's':
            printf("\n Running Mulder's algorithm (memory safe version)...");
            decomposicoes_encontradas = encontra_decomp_mulder_safe(lista_polinomios, deg(expr_simplificada), expr_simplificada, lista_literais);
            break;
            
    }
    para_cronometro(&cronometro);
    
    //verifica quantos PB's de fato fizeram parte de alguma decomposicao
    lista_polinomios = remove_polinomios_nao_pertencentes(decomposicoes_encontradas,lista_polinomios);
    
    //contagem p—s-encontrar decomposicoes
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    
    printf("\n Number of Base-Polynomials effectively used in any translinear polynomial (extremely reduced set T5) is: %d",contador);
    
    //imprime apenas os polinomios utilizados em alguma decomposicao
    printf("\n The base polynomials used in any of the translinear polynomials found are:\n");
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        imprime_lista_expr_expandida(percorre_polinomios->polinomio->P, lista_literais);
        //imprime o identificador do polinomio
        printf(" \t%d\n",percorre_polinomios->polinomio->id);
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }

    
    //elimina o vetor de polin™mios translineares
    destroi_vetor_decomp(decomposicoes_encontradas);
    destroi_lista_sementes(lista_sementes);
    decomposicoes_encontradas = NULL;
    lista_sementes = NULL;
    
    //destroi o vetor de polinomios
    while (lista_polinomios != NULL)
    {
        percorre_polinomios = lista_polinomios->proximo_polinomio;
        destroi_lista_expr_expandida(lista_polinomios->polinomio->P);
        free(lista_polinomios->polinomio);
        free(lista_polinomios);
        lista_polinomios = percorre_polinomios;
    }
    
    //destroi estruturas de dados auxiliares
    destroi_lista_expr_expandida(polinomio_base);
    destroi_tabela_literais(lista_literais);
    destroi_lista(lista_token);
    destroi_arvore_expr(arvore);
    destroi_lista_expr(expressao_RPN);
    destroi_lista_expr_expandida(expr_simplificada);
    printf("\nPress any key to continue...");
    getchar();
    return 0;
}
