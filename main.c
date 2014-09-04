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
#if defined DEBUG_POLYDIV || defined DEBUG_PARFRAC || defined (DEBUG_MEMORIA)
    char    divisor_entrada[102];
    token   *token_divisor;
    token *expressao_RPN_divisor = NULL;
    arvore_expr *arvore_divisor = NULL;
    lista_expr  *expr_expandida_divisor = NULL;
    lista_expr  *expr_simplificada_divisor = NULL;
    lista_expr  *quociente = NULL;
    lista_expr  *resto = NULL;
    lista_expr  *prova_real;

#endif
#if defined (DEBUG_PARFRAC) || defined (DEBUG_MEMORIA)
    lista_expr *divisor_P1 = NULL;
    lista_expr *divisor_P2 = NULL;
    lista_expr *numerador_a = NULL;
    lista_expr *numerador_b = NULL;
    lista_expr *quociente_parfrac = NULL;
    lista_expr *aux1;
    lista_expr *aux2;
#endif
#if defined (DEBUG_VECTOR_GEN)
    lista_expr *polinomio_base;
    vetor_polinomios *lista_polinomios = NULL;
    vetor_polinomios *percorre_polinomios;
    vetor_sementes  *lista_sementes = NULL;
    vetor_sementes  *percorre_Sementes = NULL;
    lista_expr      *prova_1, *prova_2, *prova_3;
    int             contador;
    vetor_decomp    *decomposicoes_encontradas;
    vetor_decomp    *ptr_decomposicoes;
    int             ok_count = 0;
    int             fail_count = 0;
    int             lim_sup = 1;
    int             lim_inf = -1;
    
    
#endif
    
#if !defined (DEBUG_POLYDIV) && !defined (DEBUG_PARFRAC) &&!defined (DEBUG_VECTOR_GEN) && !defined (DEBUG_MEMORIA)
    lista_expr *polinomio_base;
    vetor_polinomios *lista_polinomios = NULL;
    vetor_polinomios *percorre_polinomios;
    vetor_sementes  *lista_sementes = NULL;
    vetor_sementes  *percorre_sementes = NULL;
   // lista_expr      *prova_1, *prova_2, *prova_3;
    int             contador,contador2;
    vetor_decomp    *decomposicoes_encontradas;
//    vetor_decomp    *ptr_decomposicoes;
//    int             ok_count = 0;
//    int             fail_count = 0;
    int             lim_inferior;
    int             lim_superior;
    unsigned int    diogo_parfrac, mulder_parfrac;
//    vetor_correlacao_direta *correlacoes;
    time_t   cronometro;
#endif
    

    int i;
    token   *lista_token = NULL;
    tabela_literais *lista_literais = NULL;
    token *expressao_RPN = NULL;
    arvore_expr *arvore = NULL;
    lista_expr  *expr_expandida = NULL;
    lista_expr  *expr_simplificada = NULL;
    
 /*   //teste de tipos
    
    printf("\n normal: %lu",sizeof(int));
    printf("\n 64-bit: %lu",sizeof(int64_t));
    printf("\n 32-bit: %lu",sizeof(int32_t));
    printf("\n 16-bit: %lu",sizeof(int16_t));
    printf("\n 8-bit: %lu",sizeof(int8_t));
    printf("\n char: %lu",sizeof(char));
    
    return 0;*/
    
    
#if !defined (DEBUG_POLYDIV) && !defined (DEBUG_PARFRAC) &&!defined (DEBUG_VECTOR_GEN) && !defined (DEBUG_MEMORIA)
    //inicio
    equacao_entrada[0] = '\0';
    if (argv[1] != NULL)
    {
        i = 1;
        while (argv[i] != NULL)
        {
            strcat(equacao_entrada, (char *)argv[i]);
            i++;
        }
        lim_inferior = -3;
        lim_superior = +3;
        
    }
    else
    {
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
    }
    
    
#endif
    
#if defined (DEBUG_VECTOR_GEN)
   
     sprintf(equacao_entrada, "x^3 + x^2*y - x*z^2 + y*z^2");
    lim_inf = -1;
    lim_sup = +1;
    //sprintf(equacao_entrada, "X^3 + X^2*Y - X*Z^2 + Y*Z^2");
    

#endif

    //teste da divisao polinomial
#if defined (DEBUG_POLYDIV)
    //polinomio numerador
    sprintf(equacao_entrada, "36*(6*x +5*y +6*z )*(x +y )*(-x +5*y +z ) -36*(5*x +5*y +4*z )*(5*x +5*y +3*z)*y" );
    //sprintf(equacao_entrada, "x^3 + x^2*y - x*z^2 + y*z^2");
    
#endif
    
#if defined (DEBUG_PARFRAC)
    //teste da expansao em fracoes parciais
    //as entradas ser‹o igual no teste da divis‹o polinomial, com a diferena de que o divisor ser‡ quebrado em P1 e P2
    //A divis‹o ser‡ feita assim como no teste da divis‹o polinomial, e o resto ser‡ processado junto com P1 e P2 para 
    //expans‹o em fracoes parciais
    
    //polinomio numerador
    //sprintf(equacao_entrada, "X^3 + X^2*Y - X*Z^2 + Y*Z^2");
    //sprintf(equacao_entrada, "X");
    //sprintf(equacao_entrada, "x^5*(y-2*z+2*k)-y^2*(x+z-k)^2*(x-2*z)*k");
     sprintf(equacao_entrada, "x^3 + x^2*y - x*z^2 + y*z^2");
#endif
    
#if defined (DEBUG_MEMORIA)
    sprintf(equacao_entrada, "(2*x-y)*(-x+2*z)*(z-2*y)*(x+y) - (z+y)*(z-y)*(x-z)*(y+2*x)");
#endif
    
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
    
    //Simplifica‹o da ‡rvore de express›es
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

#if defined (DEBUG_POLYDIV) 
    //polinomio denominador
    sprintf(divisor_entrada, "(6*x +5*y +6*z )");
    //sprintf(divisor_entrada, " 2*(x +y +z )*(x )*(z ) -(x +y )*(x +z )*(x +z )");
    token_divisor = le_tokens(divisor_entrada);
    constroi_tabela_literais(&lista_literais, token_divisor);
    expressao_RPN_divisor = constroi_lista_expr(token_divisor);    
    arvore_divisor = constroi_arvore_expr(expressao_RPN_divisor);
    expr_expandida_divisor = constroi_lista_expressoes_exp(arvore_divisor);
    expr_simplificada_divisor = simplifica_expr_expandida(expr_expandida_divisor);
    destroi_lista_expr_expandida(expr_expandida_divisor);
    expr_simplificada_divisor = lexdegbubblesort(expr_simplificada_divisor);
#endif
    

    
#if defined (DEBUG_PARFRAC) ||defined(DEBUG_MEMORIA)
    //polinomios denominadores
    //P1
    //sprintf(divisor_entrada, "I0 +Iin");
    sprintf(divisor_entrada, "-x+z");
    token_divisor = le_tokens(divisor_entrada);
    constroi_tabela_literais(&lista_literais, token_divisor);
    expressao_RPN_divisor = constroi_lista_expr(token_divisor);    
    arvore_divisor = constroi_arvore_expr(expressao_RPN_divisor);
    expr_expandida_divisor = constroi_lista_expressoes_exp(arvore_divisor);
    expr_simplificada_divisor = simplifica_expr_expandida(expr_expandida_divisor);
    destroi_lista_expr_expandida(expr_expandida_divisor);
    expr_simplificada_divisor = lexdegbubblesort(expr_simplificada_divisor);
    //salva P1
    divisor_P1 = expr_simplificada_divisor;

    
    //reseta os ponteiros 
    destroi_lista(token_divisor);	
    destroi_arvore_expr(arvore_divisor);
    destroi_lista_expr(expressao_RPN_divisor);
    expr_expandida_divisor = NULL;
    expr_simplificada_divisor = NULL;

    
    //P2
   // sprintf(divisor_entrada, "Iin+I0+Iout");
    sprintf(divisor_entrada, "x");
    token_divisor = le_tokens(divisor_entrada);
    constroi_tabela_literais(&lista_literais, token_divisor);
    expressao_RPN_divisor = constroi_lista_expr(token_divisor);    
    arvore_divisor = constroi_arvore_expr(expressao_RPN_divisor);
    expr_expandida_divisor = constroi_lista_expressoes_exp(arvore_divisor);
    expr_simplificada_divisor = simplifica_expr_expandida(expr_expandida_divisor);
    destroi_lista_expr_expandida(expr_expandida_divisor);
    expr_simplificada_divisor = lexdegbubblesort(expr_simplificada_divisor);
    
    //salva P2
    divisor_P2 = expr_simplificada_divisor;
    
    //multiplica P1 e P2 e procede normalmente a divis‹o
    expr_simplificada_divisor = NULL;
    expr_expandida_divisor = NULL;
    expr_expandida_divisor = multiplica_expr(divisor_P1 , divisor_P2 );
    expr_simplificada_divisor = simplifica_expr_expandida(expr_expandida_divisor);
    destroi_lista_expr_expandida(expr_expandida_divisor);
    expr_simplificada_divisor = lexdegbubblesort(expr_simplificada_divisor);

#endif

#if defined DEBUG_POLYDIV || defined DEBUG_PARFRAC 


    //divisao polinomial
    if (polydiv(expr_simplificada, expr_simplificada_divisor, &quociente, &resto))
    {
        
        
        printf("\n divisao realizada com sucesso");
        printf("\n equacao de entrada: %s ",equacao_entrada);
        printf("\n equacao do divisor: %s ",divisor_entrada);
        printf("\n dividendo: ");
        imprime_lista_expr_expandida(expr_simplificada,lista_literais);
        printf("\n divisor:   ");
        imprime_lista_expr_expandida(expr_simplificada_divisor,lista_literais);
        printf("\n resto:     ");
        imprime_lista_expr_expandida(resto,lista_literais);
        printf("\n quociente: ");
        imprime_lista_expr_expandida(quociente,lista_literais);
#if !defined (DEBUG_PARFRAC) //preciso dos ponteiros de resto e quocientes intactos para realizar a expans‹o
        //tirar a prova real. dividendo - (Quociente * Divisor + resto)  = 0
        prova_real = multiplica_expr(quociente, expr_simplificada_divisor);
        soma_expr(prova_real,resto);
        subtrai_expr(&expr_simplificada, prova_real);
        printf("\n resultado da prova real: ");
        expr_simplificada = simplifica_expr_expandida(expr_simplificada);
        imprime_lista_expr_expandida(expr_simplificada,lista_literais);
#endif



    }
    else
        printf("\n erro na divisao polinomial");
        
#endif
    
#if defined (DEBUG_MEMORIA) 
    resto = expr_simplificada;
    while(TRUE)
    {

        
#endif
    
#if defined (DEBUG_PARFRAC) || defined (DEBUG_MEMORIA)
    if(partial_fraction_expansion(resto, divisor_P1, divisor_P2, &quociente_parfrac, &numerador_a, &numerador_b))
    {
#if !defined (DEBUG_MEMORIA) 
        printf("\n Expansao em fracoes parciais realizada com sucesso\n (");
        imprime_lista_expr_expandida(resto, lista_literais);
        printf(")/(");
        imprime_lista_expr_expandida(divisor_P1, lista_literais);
        printf(")*(");
        imprime_lista_expr_expandida(divisor_P2, lista_literais);
        printf(") = ");
        imprime_lista_expr_expandida(quociente,lista_literais);
        printf("+ (");
        imprime_lista_expr_expandida(numerador_a, lista_literais);
        printf(")/(");
        imprime_lista_expr_expandida(divisor_P1, lista_literais);
        printf(")+(");
        imprime_lista_expr_expandida(numerador_b, lista_literais);
        printf(")/(");
        imprime_lista_expr_expandida(divisor_P2, lista_literais);
        printf(")");
        
        printf("\nprova real = ");
        //tirar a prova real. 
        prova_real = multiplica_expr(numerador_a, divisor_P2);
        aux1 = multiplica_expr(numerador_b, divisor_P1);
        soma_expr(prova_real,aux1);
        
        aux1 = multiplica_expr(quociente, divisor_P1);
        aux2 = multiplica_expr(aux1, divisor_P2);
        soma_expr(prova_real, aux2);
        
        destroi_lista_expr_expandida(aux1);
        subtrai_expr(&expr_simplificada, prova_real);
        printf("\n resultado da prova real: ");
        expr_simplificada = simplifica_expr_expandida(expr_simplificada);
        imprime_lista_expr_expandida(expr_simplificada,lista_literais);
        
    }
    else
       printf("\n Nao foi possivel encontrar uma expansao satisfatoria");
#else
        aux1 = NULL;
    }
#endif



    
#endif
 
#if defined (DEBUG_MEMORIA) 
    //limpeza de ponteiros

    destroi_lista_expr_expandida(quociente_parfrac);
    quociente_parfrac = NULL;
    if (numerador_a != NULL)
        destroi_lista_expr_expandida(numerador_a);
    if (numerador_b != NULL) 
        destroi_lista_expr_expandida(numerador_b);
    numerador_a = NULL;
    numerador_b = NULL;
}
#endif


    

    //limpeza de ponteiros
    
#if defined (DEBUG_POLYDIV) || defined (DEBUG_PARFRAC) 
    destroi_lista(token_divisor);	
    destroi_arvore_expr(arvore_divisor);
    destroi_lista_expr(expressao_RPN_divisor);
    destroi_lista_expr_expandida(expr_simplificada_divisor);
    destroi_lista_expr_expandida(resto);
    destroi_lista_expr_expandida(quociente);
    resto = NULL;
    quociente = NULL;
    
#endif
    
#if defined (DEBUG_PARFRAC)
    destroi_lista_expr_expandida(divisor_P1);
    destroi_lista_expr_expandida(divisor_P2);
    if (numerador_a != NULL)
        destroi_lista_expr_expandida(numerador_a);
    if (numerador_b != NULL) 
        destroi_lista_expr_expandida(numerador_b);
    numerador_a = NULL;
    numerador_b = NULL;
#endif


    
#if defined (DEBUG_VECTOR_GEN)
    //imprimir a lista de literais e construir o polinomio base
        
    polinomio_base = gera_polinomio_base(lista_literais);
    
    printf("\n equacao:");
    imprime_lista_expr_expandida(expr_simplificada, lista_literais);
    
    //imprime o polinomio base
    printf("\n o polinomio base gerado e: \n");
    imprime_lista_expr_expandida(polinomio_base, lista_literais);
    
    //construcao dos vetores, aqui com coeficientes entre -1 e 1
    lista_polinomios = gera_vetor(lista_polinomios, polinomio_base, polinomio_base, lim_inf, lim_sup);
    
    //rebobina a lista
    while (lista_polinomios->polinomio_anterior != NULL) 
    {
        lista_polinomios = lista_polinomios->polinomio_anterior;
    }
    
    //elimina o polinomio nulo
    lista_polinomios = elimina_zero(lista_polinomios);
    
    //elimina os polinomios inteiramente negativos
    lista_polinomios = remove_polinomios_negativos(lista_polinomios);
    
    //elimina os polinomios redundantes
    lista_polinomios = remove_polinomios_redundantes(lista_polinomios);
    
    //imprime as combina›es
    printf("\nOs vetores gerados e seus identificadores respsectivos sao:\n");
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
    printf("\n numero de vetores gerados: %d",contador);
    getchar();
    
    //gerar as sementes ->aqui j‡ est‡ dando pau no windows
    lista_sementes = gera_vetor_semente(lista_polinomios, expr_simplificada);

    //aqui eu j‡ posso destruir a lista de polinomios
    //destroi o vetor de polinomios 
    while (lista_polinomios != NULL) 
    {
        percorre_polinomios = lista_polinomios->proximo_polinomio;
        destroi_lista_expr_expandida(lista_polinomios->polinomio->P);
        free(lista_polinomios->polinomio);
        free(lista_polinomios);
        lista_polinomios = percorre_polinomios;
    }
    destroi_lista_expr_expandida(polinomio_base);
    
    
    //imprime as sementes
    percorre_Sementes = lista_sementes;
    contador = 0;
    while (percorre_Sementes!= NULL) 
    {
        contador++;
        printf("\n");
        imprime_lista_expr_expandida(percorre_Sementes->P1.P, lista_literais);
        printf(" e ");
        imprime_lista_expr_expandida(percorre_Sementes->P2.P, lista_literais);
        //tira a prova real
        printf("prova = ");
        
        //faz a multiplica‹o cruzada
        prova_3 = multiplica_expr(percorre_Sementes->P1.P, percorre_Sementes->R2);
        prova_2 = multiplica_expr(percorre_Sementes->P2.P, percorre_Sementes->R1);
        
        //soma os dois produtos
        soma_expr(prova_3, prova_2);
        
        //multiplica o quociente por P1 e P2
        prova_1 = multiplica_expr(percorre_Sementes->quociente, percorre_Sementes->P1.P);
        prova_2 = multiplica_expr(prova_1, percorre_Sementes->P2.P);
        
        destroi_lista_expr_expandida(prova_1);
        
        //soma com o quociente multiplicado
        soma_expr(prova_3,prova_2);
        
        //subtrai tudo do polinomio de entrada
        prova_2 = copia_lista_expr(expr_simplificada);
        subtrai_expr(&prova_3, prova_2);
        prova_2 = simplifica_expr_expandida(prova_3);
        imprime_lista_expr_expandida(prova_2, lista_literais);
        
        printf("\n");
        destroi_lista_expr_expandida(prova_3);
        destroi_lista_expr_expandida(prova_2);
        percorre_Sementes = percorre_Sementes->conjunto_prox;
    }
    printf("\n numero de pares encontrados: %d",contador);
        getchar();

    //realiza a decomposi‹o 
    decomposicoes_encontradas = encontra_decomp(lista_sementes, deg(expr_simplificada),expr_simplificada, lista_literais);
   /* 
    //elimina as redundantes
    elimina_decomp_redundantes(decomposicoes_encontradas);
    
    
    //imprimir a lista de polinomios e fazer prova real
    ptr_decomposicoes = decomposicoes_encontradas;
    if (decomposicoes_encontradas == NULL) 
    {
        printf("\n nao foram encontradas decomposicoes validas");
    }
    else
    {
        while (ptr_decomposicoes != NULL) 
        {
            //tirar a prova real
            if (prova_real(ptr_decomposicoes, expr_simplificada)) 
            {
                imprime_decomposicao(ptr_decomposicoes,lista_literais);
                
                
                ok_count++;

            }
            else
                fail_count++;
                
            ptr_decomposicoes = ptr_decomposicoes->prox_decomp;
            
        }
        printf("\n decomposicoes validas encontradas = %d",ok_count);
        printf("\n decomposicoes invalidas encontradas = %d",fail_count);
        //imprimir o tamanho em bytes de uma decomposicao
        ok_count = 0;
        ptr_decomposicoes = decomposicoes_encontradas;
        while(ptr_decomposicoes!= NULL)
        {
            ok_count+= decomp_size(ptr_decomposicoes);
            ptr_decomposicoes = ptr_decomposicoes->prox_decomp;
        }
        printf("\n tamanho de uma decomposicao em bytes: %d", ok_count);
        
        getchar();
    }*/
    //elimina o lixo da memoria
    destroi_vetor_decomp(decomposicoes_encontradas);
    destroi_lista_sementes(lista_sementes);
printf("\n verificar memoria");
getchar();

    
#endif
    
//processamento principal do aplicativo
#if !defined (DEBUG_POLYDIV) && !defined (DEBUG_PARFRAC) &&!defined (DEBUG_VECTOR_GEN) && !defined (DEBUG_MEMORIA)

/***********************
 Base-Polynomial generation and reduction
 ***************************/

    //imprimir a lista de literais e construir o polinomio base
    
    polinomio_base = gera_polinomio_base(lista_literais);
if (argv[1] == NULL)
{
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
}
    //imprimir os parametros colhidos
    printf("\n Expanded, simplified and re-ordered Polynomial:");
    imprime_lista_expr_expandida(expr_simplificada, lista_literais);
    printf("\n coefficient inferior limit: %d", lim_inferior);
    printf("\n coefficient superior limit: %d\n", lim_superior);
    
//    printf("\n Press any key to continue...\n");
//    getchar();
    
    //construcao dos vetores, aqui com coeficientes especificados pelo usuario
    lista_polinomios = gera_vetor(lista_polinomios, polinomio_base, polinomio_base,lim_inferior, lim_superior);
    
    //rebobina a lista
    while (lista_polinomios->polinomio_anterior != NULL) 
    {
        lista_polinomios = lista_polinomios->polinomio_anterior;
    }


    //inserir contagem (antes de tudo)
    // printf("\nOs vetores gerados sao:\n");
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Total number of Base-Polynomials generated (T0 set) is: %d\n",contador);
 //   printf("\n Press any key to continue...\n");
 //   getchar();

inicia_cronometro(&cronometro);
    //elimina o polinomio nulo
    lista_polinomios = elimina_zero(lista_polinomios);
    
    //elimina os polinomios inteiramente negativos
    lista_polinomios = remove_polinomios_negativos(lista_polinomios);

    //inserir contagem (apos os estritamento negativos)
    // printf("\nOs vetores gerados sao:\n");
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Number of Base-Polynomials after removing strict negative BP's (T1 set) is: %d\n",contador);
para_cronometro(&cronometro);
//    printf("\n Press any key to continue...\n");
//    getchar();

    //elimina os polinomios redundantes
inicia_cronometro(&cronometro);
    lista_polinomios = remove_polinomios_redundantes(lista_polinomios);

    //contagem apos-redundantes
    // printf("\nOs vetores gerados sao:\n");
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Number of Base-Polynomials after redundancy removal (T2 set) is: %d\n",contador);
para_cronometro(&cronometro);
//    printf("\n Press any key to continue...\n");
//    getchar();

    //gerar as sementes
inicia_cronometro(&cronometro);
    lista_sementes = gera_vetor_semente(lista_polinomios, expr_simplificada);

    //conta quantos vetores-semente existem
    percorre_sementes = lista_sementes;
    contador2 = 0;
    while (percorre_sementes != NULL)
    {
        contador2++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_sementes = percorre_sementes->conjunto_prox;
    }
    printf("\n Number of Base-Polynomial pairs (T3 set) is: %d\n",contador2);
para_cronometro(&cronometro);
//    printf("\n Press any key to continue...\n");
//    getchar();

/*//imprime as sementes
percorre_sementes = lista_sementes;
contador = 0;
while (percorre_sementes!= NULL)
{
    contador++;
    printf("\n");
    imprime_lista_expr_expandida(percorre_sementes->P1.P, lista_literais);
    printf(" e ");
    imprime_lista_expr_expandida(percorre_sementes->P2.P, lista_literais);
    
    printf("\n");

    percorre_sementes = percorre_sementes->conjunto_prox;
}
printf("\n numero de pares encontrados: %d",contador);
getchar(); */

    //refazer T4 ap—s T3 ter sido gerado
    lista_polinomios = reconta_polinomios(lista_sementes,lista_polinomios);

    //contagem p—s-geracao de T3
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }
    printf("\n Number of Base-Polynomials after seed generation (reduced set T4) is: %d\n",contador);
//    printf("\n Press any key to continue...\n");
//    getchar();
/*
    //imprime apenas os polinomios utilizados em alguma decomposicao
    printf("\n The base polynomials from reduced set T4 are:\n");
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        imprime_lista_expr_expandida(percorre_polinomios->polinomio->P, lista_literais);
        //imprime o identificador do polinomio
        printf(" \t%d\n",percorre_polinomios->polinomio->id);
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }*/
/*

    //calculo do pior caso para o meu algoritmo
    lim_superior = contador2*contador2;
    printf("\n Worst case number of attempts using this algorithm is: %d",lim_superior);

    //calculo de pior caso para o mulden

    lim_superior = contador2*contador;
    printf("\n Worst case number of attempts using Mulder's algorithm is: %d",lim_superior);
    printf("\n Press any key to continue...\n");
//    getchar();*/

/***********************
Recursive Division By BP pairs
***********************/

    //First run mulder's algorithm
    printf("\n Performing Mulder's algorithm...");
    global_num_parfrac = 0;
inicia_cronometro(&cronometro);
decomposicoes_encontradas = NULL;
    decomposicoes_encontradas = encontra_decomp_mulder(lista_polinomios, deg(expr_simplificada), expr_simplificada, lista_literais);
para_cronometro(&cronometro);

    //stores how many partial fraction expansion were performed with mulder's algorithm
  //  mulder_parfrac = global_num_parfrac;

    //cleans-up memory to perform Diogo's algorithm
 //   destroi_vetor_decomp(decomposicoes_encontradas);
//    decomposicoes_encontradas = NULL;
/*    global_num_parfrac = 0;

    //performs diogo's decomposition
//    printf("\n Performing Diogo's algorithm...");
inicia_cronometro(&cronometro);
    decomposicoes_encontradas = encontra_decomp(lista_sementes,deg(expr_simplificada), expr_simplificada, lista_literais);
para_cronometro(&cronometro);
    diogo_parfrac = global_num_parfrac; */

    //testa quantos PB's de fato fizeram parte de alguma decomposicao
    lista_polinomios = remove_polinomios_nao_pertencentes(decomposicoes_encontradas,lista_polinomios);
/*
    // prints comparative results
    if (mulder_parfrac > diogo_parfrac)
     {
     printf("\n Diogo's algorithm performed less partial fraction expansions then Mulder's algorithm by %d%c",(100 -(diogo_parfrac*100)/(mulder_parfrac)),'%');
     }
     else if (diogo_parfrac > mulder_parfrac)
     {
     printf("\n Mulder's algorithm performed less partial fraction expansions then Diogo's algorithm by %d%c",(100 -(mulder_parfrac*100)/(diogo_parfrac)),'%');
     }
     else //they are equal
     {
     printf("\n Diogo's algorithm performed the same amount of partial fraction expansio as Mulder's");
     }
*/
    //contagem p—s-encontrar decomposicoes
    percorre_polinomios = lista_polinomios;
    contador = 0;
    while (percorre_polinomios != NULL)
    {
        contador++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_polinomios = percorre_polinomios->proximo_polinomio;
    }

    printf("\n Number of Base-Polynomials effectively used in any decompositions (extremely reduced set T1) is: %d",contador);
    printf("\n Press any key to continue...\n");
    global_num_parfrac = 0;
//    getchar();

    //elimina o lixo da memoria
    destroi_vetor_decomp(decomposicoes_encontradas);
    destroi_lista_sementes(lista_sementes);
    decomposicoes_encontradas = NULL;
    lista_sementes = NULL;

/**************
 Algorithm re-run after obtaining extremely reduced BP set
 ****************/

printf("\n Re-running the algorithms with the extremely reduced set T1");
    //gerar as sementes
//    lista_sementes = gera_vetor_semente(lista_polinomios, expr_simplificada);

    //conta quantos vetores-semente existem
 /*   percorre_sementes = lista_sementes;
    contador2 = 0;
    while (percorre_sementes != NULL)
    {
        contador2++;
        // imprime_lista_expr_expandida(percorre_polinomios->polinomio, lista_literais);
        // printf("\n");
        percorre_sementes = percorre_sementes->conjunto_prox;
    }
    printf("\n Number of Base-Polynomial pairs (refined T2 set) is: %d",contador2);
    printf("\n Press any key to continue...\n");
//    getchar();

    //calculo do pior caso para o meu algoritmo
    lim_superior = 0;
    for(lim_inferior=contador2;lim_inferior>=1;lim_inferior--)
    {
        lim_superior += lim_inferior;
    }

    printf("\n Worst case number of attempts using this alfgorithm is: %d",lim_superior);

    //calculo de pior caso para o mulden

    lim_superior = contador2*(contador - 1)*2;
    printf("\n Worst case number of attempts using Mulder's agorithm is: %d",lim_superior);
//    printf("\n Press any key to continue...\n");
//    getchar();
 */

    //zera novamente o mnumero de fracoes parciais
    global_num_parfrac = 0;

    //realiza o algoritmo do mulder
    printf("\n Performing Mulder's algorithm...");
inicia_cronometro(&cronometro);
    decomposicoes_encontradas = encontra_decomp_mulder(lista_polinomios, deg(expr_simplificada), expr_simplificada, lista_literais);
para_cronometro(&cronometro);
    mulder_parfrac = global_num_parfrac;

    //elimina o lixo da memoria
    destroi_vetor_decomp(decomposicoes_encontradas);
    decomposicoes_encontradas = NULL;
/*
    //refaz o nosso algoritmo para verificar se consegue s— com os efetivamente usados
    global_num_parfrac = 0;
    //realiza a decomposi‹o
    printf("\n Performing Diogo's algorithm...");
inicia_cronometro(&cronometro);
    decomposicoes_encontradas = encontra_decomp(lista_sementes,deg(expr_simplificada), expr_simplificada, lista_literais);
para_cronometro(&cronometro);

    diogo_parfrac = global_num_parfrac;
    //elimina o lixo da memoria
    destroi_vetor_decomp(decomposicoes_encontradas);

    //prints comparative result
    if (mulder_parfrac > diogo_parfrac)
    {
        printf("\n Diogo's algorithm performed less partial fraction expansions then Mulder's algorithm by %d%c",(100 -(diogo_parfrac*100)/(mulder_parfrac)),'%');
    }
    else if (diogo_parfrac > mulder_parfrac)
    {
        printf("\n Mulder's algorithm performed less partial fraction expansions then Diogo's algorithm by %d%c",(100 -(mulder_parfrac*100)/(diogo_parfrac)),'%');
    }
    else //they are equal
    {
        printf("\n Diogo's algorithm performed the same amount of partial fraction expansionsas Mulder's");
    }

    //imprime apenas os polinomios utilizados em alguma decomposicao
    printf("\n The base polynomials used in any of the decompositions found are:\n");
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
///imprime as sementes
 percorre_sementes = lista_sementes;
 contador = 0;
 while (percorre_sementes!= NULL)
 {
 contador++;
 printf("\n");
 imprime_lista_expr_expandida(percorre_sementes->P1.P, lista_literais);
 printf(" e ");
 imprime_lista_expr_expandida(percorre_sementes->P2.P, lista_literais);
 
 printf("\n");
 
 percorre_sementes = percorre_sementes->conjunto_prox;
 }
 printf("\n numero de pares encontrados: %d",contador);
 getchar();
 */
//fazer os testes do novo algoritmo aqui
//lista_sementes = ordena_vetor_semente(lista_sementes);
//correlacoes = constroi_correlacao_direta(lista_sementes);
//imprime_correlacoes_diretas(correlacoes);
//fim testes

    //aqui eu j‡ posso destruir a lista de polinomios
    //destroi o vetor de polinomios
    while (lista_polinomios != NULL)
    {
        percorre_polinomios = lista_polinomios->proximo_polinomio;
        destroi_lista_expr_expandida(lista_polinomios->polinomio->P);
        free(lista_polinomios->polinomio);
        free(lista_polinomios);
        lista_polinomios = percorre_polinomios;
    }
    destroi_lista_expr_expandida(polinomio_base);

    //e o vetor sementes
    destroi_lista_sementes(lista_sementes);


#endif
    destroi_tabela_literais(lista_literais);
    destroi_lista(lista_token);	
    destroi_arvore_expr(arvore);
    destroi_lista_expr(expressao_RPN);
    destroi_lista_expr_expandida(expr_simplificada);
printf("\nPress any key to continue...");
//getchar();
    return 0;
}
           
