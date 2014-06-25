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
    //Variáveis
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
#endif
    

    
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
    printf("\n\t Decomp_nparametrica: Aplicativo para realizacao de decomposicao nao-parametrica\n\r\t em equacoes a serem implementadas em circuito translinear.\n\r \
        Insira abaixo a equacao do polinomio desejado (100 caracteres no maximo). Os escalares deverao ser apenas \n\r\t numeros inteiros. Exemplo:\n\r\t \
        x^3 + x^2*y - x*z^2 + y*z^2\n\r  (tecle enter para utiliza-lo)");
   
    //leitura da string de entrada
    
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);  
    }
    
    //equacao padrao, caso o usuario aperte enter direto. pode ser igual a nova linha ou retorno de carro para tratar como o sistema operacional reconhece o ENTER em varias plataformas diferentes
    if (*equacao_entrada == '\n' || *equacao_entrada == '\r')
        sprintf(equacao_entrada, "x^3 + x^2*y - x*z^2 + y*z^2");
    
#endif
    
#if defined (DEBUG_VECTOR_GEN)
   
    sprintf(equacao_entrada, "3*a+b^2");
    lim_inf = -1;
    lim_sup = +1;
    //sprintf(equacao_entrada, "X^3 + X^2*Y - X*Z^2 + Y*Z^2");
    

#endif

    //teste da divisao polinomial
#if defined (DEBUG_POLYDIV)
    //polinomio numerador
    //sprintf(equacao_entrada, "(-2*Iin-Iout+I0)*(Iin+I0)");
    sprintf(equacao_entrada, "3*a+b^2");
    
#endif
    
#if defined (DEBUG_PARFRAC)
    //teste da expansao em fracoes parciais
    //as entradas serão igual no teste da divisão polinomial, com a diferença de que o divisor será quebrado em P1 e P2
    //A divisão será feita assim como no teste da divisão polinomial, e o resto será processado junto com P1 e P2 para 
    //expansão em fracoes parciais
    
    //polinomio numerador
    //sprintf(equacao_entrada, "X^3 + X^2*Y - X*Z^2 + Y*Z^2");
    sprintf(equacao_entrada, "X");
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
        
           
    //Criação das pilhas de avaliação de expressão
           
    expressao_RPN = constroi_lista_expr(lista_token);    
    arvore = constroi_arvore_expr(expressao_RPN);
#if defined DEBUG_EXPR
    imprime_arvore_expr(arvore);
    imprime_lista_expr(expressao_RPN);
#endif
    
    //Expansão da arvore de expressões
    expr_expandida = constroi_lista_expressoes_exp(arvore);
#if defined DEBUG_EXPAND
    printf("\n equacao expandida:");
    imprime_lista_expr_expandida(expr_expandida,lista_literais);
#endif   
    
    //Simplificação da árvore de expressões
    expr_simplificada = simplifica_expr_expandida(expr_expandida);
    destroi_lista_expr_expandida(expr_expandida);
#if defined DEBUG_SIMPLIFY
    printf("\n equacao simplificada:");
    imprime_lista_expr_expandida(expr_simplificada,lista_literais);
#endif 
    
    //ordenação lexdeg
   expr_simplificada = lexdegbubblesort(expr_simplificada);
#if defined DEBUG_SIMPLIFY
    printf("\n equacao reordenada lex:");
    imprime_lista_expr_expandida(expr_simplificada,lista_literais);
#endif 

#if defined (DEBUG_POLYDIV) 
    //polinomio denominador
    //sprintf(divisor_entrada, "2*Iin+I0");
    sprintf(divisor_entrada, " a^2+a*b");
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
    sprintf(divisor_entrada, "-X+Z");
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
    sprintf(divisor_entrada, "Z");
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
    
    //multiplica P1 e P2 e procede normalmente a divisão
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
#if !defined (DEBUG_PARFRAC) //preciso dos ponteiros de resto e quocientes intactos para realizar a expansão
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
    
    //imprime as combinações
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
    
    //gerar as sementes ->aqui já está dando pau no windows
    lista_sementes = gera_vetor_semente(lista_polinomios, expr_simplificada);

    //aqui eu já posso destruir a lista de polinomios
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
        
        //faz a multiplicação cruzada
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

    //realiza a decomposição 
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

    //imprimir a lista de literais e construir o polinomio base
    
    polinomio_base = gera_polinomio_base(lista_literais);
    
    //leitura dos limites superiores e inferiores dos coeficientes
    printf("\n Insira o limite inferior de coeficientes a serem utilizados (padrao -1)");
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);  
    }
    
    //tratamento do número inserido
    //se apertar enter direto, é -1
    if (*equacao_entrada == '\n' || *equacao_entrada == '\r')
        lim_inferior = -1;
    else
    {
        lim_inferior = atoi(equacao_entrada);
    }
    
    //limite superior
    printf("\n Insira o limite superior de coeficientes a serem utilizados (padrao +1)");
    if ( fgets (equacao_entrada, 100 , stdin) == NULL )
    {
        erro(ERRO_002);
        return(0);  
    }
    
    //tratamento do número inserido
    //se apertar enter direto, é -1
    if (*equacao_entrada == '\n' || *equacao_entrada == '\r')
        lim_superior = 1;
    else
    {
        lim_superior = atoi(equacao_entrada);
    }
    
    //imprimir os parametros colhidos
    printf("\n equacao:");
    imprime_lista_expr_expandida(expr_simplificada, lista_literais);
    printf("\n limite inferior: %d", lim_inferior);
    printf("\n limite superior: %d", lim_superior);
    
    printf("\n pressione qualquer tecla para continuar...\n");
    getchar();
    
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
    printf("\n numero de polinomios base total eh: %d",contador);
    printf("\n pressione qualquer tecla para continuar...\n");
    getchar();


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
    printf("\n numero de polinomios base apos eliminar negativos sao: %d",contador);
    printf("\n pressione qualquer tecla para continuar...\n");
    getchar();

    //elimina os polinomios redundantes
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
    printf("\n numero de polinomios base apos eliminar redundantes sao: %d",contador);
    printf("\n pressione qualquer tecla para continuar...\n");
    getchar();

//inserir aqui a remocao de polinomios que nao podem fazer parte de uma decomp
    //imprime as combinações
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
    printf("\n numero de polinomios base gerados: %d",contador);
    printf("\n pressione qualquer tecla para continuar...\n");
    getchar();
    
    //gerar as sementes 
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
    printf("\n numero de pares de polinomios base gerados: %d",contador2);
    printf("\n pressione qualquer tecla para continuar...\n");
    getchar();

//calculo do pior caso para o meu algoritmo
lim_superior = 0;
for(lim_inferior=contador2;lim_inferior>=1;lim_inferior--)
{
    lim_superior += lim_inferior;
}

printf("\n numero de pior caso de combinacoes meu algoritmo e: %d",lim_superior);

//calculo de pior caso para o mulden

lim_superior = contador2*(contador - 1)*2;
printf("\n numero de pior caso de combinacoes para mulder e: %d",lim_superior);
printf("\n pressione qualquer tecla para continuar...\n");
getchar();
    //calcula o numero total de testes a serem feitos
    //encontra_decomp_dummie(lista_sementes, deg(expr_simplificada));

    //encontra o total de testes pelo metodo mulder
    //encontra_decomp_dummie_mulder(lista_sementes,lista_polinomios, deg(expr_simplificada));


    //aqui eu já posso destruir a lista de polinomios
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

    //realiza a decomposição
    decomposicoes_encontradas = encontra_decomp(lista_sementes,deg(expr_simplificada), expr_simplificada, lista_literais);

    //elimina o lixo da memoria
    destroi_vetor_decomp(decomposicoes_encontradas);
    destroi_lista_sementes(lista_sementes);

#endif
    destroi_tabela_literais(lista_literais);
    destroi_lista(lista_token);	
    destroi_arvore_expr(arvore);
    destroi_lista_expr(expressao_RPN);
    destroi_lista_expr_expandida(expr_simplificada);
printf("\npressione qualquer tecla para continuar...");
getchar();
    return 0;
}
           
