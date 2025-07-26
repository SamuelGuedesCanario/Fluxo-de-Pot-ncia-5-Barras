1. Objetivo
Este código implementa o método Newton-Raphson para resolver o problema de fluxo de carga (Power Flow) em sistemas elétricos de potência. Ele calcula:
✔ Tensões nas barras (magnitude e ângulo)
✔ Fluxos de potência nas linhas
✔ Perdas no sistema

2. Fluxo do Código
O algoritmo segue as seguintes etapas:

Passo 1: Entrada de Dados
Barras: Slack, PV e PQ (com P, Q, V especificados).

Linhas: Impedâncias (R, X) e admitâncias shunt (Bsh).

Passo 2: Montagem da Matriz Ybus
Calcula a matriz de admitância nodal (Y = G + jB).

Considera admitâncias shunt (metade em cada extremidade da linha).

Passo 3: Inicialização das Variáveis
Barras PQ: V = 1.0 pu, θ = 0 rad.

Barras PV: V = valor especificado, θ = 0 rad.

Barra Slack: V e θ fixos.

Passo 4: Loop Newton-Raphson
Calcula Pcalc e Qcalc (potências injetadas).

Determina os mismatches (ΔP, ΔQ).

Verifica convergência (tolerância 1e-6).

Monta a matriz Jacobiana (J11, J12, J21, J22).

Resolve o sistema linear J * [Δθ; ΔV] = [ΔP; ΔQ].

Atualiza V e θ.

Passo 5: Pós-Processamento
Fluxos de potência nas linhas.

Perdas de transmissão.

Geração na barra Slack.

Passo 6: Saída dos Resultados
Tabela de tensões (V, θ).

Fluxos de potência ativa/reativa.

Perdas totais.

README.md (Para o Repositório GitHub)
markdown
# **Power Flow Solver - Newton-Raphson Method (MATLAB)**  

## **📌 Descrição**  
Implementação do **método Newton-Raphson** para solução do problema de **fluxo de carga** em sistemas elétricos de potência.  

## **📋 Funcionalidades**  
✔ Cálculo de **tensões nodais (V, θ)**  
✔ Cálculo de **fluxos de potência nas linhas**  
✔ Cálculo de **perdas de transmissão**  
✔ Suporte para barras **PQ, PV e Slack**  

## **⚙️ Como Usar**  
1. **Defina os dados do sistema** em `dados_barras` e `dados_linhas`.  
2. **Execute o script `fluxo_carga_NR.m`**.  
3. **Visualize os resultados**:  
   - Tensões nas barras (`V`, `θ`).  
   - Fluxos de potência (`P`, `Q`).  
   - Perdas totais.  

## **📊 Exemplo de Saída**  
Barra | V (pu) | Ângulo (graus)

1 | 1.0600 | 0.0000
2 | 1.0474 | -2.8063
3 | 1.0242 | -4.9970

text

## **Referências**  
- Kundur, P. (1994). *Power System Stability and Control*.  
- Saadat, H. (1999). *Power System Analysis*.  

## **Contato**  
Se tiver dúvidas, abra uma **issue** ou envie um e-mail para: [samuelgcanario1618@gmail.com].  
Observações Finais
Precisão: O método Newton-Raphson tem convergência quadrática.

Limitações: Pode divergir para sistemas mal condicionados.

Extensões Futuras:

Implementar controles automáticos de tensão.

Adicionar análise de contingências.
