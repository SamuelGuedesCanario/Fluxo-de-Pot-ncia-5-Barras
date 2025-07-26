function fluxo_carga_tensao()
clear; clc; close all;

%% Dados do sistema
Sbase = 100; % MVA

% Valores esperados de tensão (que você deseja)
V_esperado = [1.06; 1.0474; 1.0242; 1.0236; 1.0179];

% [Barra Tipo Pg Qg Pd Qd Vspec]
dados = [
    1 1 0    0   0   0   1.06;
    2 2 40   30  20  10   1.00;
    3 3 0    0   45  15   1.00;
    4 3 0    0   40  5    1.00;
    5 3 0    0   60  10   1.00
];

%% Inicialização
nb = size(dados,1);
tipo = dados(:,2);
PQ = find(tipo == 3);
PV = find(tipo == 2);
slack = find(tipo == 1);

% Conversão para pu
P = (dados(:,3) - dados(:,5)) / Sbase;
Q = (dados(:,4) - dados(:,6)) / Sbase;
Vspec = dados(:,7);

% Condições iniciais
V = ones(nb,1); 
V(slack) = Vspec(slack);
V(PV) = Vspec(PV);
theta = zeros(nb,1);

% Parâmetros do método
tol = 1e-6;
max_iter = 10;
tolerancia_tensao = 0.05; % Tolerância apenas para tensão

%% Matriz de admitâncias
Y = zeros(nb);
LT = [
    1 2 0.02 0.06 0.03;
    1 3 0.08 0.24 0.025;
    2 3 0.06 0.25 0.02;
    2 4 0.06 0.18 0.02;
    2 5 0.04 0.12 0.015;
    3 4 0.01 0.03 0.01;
    4 5 0.08 0.24 0.025
];

for k = 1:size(LT,1)
    i = LT(k,1); j = LT(k,2);
    R = LT(k,3); X = LT(k,4); Bsh = LT(k,5);
    Yk = 1/(R + 1i*X);
    Y(i,i) = Y(i,i) + Yk + 1i*Bsh/2;
    Y(j,j) = Y(j,j) + Yk + 1i*Bsh/2;
    Y(i,j) = Y(i,j) - Yk;
    Y(j,i) = Y(i,j);
end

G = real(Y); B = imag(Y);

%% Método de Newton-Raphson
for iter = 1:max_iter
    % Cálculo das potências
    Pcalc = zeros(nb,1);
    Qcalc = zeros(nb,1);
    for i = 1:nb
        for k = 1:nb
            Pcalc(i) = Pcalc(i) + V(i)*V(k)*(G(i,k)*cos(theta(i)-theta(k)) + B(i,k)*sin(theta(i)-theta(k)));
            Qcalc(i) = Qcalc(i) + V(i)*V(k)*(G(i,k)*sin(theta(i)-theta(k)) - B(i,k)*cos(theta(i)-theta(k)));
        end
    end

    % Mismatches
    dP = P - Pcalc;
    dQ = Q - Qcalc;
    dP(slack) = [];
    dQ([slack; PV]) = [];
    mismatch = [dP; dQ];

    % Verificação de convergência
    if max(abs(mismatch)) < tol
        fprintf('\nConvergência alcançada na iteração %d\n', iter);
        break;
    end

    % Construção da Jacobiana
    [J11, J12, J21, J22] = build_jacobian(nb, PQ, PV, V, theta, G, B, Pcalc, Qcalc);
    J = [J11 J12; J21 J22];
    
    % Solução do sistema
    dx = J \ mismatch;
    
    % Atualização das variáveis
    theta(2:end) = theta(2:end) + dx(1:nb-1);
    V(PQ) = V(PQ) + dx(nb:end);
    V(PV) = Vspec(PV); % Mantém tensão especificada nas PV
end

%% Comparação com valores esperados (apenas tensão)
fprintf('\nComparação de Tensões (Tolerância = %.2f pu):\n', tolerancia_tensao);
fprintf('Barra | V esperado | V calculado | Diferença | Status\n');
for i = 1:nb
    dif_V = abs(V(i) - V_esperado(i));
    status = 'Dentro';
    if dif_V > tolerancia_tensao
        status = 'Fora';
    end
    fprintf('%4d  | %9.4f  | %10.4f  | %8.4f  | %s\n',...
            i, V_esperado(i), V(i), dif_V, status);
end

%% Verificação final da tensão
erro_V = max(abs(V - V_esperado));

fprintf('\nErro máximo em V: %.4f pu\n', erro_V);

if erro_V < tolerancia_tensao
    fprintf('Todos valores de tensão dentro da tolerância de %.2f pu\n', tolerancia_tensao);
else
    fprintf('Atenção: Algumas tensões fora da tolerância de %.2f pu\n', tolerancia_tensao);
end

%% Resultados completos (para referência)
fprintf('\nResultados completos:\n');
fprintf('Barra |  V (pu)  | Ângulo (graus)\n');
for i = 1:nb
    fprintf('%4d  |  %7.4f |   %9.4f\n', i, V(i), rad2deg(theta(i)));
end

%% Função para construção da Jacobiana
function [J11, J12, J21, J22] = build_jacobian(nb, PQ, PV, V, theta, G, B, Pcalc, Qcalc)
    % Submatrizes da Jacobiana
    J11 = zeros(nb-1); % dP/dTheta
    J12 = zeros(nb-1, length(PQ)); % dP/dV
    J21 = zeros(length(PQ), nb-1); % dQ/dTheta
    J22 = zeros(length(PQ)); % dQ/dV
    
    % Barras não slack (para J11 e J12)
    buses_no_slack = setdiff(1:nb, 1);
    
    % Preenche J11 (dP/dTheta)
    for i = 1:(nb-1)
        m = buses_no_slack(i);
        for k = 1:(nb-1)
            n = buses_no_slack(k);
            if m == n
                J11(i,k) = -Qcalc(m) - B(m,m)*V(m)^2;
            else
                J11(i,k) = V(m)*V(n)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
            end
        end
    end
    
    % Preenche J12 (dP/dV)
    for i = 1:(nb-1)
        m = buses_no_slack(i);
        for k = 1:length(PQ)
            n = PQ(k);
            if m == n
                J12(i,k) = Pcalc(m)/V(m) + G(m,m)*V(m);
            else
                J12(i,k) = V(m)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
            end
        end
    end
    
    % Preenche J21 (dQ/dTheta)
    for i = 1:length(PQ)
        m = PQ(i);
        for k = 1:(nb-1)
            n = buses_no_slack(k);
            if m == n
                J21(i,k) = Pcalc(m) - G(m,m)*V(m)^2;
            else
                J21(i,k) = -V(m)*V(n)*(G(m,n)*cos(theta(m)-theta(n)) + B(m,n)*sin(theta(m)-theta(n)));
            end
        end
    end
    
    % Preenche J22 (dQ/dV)
    for i = 1:length(PQ)
        m = PQ(i);
        for k = 1:length(PQ)
            n = PQ(k);
            if m == n
                J22(i,k) = Qcalc(m)/V(m) - B(m,m)*V(m);
            else
                J22(i,k) = V(m)*(G(m,n)*sin(theta(m)-theta(n)) - B(m,n)*cos(theta(m)-theta(n)));
            end
        end
    end
end
end