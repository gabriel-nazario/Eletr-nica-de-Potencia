clear; close all; clc;
format long


%% COMPENSADOR TRIFASICO 4L

% Parametros do Circuito
Rg = 0.1;
Lg = 5e-3;
Xg = Lg*(2*pi*60);

% Relação pot da carga
P = 1000; %(W)
fp = 0.8; %(at)
V = 110; %(Vrms)
Q = (P/fp)*sin(acos(fp));
S = complex(P, Q);

Rl = real(V^2/S);
Ll = -(imag(V^2/S))/(2*pi*60);

% Parametros de Simulação Numerica
h = 1E-7;  % -> Passo de calculo
t = 0;     % -> Tempo inicial de Simulação
tf = 0.2;

% Parametros de Gravação para plot
tsave0 = 0;
tsave = tsave0;
npt = 20000;
hsave = (tf-tsave0)/npt;

if hsave < h
    hsave = h;
    npt = (tf-tsave0)/hsave;
end 

%Variavel de interação do salvamento
n = 0;

% Condições iniciais do Sistema
Vc = 600;     % Tensão maxima do barramento CC

Eg = 110*sqrt(2);

% Vg de regime
Igx    = roots([Rg -Eg P]);
Igx    = min(Igx);
Igx * sqrt(2);
faseg = atan((-Igx*Xg)/(Eg - Igx*Rg));
Vgrefx = -Igx*Xg/sin(faseg);

vg_ref = Vgrefx
vl_ref = 110*sqrt(2);

thetal_ref = 0;

f_ref = 60;
w_ref = 2*pi*f_ref;

ig1 = 0;
il1 = 0;

ig2 = 0;
il2 = 0;

ig3 = 0;
il3 = 0;

% Tensões iniciais
vg10_int = 0;
vl10_int = 0;

vg30_int = 0;
vl30_int = 0;

vg0_int = 0;
vl0_int = 0;

vg2_int = 0;
vl2_int = 0;

% Formação do sinal triangular
ttriangle = 0;           % Tempo inicial da onda triangular
ftriangle = 10E3;        % Frequência da onda triangular de 10 kHz
htriangle = 1/ftriangle; % Período da onda triangular
vtriangle = Vc/2;         % Tensão máxima da onda triangular
dtriangle = Vc/htriangle; % Derivada da onda triangular (dV/dT)
sign = -1;               % Sinal inicial da derivada (Comportamento decrescente)


% Inicio do looping de iterações
while(t < tf)
    t = t + h;

    if t >= ttriangle
        ttriangle = ttriangle + htriangle;

        if vtriangle <= 0
            vtriangle = -Vc/2;
            sign = 1;
        else
            vtriangle = Vc/2;
            sign = -1;
        end
        
        %Tensões de referencia na entrada
        Vg1_ref = vg_ref*cos(w_ref*t + faseg);
        Vg2_ref = vg_ref*cos(w_ref*t + 2*pi/3 + faseg);
        Vg3_ref = vg_ref*cos(w_ref*t - 2*pi/3 + faseg);

        %Tensões de referencia na saída
        Vl1_ref = vl_ref*cos(w_ref*t + thetal_ref);
        Vl2_ref = vl_ref*cos(w_ref*t + thetal_ref + 2*pi/3);
        Vl3_ref = vl_ref*cos(w_ref*t + thetal_ref - 2*pi/3);

        %% PWM      
        Vg10_ref = Vg1_ref - Vg2_ref; % Vg12
        Vg30_ref = Vg3_ref - Vg2_ref; % Vg32
        
        Vl10_ref = Vl1_ref - Vl2_ref; % Vl12
        Vl30_ref = Vl3_ref - Vl2_ref; % Vl32
 
        % Media dos valores a partir do valor e entrada/periodo
        vg10_med = vg10_int*ftriangle; 
        vl10_med = vl10_int*ftriangle;
        
        vg30_med = vg30_int*ftriangle; 
        vl30_med = vl30_int*ftriangle;

        vg2_med = vg2_int*ftriangle; 
        vl2_med = vl2_int*ftriangle;        
        
        vg0_med = vg0_int*ftriangle;
        vl0_med = vl0_int*ftriangle;
        
        % - Zerando a integral

        vg10_int = 0;
        vl10_int = 0;

        vg30_int = 0;
        vl30_int = 0;

        vg2_int = 0;
        vl2_int = 0;

        vg0_int = 0;
        vl0_int = 0;
    end

    % Passo da iteração
    vtriangle = vtriangle + sign*dtriangle*h;
    
    % Comparação do sinal referencia com o triangular para p PWM
    %-Entrada 10-
    if Vg10_ref >= vtriangle
        qg1 = 1;
    else
        qg1 = 0;
    end
    %-Saida 10-
    if Vl10_ref >= vtriangle
        ql1 = 1;
    else
        ql1 = 0;
    end
    %-Entrada 30-
    if Vg30_ref >= vtriangle
        qg3 = 1;
    else
        qg3 = 0;
    end
    %-Saida 30-
    if Vl30_ref >= vtriangle
        ql3 = 1;
    else
        ql3 = 0;
    end

    % Relação disposta no material 3L (2qj-1)
    vg10 = (2*qg1 - 1)*(Vc/2);
    vl10 = (2*ql1 - 1)*(Vc/2);

    vg30 = (2*qg3 - 1)*(Vc/2);
    vl30 = (2*ql3 - 1)*(Vc/2);
    
    % Apenas para mostrar o sinal senoidal da entrada
    eg1 = Eg*cos(w_ref*t);
    eg2 = Eg*cos(w_ref*t + 2*pi/3);
    eg3 = Eg*cos(w_ref*t - 2*pi/3);
    
    % Relação disposta no material 4L
    % - Entrada -
    vg0 = (vg10 + vg30)/3;
    vg1 = 2*vg10/3 - vg30/3;
    vg3 = 2*vg30/3 - vg10/3;
    vg2 = -(vg10 + vg30)/3;
    
    % - Saida -
    vl0 = (vl10 + vl30)/3;
    vl1 = 2*vl10/3 - vl30/3;
    vl3 = 2*vl30/3 - vl10/3;
    vl2 = -(vl10 + vl30)/3;
    
    % Relação disposta no material
    vg10_int = vg10_int + vg10*h;
    vl10_int = vl10_int + vl10*h;

    vg30_int = vg30_int + vg30*h;
    vl30_int = vl30_int + vl30*h;

    vg2_int = vg2_int + vg2*h;
    vl2_int = vl2_int + vl2*h;
    
    vg0_int = vg0_int + vg0*h;
    vl0_int = vl0_int + vl0*h;

    % Integração Numerica a partir das EDOs
    % - Entrada -
    ig1 = (1-(h*Rg/Lg))*ig1 + (h/Lg)*(eg1-vg1);
    ig3 = (1-(h*Rg/Lg))*ig3 + (h/Lg)*(eg3-vg3);
    ig2 = - ig1 - ig3;

    % - Saida -
    il1 = (1-(h*Rl/Ll))*il1 + (h/Ll)*(vl1);
    il3 = (1-(h*Rl/Ll))*il3 + (h/Ll)*(vl3);
    il2 = - il1 - il3;

    % Salvando as variaveis para plot 
    if tsave <= t
        tsave = tsave + hsave;
        n = n + 1;
        Ts(n) = t;                        %#ok<SAGROW>
        
        % - Fonte e barramento -
        vcs(n) = Vc;                      %#ok<SAGROW>
        eg1s(n) = eg1;                      %#ok<SAGROW>
        eg2s(n) = eg2;                      %#ok<SAGROW>
        eg3s(n) = eg3;                      %#ok<SAGROW>
        
        % - Valores de referencia (puros) - 
        Vg1_refs(n) = Vg1_ref;              %#ok<SAGROW>
        Vl_refs(n) = Vl1_ref;              %#ok<SAGROW>

        Vg2_refs(n) = - Vg1_ref - Vg3_ref;              %#ok<SAGROW>
        Vl2_refs(n) = - Vl1_ref - Vg3_ref;              %#ok<SAGROW>

        Vg3_refs(n) = Vg3_ref;              %#ok<SAGROW>
        Vl3_refs(n) = Vl3_ref;              %#ok<SAGROW>

        vg0s(n) = vg0;            %#ok<SAGROW> 
        vl0s(n) = vl0;            %#ok<SAGROW> 
        
        % - Comparados ao triangular - 
        Vg10_refs(n) = Vg10_ref;          %#ok<SAGROW>
        Vl10_refs(n) = Vl10_ref;          %#ok<SAGROW>

        Vg30_refs(n) = Vg30_ref;          %#ok<SAGROW>
        Vl30_refs(n) = Vl30_ref;          %#ok<SAGROW>

        % - Relação de chaveamento pos PWM - 
        vg10s(n) = vg10;          %#ok<SAGROW>
        vl10s(n) = vl10;          %#ok<SAGROW>

        vg30s(n) = vg30;          %#ok<SAGROW>
        vl30s(n) = vl30;          %#ok<SAGROW>
        
        % - Media das tensoes -

        vg0_meds(n) = vg0_med;    %#ok<SAGROW> 
        vl0_meds(n) = vl0_med;    %#ok<SAGROW> 

        vg10_meds(n) = vg10_med;          %#ok<SAGROW>
        vl10_meds(n) = vl10_med;          %#ok<SAGROW>

        vg2_meds(n) = vg2_med;          %#ok<SAGROW>
        vl2_meds(n) = vl2_med;          %#ok<SAGROW>

        vg30_meds(n) = vg30_med;          %#ok<SAGROW>
        vl30_meds(n) = vl30_med;          %#ok<SAGROW>
        
        % - Resultados das entradas e saidas (Tensão) -

        vg0s(n) = vg0;            %#ok<SAGROW>
        vg1s(n) = vg1;            %#ok<SAGROW>
        vl1s(n) = vl1;            %#ok<SAGROW>

        vg2s(n) = vg2;            %#ok<SAGROW>
        vl2s(n) = vl2;            %#ok<SAGROW>

        vg3s(n) = vg3;            %#ok<SAGROW>
        vl3s(n) = vl3;            %#ok<SAGROW>

        % - Resultados das entradas e saidas (Corrente) -
        ig1s(n) = ig1;            %#ok<SAGROW>
        il1s(n) = il1;            %#ok<SAGROW>

        ig2s(n) = ig2;            %#ok<SAGROW>
        il2s(n) = il2;            %#ok<SAGROW>

        ig3s(n) = ig3;            %#ok<SAGROW>
        il3s(n) = il3;            %#ok<SAGROW>
        
        % - Estado de ativação das chaves -
        qg1s(n) = qg1;            %#ok<SAGROW>
        ql1s(n) = ql1;            %#ok<SAGROW>

        qg3s(n) = qg3;            %#ok<SAGROW>
        ql3s(n) = ql3;            %#ok<SAGROW>

        vtriangles(n)= vtriangle; %#ok<SAGROW>
    
        % Como se trata de um tamanho móvel desses arrays, é necessario
        % usar esse comentario: %#ok<SAGROW> para evitar warnings
    end
end

% Plots

%---- Entradas do inversor ----
figure('Name','Tensão da Rede')
plot(Ts,eg1s,Ts,eg2s,Ts,eg3s,Ts,vcs,'r','LineWidth',1.5),zoom
title('Tensão de rede do barramento', 'FontSize', 18)
legend('eg1','eg2','eg3','E_barramento')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensões no polo G1 ----
figure('Name','Tensão do conversor lado G')
subplot(3,1,1)
plot(Ts,Vg1_refs,Ts,vg1s,Ts,vg10_meds-vg0_meds,'r','LineWidth',1.5),zoom
title('vg1')
legend('vg1_{ref}','vg1_{chaveado}','vg1_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensões no polo G2 ----
subplot(3,1,2)
plot(Ts,Vg2_refs,Ts,vg2s,Ts,vg2_meds,'r','LineWidth',1.5),zoom
title('vg2')
legend('vg2_{ref}','vg2_{chaveado}','vg2_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensões no polo G3 ----
subplot(3,1,3)
plot(Ts,Vg3_refs,Ts,vg3s,Ts,vg30_meds-vg0_meds,'r','LineWidth',1.5),zoom
title('vg3')
legend('vg3_{ref}','vg3_{chaveado}','vg3_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor


%---- Tensão no polo L1 
figure('Name','Tensão do conversor lado L')
subplot(3,1,1)
plot(Ts,vl1s,Ts,vl10_meds-vl0_meds,'r','LineWidth',1.5),zoom
title('Tensão do Conversor lado L1')
legend('vl1_{chaveado}','vl1_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensão no polo L2 
subplot(3,1,2)
plot(Ts,vl2s,Ts,vl2_meds,'r','LineWidth',1.5),zoom
title('Tensão do Conversor lado L2')
legend('vl2_{chaveado}','vl2_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensão no polo L3 
subplot(3,1,3)
plot(Ts,vl3s,Ts,vl30_meds-vl0_meds,'r','LineWidth',1.5),zoom
title('Tensão do Conversor lado L3')
legend('vl3_{chaveado}','vl3_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Vgx0
figure('Name','VgX0')
subplot(2,1,1)
plot(Ts,Vg10_refs,Ts, vg10s, Ts, vg10_meds,'r-','LineWidth',2),zoom
title('vg10')
legend('vg10_{ref}','vg10_{chaveado}','vg10_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

subplot(2,1,2)
plot(Ts,Vg30_refs,Ts, vg30s, Ts, vg30_meds,'r-','LineWidth',2),zoom
title('vg30')
legend('vg30_{ref}','vg30_{chaveado}','vg30_{med}')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Correntes
figure('Name','CorrentesX')
subplot(2,1,1)
plot(Ts,ig1s, Ts, ig2s, Ts, ig3s,'r-','LineWidth',2),zoom
title('ig')
legend('ig1','ig2','ig3')
xlabel("Tempo (s)")
ylabel("Corrente (A)")
grid minor

subplot(2,1,2)
plot(Ts,il1s, Ts, il2s, Ts, il3s,'r-','LineWidth',2),zoom
title('il')
legend('il1','il2','il3')
xlabel("Tempo (s)")
ylabel("Corrente (A)")
grid minor

%---- Tensão do sinal triangular em g1
figure('Name','Sinal triangular em g1')
subplot(2,1,1)
plot(Ts,vtriangles, Ts, Vg10_refs, Ts, qg1s.*100,'r-','LineWidth',2),zoom
title('Sinal triangular em g1')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensão do sinal triangular em g3
subplot(2,1,2)
plot(Ts,vtriangles, Ts, Vg30_refs, Ts, qg3s.*100,'r-','LineWidth',2),zoom
title('Sinal triangular em g3')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensão do sinal triangular em L1
figure('Name','Sinal triangular em L1')
subplot(2,1,1)
plot(Ts,vtriangles, Ts, Vl10_refs, Ts, ql1s.*100,'r-','LineWidth',2),zoom
title('Sinal triangular em l1')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

%---- Tensão do sinal triangular em L3
subplot(2,1,2)
plot(Ts,vtriangles, Ts, Vl30_refs, Ts, ql3s.*100,'r-','LineWidth',2),zoom
title('Sinal triangular em l3')
xlabel("Tempo (s)")
ylabel("Tensão (V)")
grid minor

