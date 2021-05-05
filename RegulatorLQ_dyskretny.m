%Zadanie nr1: Projekt regulatora optymalnego LQ_D od stanu
%Zakres zadania:
%Analiza czesci ciaglej - identyfiakcja i opis dzialania instrukcji
%Projekt czesci dyskretnej - dyskretny regulator LQ_D statyczny i
%dynamiczny
%Charakterystyki rozwiazac ciaglych i dyskretnych dla R=[0.1 1.0 10]
%Opcjonalnie: wyznaczenie czestotliwosci patologicznej
%Analiza - (dla czesci ciaglej i dyskretnej) - wnioski ilosciowe i
%jakosciowe dla UR_A i UR_D; rekomendacje
%Uwagi: wykresy z legend; krotkie jednoznaczne komentarze

%Autorzy: Badelek Piotr & Fusiara Jakub
%Grupa: WMT18AP1S1
%Data: 19.04.2021
%                         1. Modele numeryczne 
clear;
close all;
clear all;

s=tf('s');                                                                 % Korzystamy z funckji tf (Transfer function model), czyli tworzymy przestrzen s
                                                                           
disp('-------------Obiekt nominalny------------------');
G=1/s^2 ;                                                                  % Stworzenie transmitancji czlonu calkujacego 2 rzedu 
[A, B, C, D]=ssdata(G);                                                    % Utworzenie macierzy modelu stanu na podstawie transmitancji G                                                 
Gss=ss(G);                                                                 % Stworzenie modelu w przestrzeni stanu z transmitancji G
n=size(A,1);                                                               % Przypisanie do zmiennej n, ilosci wierszy w macierzy A
Q=eye(n);                                                                  % Stworzenie macierzy jednostkowej, o rozmiarze rownej wartosci zmiennej n 

i = 1;                                                                     % Zmienna pomocnicza do pozycjonowania wykresow

R=0.1;                                                                     % Macierz wag, kara na sygnal sterowania
K=lqr(A, B, Q, R);                                                         % Rozwiazanie rownania Ricattiego
Acl=A-B*K;                                                                 % Obliczenie macierzy dla calego zamknietego obiektu regulacji, z uwzglednieniem sterowania regulatorem statycznym LQR
Gclz_q=ss(Acl,B,C,0);                                                      % Utworzenie modelu w przestrzeni stanu sterowania z regulatorem statycznym bez obserwatora
F_q=1/dcgain(Gclz_q);                                                      % Obliczenie wzmocnienia prekompensatora 
Gclu_q=ss(Acl,B,-K,0);                                                     % Utworzenie modelu w przestrzeni stanu ukladu regulacji z regulatorem statycznym                                             
P=eig(Acl);                                                                % Wyodrebnienie wartosci wlasnych macierzy Acl
L=place(A' ,C',3*P);                                                       % Przypisanie do zmiennej L wzmocnienia estymatora
Css=reg(Gss,K,L');                                                         % Stworzenie modelu regulatora LQR, o wzmocnieniu K, oraz estymatora o wzmocnieniu L                                                   
Gclz_o=feedback(Gss,Css,1);                                                % Stworzenie zamknietego ukladu regulacji, z modelem LQR w petli ujemnego sprzezenia zwrotnego
F_o=1/dcgain(Gclz_o);                                                      % Obliczenie wzmocnienia prekompensatora dla ukladu regulacji z estymatorem
Gclu_o=feedback(Gss*Css,1,1);                                              % Stworzenie modelu zamknietego ukladu regulacji z regulatorem LQR i estymatorem
tfinal=15;                                                                 % Maksymalny czas cyklu

y = stepinfo(F_o*Gclz_o);                                                  % Zapisanie w zmiennej wartosci potrzebnej do obliczenia czasu probkowania
Ts = y.RiseTime/12;                                                        % Obliczenie czasu probkowania
Gss_d=c2d(Gss,Ts);                                                         % Zamiana modelu ciaglego na model dyskretny z czestotliwoscia probkowania Ts
gs = tf(Gss_d);                                                            % Tworzymy model transmitancyjny korzystajac z modelu dyskretnego
[A_d,B_d,C_d,D_d, Ts] = ssdata(Gss_d);                                     % Utworzenie macierzy modelu stanu na podstawie transmitancji G  
K_d = dlqr(A_d, B_d, Q, R);                                                % Stworzenie modelu regulatora LQR z uwzglednieniem macierzy wag Q i R 
Acl_d = A_d-B_d*K_d;                                                       % Stworzenie modelu zamknietego ukladu regulacji, z regulatorem LQR w dziedzinie dyskretnej
Gclz_qd = ss(Acl_d, B_d, C_d, 0, Ts);                                      % Utworzenie modelu w przestrzeni stanu dla sterowania z regulatorem statycznym (Bez obserwatora)
F_qd = 1/dcgain(Gclz_qd);                                                  % Wyznaczenie czlonu kompensacyjnego, potrzebnego do wzmocnienia sygnalu w stanie ustalonym (Bez obserwatora)
Gclu_qd = ss(Acl_d, B_d, -K_d, 0, Ts);                                     % Stworzenie modelu zamknietego ukladu regulacji, z regulatorem LQR
P_d= eig(Acl_d);                                                           % Wyodrebnienie wartosci wlasnych macierzy Acl_d
L_d = place(A_d', C_d',P_d);                                               % Przypisanie do zmiennej wzmocnienia estymatora
Css_d=reg(Gss_d,K_d,L_d');                                                 % Stworzenie dyskretnego modelu regulatora LQR
Gclz_od = feedback(Gss_d,Css_d,1);                                         % Zamkniecie modelu ukladu regulacji petla ujemnego sprzezenia zwrotnego
F_od = 1/dcgain(Gclz_od);                                                  % Obliczenie wartosci wzmocnienia prekompensatora dla dyskretnego ukladu regulacji z estymatorem
Gclu_od = feedback(Gss_d*Css_d,1,1);                                       % Stworzenie modelu zamknietego, gdzie szeregowo polaczony jest dyskretny obiekt regulacji i dyskretny regulator LQR 
figure(1)                                                                  % Funkcja wywolujaca okno do rysowania wykresow
subplot(2,3,i)                                                             % Wywolanie funkcji rysujacej wykres na odpowiednim miejscu. Zmienna 'i' okresla poziome polozenie konkretnego wykresu
step(F_q*Gclz_q, 'g', F_o*Gclz_o, 'r', F_qd*Gclz_qd, 'b',...               % Wyrysowanie odpowiedzi skokowych utworzonych wczesniej modeli
F_od*Gclz_od, 'y', tfinal)      
title(['Odp. skokowa ukladu zamknietego R=', num2str(R)])                  % Nadanie tytulu wykresowi
grid;                                                                      % Wlaczenie siatki	
xlim([0 15]);                                                              % Definiowanie wartosci granicznych na wykresie
ylim([0 1.1]);
legend('ciagly, bez estymatora', 'ciagly, z estymatorem',...
    'dyskretny, bez estymatora', 'dyskretny, z estymatorem');              % Dodanie legendy


subplot(2,3,i+3)                                                           % Wywolanie funkcji rysujacej wykres na odpowiednim miejscu
step(F_q*Gclu_q, 'g', F_o*Gclu_o, 'r', F_qd*Gclu_qd, 'b',...               % Wyrysowanie odpowiedzi skokowych utworzonych wczesniej modeli
F_od*Gclu_od, 'y', tfinal) 
title(['Sygnal ster. dla R=', num2str(R)])                                 % Nadanie tytulu wykresowi                                     
grid;                                                                      % Wlaczenie siatki
legend('ciagly, bez estymatora', 'ciagly, z estymatorem',...
'dyskretny, bez estymatora', 'dyskretny, z estymatorem');


R=1;                                                                       % Macierz wag (kara jaka obciazany jest sygnal sterowania)
K=lqr(A, B, Q, R);                                                         % Rozwiazanie rownania Ricattiego
Acl=A-B*K;                                                                 % Obliczenie macierzy dla calego zamknietego obiektu regulacji, z uwzglednieniem sterowania regulatorem statycznym LQR
Gclz_q=ss(Acl,B,C,0);                                                      % Utworzenie modelu w przestrzeni stanu sterowania z regulatorem statycznym bez obserwatora
F_q=1/dcgain(Gclz_q);                                                      % Obliczenie wzmocnienia prekompensatora 
Gclu_q=ss(Acl,B,-K,0);                                                     % Utworzenie modelu w przestrzeni stanu ukladu regulacji z regulatorem statycznym                                             
P=eig(Acl);                                                                % Wyodrebnienie wartosci wlasnych macierzy Acl
L=place(A' ,C',3*P);                                                       % Przypisanie do zmiennej L wzmocnienia estymatora
Css=reg(Gss,K,L');                                                         % Stworzenie modelu regulatora LQR, o wzmocnieniu K, oraz estymatora o wzmocnieniu L                                                   
Gclz_o=feedback(Gss,Css,1);                                                % Stworzenie zamknietego ukladu regulacji, z modelem LQR w petli ujemnego sprzezenia zwrotnego
F_o=1/dcgain(Gclz_o);                                                      % Obliczenie wzmocnienia prekompensatora dla ukladu regulacji z estymatorem
Gclu_o=feedback(Gss*Css,1,1);                                              % Stworzenie modelu zamknietego ukladu regulacji z regulatorem LQR i estymatorem
tfinal=15;                                                                 % Maksymalny czas cyklu
                                                                           % Wyznacznie modelu dyskretnego
                                                                           % Czas probkowania dobrany na podstawie czasu narastania dla ukladu zamknietego z obserwatorem                                                                                                                                                      
y = stepinfo(F_o*Gclz_o);                                                  % Zapisanie w zmiennej wartosci potrzebnej do obliczenia czasu probkowania
Ts = y.RiseTime/12;                                                        % Obliczenie czasu probkowania
Gss_d=c2d(Gss,Ts);                                                         % Zamiana modelu ciaglego na model dyskretny z czestotliwoscia probkowania Ts
gs = tf(Gss_d);                                                            % Tworzymy model transmitancyjny korzystajac z modelu dyskretnego
[A_d,B_d,C_d,D_d, Ts] = ssdata(Gss_d);                                     % Utworzenie macierzy modelu stanu na podstawie transmitancji G  
K_d = dlqr(A_d, B_d, Q, R);                                                % Stworzenie modelu regulatora LQR z uwzglednieniem macierzy wag Q i R 
Acl_d = A_d-B_d*K_d;                                                       % Stworzenie modelu zamknietego ukladu regulacji, z regulatorem LQR w dziedzinie dyskretnej
Gclz_qd = ss(Acl_d, B_d, C_d, 0, Ts);                                      % Utworzenie modelu w przestrzeni stanu dla sterowania z regulatorem statycznym (Bez obserwatora)
F_qd = 1/dcgain(Gclz_qd);                                                  % Wyznaczenie czlonu kompensacyjnego, potrzebnego do wzmocnienia sygnalu w stanie ustalonym (Bez obserwatora)
Gclu_qd = ss(Acl_d, B_d, -K_d, 0, Ts);                                     % Stworzenie modelu zamknietego ukladu regulacji, z regulatorem LQR
P_d= eig(Acl_d);                                                           % Wyodrebnienie wartosci wlasnych macierzy Acl_d
L_d = place(A_d', C_d',P_d);                                               % Przypisanie do zmiennej wzmocnienia estymatora
Css_d=reg(Gss_d,K_d,L_d');                                                 % Stworzenie dyskretnego modelu regulatora LQR
Gclz_od = feedback(Gss_d,Css_d,1);                                         % Zamkniecie modelu ukladu regulacji petla ujemnego sprzezenia zwrotnego
F_od = 1/dcgain(Gclz_od);                                                  % Obliczenie wartosci wzmocnienia prekompensatora dla dyskretnego ukladu regulacji z estymatorem
Gclu_od = feedback(Gss_d*Css_d,1,1);                                       % Stworzenie modelu zamknietego, gdzie szeregowo polaczony jest dyskretny obiekt regulacji i dyskretny regulator LQR 
i = i+1;
figure(1)                                                                  % Funkcja wywolujaca okno do rysowania wykresow
subplot(2,3,i)                                                             % Wywolanie funkcji rysujacej wykres na odpowiednim miejscu. Zmienna 'i' okresla poziome polozenie konkretnego wykresu
step(F_q*Gclz_q, 'g', F_o*Gclz_o, 'r', F_qd*Gclz_qd, 'b',...               % Wyrysowanie odpowiedzi skokowych utworzonych wczesniej modeli
F_od*Gclz_od, 'y', tfinal)      
title(['Odp. skokowa ukladu zamknietego R=', num2str(R)])                  % Nadanie tytulu wykresowi
grid;                                                                      % Wlaczenie siatki	
xlim([0 15]);                                                              % Definiowanie wartosci granicznych na wykresie
ylim([0 1.1]);
legend('ciagly, bez estymatora', 'ciagly, z estymatorem',...
    'dyskretny, bez estymatora', 'dyskretny, z estymatorem');              % Dodanie legendy


subplot(2,3,i+3)                                                           % Wywolanie funkcji rysujacej wykres na odpowiednim miejscu
step(F_q*Gclu_q, 'g', F_o*Gclu_o, 'r', F_qd*Gclu_qd, 'b',...               % Wyrysowanie odpowiedzi skokowych utworzonych wczesniej modeli
F_od*Gclu_od, 'y', tfinal) 
title(['Sygnal ster. dla R=', num2str(R)])                                 % Nadanie tytulu wykresowi                                     
grid;                                                                      % Wlaczenie siatki
legend('ciagly, bez estymatora', 'ciagly, z estymatorem',...
'dyskretny, bez estymatora', 'dyskretny, z estymatorem');


R=10;                                                                      % Macierz wag (kara jaka obciazany jest sygnal sterowania)
K=lqr(A, B, Q, R);                                                         % Rozwiazanie rownania Ricattiego
Acl=A-B*K;                                                                 % Obliczenie macierzy dla calego zamknietego obiektu regulacji, z uwzglednieniem sterowania regulatorem statycznym LQR
Gclz_q=ss(Acl,B,C,0);                                                      % Utworzenie modelu w przestrzeni stanu sterowania z regulatorem statycznym bez obserwatora
F_q=1/dcgain(Gclz_q);                                                      % Obliczenie wzmocnienia prekompensatora 
Gclu_q=ss(Acl,B,-K,0);                                                     % Utworzenie modelu w przestrzeni stanu ukladu regulacji z regulatorem statycznym                                             
P=eig(Acl);                                                                % Wyodrebnienie wartosci wlasnych macierzy Acl
L=place(A' ,C',3*P);                                                       % Przypisanie do zmiennej L wzmocnienia estymatora
Css=reg(Gss,K,L');                                                         % Stworzenie modelu regulatora LQR, o wzmocnieniu K, oraz estymatora o wzmocnieniu L                                                   
Gclz_o=feedback(Gss,Css,1);                                                % Stworzenie zamknietego ukladu regulacji, z modelem LQR w petli ujemnego sprzezenia zwrotnego
F_o=1/dcgain(Gclz_o);                                                      % Obliczenie wzmocnienia prekompensatora dla ukladu regulacji z estymatorem
Gclu_o=feedback(Gss*Css,1,1);                                              % Stworzenie modelu zamknietego ukladu regulacji z regulatorem LQR i estymatorem
tfinal=15;                                                                 % Maksymalny czas cyklu


                                                                           % Wyznacznie modelu dyskretnego                                                                                                                                                                                                                              
y = stepinfo(F_o*Gclz_o);                                                  % Zapisanie w zmiennej wartosci potrzebnej do obliczenia czasu probkowania
Ts = y.RiseTime/12;                                                        % Obliczenie czasu probkowania
Gss_d=c2d(Gss,Ts);                                                         % Zamiana modelu ciaglego na model dyskretny z czestotliwoscia probkowania Ts
gs = tf(Gss_d);                                                            % Tworzymy model transmitancyjny korzystajac z modelu dyskretnego
[A_d,B_d,C_d,D_d, Ts] = ssdata(Gss_d);                                     % Utworzenie macierzy modelu stanu na podstawie transmitancji G  
K_d = dlqr(A_d, B_d, Q, R);                                                % Stworzenie modelu regulatora LQR z uwzglednieniem macierzy wag Q i R 
Acl_d = A_d-B_d*K_d;                                                       % Stworzenie modelu zamknietego ukladu regulacji, z regulatorem LQR w dziedzinie dyskretnej
Gclz_qd = ss(Acl_d, B_d, C_d, 0, Ts);                                      % Utworzenie modelu w przestrzeni stanu dla sterowania z regulatorem statycznym (Bez obserwatora)
F_qd = 1/dcgain(Gclz_qd);                                                  % Wyznaczenie czlonu kompensacyjnego, potrzebnego do wzmocnienia sygnalu w stanie ustalonym (Bez obserwatora)
Gclu_qd = ss(Acl_d, B_d, -K_d, 0, Ts);                                     % Stworzenie modelu zamknietego ukladu regulacji, z regulatorem LQR
P_d= eig(Acl_d);                                                           % Wyodrebnienie wartosci wlasnych macierzy Acl_d
L_d = place(A_d', C_d',P_d);                                               % Przypisanie do zmiennej wzmocnienia estymatora
Css_d=reg(Gss_d,K_d,L_d');                                                 % Stworzenie dyskretnego modelu regulatora LQR
Gclz_od = feedback(Gss_d,Css_d,1);                                         % Zamkniecie modelu ukladu regulacji petla ujemnego sprzezenia zwrotnego
F_od = 1/dcgain(Gclz_od);                                                  % Obliczenie wartosci wzmocnienia prekompensatora dla dyskretnego ukladu regulacji z estymatorem
Gclu_od = feedback(Gss_d*Css_d,1,1);                                       % Stworzenie modelu zamknietego, gdzie szeregowo polaczony jest dyskretny obiekt regulacji i dyskretny regulator LQR 
i = i+1;
figure(1)                                                                  % Funkcja wywolujaca okno do rysowania wykresow
subplot(2,3,i)                                                             % Wywolanie funkcji rysujacej wykres na odpowiednim miejscu. Zmienna 'i' okresla poziome polozenie konkretnego wykresu
step(F_q*Gclz_q, 'g', F_o*Gclz_o, 'r', F_qd*Gclz_qd, 'b',...               % Wyrysowanie odpowiedzi skokowych utworzonych wczesniej modeli
F_od*Gclz_od, 'y', tfinal)      
title(['Odp. skokowa ukladu zamknietego R=', num2str(R)])                  % Nadanie tytulu wykresowi
grid;                                                                      % Wlaczenie siatki	
xlim([0 15]);                                                              % Definiowanie wartosci granicznych na wykresie
ylim([0 1.1]);
legend('ciagly, bez estymatora', 'ciagly, z estymatorem',...
    'dyskretny, bez estymatora', 'dyskretny, z estymatorem');              % Dodanie legendy


subplot(2,3,i+3)                                                           % Wywolanie funkcji rysujacej wykres na odpowiednim miejscu
step(F_q*Gclu_q, 'g', F_o*Gclu_o, 'r', F_qd*Gclu_qd, 'b',...               % Wyrysowanie odpowiedzi skokowych utworzonych wczesniej modeli
F_od*Gclu_od, 'y', tfinal) 
title(['Sygnal ster. dla R=', num2str(R)])                                 % Nadanie tytulu wykresowi                                     
grid;                                                                      % Wlaczenie siatki
legend('ciagly, bez estymatora', 'ciagly, z estymatorem',...
'dyskretny, bez estymatora', 'dyskretny, z estymatorem');
%{
---------------------------------------------------------------------------
                         2.1 Wnioski ilosciowe:

2.1.1  Czasy regulacji dla modeli ciaglych bez estymatora wynosily kolejno:
     - Dla R = 0.1, tr = 4.1 [s]
     - Dla R = 1, tr = 4.35 [s]
     - Dla R = 10, tr = 10 [s]
    
2.1.2  Czasy regulacji dla modeli ciaglych z esymatorem wynosily kolejno:
     - Dla R = 0.1, tr = 4.39 [s]
     - Dla R = 1, tr = 4.61 [s]
     - Dla R = 10, tr = 10.3 [s] 
    
2.1.3  Czasy regulacji dla modeli dyskretnych bez estymatora wynosily 
       kolejno:
     - Dla R = 0.1, tr = 4.11 [s]
     - Dla R = 1, tr = 4.36 [s]
     - Dla R = 10, tr = 10 [s]
    
2.1.4  Czasy regulacji dla modeli dyskretnych z estymatorem wynosily 
       kolejno:
     - Dla R = 0.1, tr = 6.03 [s]
     - Dla R = 1, tr = 6.28 [s]
     - Dla R = 10, tr = 13.6 [s] 
   
2.1.5. Przeregulowanie nie wystepuje przy R = 0.1  
    
2.1.6  Przeregulowanie dla modeli ciaglych bez estymatora wynosilo kolejno:
     - Dla R = 1, P = 0.433 [%]
     - Dla R = 10, P = 2.51 [%]
    
2.1.7  Przeregulowanie dla modeli ciaglych z estymatorem wynosilo kolejno:
     - Dla R = 1, P = 0.416 [%]
     - Dla R = 10, P = 2.52 [%]
    
2.1.8  Przeregulowanie dla modeli dyskretnych bez estymatora wynosilo 
       kolejno:
     - Dla R = 1, P = 0.442 [%]
     - Dla R = 10, P = 2.52 [%]
    
2.1.9  Przeregulowanie dla modeli dyskretnych z estymatorem wynosilo 
       kolejno:
     - Dla R = 1, P = 0.426 [%]
     - Dla R = 10, P = 3.14 [%]
    
2.1.10 Czas ustalania sygnalu sterujacego dla modeli ciaglych bez 
       estymatora wynosil kolejno:
     - Dla R = 0.1, tr = 3.12 [s]
     - Dla R = 1, tr = 5.17 [s]
     - Dla R = 10, tr = 8.86 [s] 
    
2.1.11 Czas ustalania sygnalu sterujacego dla modeli ciaglych z estymatorem 
       wynosil kolejno:
     - Dla R = 0.1, tr = 4.4 [s]
     - Dla R = 1, tr = 6.03 [s]
     - Dla R = 10, tr = 9.72 [s]
   
2.1.12 Czas ustalania sygnalu sterujacego dla modeli dyskretnych bez 
       estymatora wynosil kolejno:
     - Dla R = 0.1, tr = 3.42 [s]
     - Dla R = 1, tr = 5.24 [s]
     - Dla R = 10, tr = 8.85 [s]
 
2.1.13 Czas ustalania sygnalu sterujacego dla modeli ciaglych z estymatorem 
       wynosil kolejno:
     - Dla R = 0.1, tr = 6.89 [s]
     - Dla R = 1, tr = 8.39 [s]
     - Dla R = 10, tr = 17.2 [s]
 
2.1.14 Wartosć sygnalu sterujacego dla modeli ciaglych bez estymatora 
       wynosi kolejno:
     - Dla R = 0.1, wynosi -3.16
     - Dla R = 1, wynosi -1
     - Dla R = 10, wynosi -0.316 
    
2.1.15 Wartosć sygnalu sterujacego dla modeli ciaglych z estymatorem wynosi
       kolejno:
     - Dla R = 0.1, wynosi -1.12
     - Dla R = 1, wynosi -0.474
     - Dla R = 10, wynosi -0.168

2.1.16 Wartosć sygnalu sterujacego dla modeli dyskretnych bez estymatora 
       wynosi kolejno:
     - Dla R = 0.1, wynosi -2.11
     - Dla R = 1, wynosi -0.815
     - Dla R = 10, wynosi -0.272
    
2.1.17 Wartosć sygnalu sterujacego dla modeli dyskretnych z estymatorem
       wynosi kolejno:
     - Dla R = 0.1, wynosi -0.346
     - Dla R = 1, wynosi -0.169
     - Dla R = 10, wynosi -0.0637    
---------------------------------------------------------------------------    
    
                         2.2 Wnioski jakosciowe:
    
2.2.1 Wraz ze zwiekszaniem wartosci R, zmniejsza sie wartosć sygnalu 
      sterowania.

2.2.2 Wraz ze wzrostem wartosci R, czas narastania odpowiedzi
      skokowej Tn rosnie, a co za tym idzie, rosnie takze czas probkowania
      Ts.

2.2.3 Wraz ze wzrostem wartosci macierzy wag R, rosnie czas ustalania.
    
2.2.4 Zwiekszenie wartosci R, czyli kary nakladanej na sygnal sterowania,
      skutecznie ten sygnal zmniejszalo.
    
2.2.5 W pierwotnym programie bazowym, uzywalismy do rysowania odpowiedzi 
      skokowej dla ukladu z obserwatorem prekompensatora, ktory byl
      wyznaczony dla ukladu bez obserwatora.
      Zamienilismy wyrazenie (F_q*Gclz_o) na (F_o*Gclz_o).

2.2.6 Modele dyskretne z estymatorem znaczaco odbiegaja od pozostalych
      modeli. Najwolniej osiagaja one wartosć ustalona. 

2.2.7 Modele z esymatorem cechowaly sie gorsza dynamika przy nizszej
      wartosci kosztu J. 
---------------------------------------------------------------------------
                         3. Rekomendacja:
    
3.1   Jezeli zalezy nam na braku przeregulowania i szybkim czasie 
      regulacji, najlepiej sprawdza sie regulatory bez estymatora z R=0.1.

3.2   Regulator dyskretny z estymatorem ma niezaleznie od wartosci R
      najdluzszy czas regulacji i najwieksze przeregulowanie (jezeli
      wystepuje).

3.3   Wsrod opracowanych regulatorow, najgorzej wypadaja te z R=10, 
      posiadaja najwieksze czasy ustalania oraz najwieksze wartosci 
      przeregulowan.

%}