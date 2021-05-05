%{
'Dyskretny regulator Kesslera w serwomechanizmie DC'
Wyznaczanie transmitancji dyskretnej regulatora Kesslera oraz modelu
numerycznego dyskretnego ukladu regulacji w oparciu o:
a) kryterium optimum modulu (z tablic i na podst. filtra modelu UZ)
b) kryterium symetryczne (z tablic i na podst. filtra modelu UZ)
oraz wyznaczyć przebiegi odpowiedzi skokowych omega(kTp) ukladow regulacji
z wyznaczonymi regulatorami dla przypadkow a) i b) dla przyjetego czasu Tp.
Autorzy: Badelek Piotr & Fusiara Jakub
Grupa: WMT18AP1S1


Przyjete oznaczenia:
cefi - iloczyn stalej elektrycznej i strumienia magnetycznego 
Tm - mechaniczna stala czasowa silnika DC
Te - elektryczna stala czasowa silnika DC
Tp - czas probkowania
Ko - wzmocnienie obiektu
Go - transmitancja zastepcza toru glownego sterowania
Tsigma1 - stala czasowa czlonu PWM
Tsigma2 - stala czasowa przetwornika C/A
Tsigma3 - stala czasowa wzmacniacza
%}

%Utworzenie modelu silnika pradu stalego dla zadanych parametrow:
cefi = 2.62;                                                               %[V*s]
Tm = 0.18;                                                                 %[s]
Te = 0.03;                                                                 %[s]
Tsigma1 = 0.02;                                                            %[s]
Tsigma2 = 0.06;                                                            %[s]
Tsigma3 = 0.01;                                                            %[s]

Tsigma = Tsigma1 + Tsigma2 + Tsigma3;                                      % Zastepcza stala czasowa
Tp = 0.15;                                                                 % Okreslenie czasu probkowania
%--------------------------------------------------------------------------
%                      1. Modele numeryczne 
s=tf('s');                                                                 % Zdefiniowanie zmiennej zespolonej s 
%{
                                                                            Postać transmitancji toru glownego sterowania z czlonem o dominujacej
                                                                            stalej Tm oraz czlonem reprezentujacym dynamike zastepcza:
%}
Ko=1/cefi;                                                                 % Wzmocnienie obiektu                     
Go=Ko/((Tm*s+1)*(Tsigma*s+1));                                             % Transmitancja zastepcza toru glownego ster.

% a) Kryterium optimum modulu

% Kryterium optimum modulu z tablic
Tr1=Tm;                                                                    % Kompensacja dominujacej stalej czasowej
Kr1=Tm/(2*Ko*Tsigma);                                                      % Nastawy, Zapis wzmocnienia regulatora
Greg_M=Kr1*(1+Tr1*s)/s;                                                    % Zapis transmitancji ciaglej regulatora
%dyskretyzacja
GkM_d = c2d(Greg_M, Tp, 'zoh');                                            % Dyskretyzacja transmitancji ciaglej regulatora
Go_d = c2d(Go, Tp, 'zoh');                                                 % Dyskretyzacja transmitancji zastepczej toru glownego ster.
GUR_M1_d = GkM_d*Go_d/(1+GkM_d*Go_d);                                      % Transmitancja dyskretna ukladu regulacji

                                                                           % Zdefiniowanie czasu symulacji
% Kryterium optimum modulu na podstawie filtra UZ
GZ_M2=1/(1+2*Tsigma*s+2*Tsigma^2*s^2);                                     % Transmitancja ukladu zamknietego
Greg_M2=GZ_M2/((1-GZ_M2)*Go);                                              % Transmitancja regulatora
%dyskretyzacja
GUR_M2_d=c2d(GZ_M2,Tp,'zoh');                                              % Dyskretyzacja transmitancji ukladu regulacji
Greg_M2_d=c2d(Greg_M2,Tp,'zoh');                                           % Dyskretyzacja transmitancji regulatora


% b) Kryterium symetryczne

%Kryterium symetryczne z tablic
Tr = 4*Tsigma;                                                             % Obliczenie stalej czasowej regulatora                                                           
Kr2= Tm/(8*Ko*Tsigma^2);                                                   % Nastawy, Zapis wzmocnienia regulatora
Greg_S= Kr2*(1+Tr*s)/s;                                                    % Transmitancja ciagla regulatora
GUR_S= Greg_S*Go/(1+Greg_S*Go);                                            % Transmitancja ciagla ukladu regulacji
% Dyskretyzacja 
GkS_d= c2d(Greg_S,Tp,'zoh');                                               % Dyskretyzacja transmitancji ciaglej regulatora
GUR_S1_d=c2d(GUR_S,Tp,'zoh');                                              % Transmitancja dyskretna ukladu regulacji


% Kryterium symetryczne na podstawie filtra UZ
GUR_S2=(4*Tsigma*s+1)/(8*Tsigma^3*s^3+8*Tsigma^2*s^2+4*Tsigma*s+1);        % Transmitancja ciagla ukladu zamknietego
Greg_S2=GUR_S2/((1-GUR_S2)*Go);                                            % Transmitancja ciagla regulatora
%Dyskretyzacja
GUR_S2_d=c2d(GUR_S2,Tp,'zoh');                                             % Dyskretyzacja transmitancji ukladu zamknietego
Greg_S2_d=c2d(Greg_S2,Tp,'zoh');                                           % Dyskretyzacja transmitancji regulatora 

tfinal = 15;
figure(1),step(GUR_M1_d,GUR_S1_d,tfinal)                                   % Odpowiedz skokowa ukladu
title('Odp. skokowa omega(kTp) dla obiektu z reg. Kesslera dla metody tablicowej')
legend('Kryterium modulu', 'Kryterium symetryczne')
grid
figure(2),step(GUR_M2_d,GUR_S2_d,tfinal)                                   % Odpowiedz skokowa ukladu
title('Odp. skokowa omega(kTp) dla obiektu z reg. Kesslera na podstawie filtra modelu UZ')
legend('Kryterium modulu', 'Kryterium symetryczne')
grid

%{
---------------------------------------------------------------------------
                         2.1 Wnioski ilosciowe

 2.1.1 Kryterium modulu z tablic:
         Czas narastania Tn = 1.95 [s]
         Czas regulacji Tr = 3 [s]
         Przeregulowanie P = 0 [%]
         Uchyb ustalony E = 0 

 2.1.2 Kryterium modulu na podstawie filtra UZ:
         Czas narastania Tn = 0.317 [s]
         Czas regulacji Tr = 0.762 [s]
         Przeregulowanie P = 4.18 [%]
         Uchyb ustalony E = 0 

 2.1.3 Kryterium symetryczne z tablic:
         Czas narastania Tn = 0.845 [s]
         Czas regulacji Tr = 1.8 [s]
         Przeregulowanie P = 0 [%] 
         Uchyb ustalony E = 0

 2.1.4 Kryterium symetryczne na podstawie filtra UZ:
         Czas narastania Tn = 0.224 [s] 
         Czas regulacji Tr = 1.49 [s]
         Przeregulowanie P = 40.3 [%]
         Uchyb ustalony E = 0

---------------------------------------------------------------------------
                        2.2 Wnioski jakosciowe

 2.2.1 Najdluzszym czasem narastania sygnalu odpowiedzi charakteryzowal
       sie uklad wykorzystujacy kryterium optimum modulu, wyznaczony za
       pomoca tablic. Najkrotszy czas narastania otrzymalismy dla ukladu
       wykorzystujacego kryterium symetryczne na podstawie filtra ukladu
       zamknietego. 

 2.2.2 Najdluzszy czas regulacji otrzymalismy dla ukladu wykorzystujacego
       optimum modulu, wyznaczony za pomoca tablic. Najkrotszy czas
       regulacji zostal uzyskany przy ukladzie wykorzystujacym optimum
       modulu na podstawie filtra UZ. 

 2.2.3 Najwieksze przeregulowanie sygnalu odpowiedzi uzyskalismy dla ukladu
       wykorzystujacego kryterium symetryczne na podstawie filtru UZ.Dla
       ukladow wykorzystujacych metode z tablic nie uzyskalismy
       przeregulowan. 

 2.2.4 Dla wszystkich wyznaczonych odpowiedzi skokowych uzyskalismy
       zerowy uchyb ustalony.

---------------------------------------------------------------------------
                         3.  Rekomendacja

   3.1 Powyzej zastosowane kryteria zapewniaja dokladnosć serwomechanizmu
       dzieki uzyskanym zerowym uchybom ustalonym odpowiedzi skokowych.

   3.2 Najszybsze uzyskanie i utrzymanie zadanej predkosci obrotowej
       serwomechanizmu zapewnia zastosowanie regulatora Kesslera z
       wykorzystaniem kryterium optimum modulu na podstawie filtra ukladu
       zamknietego.
%}
%}