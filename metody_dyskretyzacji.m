clc
clear

% Dane
T0 = 5 * 0.25;
Tp = 0.1 * T0;
k = 2.5;
 
% transmitancja G(s)
mian_G = [T0 1];
Gs = tf(k, mian_G);
step(Gs)
hold on

% metoda Tustina
licznik_Tustin = [k*Tp/(2*T0+Tp) k*Tp/(2*T0+Tp)];
mian_Tustin = [1 -(2*T0-Tp)/(2*T0+Tp)];
H = tf(licznik_Tustin, mian_Tustin, Tp);
step(H)
hold on

% metoda Eulera „wstecz"
licznik_Euler_wstecz = [(k*Tp)/(Tp+T0)];
mian_Euler_wstecz = [1 -(T0/(Tp+T0))];
F = tf(licznik_Euler_wstecz, mian_Euler_wstecz, Tp);
step(F)
hold on

% metoda Eulera „wprzód"
licznik_Euler_wprzod = [(k*Tp)/T0];
mian_Euler_wprzod = [1 ((Tp-T0)/T0)];
E = tf(licznik_Euler_wprzod, mian_Euler_wprzod, Tp);
step(E)
hold on

% c2d czyli dyskretyzacja matlaba
Gd = c2d(Gs,Tp);
step(Gd)
hold on
legend('G(s)','Metoda Tustina','Metoda Eulera „wstecz"','Metoda Eulera „wprzód"','Funkcja dyskretyzacji Matlaba')

    
    

