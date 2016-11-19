Function
% Funcion que usa el metodo de Crank-Nicolson para aproximar la solucion de la ecuacion diferencial parcial:
%	du/dt (x,t) - alpha^2 d^2u/dx^2 (x,t) = 0
% Sujeta a las condicones de frontera:
%	u(0,t) = u(L,t) = 0, 0<t<T
% Sujeta a las condiciones iniciales:
%	u(x,0) = f(x), 0<=x<=L
% Los argumentos de la función son:
%	L : Limite superior de la variable espacial (x).
%	T :Limite superior de la variable temporal (t).
%	alpha :Raiz cuadrada de el coeficiente de la segunda derivada espacial.
%	m :Numero de intervalos espaciales para la creacion de la malla.
%	n :Numero de intervalos Temporales para la creacion de la malla.