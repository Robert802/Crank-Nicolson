function w = crnkNclsn(L,T,alpha,m,n,funk)
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
%	w :Aproximación a u(x,t) de manera discreta en x y t.
%	funk: Funcion que representa el valor inicial de u(x,0) = f(x).

%inicializar h, k, lambda, y w
h = L./m;
k = T./n;
lambda = (alpha.^2).*k./(h.^2);
w = sparse(m+1,n);

%inicializar primera columna de la matriz w
w(2:m,1) = funk(h.*(1:(m-1)).')

%Se procede a construir las matrices A y B con base en lambda
Adiag = repmat([-lambda./2 1+lambda -lambda./2],m-1,1);
Bdiag = repmat([lambda./2 1-lambda lambda./2],m-1,1);
A = spdiags(Adiag,[-1 0 1],m-1,m-1)
B = spdiags(Bdiag,[-1 0 1],m-1,m-1)

%Se completa la aproximacion w con un ciclo para solucionar el sistema:
%	A w(:,i+1) = B w(:,i)
for i = 1: (n-1)
	w(2:m,i+1) = A\(B*w(2:m,i));
end

%Para presentar la matriz la pasamos de sparse a su version completa.
w = full(w);
end


