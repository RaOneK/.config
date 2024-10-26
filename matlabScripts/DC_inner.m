% OPTIMIZE
% precompute 1/dx 1/dy 1/dxMID, 1/dyMID. multiplication is faster than
% division.
clc; clear; close all;

m = 40; 
n = 40;
for i = 1:m-1
    blkDx(i,i)=1;
    blkDx(i+1,i)=-1;
end

for i = 1:n
    Dx(  m*(i-1)+1:m*(i-1)+m  ,  (i-1)*(m-1)+1:(i-1)*(m-1)+m-1  )= blkDx;
end

for i = 1:(n-1)*m
    Dy (i,i)= 1;
    Dy (i+m,i)= -1;
end

D( 1:m*n , 1: (m-1)*n) = Dx;
D( 1:m*n , (m-1)*n+1:(m-1)*n +(n-1)*m) = Dy;


for i =1:(n-1)*(m-1)
    Cy(i,i) = 1;
    Cy(i+(m-1),i)=-1;
end

for i=1:m-1
    Ccx(i,i)=-1;
    Ccx(i+1,i)=1;
end

for i = 1:n-1
    Cx((i-1)*m+1:(i-1)*m+m,(i-1)*(m-1)+1:(i-1)*(m-1)+(m-1))=Ccx;
end

C (1:n*(m-1) , 1:(n-1)*(m-1)) = Cy;
C (n*(m-1)+1:n*(m-1)+(n-1)*m , 1:(n-1)*(m-1)) = Cx;

G=-D';

%bdry nodes

figure();

x0=50;
y0=550;
width=500;
height=700;
set(gcf,'position',[x0,y0,width,height])
x = [0 24];
y = [0 9];
clim([-1 1])
imagesc((C))
title('C');
colorbar
figure();
x0=50;
y0=550;
width=700;
height=500;
set(gcf,'position',[x0,y0,width,height])
x = [0 16];
y = [0 24];
clim([-1 1])
imagesc((D))
title('D');
colorbar
figure();
x0=50;
y0=550;
width=700;
height=500;
set(gcf,'position',[x0,y0,width,height])
x = [0 16];
y = [0 24];
clim([-1 1])
imagesc(C')
title('CT');
colorbar

% close all;
% params
% m = 10; 
% n = 10;
kx = 1.2;
ky = 1.3;
lx = 1.0;
ly = 1.0;
U = 1;
Re = 100.0;
dt = 0.01;
nu = (m-1)*n;
nv = m*(n-1);
nq = nu + nv;

% x,y,dx,dy
x = zeros(1,m+1)';
y = zeros(1,n+1)';
dx = zeros(1,m)';
dy = zeros(1,n)';
% X, DX gem
dx1 = lx  * (1-kx) /(1-kx^m);
if (abs(kx - 1) < 1e-9) 
    dx1= lx/m;
end
x(1) = 0.0;
x(m+1) = lx;
dx(1) = dx1;
for i = 2:m
    dx(i) = kx*dx(i-1);
    x(i) = x(i-1) + dx(i-1);
end

% Y, DY gen
dy1 = ly  * (1-ky) /(1-ky^n);
if (abs(ky - 1) < 1e-9) 
    dy1= ly/n;
end
y(1) = 0.0;
y(n+1) = ly;
dy(1) = dy1;
for i = 2:n
    dy(i) = ky*dy(i-1);
    y(i) = y(i-1) + dy(i-1);
end

% R GEN
dyq = zeros((m-1)*n, 1);
for i = 1:n
    dyq((i-1)*(m-1)+1:i*(m-1)) = repmat(dy(i), m-1, 1);
end
dxq = zeros((n-1)*m, 1); % Preallocate dxq with zeros
index = 1; % Initialize an index variable
for i = 1:n-1
    dxq(index:index+length(dx)-1) = dx; % Fill dxq
    index = index + length(dx); % Update the index
end
R = [diag([dyq; dxq])];


% M GEN
M = zeros((m-1)*n+m*(n-1),1);
size(M)
ctr = 0;
rowShift = 0;
for j = 1:n
    for i = 2:m
        M((i-1)+rowShift) = 0.5*(dx(i)+dx(i-1));
        ctr = ctr+1;
    end
    % ctr = ctr + 1; % SHIFT BY ONE ROW
    rowShift = ctr;
end

% for i = 1:n % THATS EXTRA UNKNOWN ROW AT OUTLET LOOKS LIKE
%     M(i*m) = dx(m);
% end


rowShift = (m-1)*n;
ctr = 0;
for i = 2:n
    for j = 1:m
        M(rowShift+(i-2)*(m)+j) = 0.5*(dy(i)+dy(i-1));
    end
end
M = diag(M);

    nu = (m - 1) * n;
    nv = m * (n - 1);
    %!!!!!!!!!!!!!!!!!!!! Laplace !!!!!!!!!!!!!!!!!!!!!!!!!!!
    newL = zeros(m*(n-1)+(m-1)*n,m*(n-1)+(m-1)*n);
    Lxu = zeros((m-1)*n,(m-1)*n);
    Lyu = zeros((m-1)*n,(m-1)*n);
    Lxv = zeros(m*(n-1),m*(n-1));
    Lyv = zeros(m*(n-1),m*(n-1));

    %!!!!!!!!!!!!!!!! Lyu : d^2u/dy^2 !!!!!!!!!!!!!!!!!!!!!!!

    Imax= (m-1)*n;
    BC = zeros(m*n+m*(n-1),1);
    for I = 1:m-1
        % bottom
        hc = dyq(I);
        hn = 0.5*(dyq(I+(m-1)) + dyq(I)); 
        Lyu(I,I) = -2*(hn+hc)/(hc*hc*hn);
        Lyu(I,I + m-1) = 1/(hn*hc);

        % % top
        J = Imax + 1 - I;
        hc = dyq(J);
        hs = 0.5*(dyq(J-(m-1))+dyq(J));
        Lyu(J,J) = -(2*hc+hs)/(hc*hc*hs);
        Lyu(J,J - (m-1)) = 1/(hs*hc);

        % BCtop
        BC(J) = BC(J) + 2*U/hc/hc/Re;
    end

    % INNER PART
    for I = m:Imax-(m-1)
        hc = dyq(I);
        hn = 0.5*(dyq(I)+dyq(I+m-1));
        hs = 0.5*(dyq(I)+dyq(I-(m-1)));

        Lyu(I,I) = - (hs+hn) /(hs*hc*hn); % u  times 
        Lyu(I,I + m-1) = 1/(hn*hc); % u east times
        Lyu(I,I -(m-1))= 1/(hs*hc); % u west times
    end

    %!!!!!!!!!!!!!!!! Lxu : d^2u/dx^2 !!!!!!!!!!!!!!!!!!!!!!!

    for I = 1:n
        dxLu((I-1)*(m-1)+1:I*(m-1)) = dx(2:m);
    end

    ctr = 2; % row counter for inner part
    for I = 0:n-1
        %!left
        hc = 0.5*( dxLu(I*(m-1)+1)+dx(1) );
        hw = dx(1); %! SHOULD BE dx(1)/2 i think? NO dx is distance between two velocities on edges
        he = dxLu(I*(m-1)+1);
        %! VALUE AT THE FACE ON LEFT EDGE Lxu(I*(m-1)+1 , I*(m-1)+1 - 1) -> moved to RHS in BC subroutine
        Lxu(I*(m-1)+1 , I*(m-1)+1 ) = -2/(hw*he);
        Lxu(I*(m-1)+1 , I*(m-1)+1 + 1 ) = 1/(he*hc);

        % ! right
        hc = dxLu((I+1) * (m-1));
        he = dxLu((I+1) * (m-1));
        hw = dxLu((I+1) * (m-1));
        % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW % NEW 
        hc = 0.5*( dxLu((I+1)* (m-1)) + dxLu((I+1)*(m-1) -1));
        he = dxLu((I+1)* (m-1));
        hw = dxLu((I+1)*(m-1) -1 );

        % ! TWO POINT SCHEME
        Lxu((I+1)* (m-1) , (I+1)*(m-1) )    = -2/(hw*he) + 1/(hc*he);
        Lxu((I+1)*(m-1) , (I+1)*(m-1) -1 )  = 1/(hc*hw);

        % ! INTERIOR
        for J = I*(m-1)+2:(I+1)*(m-1)-1
            hc = 0.5*(dxLu(J)+dxLu(J-1));
            he = dxLu(J);
            hw = dxLu(J-1);
            Lxu( J , J ) = - (hw+he)/(hw*he*hc);
            Lxu( J , J+1 ) = 1/(he*hc);
            Lxu( J , J-1 ) = 1/(hw*hc);
            % ctr = ctr + 1;
        end
        % ctr = ctr + 2;
    end

    %%!!!!!!!!!!!!!!!!!! Lyv : d^2v/dy^2 !!!!!!!!!!!!!!!!
    for I=0:n-2
        dyLv(I*m+1:I*m+m)=dy(I+2);
    end
    Imax= m*(n-1);
    for I = 1:m
        %!bottom
        hc = 0.5*(dyLv(I)+dy(1));
        hn = dyLv(I);
        hs = dy(1);
        Lyv(I,I) = -(hs+hn)/(hs*hc*hn);
        Lyv(I,I + m) = 1/(hn*hc);

        %!top
        J = Imax + 1 - I;
        hc = 0.5*(dyLv(J)+dyLv(J-m));
        hn = dyLv(J);
        hs = dyLv(J-m);
        Lyv(J,J) = -(hs+hn)/(hs*hc*hn);
        Lyv(J,J - m ) = 1/(hs*hc);
    end
    % ! INNER
    for I = m + 1:Imax-(m);
        hc = 0.5*(dyLv(I)+dyLv(I-m));
        hn = dyLv(I);
        hs = dyLv(I-m);
        Lyv(I,I) = -(hs+hn)/(hs*hc*hn);
        Lyv(I,I + m) = 1/(hn*hc);
        Lyv(I,I - m)= 1/(hs*hc);
    end


    % !!!!!!!!!!!!!!!!!! Lxv : d^2v/dx^2 !!!!!!!!!!!!!!!!
    for I = 0:n-2
        % !left
        hc = dxq(I*m+1);
        he = 0.5*( dxq(I*m+1)+dxq(I*m+1+1)); %! was I+m+1+1 typo
        Lxv(I*m + 1 , I*m + 1 ) = - (2*he+hc)/(hc*hc*he);
        Lxv(I*m + 1 , I*m + 1 + 1 ) = 1/(hc*he);

        % !right
        hc = dxq((I+1)*m);
        hw = 0.5*( dxq((I+1)* m) + dxq((I+1)*m -1) );
        Lxv((I+1)* m , (I+1)*m ) = - (2*hw+hc)/(hc*hc*hw);
        Lxv((I+1)* m , (I+1)*m -1 ) = 1/(hc*hw);

        % ! INNER
        for J = I* m + 2:(I+1)* m -1
            hc = dxq(J);
            he = 0.5*(dxq(J)+dxq(J+1));
            hw = 0.5*(dxq(J)+dxq(J-1)); 
            Lxv( J , J ) = - (he+hw)/(he*hc*hw);
            Lxv( J , J+1 ) = 1/(he*hc);
            Lxv( J , J-1 ) = 1/(hw*hc);
        end
    end
    newL(1:nu,1:nu) = Lxu + Lyu;
    newL(nu+1:nu+nv,nu+1:nu+nv) = Lxv + Lyv;

newA = (M*newL*inv(R));
A = newA;
figure();
symmCheck = round((A-A')*1e12)/1e12;
spy(symmCheck); % precise up to 14 digits after comma