clc; clear; close all;
figure('Renderer', 'painters', 'Position', [1900 1200 600 400])
% params
m = 4; 
n = 4;
kx = 1.02;
ky = 0.98;
lx = 1.0;
ly = 1.0;
U = 1;
Re = 100.0;
dt = 0.01;
nu = m*n;
nv = m*(n-1);

u = zeros(1,m*n)';
v = zeros(1,m*(n-1))';

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
for I = 2:m
    dx(I) = kx*dx(I-1);
    x(I) = x(I-1) + dx(I-1);
end

% Y, DY gen
dy1 = ly  * (1-ky) /(1-ky^n);
if (abs(ky - 1) < 1e-9) 
    dy1= ly/n;
end
y(1) = 0.0;
y(n+1) = ly;
dy(1) = dy1;
for I = 2:n
    dy(I) = ky*dy(I-1);
    y(I) = y(I-1) + dy(I-1);
end

% R GEN
dyq = zeros(m*n, 1);
for I = 1:n
    dyq((I-1)*(m)+1:I*(m)) = repmat(dy(I), m, 1);
end
dxq = zeros((n-1)*m, 1); % Preallocate dxq with zeros
index = 1; % Initialize an index variable
for I = 1:n-1
    dxq(index:index+length(dx)-1) = dx; % Fill dxq
    index = index + length(dx); % Update the index
end
R = [diag([dyq; dxq])];


% M GEN
M = zeros((m-1)*n+m*(n-1),1);
ctr = 0;
rowShift = 0;
for j = 1:n
    for I = 2:m
        M((I-1)+rowShift) = 0.5*(dx(I)+dx(I-1));
        ctr = ctr+1;
    end
    ctr = ctr + 1;
    rowShift = ctr;
end
for I = 1:n
    M(I*m) = dx(m);
end
rowShift = (m)*n;
ctr = 0;
for I = 2:n
    for j = 1:m
        M(rowShift+(I-2)*(m)+j) = 0.5*(dy(I)+dy(I-1));
    end
end
M = diag(M);



%!!!!!!!!!!!!!!!!!!!! Laplace !!!!!!!!!!!!!!!!!!!!!!!!!!!
L = zeros(m*(n-1)+m*n,m*(n-1)+m*n);
Lxu = zeros(m*n,m*n);
Lyu = zeros(m*n,m*n);
Lxv = zeros(m*(n-1),m*(n-1));
Lyv = zeros(m*(n-1),m*(n-1));

%!!!!!!!!!!!!!!!! Lyu : d^2u/dy^2 !!!!!!!!!!!!!!!!!!!!!!!

Imax= (m)*n;
BC = zeros(m*n+m*(n-1),1);
for I = 1:m       
    % bottom
    hc = dyq(I);
    hn = 0.5*(dyq(I+m) + dyq(I)); 
    Lyu(I,I) = -2*(hn+hc)/(hc*hc*hn);
    Lyu(I,I + m) = 1/(hn*hc);
    
    % % top
    J = Imax + 1 - I;
    hc = dyq(J);
    hs = 0.5*(dyq(J-(m))+dyq(J));
    Lyu(J,J) = -(2*hc+hs)/(hc*hc*hs);
    Lyu(J,J - (m)) = 1/(hs*hc);

    % BCtop
    BC(J) = BC(J) + 2*U/hc/hc/Re;
end

% INNER PART
for I = m+1:Imax-(m)
    hc = dyq(I);
    hn = 0.5*(dyq(I)+dyq(I+m));
    hs = 0.5*(dyq(I)+dyq(I-(m)));

    Lyu(I,I) = - (hs+hn) /(hs*hc*hn); % u  times 
    Lyu(I,I + m) = 1/(hn*hc); % u east times
    Lyu(I,I -(m))= 1/(hs*hc); % u west times
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
    Lxu(I*(m)+1 , I*(m)+1 ) = -2/(hw*he);
    Lxu(I*(m)+1 , I*(m)+1 + 1 ) = 1/(he*hc);
   
    % ! right
    hc = dxLu((I+1) * (m-1));
    he = dxLu((I+1) * (m-1));
    hw = dxLu((I+1) * (m-1));

    % ! TWO POINT SCHEME
    Lxu((I+1)* (m) , (I+1)*(m) )    = -2/(hw*he) + 1/(hc*he);
    Lxu((I+1)*(m) , (I+1)*(m) -1 )  = 1/(hc*hw);

    % ! INTERIOR
    for J = I*(m-1)+2:(I+1)*(m-1)
        hc = 0.5*(dxLu(J)+dxLu(J-1));
        he = dxLu(J);
        hw = dxLu(J-1);
        Lxu( ctr , ctr ) = - (hw+he)/(hw*he*hc);
        Lxu( ctr , ctr+1 ) = 1/(he*hc);
        Lxu( ctr , ctr-1 ) = 1/(hw*hc);
        ctr = ctr + 1;
    end
    ctr = ctr + 2;
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
L(1:nu,1:nu) = Lxu + Lyu;
L(nu+1:nu+nv,nu+1:nu+nv) = Lxv + Lyv;

% A = (M*L*inv(R));
% spy(round((A-A')*100)/100);





% !! ADVECT !!
uu = zeros(1,m*n);
vv = zeros(1,m*n);
uv = zeros(1,(m+1)*(n+1));

u = u;
v = v;
uleftBC = ones(1,m);

ctr = 0;
% !!!!!!!!!!!!!!!!!!! computing uu, vv & uv !!!!!!!!!!!!!!!!!!!!!!!!!
for I = 1:m*n
    ctr = ctr+1;
    % ! indices in the form of matrix uu_ii,jj
    arg=floor(((I-1)/m)); %!arg=floor(float((I-1)/m));
    jj = round(arg) + 1;
    ii = mod ( (I-1),m) + 1;
    
    
    J1= I- 1;      %! index of u_i,j in the u vector
    Jr= I;  %! index of u_i+1,j in the u vector
    
    if (ii==1)
        uu(I) = ((uleftBC(jj)+u(Jr)) /2 )^2;
    else
        if (ii==m)
            uu(I) = ( (u(J1)+u(Jr)) /2 )^2;
        else
             uu(I) = ( (u(J1) + u(Jr)) /2 )^2;
        end
    end
        
    J2= (jj - 2)*(m) + ii;       %!index of v_i,j in the v vector
    Ju= (jj + 1 - 2)*(m) + ii;   %!index of v_i,j+1 in the v vector
    
    if (jj==1) 
        vv(I) = ( v(Ju) /2 )^2;
        else
        if (jj==n) 
            vv(I) = ( v(J2) /2 )^2;
        else
            vv(I) = ( (v(J2) + v(Ju)) /2 )^2;
        end
    end
    
    
%     Jlw = (jj-1 -1)*(m-1) + ii - 1;  %!index of u_i,j-1 in the u-vector
%     Jlf = (jj - 2)*(m) + ii - 1;     %!index of v_i-1,j in the v-vector
%     
%     if  ((ii>1) && (jj>1)) 
%         uv(I) =  ((dy(jj-1)*u(J1) + dy(jj)*u(Jlw))/(dy(jj)+dy(jj-1) )) * (((dx(ii-1)*v(J2) + dx(ii)*v(Jlf))) / (dx(ii)+dx(ii-1)));
%     else
%         if (jj == 1) 
%            uv(I) = 0;
%         else
%             if (ii == 1) 
%                uv(I) = 1;
%             end
%         end
%     end
%     fprintf("I = %3d(%1d,%d); J1 = %3d; Jr = %3d; J2 = %3d; Ju = %3d; Jlw = %3d; Jlf = %3d;\n",I,ii,jj,J1,Jr,J2,Ju,Jlw,Jlf);
end
ari = [];
arj = [];
uleftBCcorner = ones(1,m+1); % values of uv at the corner
% botBCcorners
for ii = 1:m+1
    uv(ii) = 0;
end
% leftBCcorners
for jj=2:n
    I = (jj-1)*(m+1)+1;
    ui = (uleftBC(jj-1) * dy(jj) + uleftBC(jj) * dy(jj-1))/(dy(jj) + dy(jj-1));
    vi = v(1 + (jj-2)*(m))/2; % inner nonzero, outerGhost is zero
    uv(I) = ui*vi;
    fprintf('I = %3d. ui=%3.1f;vi=%4.1f\n',I,ui,vi);
end
% topBCcorners
for ii = 1:m+1
    I = ii + (m+1)*n;
    uv(I) = uleftBC(m)/2;
    fprintf('I = %3d. ui=%3.1f;vi=%4.1f\n',I,ui,vi);
end
%rightBCcorners
for jj=2:n
    I = (jj)*(m+1);
    II = jj*m;
    JJ = II - m;
    ui = (u(II)*dy(jj-1) + u(II-m)*dy(jj))/(dy(jj)+dy(jj-1));
    vi = v(JJ)/2;
    uv(I)=ui*vi;
%     fprintf('I=%d; II = %d; JJ = %d\n',I,II,JJ)
fprintf('I = %3d. ui=%3.1f;vi=%4.1f\n',I,ui,vi);
end

for I = 1:(m+1)*(n+1)
    jj = floor((I-1)/(m+1))+1;
    ii = mod(I-1,(m+1))+1;

    if jj == 1
        uv(I) = 0;
    else
        if ii==1 && jj>1 && jj < n+1 % left side except lower corner
            
        end
    end

    ari = [ari ii];
    arj = [arj jj];
end

% ari
% arj
% plot(u,'x'); hold on;
% % plot(v,'o'); hold on;
% plot(uu,'o');hold on;
% % plot(vv);hold on;
% % plot(uv);hold on;
% legend('u','uu');