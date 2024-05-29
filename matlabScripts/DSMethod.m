clc; clear; close all;

M = 3; 
N = 4;
nt = 3;
dt = 1e-5;

kx = 1.05;
ky = 0.95;
lx = 1;
ly = 1;


% function names are self-explanatory
[dx, dy, x, y] = gridGen(M,N,kx,ky,lx,ly);
[dxq, dyq] = delQ(M,N,dx,dy);
[u,v,q,psi] = initialConditions(M,N,nt,dxq,dyq);
[D,G,C] = divGradCurlGen(M,N);
[L] = laplacianGen(M,N,dx,dy,dxq,dyq);
[rnadvect] = advection(M,N,dx,dy,u,v,dxq,dyq,nt,dt);

function [dx, dy, x, y] = gridGen(M,N,kx,ky,lx,ly)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           GRID GENERATION         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dx = zeros(M,1);
    dy = zeros(N,1);
    x = zeros(M+1,1);
    y = zeros(N+1,1);

    % dx, x
    dx1 = lx *(1-kx)/(1-kx^M);
    if abs(kx-1)<1e-10
        dx1=lx/M;
    end
    x(1) = 0.0;
    x(M+1) = lx;
    dx(1) = dx1;
    for i = 2:M
        dx(i) = kx*dx(i-1);
        x(i) = x(i-1)+dx(i-1);
    end
    
    % dy, y
    dy1 = ly  * (1-ky) /(1-ky^N);
    
    if abs(ky - 1)<1e-10
        dy1= ly/N;
    end
    
    y(1)= 0.0;
    y(N+1) = ly;
    dy(1) = dy1;
    
    for i = 2:N
        dy(i) = ky*dy(i-1);
        y(i) = y(i-1)+dy(i-1);
    end
end
function [dxq, dyq] = delQ(M,N,dx,dy)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       delQ has two arrays: dyq, dxq.                              %
    %       These are modified arrays of dx, dy to multiply q by.       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dyq = zeros((M+1)*N,1);     % (M+1)N dyq, u components times dy
    dxq = zeros(M*(N+1),1);     % M(N+1) dxq, v components times dx
    for i = 1:N                 % all cell rows, 
                                % M+1 element in each cell row (incl. bdry)
        jmin = (i-1)*(M+1)+1;   % start index
        jmax = i* (M+1);        % end index
        dyq(jmin:jmax)=dy(i);   % fill deltaY of this cell row
                                % repeats for each cell row
    end
    
    for i = 1:N+1               % all cell bot faces, total of N, 
                                % + one extra on top bdry
        jmin = (i-1)*M+1;       % start index
        jmax =  i*M;            % end index
        dxq(jmin:jmax) = dx;    % fill deltaX of faces, 
                                % row has distinct dx
    end
end
function [u,v,q,psi] = initialConditions(M,N,nt,dxq,dyq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Changes first columns of u,v,q,psi to initial conditions    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    nu = (M+1)*(N);
    nv = (M)*(N+1);
    nq = nu + nv;
    npsi = (M+1)*(N+1);

    u = zeros(nu,nt);
    v = zeros(nv,nt);
    q = zeros(nq,nt);
    psi = zeros(npsi,nt); % we do not need initial psi values

    u(:,1) = 1;     % uniform unitary horizontal u component
    v(:,1) = 1e-5;  % slightly positive uniform vertical v component
    
    q(:,1) = [u(:,1);v(:,1)].*[dyq; dxq]; %elementwise .* multiplication
end
function [D,G,C] = divGradCurlGen(M,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DIVERGENCE          %
%           & GRADIENT          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu = (M+1)*N;
nv = M*(N+1);
nc = M*N;
nn = (M+1)*(N+1);

D = zeros(nc,nu+nv);
Du = zeros(nc,nu);
Dv = zeros(nc,nv);

colShift = 0;
rowShift = 0; 
for j = 1:N % row iter
    for i = 1:M % fix row, run column
        Du(i+rowShift,i+colShift) = -1;
        Du(i+rowShift,i+1+colShift) = 1;
    end
    colShift = colShift + M+1;
    rowShift = rowShift + M;
end
% right part
for i = 1:nc
    Dv(i,i)=-1;
    Dv(i,i+M) = 1;
end
% concat
D = [Du Dv];

% GRAD
G = -D';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               CURL            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = zeros(nu+nv,nn);
Cu = zeros(nu,nn);
Cv = zeros(nv,nn);
% upper part
for i = 1:nu
    Cu(i,i)=-1;
    Cu(i,i+M+1) = 1;
end
% lower part
colShift = 0;
rowShift = 0; 
for j = 1:N+1
    for i = 1:M % run over each row
        Cv(i+rowShift,i+colShift) = 1;
        Cv(i+rowShift,i+1+colShift) = -1;
    end
    colShift = colShift + M+1;
    rowShift = rowShift + M;
end
%concat
C = [Cu;Cv];

end
% function [bcun,bcvn] = boundaryCond(M,N,u,v,dx,dy,dxq,dyq) % BC in Laplacian 
%     bcn = zeros()
% end



function [uLeftBC,uLeftBCVirtual] = LeftBC(M,N,nt,dt)      % hardcoded BC values at the left bdry
    uLeftBC = zeros((M+1)*N,1);             % only few will be filled, rest are 0
    uLeftBCVirtual = zeros(N,1);            % single column non-zero only at left bc
    for i = 1:N                             % go over all rows in grid
        iGlobal = (i-1)*(M+1)+1;            % get global index
        uLeftBC(iGlobal) = 0.2;             % fill bc value from left bdry
        uLeftBCVirtual(i) = 0.1;            % value outside the boundary at the next u
%         fprintf("%d\n",iGlobal);            % global index test print
    end
end
function [vBotBC] = BotBC(M,N,nt,dt)        % hardcoded BC values at the bot bdry
    vBotBC = zeros((N+1)*M,1);              % only few will be filled, rest are 0
    for i = 1:M                             % go over all cols in grid
        iGlobal = i;                        % get global index
        vBotBC(iGlobal) = 0;                % fill bc value from bot bdry
%         fprintf("%d\n",iGlobal);            % global index test print
    end
end
function [uRightBC,vRightBC] = RightBC(M,N,u,v)
    
end
% RIGHT AND TOP OPEN BC ARE HARDCODED INTO CORRESPONDING ADVECT/LAPLACIAN
% function [uuRightBC,vvRightBC] = RightBC(M,N,nt,dt,u,v) %
%     uuRightBC = zeros(N,1);     % all be filled
%     vvRightBC = zeros(N,1);     % all be filled
%         
%     
% end

function [rnadvect] = advection(M,N,dx,dy,u,v,dxq,dyq,nt,dt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           ADVECTION           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:(M+1)*N % TEST LOOP TO VERIFY ALGEBRA
        u(i) = i/((M+1)*N);
    end

    for i = 1:(N+1)*M % TEST LOOP TO VERIFY ALGEBRA
        v(i) = i/((N+1)*M);
    end
    nu = (M+1)*N;
    nv = M*(N+1);
    nq = nu+nv;
    nc = (M)*(N);
    nn = (M+1)*(N+1);
    rnadvect = zeros(nq,1);
    uu = zeros(nc,1);
    vv = zeros(nc,1);
    uv = zeros(nn,1);
    [uLeftBC,uLeftBCVirtual] = LeftBC(M,N,nt,dt);
    [vBotBC] =  BotBC(M,N,nt,dt);
    duudx = zeros(nu,1);
    duvdy = zeros(nu,1);
    advU = zeros(nu,1); % dummy, to keep in mind preallocated memory, array is used later
    dvvdy = zeros(nv,1);
    dvudx = zeros(nv,1);
    advV = zeros(nv,1); % dummy, to keep in mind preallocated memory, array is used later
    uuLeftVirtual = 0; % dummy, to keep in mind preallocated memory, value is reassigned later
    uRightVirtual = zeros(N,1);
    vTopVirtual = zeros(M,1);
    % go over all cells to precompute uu,vv at the cell cetres
    for iCell = 1:M*N
        col = mod(iCell-1, M) + 1;          % column index 
        row = ceil(iCell/M);                % row index 
        % indices computaiton might be different in other 
        % programming languages
        ijleft = (row-1)*(M+1) + col;
        ijright = ijleft+1;
        ijtop = (row)*M+col;
        ijbot = (row-1)*M+col;
        
        %RIGHT and TOP uu vv's values outside the boundary are computed
        %when evaluating duudx dvvdy
        if col==1       % left bdry inside
            uu(iCell) = 0.25*(uLeftBC(ijleft)+u(ijright))^2; % [(u_{i,j}+u_{i+1,j})/2]^2
        else            % inner part
            uu(iCell) = 0.25*(u(ijleft)+u(ijright))^2; % [(u_{i,j}+u_{i+1,j})/2]^2
        end
        if row == 1     % bot bdry inside
            vv(iCell) = 0.25*(vBotBC(ijbot)+v(ijtop))^2; % [(v_{i,j}+v_{i,j+1})/2]^2
        else            % inner part
            vv(iCell) = 0.25*(v(ijbot)+v(ijtop))^2; % [(v_{i,j}+v_{i,j+1})/2]^2
        end

%         fprintf("row:%3d\t" + ...
%             "col:%3d\t" + ...
%             "L Indx:%3d\t" + ...
%             "R Indx:%3d\t" + ...
%             "B Indx:%3d\t" + ...
%             "T Indx:%3d\n", ...
%             row, col,ijleft,ijright,ijbot,ijtop);
    end

    % compute uv
    for iuvGlobal = 1:(M+1)*(N+1)
        row = floor((iuvGlobal-1)/(M+1))+1;
        col = mod(iuvGlobal-1,M+1)+1;

        iTopU = iuvGlobal; % same as iuvGlobal
        iBotU = iuvGlobal - (M+1);
        iLeftV = col - 1 + (row-1)*(M);
        iRightV = iLeftV + 1;
        if (row == 1) || (col == 1) % BOT or LEFT
            uv(iuvGlobal) = 0;
            % we will use equality of duv/dx duv/dy of inner part
            % due to linear interpolation when computing central derivatives
        elseif (col == M+1) && (row ~= N+1) % RIGHT
%             uv(iuvGlobal) = -1e-4; % fill with distinct values for
%             % DEBUGGING
            % extrapolate V from inside to boundary FIXED FROM UV EXTRAP
            vExtrap = v(iLeftV-1)+(1.5*dx(col-1)+dx(col-2)/2)*(v(iLeftV)-v(iLeftV-1))/(dx(col-1)/2+dx(col-2)/2);
            uv(iuvGlobal) = ((u(iTopU)*dy(row-1)+u(iBotU)*dy(row))/(dy(row-1)+dy(row)))...
                *((vExtrap*dx(col-1)+v(iLeftV)*dx(col-1))/(dx(col-1)+dx(col-1)));
        elseif (row == N+1) && (col ~= M+1) % TOP
%             uv(iuvGlobal) = -1e-4; % fill with distinct values for
            % extrapolate from inside to boundary
            uExtrap = u(iBotU-(M+1))+(1.5*dy(row-1)+dy(row-2)/2)*(u(iBotU)-u(iBotU-(M+1)))/(dy(row-1)/2+dy(row-2)/2);
%             fprintf('uExtrap: %18.15f\tu(iBotU-(M+1)): %5.2f\t u(iBotU): %5.2f\n',uExtrap,u(iBotU-(M+1)),u(iBotU));
            uv(iuvGlobal) = ((uExtrap*dy(row-1)+u(iBotU)*dy(row-1))/(dy(row-1)+dy(row-1)))...
                *((v(iRightV)*dx(col-1)+v(iLeftV)*dx(col))/(dx(col-1)+dx(col)));
        elseif (col == M+1) && (row == N+1)
            % three point interpolation on rectangle 
            % half sum of opposite corner values are equal to each other
            % fLB + fTR = fLT + fRB
            vExtrap = v(iLeftV-1)+(1.5*dx(col-1)+dx(col-2)/2)*(v(iLeftV)-v(iLeftV-1))/(dx(col-1)/2+dx(col-2)/2);
            uExtrap = u(iBotU-(M+1))+(1.5*dy(row-1)+dy(row-2)/2)*(u(iBotU)-u(iBotU-(M+1)))/(dy(row-1)/2+dy(row-2)/2);
%             fprintf('uExtrap: %18.15f\t vExtrap: %18.15f\t u(iBotU-(M+1)): %5.2f\t u(iBotU): %5.2f\n',uExtrap,vExtrap,u(iBotU-(M+1)),u(iBotU));
%             uv(iuvGlobal) = uv(iuvGlobal-1) + uv(iuvGlobal-(M+1)) - uv(iuvGlobal-(M+1)-1);
            uv(iuvGlobal) = ((uExtrap*dy(row-1)+u(iBotU)*dy(row-1))/(dy(row-1)+dy(row-1)))...
                *((vExtrap*dx(col-1)+v(iLeftV)*dx(col-1))/(dx(col-1)+dx(col-1)));
        else
            uv(iuvGlobal) = ((u(iTopU)*dy(row-1)+u(iBotU)*dy(row))/(dy(row-1)+dy(row)))...
                *((v(iRightV)*dx(col-1)+v(iLeftV)*dx(col))/(dx(col-1)+dx(col)));
        end
%         fprintf("uv. col: %3d\t row: %3d\t iLeftV: %3d\t iRightV: %3d\t iTopU: %3d\t iBotU: %3d \n", col,row,iLeftV,iRightV,iTopU,iBotU);
    end

    % go over all velocities to compute central derivatives
    % HARDCODED dirichlet BC's LEFT BOT
    % SOFTCODED extrapolation values at outlet/open bdry
    % u components
    cellCtr = 0;
    for iuGlobal=1:(M+1)*N
        col = mod(iuGlobal-1,M+1)+1;
        row = floor((iuGlobal-1)/(M+1))+1;
        if col==1 % LEFT bdry
            cellCtr = cellCtr+1; % cell indices change
            % [(uleftBC + uleftBCVirtual)/2]^2
            uuLeftVirtual = 0.25*(uLeftBC(iuGlobal)+uLeftBCVirtual(row))^2; 
            % virtual dx outside is the same as inside, but values of uu
            % are not
            duudx(iuGlobal) = (uu(cellCtr)-uuLeftVirtual)/(dx(col)); 
            % since v = 0 at left bdry uv=0 as well 
            duvdy(iuGlobal) = 0; 
        elseif col==M+1 % RIGHT bdry
            % linear extrapolation from inside to outside means
            % derivatives in x direction are equal
%             duudx(iuGlobal) = duudx(iuGlobal-1); % due to linear
%             extrapolation above WAS WRONG, FIXED
            % extrapolate u outside, then use value to compute uu
            iMidU = iuGlobal;
            iLeftU = iMidU-1;
            uRightVirtual(row) = (1.5)*(u(iMidU)-u(iLeftU))+u(iLeftU); 
%             fprintf("col-1: %3d\t iMidU: %3d\t iLeftU: %3d\t uVirt:%18.15f\n",col-1,iMidU,iLeftU,uRightVirtual(row));
            duudx(iuGlobal) = 2*(uRightVirtual(row)*uRightVirtual(row) - uu(cellCtr-1))/(dx(col-1)+dx(col-1));
            % extrapolate uv from inside to the boundary manually
            % ACTUALLY is WRONG, function is quadratic (uv), need square
            % root extrapolation, not linear. 
%             uvBot = uv(iuGlobal-2)+(dx(col-1)+dx(col-2))*(uv(iuGlobal-1)-uv(iuGlobal-2))/dx(col-2);
%             uvTop = uv(iuGlobal-2+M+1)+(dx(col-1)+dx(col-2))*(uv(iuGlobal-1+M+1)-uv(iuGlobal-2+M+1))/dx(col-2);
%             duvdy(iuGlobal) = (uvTop-uvBot)/dy(row);
            % extrapolation in precomputed uv, use standard central
            % difference formula
            iuvBot = iuGlobal;
            iuvTop = iuGlobal + M+1;
            duvdy(iuGlobal) = (uv(iuvTop) - uv(iuvBot))/dy(row);
        else % INNER part
            cellCtr = cellCtr+1;
            duudx(iuGlobal) = 2*(uu(cellCtr) - uu(cellCtr-1))/(dx(col)+dx(col-1));
            iuvBot = iuGlobal;
            iuvTop = iuvBot + M+1;
            duvdy(iuGlobal) = (uv(iuvTop) - uv(iuvBot))/dy(row);
        end
%         fprintf("duu duv. col: %3d\t row: %3d\t \n", col,row);
    end
    advU = duudx + duvdy;
    
    % v components
    for ivGlobal=1:(N+1)*M
        col = mod(ivGlobal-1,M)+1;
        row = floor((ivGlobal-1)/(M))+1; 
        if row == 1 % bottom bdry
            vvVirtualBot = -vv(ivGlobal);   % use symmetric flow at wall
            uvLeft = 0;                     % wall bc v=0
            uvRight = 0;                    % wall bc v=0
            dvvdy(ivGlobal) = (vv(ivGlobal) - vvVirtualBot)/(dy(row));
            dvudx(ivGlobal) = (uvRight-uvLeft)/dx(col); % should be zero
        elseif row == N+1 % top bdry
            iuvLeft = (row-1)*(M+1)+col; 
            iuvRight = iuvLeft + 1;
            % linear extrapolation from inside to bdry
            % means derivatives in y direction are equal
%             dvvdy(ivGlobal) = dvvdy(ivGlobal-(M));
            % extrapolation above WAS WRONG, FIXED
            iBotV = ivGlobal-(M);
            iMidV = ivGlobal;
            vTopVirtual(col) = 1.5*(v(iMidV)-v(iBotV))+v(iBotV);
            % ivGlobal index is same as vv index at cell above
            dvvdy(ivGlobal) = (vTopVirtual(col)*vTopVirtual(col) - vv(ivGlobal-M))/dy(row-1);
            % extrapolation does not affect dvudx
            dvudx(ivGlobal) = (uv(iuvRight)-uv(iuvLeft))/dx(col);
%             fprintf('vTopVirtual(%3d): %4.2f \t dvvdy(%3d): %18.15f\t dvudx(%3d): %18.15f iBotV: %3d \t iMidV: %3d \n',col,vTopVirtual(col),ivGlobal,dvvdy(ivGlobal),ivGlobal,dvudx(ivGlobal),ivGlobal-(M),iMidV);
        else % inner part
            iuvLeft = (row-1)*(M+1)+col; 
            iuvRight = iuvLeft + 1;
            dvvdy(ivGlobal) = 2*(vv(ivGlobal) - vv(ivGlobal - M))/(dy(row)+dy(row-1));
            dvudx(ivGlobal) = (uv(iuvRight)-uv(iuvLeft))/dx(col);
        end
%         fprintf("dvvdy(%2d)=%19.15f\t dvudx(%2d)=%19.15f\t col: %3d\t row: %3d\t \n", ivGlobal, dvvdy(ivGlobal), ivGlobal,dvudx(ivGlobal), col,row); 
    end
    advV = dvvdy+dvudx;
    rnadvect = [advU; advV];
end
function [L] = laplacianGen(M,N,dx,dy,dxq,dyq)
    % Initialize dimensions
    nu = (M+1)*N;
    nv = M*(N+1);
    nq = nu+nv;
    L = zeros(nq,nq);

    % Initialize submatrices
    Lyu = zeros(nu,nu); % !!!!!!!!!!!!!!!!! Lyu : d^2u/dy^2 !!!!!!!!!!!!!!!!!!
    Lxu = zeros(nu,nu); % !!!!!!!!!!!!!!!!! Lxu : d^2u/dx^2 !!!!!!!!!!!!!!!!!!
    Lxv = zeros(nv,nv); % !!!!!!!!!!!!!!!!! Lxv : d^2v/dx^2 !!!!!!!!!!!!!!!!!!
    Lyv = zeros(nv,nv); % !!!!!!!!!!!!!!!!! Lyv : d^2v/dy^2 !!!!!!!!!!!!!!!!!!
    dxLu = zeros(1,nu); % Grid spacing for u in x-direction
    dyLv = zeros(1,nv); % Grid spacing for v in y-direction

    iUGlobalMax = (M+1)*N;

    % Compute Laplacian (Luyy) for u in y-direction (vertical)
    for iUGlobal = 1:iUGlobalMax
        row = ceil(iUGlobal/(M+1));
        col = mod(iUGlobal-1,M+1)+1;
        hc = dyq(iUGlobal);
        if row==1 % only row affects Luyy because of derivative direction 
            hn = 0.5*(dyq(iUGlobal)+dyq(iUGlobal+(M+1)));
            hs = hc;
            
            Lyu(iUGlobal,iUGlobal) = -2/(hs*hn)-1/(hs*hc); % main diag
            Lyu(iUGlobal,iUGlobal+(M+1)) = 1/(hn*hc); % top neighbour
            % bot neighbour is across the wall, hence interpolation
            
        elseif row==(N) % top row, affected by Dong's BC
%             fprintf('iUGlobal:%3d\t row: %3d\t col: %3d\n',iUGlobal,row, col);
        else % inner part
            hn = 0.5*(dyq(iUGlobal)+dyq(iUGlobal+(M+1)));
            hs = 0.5*(dyq(iUGlobal)+dyq(iUGlobal-(M+1)));
            Lyu(iUGlobal,iUGlobal) = -2/(hs*hn); % main diag
            Lyu(iUGlobal,iUGlobal+(M+1)) = 1/(hn*hc); % top neighbour
            Lyu(iUGlobal,iUGlobal-(M+1)) = 1/(hs*hc); % bot neighbour
        end
    end

    % Compute Laplacian for u in x-direction (horizontal)
    % can be combined with upper loop, but L is precomputed only once,
    % hence can afford
    idxqRight = 1; % index is changed everywhere except very last column
    for iUGlobal = 1:iUGlobalMax
        row = ceil(iUGlobal/(M+1));
        col = mod(iUGlobal-1,M+1)+1;
        
        if (col == 1)   %leftmost boundary part
            % no need to fill since it will be removed later. 
            % part kept for dimension consistency with other operators: i.e
            % curl div grad
            idxqRight = idxqRight + 1;
        elseif (col==2) % leftmost inner part
            he = dxq(idxqRight);
            hw = dxq(idxqLeft);
            hc = 0.5*(he+hw);
            
            Lxu(iUGlobal,iUGlobal-1) = 1/(hw*hc);   % left
            % No left element actually, since it is at the boundary. 
            % Corresponding columns and rows will be removeggjd later.
            % Explicit RHS is in separate module generating Laplacian BC.
            Lxu(iUGlobal,iUGlobal) = -2/(hw*he);    % main diag
            Lxu(iUGlobal,iUGlobal+1) = 1/(hc*he);   % right

            idxqRight = idxqRight + 1; 
            fprintf('iUGlobal:%3d\t row: %3d\t col: %3d\t idxRight: %3d\t idxLeft: %3d \n',iUGlobal,row, col,idxqRight,idxqLeft);
        elseif col == (M+1) % rightmost boundary DONG BC
            
        else % inner part
            he = dxq(idxqRight);
            hw = dxq(idxqLeft);
            hc = 0.5*(he+hw);
            Lxu(iUGlobal,iUGlobal-1) = 1/(hw*hc);   % left
            Lxu(iUGlobal,iUGlobal) = -2/(hw*he);    % main diag
            Lxu(iUGlobal,iUGlobal+1) = 1/(hc*he);   % right
            idxqRight = idxqRight + 1;
        end
        idxqLeft = idxqRight - 1;
    end



    % Compute Laplacian for v in y-direction (vertical)
    for i=1:N-1
%         dyLv((i-1)*M+1:i*M)=dy(i+1);
    end

    iUGlobalMax= M*(N-1);

    for I = 1:M
%         % Bottom boundary
%         hc = 0.5*(dyLv(I)+dy(1));
%         hn = dyLv(I);
%         hs = dy(1);
%         Lyv(I,I) = -(hs+hn)/(hs*hc*hn);
%         Lyv(I,I + M) = 1/(hn*hc);
% 
%         % Top boundary
%         J = iUGlobalMax + 1 - I;
%         hc = 0.5*(dyLv(J)+dyLv(J-M));
%         hn = dyLv(J);
%         hs = dyLv(J-M);
%         Lyv(J,J) = -(hs+hn)/(hs*hc*hn);
%         Lyv(J,J - M ) = 1/(hs*hc);
    end

    % Compute Laplacian for v in y-direction (inner part)
    for I = M + 1:iUGlobalMax-(M)
%         hc = 0.5*(dyLv(I)+dyLv(I-M));
%         hn = dyLv(I);
%         hs = dyLv(I-M);
%         Lyv(I,I) = -(hs+hn)/(hs*hc*hn);
%         Lyv(I,I + M) = 1/(hn*hc);
%         Lyv(I,I - M)= 1/(hs*hc);
    end

    % Compute Laplacian for v in x-direction (horizontal)
    for I = 1:N-1
%         % Left boundary
%         hc = dxq((I-1)*M+1);
%         he = 0.5*( dxq((I-1)*M+1)+dxq(I*M+1+1) );
%         Lxv((I-1)*M + 1 , (I-1)*M + 1 ) = - (2*he+hc)/(hc*hc*he); 
%         Lxv((I-1)*M + 1 , (I-1)*M + 1 + 1 ) = 1/(hc*he);
% 
%         % Right boundary
%         hc = dxq(I*M);
%         hw = 0.5*( dxq(I*M) + dxq(I*M -1) ); 
%         Lxv(I*M , I*M ) = - (2*hw+hc)/(hc*hc*hw); 
%         Lxv(I*M , I*M -1 ) = 1/(hc*hw);
% 
%         for J = (I-1)*M + 2:I*M -1
%             hc = dxq(J);
%             he = 0.5*(dxq(J)+dxq(J+1));
%             hw = 0.5*(dxq(J)+dxq(J-1));
% 
%             Lxv( J , J ) = - (he+hw)/(he*hc*hw);
%             Lxv( J , J+1 ) = 1/(he*hc);
%             Lxv( J , J-1 ) = 1/(hw*hc);
%         end
    end

    % Assemble Laplacian matrix
    Lxu = Lxu + Lyu;
    Lxv = Lxv + Lyv;

    L(1:nu,1:nu) = Lxu;
    L(nu+1:nq,nu+1:nq) = Lxv;
end




