clc; clear; close all;

M = 40;      % grid size in x-dir
N = 40;      % grid size in y-dir
nt = 90000;     % number of timesteps
dt = 1e-5;  % timestep size
Re = 500;   % Reynolds number
U0 = 1;     % characteristic velocity scale
delta = 1e-5; % chosen non-dimensional positive constant that is suﬃciently small
filename = ['streamfunction_dim' num2str(sqrt(M*N)) '_it' num2str(nt) '_colored.gif'];

kx = 1.05;  % grid ratio in x-dir
ky = 0.95;  % grid ratio in y-dir
kx = 1;
ky = 1;
lx = 1;     % x dimension scale
ly = 1;     % y dimension scale

nu = (M+1)*(N);
nv = (M)*(N+1);
nq = nu + nv;

% function names are self-explanatory, comments are in corresponding
% functions descriptions
fprintf("Generating grid:");
[dx, dy, x, y] = gridGen(M,N,kx,ky,lx,ly);
fprintf(" done\n");
fprintf("delQ:");
[dxq, dyq] = delQ(M,N,dx,dy);
fprintf(" done\n");
fprintf("Initial Conditions:");
[u,v,q,psi] = initialConditions(M,N,nt,dxq,dyq);
fprintf(" done\n");
fprintf("Generating G D C:");
[D,G,C] = divGradCurlGen(M,N);
fprintf(" done\n");
fprintf("Generating L:");
[L] = laplacianGen(M,N,dx,dy,dxq,dyq);
fprintf(" done\n");
fprintf("Generating M R:");
[Mass,R] = scalingMatrices(M,N,dx,dy,dxq,dyq);
fprintf(" done\n");
fprintf("Computing inverse of M, R:");
Mass = sparse(Mass);
R = sparse(R);
Rinv = R^(-1);
Minv = Mass^(-1);
fprintf(" done\n");
% timeIndex = 2;
sparseC = sparse(C);
% sparseRinv = sparse(Rinv);
fprintf("Computing LHS of system:");
% LHS = ;
fprintf(" done\n");
fprintf("Converting LHS to sparse format:");
LAP = sparse(laplacianGen(M,N,dx,dy,dxq,dyq));
hatA = sparse(eye(nq,nq) - (dt/Re) * 0.5*LAP); % I - dt/(2Re)*hatL
A = sparseC'*M*hatA*sparseC;
fprintf(" done\n");
flagDiverged = 0;
fprintf("Initializing time loop: ");
stepPlot = 100;
for timeIndex = 3:nt
    [uLeftBC,uLeftBCVirtual] = LeftBC(M,N,timeIndex,dt);
    [rnadvect] = advection(M,N,dx,dy,u,v,dxq,dyq,nt,dt);
    rnADV = (1.5*advection(M,N,dx,dy,u,v,dxq,dyq,timeIndex-1,dt) - 0.5*advection(M,N,dx,dy,u,v,dxq,dyq,timeIndex-2,dt));
    qn = Rinv*([u(:,timeIndex-1);v(:,timeIndex-1)]);
    explicitLaplacian = (1/Re)*0.5*LAP*qn;
    [LbcLHS] = (1/Re) *0.5* ...
        (laplacianRHS(M,N,dx,dy,dxq,dyq,u,v,uLeftBC,timeIndex) ...
        +laplacianRHS(M,N,dx,dy,dxq,dyq,u,v,uLeftBC,timeIndex-1));
    [pRHS] = pressureBC(M,N,dx,dy,dxq,dyq,u,v,timeIndex,Re,U0,delta);
    B = dt*sparseC'*(Mass*(qn/dt + explicitLaplacian-rnADV-LbcLHS)-pRHS);
%     denseB = C'*(Mass*bcLaplacian-Mass*rnADV-pRHS);
    psi = mldivide(A,B);
%     psi = linsolve(denseA,denseB);
    uvnew = Rinv*sparseC*psi;
    u(:,timeIndex) = uvnew(1:nu);
    v(:,timeIndex) = uvnew(nu+1:nq);
    fprintf('\n%5d. max u = %7.5f\t min u = %7.5f\t max v = %7.5f\t min v = %7.5f',timeIndex,max((u(:,timeIndex))),min((u(:,timeIndex))),max((v(:,timeIndex))),min((v(:,timeIndex))));
    if (mod(timeIndex,stepPlot)==0) || (timeIndex==3)
        subplot(1,2,1);
        title('Streamfunction contour');
        psiPlot = reshape(psi,M+1,N+1);
        contour(flipud(psiPlot'), 'LineColor','k', 'LevelStep',0.05); 
        subplot(1,2,2);
        title('minmax u minmax v');
        plot(timeIndex,max((u(:,timeIndex))),'or'); hold on;
        plot(timeIndex,min((u(:,timeIndex))),'xr'); hold on;
        plot(timeIndex,max((v(:,timeIndex))),'og'); hold on;
        plot(timeIndex,min((v(:,timeIndex))),'xg'); hold on;
%         subplot(3,1,3)
%         quiver(X,Y,u(:,timeIndex),v(:,timeIndex));
        drawnow;
            % gif creator
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if timeIndex == 3
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
    end
    if (max(abs(u(:,timeIndex)))) > 100
        flagDiverged = 1;
        break;
    end
    if (max(abs(v(:,timeIndex)))) > 100
        flagDiverged = 1;
        break;
    end
    
end
if flagDiverged == 1
    fprintf("\n solutions diverged\n");
else
    fprintf("\n done\n");
end
% spy(round((A-A')*1e14)/1e14)
% size(LHS)
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
    v(:,1) = 0;  % slightly positive uniform vertical v component
    
    q(:,1) = [u(:,1);v(:,1)].*[dyq; dxq]; %elementwise .* multiplication
    
    u(:,2) = u(:,1);
    v(:,2) = v(:,1);
%     for i = 1:(M+1)*N % TEST LOOP TO VERIFY ALGEBRA
%         u(i) = i/((M+1)*N);
%     end
% 
%     for i = 1:(N+1)*M % TEST LOOP TO VERIFY ALGEBRA
%         v(i) = i/((N+1)*M);
%     end
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
function [uLeftBC,uLeftBCVirtual] = LeftBC(M,N,timeIndex,dt)      % BC values at the left bdry
    uLeftBC = zeros((M+1)*N,1);             % only few will be filled, rest are 0
    uLeftBCVirtual = zeros(N,1);            % single column non-zero only at left bc
    for i = 1:N                             % go over all rows in grid
        iGlobal = (i-1)*(M+1)+1;            % get global index
        uLeftBC(iGlobal) = 0.2;             % fill bc value from left bdry
        uLeftBCVirtual(i) = 0.1;            % value outside the boundary at the next u
%         fprintf("%d\n",iGlobal);            % global index test print
    end
end
function [vBotBC] = BotBC(M,N,nt,dt)        % BC values at the bot bdry
    vBotBC = zeros((N+1)*M,1);              % only few will be filled, rest are 0
    for i = 1:M                             % go over all cols in grid
        iGlobal = i;                        % get global index
        vBotBC(iGlobal) = 0;                % fill bc value from bot bdry
%         fprintf("%d\n",iGlobal);            % global index test print
    end
end
% function [uRightBC,vRightBC] = RightBC(M,N,u,v)
%     
% end
% RIGHT AND TOP OPEN BC MATRIX COEFFS ARE HARDCODED INTO CORRESPONDING ADVECT/LAPLACIAN
% RHS LAPLACIAN EXPLICIT TERMS ARE IN laplacianRHS function
function [rnadvect] = advection(M,N,dx,dy,u,v,dxq,dyq,timeIndex,dt)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           ADVECTION           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nu = (M+1)*N;
    nv = M*(N+1);
    nq = nu+nv;
    nc = (M)*(N);
    nn = (M+1)*(N+1);
    rnadvect = zeros(nq,1);
    uu = zeros(nc,1);
    vv = zeros(nc,1);
    uv = zeros(nn,1);
    [uLeftBC,uLeftBCVirtual] = LeftBC(M,N,timeIndex,dt);
    [vBotBC] =  BotBC(M,N,timeIndex,dt);
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

    
    % Compute Laplacian (Luyy) for u in y-direction (vertical)
    iUGlobalMax = (M+1)*N;
    for iUGlobal = 1:iUGlobalMax
        row = ceil(iUGlobal/(M+1));
        col = mod(iUGlobal-1,M+1)+1;
        hc = dyq(iUGlobal);
%         fprintf('\nLrow: %3d',iUGlobal);
        if row==1 % only row affects Luyy because of derivative direction 
            hn = 0.5*(dyq(iUGlobal)+dyq(iUGlobal+(M+1)));
            hs = hc;
            Lyu(iUGlobal,iUGlobal) = -2/(hs*hn)-1/(hs*hc); % main diag
            Lyu(iUGlobal,iUGlobal+(M+1)) = 1/(hn*hc); % top neighbour
            % bot neighbour is across the wall, hence interpolation
        elseif row==(N) % top row, affected by Dong's BC
% TO BE FILLED WITH DONG BC
            hn = hc;
            hs = 0.5*(dyq(iUGlobal)+dyq(iUGlobal-(M+1)));
            % no top neighbour, central element is uppermost unknown
            Lyu(iUGlobal,iUGlobal) = (-2*hc+hn)/(hs*hn); % main diag
            Lyu(iUGlobal,iUGlobal-(M+1)) = 1/(hs*hc); % bot neighbour
            % fprintf('iUGlobal:%3d\t row: %3d\t col: %3d\n',iUGlobal,row, col);
        else % inner part
            hn = 0.5*(dyq(iUGlobal)+dyq(iUGlobal+(M+1)));
            hs = 0.5*(dyq(iUGlobal)+dyq(iUGlobal-(M+1)));
            Lyu(iUGlobal,iUGlobal) = -2/(hs*hn); % main diag
            Lyu(iUGlobal,iUGlobal+(M+1)) = 1/(hn*hc); % top neighbour
            Lyu(iUGlobal,iUGlobal-(M+1)) = 1/(hs*hc); % bot neighbour
        end
    end

    % Compute Laplacian (Luxx) for u in x-direction (horizontal)
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
            % No left element since it is at the boundary. 
            % Corresponding columns and rows will be removed later.
            % Explicit RHS is in separate module generating Laplacian BC.
            Lxu(iUGlobal,iUGlobal) = -2/(hw*he);    % main diag
            Lxu(iUGlobal,iUGlobal+1) = 1/(hc*he);   % right
            idxqRight = idxqRight + 1; % increment dx index
            % fprintf('iUGlobal:%3d\t row: %3d\t col: %3d\t idxRight: %3d\t idxLeft: %3d \n',iUGlobal,row, col,idxqRight,idxqLeft);
        elseif col == (M+1) % rightmost boundary DONG BC
% TO BE FILLED WITH DONG BC% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG 
% FILLED DONG BC
% p.27 first equation check. CHECKED
            hw = dxq(idxqLeft);
            hc = 0.5*(he+hw);
            he = hc; 
            Lxu(iUGlobal,iUGlobal-1) = 1/(hw*hc);   % left
            Lxu(iUGlobal,iUGlobal) = 1/(hw*hc);    % main diag
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


    hn=0;
    hs=0;
    hc=0;
    % Compute Laplacian (Lvyy) for v in y-direction (vertical)
    iVGlobalMax = (N+1)*M;
    for iVGlobal = 1:iVGlobalMax
        row = ceil(iVGlobal/(M));
        col = mod(iVGlobal-1,M)+1;
        if row == 1 % bottommost boundary part
            % no need to fill since it will be removed later. 
            % part kept for dimension consistency with other operators: i.e
            % curl div grad 
        elseif row == 2  % bottommost inner part
            hn = dy(row);
            hs = dy(row-1);
            hc = 0.5*(hn+hs);
            % south neighbour at wall is zero 
            % Lyv(iVGlobal,iVGlobal-(M)) = 1/(hs*hc); % south
            Lyv(iVGlobal,iVGlobal) = -2/(hn*hs);    % center
            Lyv(iVGlobal,iVGlobal+(M)) = 1/(hc*hn); % north
            
        elseif row == N+1 % DONG BC TOP
% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG % DONG% DONG 
% DONG BC FILLED
% p.27 first equation (change w->s, e->n, u->v)
            hs = dy(row-1);
            hc = 0.5*(hn+hs);
            hn = hc; 
            Lyv(iVGlobal,iVGlobal-(M)) = 1/(hs*hc); % south
            Lyv(iVGlobal,iVGlobal) = 1/(hc*hs);    % center
%             Lyv(iVGlobal,iVGlobal+(M)) = 1/(hc*hn); % north
        else
            hn = dy(row);
            hs = dy(row-1);
            hc = 0.5*(hn+hs);
            Lyv(iVGlobal,iVGlobal-(M)) = 1/(hs*hc); % south
            Lyv(iVGlobal,iVGlobal) = -2/(hn*hs);    % center
            Lyv(iVGlobal,iVGlobal+(M)) = 1/(hc*hn); % north
        end
%         fprintf('i:%3d row:%3d col:%3d hn: %18.15f hc: %18.15f hs: %18.15f\n',iVGlobal,row,col,hn,hc,hs);
    end
    % Compute Laplacian (Lvxx) for v in x-direction (horizontal)
    for iVGlobal = 1:iVGlobalMax
        row = ceil(iVGlobal/(M));
        col = mod(iVGlobal-1,M)+1;
        
        if col==1
            hc = dxq(iVGlobal);
            hw = hc;
            he = 0.5*(dxq(iVGlobal+1)+hc);
            % no left element, since v at (and outside) left bdry = 0
%             Lxv(iVGlobal,iVGlobal-1) = 1/(hw*hc); % left
            Lxv(iVGlobal,iVGlobal) = -2/(hw*he);    % diag
            Lxv(iVGlobal,iVGlobal+1) = 1/(hc*he);   % right
        elseif col==M
% DONG BC FILL
% FIX TYPO IN NOTES, hw changed to he in main diag element. Reasoning: hw+he=2hc => he = 2hc-hw => -2hc+hw = -he
            hc = dxq(iVGlobal);
            hw = 0.5*(dxq(iVGlobal-1)+hc);
            he = hc;
            Lxv(iVGlobal,iVGlobal-1) = 1/(hw*hc);   % left
            Lxv(iVGlobal,iVGlobal) = (-2*hc+hw)/(hw*he);    % diag
        else
            hc = dxq(iVGlobal);
            he = 0.5*(dxq(iVGlobal+1)+hc);
            hw = 0.5*(dxq(iVGlobal-1)+hc);
            Lxv(iVGlobal,iVGlobal-1) = 1/(hw*hc);   % left
            Lxv(iVGlobal,iVGlobal) = -2/(hw*he);    % diag
            Lxv(iVGlobal,iVGlobal+1) = 1/(hc*he);   % right
        end
        
%         fprintf('i:%3d row:%3d col:%3d hw: %18.15f hc: %18.15f hw: %18.15f\n',iVGlobal,row,col,hw,hc,he);
    end
    
    % Assemble Laplacian matrix
    Lxu = Lxu + Lyu;
    Lxv = Lxv + Lyv;
    
    
    L(1:nu,1:nu) = Lxu;
    L(nu+1:nq,nu+1:nq) = Lxv;
    % remove left boundary velocities from matrix
    for iUGlobal = 1:iUGlobalMax    
        row = ceil(iUGlobal/(M+1));
        col = mod(iUGlobal-1,M+1)+1;
        if col == 1 
            Lxu(iUGlobal,:)=0;
            Lxu(:,iUGlobal)=0;
        end
    end
    % remove bottom boundary velocities from matrix
    for iVGlobal = 1:iVGlobalMax
        row = ceil(iVGlobal/(M));
        col = mod(iVGlobal-1,M)+1;
%         fprintf('\nrrrow %3d',row);
        if row == 1
            Lxv(iVGlobal,:) = 0;
            Lxv(:,iVGlobal) = 0;
        end
    end
    L(1:nu,1:nu) = Lxu;
    L(nu+1:nq,nu+1:nq) = Lxv;
end
function [LbcLHS] = laplacianRHS(M,N,dx,dy,dxq,dyq,u,v,uLeftBC,timeIndex)
    nu = (M+1)*N;
    nv = M*(N+1);
    nq = nu+nv;
    LbcLHS = zeros(nq,1);
    vprev = v(:,timeIndex-1);   % prev timestep v
    uprev = u(:,timeIndex-1);   % prev timestep u
    % left BC u
    iUGlobal = 2;
    for iLeftU = 1:N
        hw = dx(1);
        he = dx(2);
        hc = 0.5*(he+hw);
        LbcLHS(iUGlobal) = uLeftBC(iUGlobal-1) / (hc*hw);
%         fprintf('LbcLHS(%3d): %18.15f\n',iUGlobal,LbcLHS(iUGlobal));
        iUGlobal = iUGlobal + M+1;
    end
    
    iUGlobal = M+1;
    iBotV = M;
    for iRightU = 1:N
        deltay = dy(iRightU);
        iTopV = iBotV + M;
        iBotVLeftNeighbour = iBotV - 1;
        iTopVLeftNeighbour = iTopV - 1;
        hw = 0.5*(dx(M)+dx(M-1));
        he = dx(M);
        hc = 0.5*(he+hw);
        vTop = (0.5*he + hw)*(vprev(iTopV) - vprev(iTopVLeftNeighbour))/hw + vprev(iTopVLeftNeighbour);
        vBot = (0.5*he + hw)*(vprev(iBotV) - vprev(iBotVLeftNeighbour))/hw + vprev(iBotVLeftNeighbour);
        LbcLHS(iUGlobal) = - deltay * (vTop - vBot) / hc;
%         fprintf('LbcLHS(%3d): %18.15f\n',iUGlobal,LbcLHS(iUGlobal));
        % shift indices
        iUGlobal = iUGlobal + M+1;
        iBotV = iTopV;
    end


    iVGlobal = nq-M+1;
    iLeftU = (M+1)*(N-1)+1;
    for iCol = 1:M
        deltax = dx(iCol);
        iRightU = iLeftU + 1;
        iLeftUBotNeighbour = iLeftU - (M+1);
        iRightUBotNeighbour = iRightU - (M+1);
        hs = 0.5*(dy(N)+dy(N-1));
        hn = dy(N);
        hc = 0.5*(hn+hs);
        uLeft = (0.5*hn + hs) * (uprev(iLeftU) - uprev(iLeftUBotNeighbour))/hs + uprev(iLeftUBotNeighbour);
        uRight = (0.5*hn + hs) * (uprev(iRightU) - uprev(iRightUBotNeighbour))/hs + uprev(iRightUBotNeighbour);
        LbcLHS(iVGlobal) = - deltax * (uRight - uLeft) / hc;
%         fprintf('LbcLHS(%3d): %18.15f\n',iVGlobal,LbcLHS(iVGlobal));
        % shift indices
        iVGlobal = iVGlobal + 1;
        iLeftU = iRightU;
    end
end
function [pRHS] = pressureBC(M,N,dx,dy,dxq,dyq,uInterpW,v,timeIndex,Re,U0,delta)
    nu = (M+1)*N;
    nv = M*(N+1);
    nq = nu+nv;
    vprev = v(:,timeIndex-1);   % prev timestep v
    uprev = uInterpW(:,timeIndex-1);   % prev timestep u
    pRHS  = zeros(nq,1);
    
    iUstart = M+1;
    iVstart = nv-(M)+1;
    iUSW = nu-(M+1)+1;
    iVSW = M; % indexation in vprev starts from 1
    for iU = iUstart:M+1:nu
        
        iVNW = iVSW + M;
        iVNWW = iVNW - 1;
        iVSWW = iVSW - 1;
        hc = dx(M);
        hw = 0.5*(dx(M-1)+hc);
        he = 0.5*hc; % v located at bdry
        vInterpS = (hw+he)*(vprev(iVSW) - vprev(iVSWW))/hw + vprev(iVSWW);
        vInterpN = (hw+he)*(vprev(iVNW) - vprev(iVNWW))/hw + vprev(iVNWW);
        vInterp = 0.5*(vInterpS + vInterpN);
        kineticEnergy = 0.5 * (uprev(iU)*uprev(iU)+vInterp*vInterp);
        Redudx = Re*(uprev(iU)-uprev(iU-1))/dx(M);
        stepfun = 0.5 * (1-tanh(uprev(iU)/(U0 * delta)));
        pRHS(iU) = Redudx - (kineticEnergy * stepfun);
%         fprintf('vExtrap: %18.15f; pRHS(%3d)=%18.15f; VNE: %3d; VNEE: %3d; VSE: %3d; VSEE: %3d\n', vExtrap, iU,pRHS(iU), iVNW,iVNWW, iVSW,iVSWW);
        iVSW = iVSW + M;
    end

    for iV = iVstart:nv
        iUSSW = iUSW - (M+1);
        iUSE = iUSW + 1;
        iUSSE = iUSE - (M+1);
        hn = 0.5*dy(N);
        hc = dy(N);
        hs = 0.5*(hc + dy(N-1));
        uInterpW = (hs+hn)*(uprev(iUSW)-uprev(iUSSW))/hs + uprev(iUSSW);
        uInterpE = (hs+hn)*(uprev(iUSE)-uprev(iUSSE))/hs + uprev(iUSSE);
        uInterp = 0.5*(uInterpW + uInterpE);
        kineticEnergy = 0.5 * (vprev(iV)*vprev(iV) + uInterp*uInterp);
        Redvdy = Re * (vprev(iV) - vprev(iV - M))/dy(N);
        stepfun = 0.5 * (1 - tanh(vprev(iV)/(U0*delta)));
        pRHS(iV+nu) = Redvdy - (kineticEnergy * stepfun);
%         fprintf('%3d %3d %3d %3d | iV: %3d | Redvdy: %18.15f | velocitySquare: %18.15f | stepfun: %18.15f \n', iUSW,iUSSW,iUSE,iUSSE, iV,Redvdy,kineticEnergy,stepfun);
        iUSW = iUSW + 1;
    end
end
function [Mass,R] = scalingMatrices(M,N,dx,dy,dxq,dyq)
    R = diag([dyq; dxq]);
    Mass = zeros((M+1)*N+M*(N+1),1);
    rowShift = 0;
    for row = 1:N
        for col = 1:M+1
            if col == 1
                Mass(col+rowShift) = dx(1);
            elseif col == M+1
                Mass(col+rowShift) = dx(M);
            else
                Mass(col+rowShift) = 0.5*(dx(col)+dx(col-1));
            end
%             fprintf('%3d. M=%18.15f\n',col+rowShift,Mass(col+rowShift))
        end
        rowShift = rowShift + (M+1);
    end
    rowShift = (M+1)*N;
    for row = 1:N+1
        for col = 1:M
            if row == 1
                Mass(rowShift+col) = dy(row);
            elseif row == N+1
                Mass(rowShift+col) = dy(row-1);
            else
                Mass(rowShift+col) = 0.5*(dy(row)+dy(row-1));
            end
%             fprintf('%3d. M=%18.15f\n',col+rowShift,Mass(rowShift+col))
        end
        rowShift = rowShift + (M);
    end
    Mass = diag(Mass);
end



