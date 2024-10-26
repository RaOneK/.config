clc; clear; close all
omega = -20;
k = (2*pi)/4;
Amplitude = 0.0;

stepPlot = 500;
M = 60;            % grid size in x-dir
N = 60;            % grid size in y-dir
nt =4;              % timesteps used for storage
timeMaxIter = 1e5;  % number of timesteps
dt = 1e-4;          % timestep size
Re = 8500;           % Reynolds number
U0 = 1;             % characteristic velocity scale
delta = 1e-5;       % chosen non-dimensional positive constant that is suﬃciently small
filename = ['PRESSURE_INDEPENDENT_streamfunction_dim' num2str(sqrt(M*N)) '_it' num2str(nt) '_colored_Re' num2str(Re) '_Amplitude' num2str(Amplitude) '_k' num2str(k) '_omega' num2str(omega) '.gif'];
set(gcf, 'Position', get(0, 'Screensize'));
check_flag = 0;
kx = 1.02;  % grid ratio in x-dir
ky = 1.02;  % grid ratio in y-dir
% kx = 1;
% ky = 1;
lx = 1;     % x dimension scale
ly = 1;     % y dimension scale

nu = (M)*(N);
nv = (M)*(N);
nq = nu + nv;
nn = (M)*(N);

[dx, dy, x, y] = gridGen(M,N,kx,ky,lx,ly);
X = [];
Y = [];
for i = 1:N
    X = [x'; X];
end
for i = 1:M
    Y = [flipud(y) Y];
end
X = X(1:N,1:M);
Y = Y(1:N,1:M);

fprintf('Precomputing matrices: \n'); tic;
fprintf('\tDiv, Grad, Curl: ');
[D,G,C] = divGradCurlGen(M,N);
C = sparse(C);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tdxq dyq: '); tic;
[dxq, dyq] = delQ(M,N,dx,dy);
fprintf('\t\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tInitial conditions: '); tic;
[u,v,q] = initialConditions(M,N,nt,dyq,dxq);
qprev(1:nu) = dyq.*0.99;
qprev(nu+1:nq) = dxq.*(0.0);
% psi = pinv(full(C))*qprev';
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tScaling matrices M, R: '); tic;
[MRinv,R] = scalingMatrices(M,N,dx,dy);
MRinv = sparse(MRinv);
R = sparse(R);
Rinv = R^(-1);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tLHS of system: '); tic
[L] = laplacianGen(M,N,dx,dy);
L = sparse(L);
hatA = speye(nq,nq) - dt/Re * L;
A = C'*MRinv*hatA*C;
fprintf('\t\t\tcompleted in %6.4f sec!\n', toc);

% % Check wether system is symmetric
% symmCheck = round((A-A')*1e12)/1e12;
% spy(symmCheck); % precise up to 14 digits after comma
% % central scheme with ghost/virtual nodes is not



%precompute matrix products. 
fprintf('\tPrecomputing RHS matrices: '); tic
MM = MRinv*R;
CTMM = C'*(MM);
CTMRinv = C'*MRinv;
CTMMhatARinv = CTMM*(hatA)*Rinv;
fprintf('\tcompleted in %6.4f sec!\n', toc);

% Cholesky Factorization
fprintf('LU decomposition: '); tic
[Lfact, Ufact, Pfact, Qfact] = lu(A, 'vector');
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

% initial conditions
fprintf('Boundary conditions: '); tic
uInlet=zeros(N,nt);
uFreestream=zeros(M,nt);
timeIndex = 3;
time = 3*dt;
uInlet(:,1:timeIndex)=1+Amplitude*cos(omega*time);
for t=1:3
    uFreestream(:,t) = 1+Amplitude*cos(omega*t+k*x(2:end));
    duFreestreamdx(:,t) = -k*Amplitude*sin(omega*t + k*0.5*(x(2:end)+x(1:end-1))); % actually its -du/dx
end

% particular solution
% using Pinv
% qp = pinv(D)*boundary2(M, N, dx, dy, q, u, v,dt,1,x,Amplitude,k,omega,uInlet);    % using PINV % do not uncomment as this specifies proper dims for qp
qp = zeros(nq,1);
% uniform horizontal velocity, no vertical
ctr = 1;
uparticular  = zeros(1,nu);
for i = 1:N
    uparticular(ctr:1:ctr+M-1) = uInlet(i,1);
    ctr = ctr + M;
end
qp(1:nu) = dyq.*uparticular;    
qp(nu+1:nq) = dxq.*(0.0);

q(:,1) = qprev;
q(:,2) = qprev;
qprev = qprev';
bc2 = zeros(M*N,nt);
qpadd=qp;
fprintf('\t\t\tcompleted in %6.4f sec!\n', toc);


% DSIPLAY MEMORY USAGE
% Get information about variables in the current workspace
varsInfo = whos;
% Compute total bytes used by all variables
totalBytes = sum([varsInfo.bytes]);
% Convert bytes to appropriate unit
units = {'Bytes', 'KB', 'MB', 'GB', 'TB'};
scale = floor(log(totalBytes)/log(1024)); % Determine the scale
scale = min(scale, length(units)-1);      % Ensure scale does not exceed units
memoryUsage = totalBytes / 1024^scale;    % Convert to the appropriate unit
% Display memory usage
fprintf('Current memory usage: \t\t\t%.2f %s\n', memoryUsage, units{scale+1});

uminGlob = 0.3;
uminGlobArr = uminGlob;
advPrevPrev = advection(M,N,dx,dy,u,v,2,dt,uInlet,uFreestream,duFreestreamdx);
time=time+dt; 
for timeI = 4:timeMaxIter
    timeNow = 4;
    time = time+dt;

%     compute inlet 
    uInlet(:,timeNow) = 1+Amplitude*cos(omega*time);
    uFreestream(:,timeNow) = 1+Amplitude*cos(omega*timeI+k*x(2:end));
    duFreestreamdx(:,timeNow) = -k*Amplitude*sin(omega*timeI + k*x(2:end));

%     uInlet(1,timeNow) = 0.5*uInlet(1,timeNow);

%     recalculate particular solution, qp is particular, qpadd is same but used in computation of B 
    ctr = 1;
    uparticular  = zeros(1,nu);
    for i = 1:N
        uparticular(ctr:1:ctr+M-1) = uInlet(i,timeNow);
        ctr = ctr + M;
    end
    qp(1:nu) = dyq.*uparticular;    
    qp(nu+1:nq) = dxq.*(0.0);
    qpadd=qp;
    

    [LBC] = laplacianRHS(M,N,dx,dy,u,v,q,dt,timeNow,x,Amplitude,k,omega,uFreestream,duFreestreamdx);
    [pbc] = pressure(M,N, dx, dy, q, u, v, dt, timeNow, delta, U0, Re, Amplitude, k, omega,time);

%     % zero laplacian
%     hatA = eye(nq,nq);
%     A = C'*MRinv*hatA*C;
%     LBC = zeros(nq,1);
%     pbc = LBC;

    advPrev = advection(M,N,dx,dy,u,v,timeNow-1,dt,uInlet,uFreestream,duFreestreamdx);
    [ADV] = 1.5*advPrev-0.5*advPrevPrev;
    advPrevPrev = advPrev;

%     [ADV] = zeros(nq,1);
    
    B = dt*CTMM*(1/Re*LBC-ADV-pbc) + CTMRinv*qprev - CTMMhatARinv * qpadd;

    % using precomputed decomposition (as per recommendation)
    % Using precomputed LU factorization with permutation vectors
    B_permuted = B(Pfact);                     % Apply row permutation to B
    intermediate = Lfact \ B_permuted;         % Solve L * intermediate = B_permuted
    psi_permuted = Ufact \ intermediate;       % Solve U * psi_permuted = intermediate
    
    % Apply column permutation to get the final solution
    psi = zeros(size(B));
    psi(Qfact) = psi_permuted;                 % Rearrange psi_permuted into psi according to Qfact

    % store data
    q(:,timeNow) = C*psi+qp;
    qnow = q(:,timeNow);
    qprev = q(:,timeNow);
    u(:,timeNow) = qnow(1:nu)./dyq';
    v(:,timeNow) = qnow(nu+1:nq)./dxq';

    q(:,1:timeNow-1) = q(:,2:timeNow);
    u(:,1:timeNow-1) = u(:,2:timeNow);
    v(:,1:timeNow-1) = v(:,2:timeNow);
    
    uInlet(:,1:timeNow-1) = uInlet(:,2:timeNow);
    uFreestream(:,1:timeNow-1) = uFreestream(:,2:timeNow);
    duFreestreamdx(:,1:timeNow-1) = duFreestreamdx(:,2:timeNow);
    
    if (mod(timeI,stepPlot)==0) || (timeI==4)
        fprintf('\n%5d. ',timeI);
        qh = C*psi;
        uh = qh(1:nu);
        vh = qh(nu+1:nq);
        uploth = reshape(uh,M,N);
        vploth = reshape(vh,M,N);
        up = qp(1:nu)./dyq';
        vp = qp(nu+1:nq)./dxq';
        uplotp = reshape(up,M,N);
        vplotp = reshape(vp,M,N);
        uplot = reshape(u(:,timeNow),M,N);
        vplot = reshape(v(:,timeNow),M,N);
        psiPlot = reshape(psi,M,N);
        
        sumleft = -sum((boundary2(M, N, dx, dy, q, u, v,dt,timeNow,x,Amplitude,k,omega,uInlet)));
        sumright = sum((q(M:M:nu,timeNow)));
        sumtop = sum((q(nq-M+1:1:nq,timeNow)));
        fprintf('massFluxBdry: sum=%+.6e ',sumleft+sumright+sumtop);

        % find min U
        uminNow = min(u(:,timeNow));
        if uminNow < uminGlob
            uminGlob = uminNow;
            uminGlobArr = [uminGlobArr uminGlob];
        end
        fprintf('umin=%8.5f ',uminNow);
        

        % plotting
        tiledlayout(4,6); 
        nexttile(3,[1 1]);
%         contourf(X,flipud(Y),((vploth).^2 + (uploth).^2)',50,'LineColor','none'); hold on; colorbar;
%         sih = streamslice(X,flipud(Y),(uploth'),(vploth')); 
        contourf(((vploth).^2 + (uploth).^2)',50,'LineColor','none'); hold on; colorbar;
        sih = streamslice((uploth'),(vploth')); 
%         sih = quiver(X,flipud(Y),(uploth)',(vploth)',0);
        set(sih,'Color','red');
        colorbar
        title (sprintf('psi contour (homogeneous solution)'))
        
  
        % mass flux
        nexttile(6,[1 1]);
        massSum = zeros(M*N,1);
        for iCell=1:M*N
            col = mod(iCell-1, M) + 1;          % column index 
            row = ceil(iCell/M); 
            ijleft = (row-1)*(M) + col - 1;
            ijright = (row-1)*(M) + col - 1+1;
            ijbot= nu+(row - 2) * M + col;
            ijtop= nu+(row - 2+1) * M + col;
            if col==1 && row ==1
                massSum(iCell) = -uInlet(row,timeNow)*dy(row)+q(ijright,timeNow)+q(ijtop,timeNow);
            elseif col == 1
                massSum(iCell) = -uInlet(row,timeNow)*dy(row)+q(ijright,timeNow)-q(ijbot,timeNow)+q(ijtop,timeNow);
            elseif row == 1
                massSum(iCell) = -q(ijleft,timeNow)+q(ijright,timeNow)+q(ijtop,timeNow);
            else
                massSum(iCell) = -q(ijleft,timeNow)+q(ijright,timeNow)-q(ijbot,timeNow)+q(ijtop,timeNow);
            end
        end
        massPlot = reshape(massSum,M,N);
        contourf(X,flipud(Y),massPlot',10);
%         contourf(X,flipud(Y),((vplotp).^2 + (uplotp).^2)',50,'LineColor','none'); hold on; colorbar;
%         sip = streamslice(X,flipud(Y),uplotp',vplotp');
%         sip = quiver(X,flipud(Y),uplotp',vplotp',1.5);
        axis tight
%         set(sip,'Color','red');
        colorbar
        title(sprintf('sum of mass fluxes in each cell'));
        nexttile(10,[3 2]);

        
        contourf(X,flipud(Y),(uplot.^2+vplot.^2)',50,'LineColor','none'); hold on;
%         VERT = plot(X,Y,'LineWidth',0.0625,"Color", [0.4, 0.4, 0.4, 1]); hold on; 
%         HORIZ = plot(X.',Y.','LineWidth',0.0625,"Color", [0.4, 0.4, 0.4, 1]); hold on; 
        colorbar;
        si = streamslice(X,flipud(Y),uplot',vplot',10);
        set(si,'LineWidth',0.1)
        axis tight
        set(si,'Color','red');
        title(sprintf('total velocity field at %8d',timeI));

        nexttile(1,[2 2]);
        contourf(X,flipud(Y),(uplot)',50,'LineColor','none'); hold on;
        colorbar
        title(sprintf('u velocity contour'));

        nexttile(13,[2 2]);
        contourf(X,flipud(Y),(vplot)',50,'LineColor','none'); hold on;
        colorbar
        title(sprintf('v velocity contour'));
        
        % plot velocity top
        nexttile(4,[1 2]);
        plot(0.5*(x(1:end-1)+x(2:end)),v(nv-M+1:1:nv,timeNow),'red');
        title(sprintf('top velocity BC'));

        % plot bc left
        nexttile(9,[3,1]);
        plot([0; uInlet(:,timeNow)]',[0; 0.5*(y(1:end-1)+y(2:end))]);
        title(sprintf('inlet u component BC=%8.5f',1+Amplitude*sin(omega*time)));

        % plot bc right
        nexttile(12,[3,1]);
        plot([0; u(M:M:nu,timeNow)]',[0; 0.5*(y(1:end-1)+y(2:end))]);
        title(sprintf('outlet u component BC'));

        drawnow;
                    % gif creator
                    frame = getframe(1);
                    im = frame2im(frame);
                    [imind,cm] = rgb2ind(im,256);
                    if check_flag == 0
                        imwrite(imind,cm,filename,'gif','LoopCount',inf,'DelayTime',0.2);
                        check_flag = 1;
                    else
                        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
                    end
    end
end

bc2 = boundary2(M, N, dx, dy, q, u, v,dt,timeNow,x,A,k,omega,uInlet);


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
function [u,v,q] = initialConditions(M,N,nt,dyq,dxq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Changes first columns of u,v,q,psi to initial conditions    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    nu = (M)*(N);
    nv = (M)*(N);
    nq = nu + nv;
    u = zeros(nu,nt);
    v = zeros(nv,nt);
    q = zeros(nq,nt);
    col = 1; 
    row = 1;
    velocity = 1;
    for i = 1:nu
        if col > (M)
            col = 1; 
            row = row+1;
            velocity = 1;
        end
        qp(i) = velocity*dyq(i);
        col = col + 1;
    end
    for i = nu+1:nq
        qp(i) = 0*dxq(i-nu);
    end
    qp = qp';
    q(:,1:3) = [qp qp qp]; % init cond
    u(:,1:3) = q(1:nu,1:3)./[dyq; dyq; dyq]';
    v(:,1:3) = q(nu+1:nu+nv,1:3)./[dxq; dxq; dxq]';
    
end
function [D, G, C] = divGradCurlGen(M, N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %       DIVERGENCE & GRADIENT    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nu = M * N;
    nv = M * N;
    nq = nu + nv;
    nodes = M * N;
    links = nu + nv;
    
    % Initialize index and value arrays for sparse matrix D
    I = [];
    J = [];
    S = [];
    
    % First set: D(i,i) = 1 for i = 1:nodes
    I1 = (1:nodes)';
    J1 = (1:nodes)';
    S1 = ones(nodes, 1);
    
    % Second set: D(ctr, ctr - 1) = -1
    ctr = 0;
    I2 = [];
    J2 = [];
    S2 = [];
    for i = 1:N % rows
        ctr = ctr + 2;
        if ctr - 1 >= 1
            I2 = [I2; ctr];
            J2 = [J2; ctr - 1];
            S2 = [S2; -1];
        end
        for j = 1:M - 2 % columns
            ctr = ctr + 1;
            I2 = [I2; ctr];
            J2 = [J2; ctr - 1];
            S2 = [S2; -1];
        end
    end
    
    % Third set: D(i, M*N + i) = 1
    I3 = (1:nodes)';
    J3 = M * N + (1:nodes)';
    S3 = ones(nodes, 1);
    
    % Fourth set: D(i, M*(N-1) + i) = -1 for i = M+1:nodes
    I4 = ((M + 1):nodes)';
    J4 = M * (N - 1) + I4;
    S4 = -ones(nodes - M, 1);
    
    % Combine all indices and values
    I = [I1; I2; I3; I4];
    J = [J1; J2; J3; J4];
    S = [S1; S2; S3; S4];
    
    % Create sparse matrix D
    D = sparse(I, J, S, nodes, links);
    
    % Compute G as the negative transpose of D
    G = -D';
    
    % Construct sparse matrix C
    D_upper = D(:, nu + 1:nq);
    D_lower = -D(:, 1:nu);
    C = [D_upper; D_lower];
    
end
function [dxq, dyq] = delQ(m,n,dx,dy)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       delQ has two arrays: dyq, dxq.                              %
    %       These are modified arrays of dx, dy to multiply q by.       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nu = m*n;
    nv = m*n;
    dxq = zeros(nv,1)';
    dyq = zeros(nu,1)';
    
    for i = 1:n
        jmin = (i-1)*(m)+1;
        jmax = i* (m);
        
        dyq(jmin:jmax)=dy(i);
    end
    
    for i = 1:n
        jmin = (i-1)*m+1;
        jmax =  i*m;
        
        dxq(jmin:jmax) = dx;
    end
end
function [MRinv, R] = scalingMatrices(M, N, dx, dy)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Generate Scaling Matrices MRinv and R        %
    %                 Using Sparse Matrices                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Ensure dx and dy are column vectors
    dx = dx(:);
    dy = dy(:);

    m = M;
    n = N;
    nu = m * n;
    nv = m * n;
    nq = nu + nv;

    % Extend dx and dy to include ghost points, matching your original code
    dx_ext = [dx; dx(M)]; % dx(M+1) = dx(M)
    dy_ext = [dy; dy(N)]; % dy(N+1) = dy(N)

    % Preallocate arrays for sparse matrix indices and values
    MRinv_indices = (1:nq)';
    MRinv_values = zeros(nq, 1);

    % Compute M1 values (for u components)
    M1_values = zeros(nu, 1);
    for I = 1:nu
        arg = floor((I - 1) / m);
        JJ = arg + 1;
        II = mod(I - 1, m) + 1;
        M1_values(I) = (dx_ext(II) + dx_ext(II + 1)) / 2 / dy_ext(JJ);
    end

    % Compute M2 values (for v components)
    M2_values = zeros(nv, 1);
    for I = 1:nv
        arg = floor((I - 1) / m);
        jj = arg + 2;
        ii = mod(I - 1, m) + 1;
        M2_values(I) = (dy_ext(jj) + dy_ext(jj - 1)) / 2 / dx_ext(ii);
    end

    % Combine M1 and M2 values
    MRinv_values = [M1_values; M2_values];

    % Build MRinv as a sparse diagonal matrix
    MRinv = sparse(MRinv_indices, MRinv_indices, MRinv_values, nq, nq);

    % Compute dxq and dyq using your original delQ function
    [dxq, dyq] = delQ(m, n, dx, dy);

    % Build R as a sparse diagonal matrix
    R_values = [dyq'; dxq'];
    R = sparse(MRinv_indices, MRinv_indices, R_values, nq, nq);
end

% ^ DONE
function [ADV] = advection(M,N,dx,dy,u,v,timeIndex,dt,uInlet,uFreestream,duFreestreamdx) % inlet at left needs to be hardcoded (uf)
    nu = (M)*N;
    nv = M*(N);
    uprev = u(:,timeIndex);
    vprev = v(:,timeIndex);
    nq = nu+nv;
    nc = (M)*(N);
    nn = (M+1)*(N+1);
    ADV    = zeros(nq,1);
    uuLeftVirtual   = zeros(N,1);
    uuRightVirtual  = zeros(N,1);
    vvBotVirtual    = zeros(M,1);
    vvTopVirtual    = zeros(M,1);
    uu = zeros(nc,1);
    vv = zeros(nc,1);
    uv = zeros(nn,1);
    uVerAvg_vHorDiff = zeros(nn,1);
    vVerAvg_vVerDiff = zeros(nc,1);
    uHorAvg_uHorDiff = zeros(nc,1);
    vHorAvg_uVerDiff = zeros(nn,1);
    duudx = zeros(nu,1);
    duvdy = duudx;
    dvudx = zeros(nv,1);
    dvvdy = dvudx;
    % cell centre uu vv precompute
    for iCell = 1:nc
        col = mod(iCell-1, M) + 1;          % column index 
        row = ceil(iCell/M); 
        ijleft = (row-1)*(M) + col - 1;
        ijright = (row-1)*(M) + col - 1+1;
        ijbot= (row - 2) * M + col;
        ijtop= (row - 2+1) * M + col;
%         fprintf('row: %3d, col: %3d, L:%4d, R:%4d, B:%4d, T:%4d\n',row,col,ijleft,ijright,ijbot,ijtop);
        % compute uu 
        if col == 1 % INLET
            uu(iCell) = 0.25* (uInlet(row,timeIndex)+uprev(ijright))^2;
            uHorAvg_uHorDiff(iCell) = 0.25* abs(uInlet(row,timeIndex)+uprev(ijright))*(uprev(ijright)-uInlet(row,timeIndex));
        else
            uu(iCell) = 0.25* (uprev(ijleft)+uprev(ijright))^2; % every cell centre
            uHorAvg_uHorDiff(iCell) = 0.25* (uprev(ijleft)+uprev(ijright))*(uprev(ijright)-uprev(ijleft));
        end
        
        % compute vv
        if row == 1 % BOT
            vv(iCell) = 0.25*(vprev(ijtop)+0)^2;
            vVerAvg_vVerDiff(iCell) = 0.25*(vprev(ijtop)+0)*(vprev(ijtop)-0);
        else
            vv(iCell) = 0.25*(vprev(ijtop)+vprev(ijbot))^2;
            vVerAvg_vVerDiff(iCell) = 0.25*(vprev(ijtop)+vprev(ijbot))*(vprev(ijtop)-vprev(ijbot));
        end
%         fprintf('r:%3d\t c:%3d\t L:%3d\t R:%3d\t B:%3d\t T:%3d\t uu=%18.15f\tvv=%18.15f\n',row,col,ijleft,ijright,ijbot,ijtop,uu(iCell),vv(iCell));
    end
    
    for iuv = 1:(M+1)*(N+1)
        row = floor((iuv-1)/(M+1))+1;
        col = mod(iuv-1,M+1)+1;
        ijtop = (M-1)*(row-1)+col-1;
        ijbot = ijtop - (M-1);
        ijright = (M)*(row-1)+col-M;
        ijleft = ijright -1;
        
%         fprintf('row: %4d, node: %4d\n',row,iuv);
        dy(N+1)=dy(N);
        if (col==1) 
            uv(iuv) = 0;
            uVerAvg_vHorDiff(iuv) = 0;
            vHorAvg_uVerDiff(iuv) = 0;
        elseif (row==1)
            uv(iuv) = 0;
            uVerAvg_vHorDiff(iuv) = 0;
            vHorAvg_uVerDiff(iuv) = 0;
        elseif (col==M+1) && (row~=N+1) % use Dong bc dvtangential/dnormal = 0
            hs = dy(row-1);
            hn = dy(row);
            he = dx(col-1);
            uv(iuv) = vprev(ijleft)*(hs*uprev(ijtop) + hn*uprev(ijbot))/(hn+hs);
            uVerAvg_vHorDiff(iuv) = 0; % from Dong
            vHorAvg_uVerDiff(iuv) = 0; % from Dong
        elseif (row==N+1) && (col~=M+1) % du/dy=0
            hw = dx(col-1); % top bdry case
            he = dx(col);
            uv(iuv) = uFreestream(col-1,timeIndex)*(he*vprev(ijleft) + hw*vprev(ijright))/(he+hw); % freestream at top
            uVerAvg_vHorDiff(iuv) = abs(uFreestream(col-1,timeIndex))*(hw*vprev(ijright)-he*vprev(ijleft))/(he+hw); % gamma
            vHorAvg_uVerDiff(iuv) = abs((hw*vprev(ijright)+he*vprev(ijleft))/(he+hw))*(2*uFreestream(col-1,timeIndex) - 2*uprev(ijbot)); % ug+ub=2uf <=> ug=2uf-ub => ug-ub=2uf-2ub
        elseif (col==M+1) && (row==N+1) % du/dy=0 dv/dx=0
            hs = dy(row-1);
            hn = dy(row);
            uv(iuv) = vprev(ijleft)*uFreestream(col-1,timeIndex);
            uVerAvg_vHorDiff(iuv) = 0; % dv/dx=0 then vdiff=0
            vHorAvg_uVerDiff(iuv) = abs(vprev(ijleft)) * (2*uFreestream(col-1,timeIndex) - 2*uprev(ijbot)); % gamma
        else
            hs = dy(row-1);
            hn = dy(row);
            hw = dx(col-1);
            he = dx(col);
            uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijright))/(he+hw)*(hs*uprev(ijtop) + hn*uprev(ijbot))/(hn+hs);
            uVerAvg_vHorDiff(iuv) = abs(hs*uprev(ijtop) + hn*uprev(ijbot))/(hn+hs)*(hw*vprev(ijright)-he*vprev(ijleft))/(he+hw);
            vHorAvg_uVerDiff(iuv) = abs(hw*vprev(ijright)+he*vprev(ijleft))/(he+hw)*((hs*uprev(ijtop) - hn*uprev(ijbot))/(hn+hs));
        end
%         fprintf('row:%3d, col:%3d, topU:%3d, botU:%3d\t leftV:%3d, rightV=%3d, uv(%3d)=%18.15f\n',row,col,ijtop,ijbot,ijleft,ijright,iuv,uv(iuv));
    end
%     figure;
%     contourf(reshape(uv,M+1,N+1)'); drawnow;
%     figure;
    
    % gamma=1 <=> UDS; gamma=0 <=> CDS.
%     gamma = min(1.2*dt*max(max(abs(u(:,timeIndex))),max(abs(v(:,timeIndex)))),1);
    
    %from MIT paper
%     gamma = min(1.2*dt*max(max(abs(u(:,timeIndex))),max(abs(v(:,timeIndex)))),1);
    % from textbook Chap 3 (eq 3.20) The numerical treatment of the navier-stokes
    % equations
    gamma = max(max(abs(uprev))*dt/min(dx),max(abs(vprev))*dt/min(dy));
    
    % compute x momentum central derivatives
    dy(N+1) = dy(N);
    dx(M+1) = dx(M);
    for iu = 1:nu
        row = floor((iu-1)/(M))+1;
        col = mod(iu-1,M)+1;
        leftcell = (row-1)*(M)+col;
        rightcell = leftcell + 1;
        botnode = (row-1)*(M+1)+col+1;
        topnode = botnode+M+1;
        
        deltax = 0.5*(dx(col)+dx(col+1));
        deltay = dy(row);
        duvdy(iu) = (uv(topnode) - uv(botnode) - gamma*(vHorAvg_uVerDiff(topnode) - vHorAvg_uVerDiff(botnode)))/(deltay);
        if col==M
            uuright = uu(leftcell);
            uuright = uprev(iu)*uprev(iu-1);
            vTprevprev = v(iu,timeIndex-1);
            if row>1
                vBprevprev = v(iu-M,timeIndex-1);
            else
                vBprevprev = 0;
            end
            ug = -2*deltax*(vTprevprev - vBprevprev)/deltay + uprev(iu-1);
            uuright = 0.25*(ug+uprev(iu))^2; 
            uHorAvg_uHorDiff_GHOST = 0.25*abs(ug+uprev(iu))*(ug-uprev(iu)); 
            duudx(iu) = (uuright - uu(leftcell) - 1*(uHorAvg_uHorDiff_GHOST - uHorAvg_uHorDiff(leftcell)))/(deltax);
        else
            duudx(iu) = (uu(rightcell) - uu(leftcell) - gamma * (uHorAvg_uHorDiff(rightcell) - uHorAvg_uHorDiff(leftcell)))/(deltax);
%             fprintf('iu=%3d, row:%3d, col:%3d, leftcell=%3d, rightcell=%3d, botnode=%3d, topnode=%3d, uuright=%8.4f, uuleft=%8.4f\n', ...
%             iu,row,col,leftcell,rightcell,botnode,topnode, uu(rightcell), uu(leftcell));
        end
        
    end
    
    % compute y momentum central derivatives
    
    for iv = 1:nv
        row = floor((iv-1)/(M))+1;
        col = mod(iv-1,M)+1;
        botcell = (row-1)*(M)+col;
        topcell = botcell+M;
        leftnode = (row-1)*(M+1)+col+1+M;
        rightnode = leftnode + 1;
%         fprintf('r:%3d, c:%3d, botcell=%3d, topcell=%3d, leftnode=%3d, rightnode=%3d',row,col,botcell,topcell,leftnode,rightnode);
        deltax = dx(col);
        deltay = 0.5*(dy(row)+dy(row+1));
        dvudx(iv) = (uv(rightnode) - uv(leftnode) - gamma*(uVerAvg_vHorDiff(rightnode) - uVerAvg_vHorDiff(leftnode)))/deltax;
        if row==N
            vvtop = vv(botcell); % fprintf('%18.15f\t',vvtop); % central difference
            % this was wrong: vvtop = vprev(iv)*vprev(iv-M); % fprintf('%18.15f\n',vvtop); % backward upwind difference
            
            % apply bc (vg-vdry)^(n)/dy=-(uR-uL)^n/dx, also uGR-uR=0 and
            % uGL - uL = 0 from du/dy=0 
            % OLD DONG BC below
%             uRprevprev = u(iv,timeIndex-1);
%             if col>1
%                 uLprevprev = u(iv-1,timeIndex-1);
%             else
%                 uLprevprev = uInlet(N,timeIndex-1);
%             end
%             vg = -2*deltay*(uRprevprev - uLprevprev)/deltax + vprev(iv-M);
            vg = -2*deltay*duFreestreamdx(col,timeIndex) + vprev(iv-M);
            vvtop = 0.25*(vg+vprev(iv))^2; 
            vvTopDiff = 0.25*abs(vg+vprev(iv))*(vg-vprev(iv)); 
            dvvdy(iv) = (vvtop - vv(botcell) - gamma*(vvTopDiff - vVerAvg_vVerDiff(botcell)))/deltay;
%             fprintf('\t %3d. vvtop = %18.15f, vvbot = %18.15f, dy=%18.15f, vg = %18.15f\n', iv, vvtop,vv(botcell),deltay,vg);
        else
            dvvdy(iv) = (vv(topcell) - vv(botcell) - gamma*(vVerAvg_vVerDiff(topcell) - vVerAvg_vVerDiff(botcell)))/deltay;
%             fprintf('\t vvtop = %18.15f, vvbot = %18.15f, dy=%18.15f\n', vv(topcell),vv(botcell),deltay);
        end
%         fprintf('\t dvvdy(%3d)=%18.15f\n',iv,dvvdy(iv));
    end

%     flipud(reshape(uv,M+1,N+1)')
%     [flipud(reshape(uu,M,N)') flipud(reshape(vv,M,N)')]
% 
%     [flipud(reshape(duudx,M,N)') flipud(reshape(duvdy,M,N)') ]
%     [flipud(reshape(dvudx,M,N)') flipud(reshape(dvvdy,M,N)') ]

    advU = duudx+duvdy;
    advV = dvudx+dvvdy;
    ADV = [advU; advV];

end
function [L] = laplacianGen(M,N,dx,dy)

    % Initialize dimensions
    nu = (M)*N;
    nv = M*(N);
    nq = nu+nv;
    L = zeros(nq,nq);

    % Initialize submatrices
    Luxx = zeros(nu,nu); % !!!!!!!!!!!!!!!!! Lxu : d^2u/dx^2 !!!!!!!!!!!!!!!!!!
    Luyy = zeros(nu,nu); % !!!!!!!!!!!!!!!!! Lyu : d^2u/dy^2 !!!!!!!!!!!!!!!!!!
    Lvxx = zeros(nv,nv); % !!!!!!!!!!!!!!!!! Lxv : d^2v/dx^2 !!!!!!!!!!!!!!!!!!
    Lvyy = zeros(nv,nv); % !!!!!!!!!!!!!!!!! Lyv : d^2v/dy^2 !!!!!!!!!!!!!!!!!!
    
    % Luxx
    dx(M+1) = dx(M);
    iumax = nu;
    for iu = 1:nu
        row = ceil(iu/(M));
        col = mod(iu-1,M)+1;
%         fprintf("R: %3d C: %3d\n",row,col);
        if col == 1
            hw = dx(1);
            he = dx(2);
            hc = 0.5*(he+hw);
            Luxx(iu,iu) = -2/(hw*he); 
            Luxx(iu,iu+1) = 1/(he*hc);
        elseif col == M
            hw = dx(M-1);
            he = dx(M);
            hc = 0.5*(he+hw);
            Luxx(iu,iu) = -2/(hw*he); %du/dx=-dv/dy
            Luxx(iu,iu-1) = 1/(hw*hc)+1/(he*hc); % west coeff is 1/hw / hc
        else
            he = dx(col+1);
            hw = dx(col);
            hc = 0.5*(he+hw);
            Luxx(iu,iu) = -2/(hw*he);
            Luxx(iu,iu-1) = 1/(hw*hc);
            Luxx(iu,iu+1) =  1/(he*hc);
        end
    end
    % Luyy
    dy(N+1)=dy(N);
    iumax = nu;
    for iu = 1:nu
        row = ceil(iu/(M));
        col = mod(iu-1,M)+1;
%         fprintf("R: %3d C: %3d\n",row,col);
        if row == 1
            hn = 0.5*(dy(row)+dy(row+1));
            hs = dy(row);
            hc = dy(row);
            Luyy(iu,iu) = -2/(hs*hn)-1/(hs*hc);
            Luyy(iu,iu+(M)) = 1/(hn*hc);
        elseif row == N
            hn = dy(row);
            hc = dy(row);
            hs = 0.5*(dy(row)+dy(row-1));
            Luyy(iu,iu) = -2/(hs*hn)-1/(hn*hc); % ug+uc=2uf <=> ug=2uf-uc
            Luyy(iu,iu)=-(2*hs+hc)/(hc*hc*hs); % from code
            Luyy(iu,iu-(M)) = 1/(hs*hc);
        else
            hc = dy(row);
            hs = 0.5*(dy(row)+dy(row-1));
            hn = 0.5*(dy(row)+dy(row+1));
            Luyy(iu,iu+(M)) = 1/(hn*hc);
            Luyy(iu,iu) = -2/(hs*hn);
            Luyy(iu,iu-(M)) = 1/(hs*hc); % south coeff is 1/hs/hc
        end

    end
    
    % Lvxx
    for iv = 1:nv
        row = ceil(iv/(M));
        col = mod(iv-1,M)+1;
%         fprintf("R: %3d C: %3d\n",row,col);
        if col == 1
            he = 0.5*(dx(1)+dx(2));
            hw = dx(1);
            hc = dx(1);
            Lvxx(iv,iv) = -2/(hw*he)-1/(hw*hc);
            Lvxx(iv,iv+1) = 1/(he*hc);
        elseif col == M
            he = dx(M);
            hc = dx(M);
            hw = 0.5*(dx(M)+dx(M-1));
            Lvxx(iv,iv) = -1/(hw*hc); % DONG Lvxx (used to be -2/(hw*he))
            Lvxx(iv,iv-1) = 1/(hw*hc);
        else
            hw = 0.5*(dx(col)+dx(col-1));
            he = 0.5*(dx(col)+dx(col+1));
            hc = dx(col);
            Lvxx(iv,iv-1) = 1/(hw*hc);
            Lvxx(iv,iv+1) = 1/(he*hc);
            Lvxx(iv,iv) = -2/(hw*he);
        end
    end
    
    % Lvyy
    ivmax = nv;
    for iv = 1:ivmax
        row = ceil(iv/(M));
        col = mod(iv-1,M)+1;
%         fprintf("R: %3d C: %3d\n",row,col);
        if row == 1
            hn = dy(row+1);
            hs = dy(row);
            hc = 0.5*(hn+hs);
            Lvyy(iv,iv) = -2/(hn*hs);
            Lvyy(iv,iv+M) = 1/(hc*hn);
        elseif row == N
            hs = dy(row);
            hn = dy(row+1);
            hc = 0.5*(hn+hs);
            Lvyy(iv,iv) = -2/(hn*hs); % dv/dy=0 vg = vc (at bdry)
            Lvyy(iv,iv-M) = 1/(hc*hs)+1/(hn*hc);


        else
            hs = dy(row);
            hn = dy(row+1);
            hc = 0.5*(hn+hs);
            Lvyy(iv,iv) = -2/(hn*hs);
            Lvyy(iv,iv-M) = 1/(hc*hs);
            Lvyy(iv,iv+M) = 1/(hc*hn);
        end
%         fprintf("R: %3d C: %3d\n",row,col);
    end
    
    % clear row/col from known values
    
    Luxx = Luxx + Luyy;
    Lvxx = Lvxx + Lvyy;
    
    L(1:nu,1:nu) = Luxx;
    L(nu+1:nq,nu+1:nq) = Lvxx;
end
function [LBC] = laplacianRHS(M,N,dx,dy,u,v,q,dt,timeIndex,x,A,k,omega,uFreestream,duFreestreamdx)
    nu = (M)*N;
    nv = M*(N);
    nq = nu+nv;
    LBC = zeros(nq,1);
    for iu = 1:M % bot bc
        hs = dy(1);
        hn = 0.5*(dy(1)+dy(2));
        hc = dy(1);
        LBC(iu) = LBC(iu) + 0;
%         fprintf('LBC(%3d)=%18.15f\n',iu,LBC(iu));
    end
    x = 1.0;
    for iu = nu:-1:nu-(M)+1 % TOP v component
        % tangential component
        idx = iu-nu+(M);
        x = x - dx(idx);
%         fprintf('%8.5f %3d %3d\n',x,idx,iu);
        hn = dy(N);
        hs = 0.5*(dy(N)+dy(N-1));
        hc = dy(N);
        LBC(iu) = LBC(iu) + 2*uFreestream(idx,timeIndex)/(hn*hc);

        % normal component
        hc = dy(N);
        LBC(iu+nu) = LBC(iu+nu) + duFreestreamdx(idx,timeIndex) *(hn+hs)/(hn*hc);
    end

    for iu=M:M:nu % right bc u
        vtprev = v(iu,timeIndex-1);
        
        if iu==M
            vbprev=0;
        else
            vbprev = v(iu-M,timeIndex-1);
        end
        hc = dx(M);
        he = dx(M);
        hw = dx(M);
        row = floor(iu/M);
        deltay = dy(row);
        LBC(iu) = LBC(iu) - 2*hw*(vtprev - vbprev)/deltay/(he*hc);
    end


    
end
% Add this function to compute bc2
function [bc2] = boundary2(M, N, dx, dy, q, u, v,dt,timeIndex,x,A,k,omega, uInlet)
    nu = (M-1)*N;
    nv = M*(N-1);
    bc2 = zeros(M*N,1);
    
    row = 1;
    cell = 1;
    for i = 1:M:M*N
        bc2(cell) = uInlet(row,timeIndex)*dy(row);
        cell = cell+M;
        row = row + 1;
    end
%     fprintf('SUM BC2 %18.15f\n',sum);
end
function [pRHS] = pressure(M,N, dx, dy, q, u, v, dt, timeIndex, delta, U0, Re, A, k, w,t) % add uf freestream for left inlet
    
    nu = (M)*N;
    nv = M*(N);
    nq = nu+nv;
    vprev = v(:,timeIndex-1);   % prev timestep v
    uprev = u(:,timeIndex-1);   % prev timestep u
    pRHS  = zeros(nq,1);
    
    iUstart = M;
    iVNW = M; % indexation in vprev starts from 1
    for iU = iUstart:M:nu
        
        iVSW = iVNW-M;
        iVNWW = iVNW - 1;
        iVSWW = iVSW - 1;
        hc = dx(M);
        hw = 0.5*(dx(M-1)+hc);
        he = 0.5*hc; % v located at bdry
        if iVSW==0
            vInterpS = 0.0;
        else
            vInterpS = (hw+he)*(vprev(iVSW) - vprev(iVSWW))/hw + vprev(iVSWW);
        end
        vInterpN = (hw+he)*(vprev(iVNW) - vprev(iVNWW))/hw + vprev(iVNWW);
        vInterp = 0.5*(vInterpS + vInterpN);
        kineticEnergy = 0.5 * (uprev(iU)*uprev(iU)+vInterp*vInterp);
        Redudx = (1/Re)*(uprev(iU)-uprev(iU-1))/dx(M);
        stepfun = 0.5 * (1-tanh(uprev(iU)/(U0 * delta)));
        pRHS(iU) = Redudx - (kineticEnergy * stepfun);
%         fprintf(' pRHS(%3d)=%18.15f; VNW: %3d; VNWW: %3d; VSW: %3d; VSWW: %3d\n', iU,pRHS(iU), iVNW,iVNWW, iVSW,iVSWW);
        iVNW = iVNW + M;
    end
    
    % New code for the top boundary % take uf -> x-mom -> integrate over dx
    x_coord = 0;
    ctr = 1;
    for iV = nv - M + 1 : nv  % Indices for the top boundary in v
        % Determine grid indices
        row = ceil(iV / M);   % Should be N
        col = mod(iV - 1, M) + 1;
        x_coord = x_coord + dx(ctr)*0.5;
        pRHS(nu + iV) = -(2*A*(Re*k + Re*w) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2)^3 + 2*A*tan((k*x_coord)/2 + (t*w)/2)^2*(Re*k + Re*w - A*Re*k) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2))/(Re*k*(tan((k*x_coord)/2 + (t*w)/2)^2 + 1)^2);
        x_coord = x_coord + dx(ctr)*0.5;
        ctr = ctr+1;
    end

    % extrapolate two pressure values from right boundary to top
    deltay1=(dy(N)+dy(M))*0.5;
    deltay2 = dy(N);
    pCornerTopFromRightBC = deltay2/deltay1 * (pRHS(nu) - pRHS(nu-M)) + pRHS(nu);
    pCornerRightFromTopBC = -(2*A*(Re*k + Re*w) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2)^3 + 2*A*tan((k*x_coord)/2 + (t*w)/2)^2*(Re*k + Re*w - A*Re*k) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2))/(Re*k*(tan((k*x_coord)/2 + (t*w)/2)^2 + 1)^2);
    deltap = pCornerTopFromRightBC - pCornerRightFromTopBC;
    % shift the top values to match the right (solve up to a constant) 
    pRHS(nq-M+1:1:nq) = pRHS(nq-M+1:1:nq) + deltap;
end