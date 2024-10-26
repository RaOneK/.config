clc; clear; close all
omega = -20;
k = 10;
Amplitude = -0.99;
M = 40;            % grid size in x-dir
N = 40;            % grid size in y-dir
smootheningFactor = 1e9; % bottom left corner, higher -> sharper

stepPlot = 500;

nt =4;              % timesteps used for storage
timeMaxIter = 1e4;  % number of timesteps
dt = 5e-5;          % desired timestep size
Re = 2500;           % Reynolds number
U0 = 1;             % characteristic velocity scale
delta = 1e-6;       % chosen non-dimensional positive constant that is suﬃciently small
kx = 1.01;  % grid ratio in x-dir
ky = 1.01;  % grid ratio in y-dir
kx = 1;
ky = 1;
lx = 1.0;     % x dimension scale
ly = 1.0;     % y dimension scale

nu = (M)*(N);
nv = (M)*(N);
nq = nu + nv;
nn = (M)*(N);
% generate grid
[dx, dy, x, y] = gridGen(M,N,kx,ky,lx,ly);

% Assuming dx and dy are arrays representing grid spacing in x and y directions
dx_min = min(dx);
dy_min = min(dy);

% % Calculate time step based on CFL condition (AB2 scheme)
% U_max = 1; % Maximum velocity in x-direction
% V_max = 1; % Maximum velocity in y-direction
% dt_cfl = 2 * min(dx_min, dy_min) / max(U_max, V_max);
% 
% % Optimal time step for the combined AB2 (advection) and Crank-Nicolson (diffusion) scheme
% safetyFactor = 0.01; % Safety factor to reduce dt
% dt = safetyFactor * dt_cfl;
% 
% fprintf('\nnew dt=%+.3e\n',dt)

% check for Diffusion and Courant number stability
D = dt/min(min(dx),min(dy))^2/Re;
C = 1.2*dt/min(min(dx),min(dy));
sqrt2d=sqrt(2*D);
fprintf('Diffusion number D = %8.3f ',D);
if D>0.5
    fprintf('> 0.5, must be <0.5 for explicit schemes (current scheme is implicit)\n');
%     return
end
fprintf('\nCourant Number C = %6.3f ',C);
if  C>sqrt2d
    fprintf('> %8.3f, must be <sqrt(2D)\n',sqrt2d);
%     return
end
if  C>1
    dtold = dt;
    dt = min(min(dx),min(dy))/5;
    D = dt/min(min(dx),min(dy))^2/Re;
    C = 1.2*dt/2/min(min(dx),min(dy));
    fprintf('\n\tNew diffusion number D = %8.3f ',D);
    fprintf('\n\tNew Courant number C = %8.3f ',C);
    fprintf('\n\tchanging time step from dt=%+6.4e to dt=%+6.4e to satisfy stability\n',dtold, dt);
%     return
end

aspectRatio = max(max(dy)/min(dx),max(dx)/min(dy));
fprintf('\nMaximum aspect ratio dx/dy=%6.4f (or dy/dx)\n',aspectRatio);
if aspectRatio>2000
    fprintf('Aspect ratio dx/dy=%6.4f\n is too high',aspectRatio);
    return
end

% create file
filename = ['PRESSURE_INDEPENDENT_streamfunction_dim' num2str(sqrt(M*N)) '_it' num2str(nt) '_colored_Re' num2str(Re) '_Amplitude' num2str(Amplitude) '_k' num2str(k) '_omega' num2str(omega) '.gif'];
set(gcf, 'Position', get(0, 'Screensize'));
check_flag = 0;

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

fprintf('\nPrecomputing matrices: \n'); tic;
fprintf('\tDiv, Grad, Curl: ');
[D,G,C] = divGradCurlGen(M,N);
% C = sparse(C);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tdxq dyq: '); tic;
[dxq, dyq] = delQ(M,N,dx,dy);
fprintf('\t\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tInitial conditions: '); tic;
[u,v,q] = initialConditions(M,N,nt,dyq,dxq,omega,Amplitude,k,ly,y);
% plot(u,'x'); hold on; 
% plot(v,'o'); hold off;
% plot(q);
% return
% plot(q(:,3)); return
qprev = q(:,3);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tScaling matrices M, R: '); tic;
[MRinv,R] = scalingMatrices(M,N,dx,dy);
Rinv = R^(-1);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tLHS of system: '); tic
[L] = laplacianGen(M,N,dx,dy);
hatA = speye(nq,nq) - dt/Re * L;
A = C'*MRinv*hatA*C;
% zero laplacian
hatA = speye(nq,nq);
A = C'*MRinv*hatA*C;
LBC = zeros(nq,1);
pbc = LBC;
% zero laplacian
% spdiags(C'*MRinv*(speye(nq,nq) - dt/Re * L)*C)
% return
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
duFreestreamdx = zeros(M,nt);
timeIndex = 3;
time = 3*dt;
% =(1+Amplitude*cos(omega*time))*[tanh(10*(0.5*(x(1:end-1)+x(2:end)))) tanh(10*(0.5*(x(1:end-1)+x(2:end)))) tanh(10*(0.5*(x(1:end-1)+x(2:end))))];
for t=1:3
    uInlet(:,t) = inletProfile(y,Amplitude,omega,t*dt,N,ly);
    uFreestream(:,t) = 1+Amplitude*cos(omega*t*dt+k*x(2:end));
    duFreestreamdx(:,t) = -k*Amplitude*sin(omega*t*dt + k*0.5*(x(2:end)+x(1:end-1))); % actually its -du/dx
end

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
fprintf('Current memory usage: \t\t\t%.2f %s\n        ', memoryUsage, units{scale+1});

uminGlob = 0.3;
uminGlobArr = uminGlob;
advPrevPrev = advection(M,N,dx,dy,u,v,2,dt,uInlet,uFreestream,duFreestreamdx);
time=time+dt; 
intermediateSteps = 0;
for timeI = 4:timeMaxIter
    timeNow = 4;
    time = time+dt;

%     compute inlet 
    uInlet(:,timeNow) = inletProfile(y,Amplitude,omega,time,N,ly);
    uFreestream(:,timeNow) = 1+Amplitude*cos(omega*timeI+k*x(2:end));
    duFreestreamdx(:,timeNow) = -k*Amplitude*sin(omega*timeI + k*x(2:end));


%     recalculate particular solution, qp is particular, qpadd is same but used in computation of B 
    ctr = 1;
    uparticular  = zeros(1,nu);
    for i = 1:N
        uparticular(ctr:1:ctr+M-1) = uInlet(i,timeNow);
        ctr = ctr + M;
    end
    qpadd=[(dyq.*uparticular) (dxq.*(0.0))]';
%     return

%     [LBC] = laplacianRHS(M,N,dx,dy,u,v,q,dt,timeNow,x,Amplitude,k,omega,uFreestream,duFreestreamdx);
%     [pbc] = pressure(M,N, dx, dy, q, u, v, dt, timeNow, delta, U0, Re, Amplitude, k, omega,time);

%     % zero laplacian
%     hatA = speye(nq,nq);
%     A = C'*MRinv*hatA*C;
%     LBC = zeros(nq,1);
%     pbc = LBC;

    advPrev = advection(M,N,dx,dy,u,v,timeNow-1,dt,uInlet,uFreestream,duFreestreamdx);
    [ADV] = 1.5*advPrev-0.5*advPrevPrev;
    advPrevPrev = advPrev;

%     [ADV] = zeros(nq,1);
    
    B = dt*CTMM*(1/Re*LBC-ADV-pbc) + CTMRinv*(qprev) - CTMMhatARinv * qpadd;
    
    % using precomputed decomposition (as per recommendation)
    % Using precomputed LU factorization with permutation vectors
    B_permuted = B(Pfact);                     % Apply row permutation to B
    intermediate = Lfact \ B_permuted;         % Solve L * intermediate = B_permuted
    psi_permuted = Ufact \ intermediate;       % Solve U * psi_permuted = intermediate
    
    % Apply column permutation to get the final solution
    psi = zeros(size(B));
    psi(Qfact) = psi_permuted;                 % Rearrange psi_permuted into psi according to Qfact

    % store data
    q(:,timeNow) = C*psi+qpadd;
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
    
    
    if (mod(timeI,floor(stepPlot/20))==0)
        fprintf('\b\b\b\b\b\b\b');
        fprintf('\n%5d.',intermediateSteps);
        intermediateSteps = intermediateSteps+floor(stepPlot/20);
    end

    if (mod(timeI,stepPlot)==0) || (timeI==4)
%         intermediateSteps = 0;
        fprintf('\b\b\b\b\b\b\b');
        fprintf('\n%5d. ',timeI);
        qh = C*psi;
        uh = qh(1:nu);
        vh = qh(nu+1:nq);
        uploth = reshape(uh,M,N);
        vploth = reshape(vh,M,N);
        up = qpadd(1:nu)./dyq';
        vp = qpadd(nu+1:nq)./dxq';
        uplotp = reshape(up,M,N);
        vplotp = reshape(vp,M,N);
        uplot = reshape(u(:,timeNow),M,N);
        vplot = reshape(v(:,timeNow),M,N);
        psiPlot = reshape(psi,M,N);
        vorticity = computeVorticity(u, v, dx, dy, M, N);
        
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
        fprintf('umin=%8.5f       ',uminNow);

        if (any(isnan(uminNow))) || (max(abs(u(:,timeNow)))>50)
            fprintf('\n\n Solution blew up \n\n');
            close all;
            return
        end
        

        % plotting
        tiledlayout(3,2); 
        nexttile;
        contourf(X,flipud(Y),(vploth)','LineColor','none'); hold on;
        colorbar
        title(sprintf('vh velocity contour'));

        nexttile;
        contourf(X,flipud(Y),(vplotp)','LineColor','none'); hold on;
        colorbar
        title(sprintf('vp velocity contour'));

        nexttile;
        contourf(X,flipud(Y),(uploth)','LineColor','none'); hold on;
        colorbar
        title(sprintf('uh velocity contour'));

        nexttile;
        contourf(X,flipud(Y),(uplotp)','LineColor','none'); hold on;
        colorbar
        caxis([0 2])
        title(sprintf('up velocity contour'));

        nexttile(5, [1 2]);
        si = streamslice(X,flipud(Y),uplot',vplot',10,'noarrows');
        set(si,'LineWidth',0.1)
        axis tight
        

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
function [u,v,q] = initialConditions(M,N,nt,dyq,dxq,omega,Amplitude,k,ly,y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Changes first columns of u,v,q,psi to initial conditions    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     uInlet=inletProfile(y,Amplitude,omega,dt,N,ly);
    nu = (M)*(N);
    nv = (M)*(N);
    nq = nu + nv;
    u = zeros(nu,nt);
    v = zeros(nv,nt);
    q = zeros(nq,nt);
%     col = 1; 
%     row = 1;
%     velocity = 1;
%     for i = 1:nu
%         if col > (M)
%             col = 1; 
%             row = row+1;
%             velocity = 1;
%         end
%         qp(i) = velocity*dyq(i);
%         col = col + 1;
%     end
    qp(1:nu) = 0;
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
    % old value
    % gamma = max(max(abs(uprev))*dt/min(dx),max(abs(vprev))*dt/min(dy));
    gamma = 40*dt/min(dx(1),dy(1));
    % TEST increase by 10
%     gamma = 20*gamma;
    
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
            duudx(iu) = (uuright - uu(leftcell) - gamma*(uHorAvg_uHorDiff_GHOST - uHorAvg_uHorDiff(leftcell)))/(deltax);
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
function [L] = laplacianGen(M, N, dx, dy)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %          Generate Laplacian Matrix L using            %
    %                 Sparse Matrices                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Ensure dx and dy are column vectors
    dx = dx(:);
    dy = dy(:);

    % Initialize dimensions
    nu = M * N;
    nv = M * N;
    nq = nu + nv;

    % Extend dx and dy to include ghost points
    dx(M + 1) = dx(M);
    dy(N + 1) = dy(N);

    % Initialize index and value arrays for Luxx components
    I_Luxx = [];
    J_Luxx = [];
    S_Luxx = [];

    % Build Luxx
    for iu = 1:nu
        row = ceil(iu / M);
        col = mod(iu - 1, M) + 1;

        if col == 1
            hw = dx(1);
            he = dx(2);
            hc = 0.5 * (he + hw);

            % Coefficients
            val_center = -2 / (hw * he);
            val_east = 1 / (he * hc);

            % Append indices and values
            I_Luxx = [I_Luxx; iu; iu];
            J_Luxx = [J_Luxx; iu; iu + 1];
            S_Luxx = [S_Luxx; val_center; val_east];
        elseif col == M
            hw = dx(M - 1);
            he = dx(M);
            hc = 0.5 * (he + hw);

            % Coefficients
            val_center = -2 / (hw * he);
            val_west = 1 / (hw * hc) + 1 / (he * hc);

            % Append indices and values
            I_Luxx = [I_Luxx; iu; iu];
            J_Luxx = [J_Luxx; iu; iu - 1];
            S_Luxx = [S_Luxx; val_center; val_west];
        else
            he = dx(col + 1);
            hw = dx(col);
            hc = 0.5 * (he + hw);

            % Coefficients
            val_center = -2 / (hw * he);
            val_west = 1 / (hw * hc);
            val_east = 1 / (he * hc);

            % Append indices and values
            I_Luxx = [I_Luxx; iu; iu; iu];
            J_Luxx = [J_Luxx; iu; iu - 1; iu + 1];
            S_Luxx = [S_Luxx; val_center; val_west; val_east];
        end
    end

    % Build Luyy
    I_Luyy = [];
    J_Luyy = [];
    S_Luyy = [];

    for iu = 1:nu
        row = ceil(iu / M);
        col = mod(iu - 1, M) + 1;

        if row == 1
            hn = 0.5 * (dy(row) + dy(row + 1));
            hs = dy(row);
            hc = dy(row);

            % Coefficients
            val_center = -2 / (hs * hn) - 1 / (hs * hc);
            val_north = 1 / (hn * hc);

            % Append indices and values
            I_Luyy = [I_Luyy; iu; iu];
            J_Luyy = [J_Luyy; iu; iu + M];
            S_Luyy = [S_Luyy; val_center; val_north];
        elseif row == N
            hn = dy(row);
            hc = dy(row);
            hs = 0.5 * (dy(row) + dy(row - 1));

            % Coefficients
            val_center = -(2 * hs + hc) / (hc * hc * hs);
            val_south = 1 / (hs * hc);

            % Append indices and values
            I_Luyy = [I_Luyy; iu; iu];
            J_Luyy = [J_Luyy; iu; iu - M];
            S_Luyy = [S_Luyy; val_center; val_south];
        else
            hc = dy(row);
            hs = 0.5 * (dy(row) + dy(row - 1));
            hn = 0.5 * (dy(row) + dy(row + 1));

            % Coefficients
            val_center = -2 / (hs * hn);
            val_south = 1 / (hs * hc);
            val_north = 1 / (hn * hc);

            % Append indices and values
            I_Luyy = [I_Luyy; iu; iu; iu];
            J_Luyy = [J_Luyy; iu; iu - M; iu + M];
            S_Luyy = [S_Luyy; val_center; val_south; val_north];
        end
    end

    % Build Lu as sparse matrix
    Luxx = sparse(I_Luxx, J_Luxx, S_Luxx, nu, nu);
    Luyy = sparse(I_Luyy, J_Luyy, S_Luyy, nu, nu);
    Lu = Luxx + Luyy;

    % Build Lvxx
    I_Lvxx = [];
    J_Lvxx = [];
    S_Lvxx = [];

    for iv = 1:nv
        row = ceil(iv / M);
        col = mod(iv - 1, M) + 1;

        if col == 1
            he = 0.5 * (dx(1) + dx(2));
            hw = dx(1);
            hc = dx(1);

            % Coefficients
            val_center = -2 / (hw * he) - 1 / (hw * hc);
            val_east = 1 / (he * hc);

            % Append indices and values
            I_Lvxx = [I_Lvxx; iv; iv];
            J_Lvxx = [J_Lvxx; iv; iv + 1];
            S_Lvxx = [S_Lvxx; val_center; val_east];
        elseif col == M
            he = dx(M);
            hc = dx(M);
            hw = 0.5 * (dx(M) + dx(M - 1));

            % Coefficients
            val_center = -1 / (hw * hc);
            val_west = 1 / (hw * hc);

            % Append indices and values
            I_Lvxx = [I_Lvxx; iv; iv];
            J_Lvxx = [J_Lvxx; iv; iv - 1];
            S_Lvxx = [S_Lvxx; val_center; val_west];
        else
            hw = 0.5 * (dx(col) + dx(col - 1));
            he = 0.5 * (dx(col) + dx(col + 1));
            hc = dx(col);

            % Coefficients
            val_center = -2 / (hw * he);
            val_west = 1 / (hw * hc);
            val_east = 1 / (he * hc);

            % Append indices and values
            I_Lvxx = [I_Lvxx; iv; iv; iv];
            J_Lvxx = [J_Lvxx; iv; iv - 1; iv + 1];
            S_Lvxx = [S_Lvxx; val_center; val_west; val_east];
        end
    end

    % Build Lvyy
    I_Lvyy = [];
    J_Lvyy = [];
    S_Lvyy = [];

    for iv = 1:nv
        row = ceil(iv / M);
        col = mod(iv - 1, M) + 1;

        if row == 1
            hn = dy(row + 1);
            hs = dy(row);
            hc = 0.5 * (hn + hs);

            % Coefficients
            val_center = -2 / (hn * hs);
            val_north = 1 / (hc * hn);

            % Append indices and values
            I_Lvyy = [I_Lvyy; iv; iv];
            J_Lvyy = [J_Lvyy; iv; iv + M];
            S_Lvyy = [S_Lvyy; val_center; val_north];
        elseif row == N
            hs = dy(row);
            hn = dy(row + 1);
            hc = 0.5 * (hn + hs);

            % Coefficients
            val_center = -2 / (hn * hs);
            val_south = 1 / (hc * hs) + 1 / (hn * hc);

            % Append indices and values
            I_Lvyy = [I_Lvyy; iv; iv];
            J_Lvyy = [J_Lvyy; iv; iv - M];
            S_Lvyy = [S_Lvyy; val_center; val_south];
        else
            hs = dy(row);
            hn = dy(row + 1);
            hc = 0.5 * (hn + hs);

            % Coefficients
            val_center = -2 / (hn * hs);
            val_south = 1 / (hc * hs);
            val_north = 1 / (hc * hn);

            % Append indices and values
            I_Lvyy = [I_Lvyy; iv; iv; iv];
            J_Lvyy = [J_Lvyy; iv; iv - M; iv + M];
            S_Lvyy = [S_Lvyy; val_center; val_south; val_north];
        end
    end

    % Build Lv as sparse matrix
    Lvxx = sparse(I_Lvxx, J_Lvxx, S_Lvxx, nv, nv);
    Lvyy = sparse(I_Lvyy, J_Lvyy, S_Lvyy, nv, nv);
    Lv = Lvxx + Lvyy;

    % Shift indices for v components
    [I_Lv, J_Lv, S_Lv] = find(Lv);
    I_Lv_shifted = I_Lv + nu;
    J_Lv_shifted = J_Lv + nu;

    % --------------------
    % Assemble L
    % --------------------

    % Get indices and values from Lu
    [I_Lu, J_Lu, S_Lu] = find(Lu);

    % Combine indices and values
    I_L = [I_Lu; I_Lv_shifted];
    J_L = [J_Lu; J_Lv_shifted];
    S_L = [S_Lu; S_Lv];

    % Build L as sparse matrix
    L = sparse(I_L, J_L, S_L, nq, nq);
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
    deltay1=(dy(N)+dy(N-1))*0.5;
    deltay2 = dy(N);
    pCornerTopFromRightBC = deltay2/deltay1 * (pRHS(nu) - pRHS(nu-M)) + pRHS(nu);
    pCornerRightFromTopBC = -(2*A*(Re*k + Re*w) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2)^3 + 2*A*tan((k*x_coord)/2 + (t*w)/2)^2*(Re*k + Re*w - A*Re*k) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2))/(Re*k*(tan((k*x_coord)/2 + (t*w)/2)^2 + 1)^2);
    deltap = pCornerTopFromRightBC - pCornerRightFromTopBC;
    % shift the top values to match the right (solve up to a constant) 
    % pRHS(nq-M+1:1:nq) = pRHS(nq-M+1:1:nq) + deltap;
    pRHS(M:M:nu) = pRHS(M:M:nu) - deltap;
end
function [uInlet] = inletProfile(y,A,w,time,N,ly)
    const=((1) + 2e-18 * (0.5 - rand(1, N)));
    smooth = (A*cos(w*time)*(0.5+0.5*tanh(50*(0.5*(y(2:end)+y(1:end-1))-(N-5)*ly/N)))');
%     uInlet=(1+A*cos(w*time)) + 2e-8 * (0.5 - rand(1, N)) + A*cos(w*time)*(0.5+0.5*tanh(1000*(y-9*ly/10)))';%e*[tanh(smootheningFactor*(0.5*(y(1:end-1)+y(2:end))))];
    uInlet = const + smooth;
    %     uInlet(1) = uInlet(1)/2;
end
function vorticity = computeVorticity(u, v, dx, dy, M, N)
    % computeVorticity: Computes the vorticity field from velocity components
    % u, v: velocity components (nu x nt matrices)
    % dx, dy: grid spacings in x and y directions (vectors)
    % M, N: grid dimensions
    % vorticity: vorticity field ((M-1)x(N-1) matrix)

    % Reshape u and v to M x N matrices (latest time step)
    u_latest = reshape(u(:, end), M, N);
    v_latest = reshape(v(:, end), M, N);

    % Compute du/dy at cell centers
    dy_mat = repmat(dy(1:N-1)', M-1, 1);
    du_dy = (u_latest(2:M, 2:N) - u_latest(2:M, 1:N-1)) ./ dy_mat;

    % Compute dv/dx at cell centers
    dx_mat = repmat(dx(1:M-1), 1, N-1);
    dv_dx = (v_latest(2:M, 2:N) - v_latest(1:M-1, 2:N)) ./ dx_mat;

    % Compute vorticity at cell centers
    vorticity = dv_dx - du_dy;
end
