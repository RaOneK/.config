clc; clear; close all
omega = -2*pi*4;
k = 2*pi*4;
phi = pi/2;
Amplitude = 1.0;
targetValue = Amplitude;       % Desired target Amplitude value
totalIterations = 3e1;  % Total number of iterations for the smooth increase
M = 80;            % grid size in x-dir
N = 80;            % grid size in y-dir
smootheningFactor = 1e9; % bottom left corner, higher -> sharper

stepPlot = 150;

nt =4;              % timesteps used for storage
timeMaxIter = 1e6;  % number of timesteps
dt = 1e-4;          % desired timestep size
Re = 122;           % Reynolds number
U0 = 1+targetValue;             % characteristic velocity scale
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


% check for Diffusion and Courant number stability
D = dt/min(min(dx),min(dy))^2/Re;
C = 1.2*dt/min(min(dx),min(dy));
sqrt2d=sqrt(2*D);
fprintf('Diffusion number D = %8.3f ',D);
if D>0.5
    fprintf('> 0.5, must be <0.5 for explicit schemes (current scheme is implicit)\n');
end
fprintf('\nCourant Number C = %6.3f ',C);
if  C>sqrt2d
    fprintf('> %8.3f, must be <sqrt(2D)\n',sqrt2d);
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
[u,v,q] = initialConditions(M,N,nt,dyq,dxq);
qprev = q(:,3);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);


fprintf('\tScaling matrices M, R: '); tic;
[MRinv,R] = scalingMatrices(M,N,dx,dy);
Rinv = R^(-1);
fprintf('\t\tcompleted in %6.4f sec!\n', toc);

fprintf('\tLHS of system: '); tic
% [L] = laplacianGen(M,N,dx,dy);
hatA = speye(nq,nq);
A = C'*MRinv*hatA*C;
LBC = zeros(nq,1);
pbc = LBC;
fprintf('\t\t\tcompleted in %6.4f sec!\n', toc);

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
for t=1:3
    time = dt*t;
    uInlet(:,t) = inletProfile(y,Amplitude,omega,time,N,ly,phi);
    uFreestream(:,t) = uFree(x,Amplitude,k,omega,time,phi);
    duFreestreamdx(:,t) = duFree(Amplitude,k,omega,time,x,phi);
end


% particular solution
qp = zeros(nq,1);
ctr = 1;
uparticular  = zeros(1,nu);
for i = 1:N
    uparticular(ctr:1:ctr+M-1) = uInlet(i,4);
    ctr = ctr + M;
end
qp(1:nu) = dyq.*uparticular;    
qp(nu+1:nq) = dxq.*(0.0);

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
fprintf('Current memory usage: \t\t\t%.2f %s\n        ', memoryUsage, units{scale+1});

advPrevPrev = advection(M,N,dx,dy,u,v,2,dt,uInlet,uFreestream,duFreestreamdx);
time=time+dt; 
intermediateSteps = 0;
uminGlob = 1;
uminGlobArr = [];

% Initial positions of particles at the inlet (x = 0)
[py] = linspace(0, ly, 75)';  % Initial vertical positions along the inlet
px = zeros(size(py));  % All particles start at x = 0


timeNow = 1;
qp = generateParticularSolutionFlatPlate(M, N, dxq, dyq, uInlet, timeNow);
streamslice(X, flipud(Y), reshape(qp(1:nu), M, N)', reshape(qp(nu+1:end), M, N)');
for timeI = 4:timeMaxIter
    timeNow = 4;
    time = time+dt;
    Amplitude = smoothAmplitude(timeI, targetValue, totalIterations);
    
%     compute inlet 
    uInlet(:,timeNow) = inletProfile(y,Amplitude,omega,time,N,ly,phi);
    uFreestream(:,timeNow) = uFree(x,Amplitude,k,omega,time,phi); 
    duFreestreamdx(:,timeNow) = duFree(Amplitude,k,omega,time,x,phi);

%     uInlet(1,timeNow) = 0.5*uInlet(1,timeNow);

%     recalculate particular solution, qp is particular, qpadd is same but used in computation of B 
    ctr = 1;
    qp = generateParticularSolutionFlatPlate(M, N, dxq, dyq, uInlet, timeNow);    

    advPrev = advection(M,N,dx,dy,u,v,timeNow-1,dt,uInlet,uFreestream,duFreestreamdx);
    [ADV] = 1.5*advPrev-0.5*advPrevPrev;
    advPrevPrev = advPrev;

    [ADV] = zeros(nq,1);
   
    qprev(nq-M+1:1:nq) = 0;
    qp(nq-M+1:1:nq) = 0;
    LBC(nq-M+1:1:nq) = 0;
    ADV(nq-M+1:1:nq) = 0;
    pbc(nq-M+1:1:nq) = 0;
    v(nv-M+1:1:nv,:) = 0;
    q(nq-M+1:1:nq,:) = 0;

    
    B = dt*CTMM*(1/Re*LBC-ADV-pbc) + CTMRinv*qprev - CTMMhatARinv * qp;

    % using precomputed decomposition (as per recommendation)
    % Using precomputed LU factorization with permutation vectors
    B_permuted = B(Pfact);                     % Apply row permutation to B
    intermediate = Lfact \ B_permuted;         % Solve L * intermediate = B_permuted
    psi_permuted = Ufact \ intermediate;       % Solve U * psi_permuted = intermediate
    
    % Apply column permutation to get the final solution
    psi = zeros(size(B));
    psi(Qfact) = psi_permuted;                 % Rearrange psi_permuted into psi according to Qfact

    % Compute particle displacement in every timestep
    [px, py] = advectParticles(px, py, u, v, M, N, dt, X, Y);

    % Add new particles from the inlet at each time step (particles are injected at x = 0)
    if mod(timeI,stepPlot/20)==0
        new_py = linspace(0, max(Y(:)), 75)'; % Adjust the number of particles injected
        new_px = zeros(size(new_py));
        px = [px; new_px];
        py = [py; new_py];
    end 

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
    
    
    if (mod(timeI,floor(stepPlot/20))==0)
        fprintf('\b\b\b\b\b\b\b\b');
        fprintf('\n%6d.',intermediateSteps);
        intermediateSteps = intermediateSteps+floor(stepPlot/20);
    end

    if (mod(timeI,stepPlot)==0) || (timeI==4)
        fprintf('\b\b\b\b\b\b\b\b');
        fprintf('\n%6d. ',timeI);
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
        [max_velocity, max_coord] = findMaxVelocity(u(:, timeNow), v(:, timeNow), M, N, dx, dy);
        
        % Print the maximum velocity and its scaled coordinate only at specific intervals
        fprintf('Max abs velocity: %8.5f at (%8.5f, %8.5f) scaled coordinate \n', max_velocity, max_coord(1), max_coord(2));



        if (any(isnan(max_velocity))) || (max_velocity>50)
            fprintf('\n\n Solution blew up \n\n');
            close all;
            return
        end
        

        % plotting
        tiledlayout(4,6); 
        nexttile(3,[1 1]);
        contourf(X, (Y), ((vploth).^2 + (uploth).^2)', 50, 'LineColor', 'none'); hold on; colorbar;
        sih = streamslice(X,(Y),(uploth'),(vploth'),'noarrows'); 
        set(sih,'Color','red');
        colorbar
        set(gca, 'YDir', 'reverse'); % Flips the Y-axis direction
        title (sprintf('psi contour (homogeneous solution)'))
       
        nexttile(6,[1 1]);
        % Reshape qp for u and v components
        u_particular = reshape(qp(1:nu), M, N);
        v_particular = reshape(qp(nu+1:end), M, N);
        
        % Overlay streamlines to visualize flow
        contourf(X, (Y), (uplotp.^2)', 'LineColor', 'none'); hold on;
        sih = streamslice(X,(Y),(uplotp'),(vplotp'),'noarrows'); 
        set(gca, 'YDir', 'reverse'); % Flips the Y-axis direction
        set(sih, 'Color', 'red');
        
        
        % Set plot properties
        colorbar;
        title(sprintf('Particular Solution at Time %8d', timeI));
        
        % Plotting velocity field and streamlines
        nexttile(10, [3 2]);
        % Plotting velocity field and streamlines
        nexttile(10, [3 2]);
        contourf(X, (Y), (uplot.^2 + vplot.^2)', 'LineColor', 'none'); hold on;
        si = streamslice(X, (Y), uplot', vplot', 10, 'noarrows');
        set(si, 'LineWidth', 0.1, 'Color', 'red');
        axis tight;
        colorbar;
        set(gca, 'YDir', 'reverse'); % Flips the Y-axis direction
        title(sprintf('Total velocity field at %8d', timeI));

        % Plot particles
        plotParticles(px, py, X, Y)

        nexttile(1,[2 2]);
        s = surf(flipud(Y), flipud(X), flipud(uplot'));
        s.EdgeColor = 'none';
        set(si,'Color','red');
        title(sprintf('u surf. front left = freestream, front right = inlet'));

        nexttile(13,[2 2]);
        s=surf(flipud(Y), flipud(X), flipud(vplot')); hold on;
        s.EdgeColor = 'none';
        colorbar
        title(sprintf('v surf. front left = freestream, front right = inlet'));
        
        % plot velocity top
        nexttile(4, [1 2]);
        plot(0.5 * (x(1:end-1) + x(2:end)), v(nv-M+1:1:nv, timeNow), 'r', 'LineWidth', 1.9); % Existing top velocity BC in red
        hold on;
        % Adding freestream velocity
        plot(0.5 * (x(1:end-1) + x(2:end)), uFreestream(:, timeNow), '--b', 'LineWidth', 1.2); % Freestream velocity in blue (dashed line)
        
        % Adding legend for better clarity
        legend({'v component (top BC)', 'Freestream u component'}, 'Location', 'best');
        
        hold off;
        title('top velocity BC (including freestream)');
        xlabel('x');
        ylabel('Velocity');

        % plot bc left
        nexttile(9,[3,1]);
        plot([0; uInlet(:,timeNow)]',[0; 0.5*(y(1:end-1)+y(2:end))]);
        title(sprintf('inlet u component BC=%8.5f',(uInlet(N,timeNow))));

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
    
    % Zero out entries corresponding to top boundary v-components
    topBoundaryIndices = find(J3 >= nq-M+1);
    S3(topBoundaryIndices) = 0;
    
    % Fourth set: D(i, M*(N-1) + i) = -1 for i = M+1:nodes
    I4 = ((M + 1):nodes)';
    J4 = M * (N - 1) + I4;
    S4 = -ones(nodes - M, 1);
    
    % Zero out entries corresponding to top boundary v-components
    topBoundaryIndices = find(J4 >= nq-M+1);
    S4(topBoundaryIndices) = 0;
    
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
function [ADV] = advection(M,N,dx,dy,u,v,timeIndex,dt,uInlet,uFreestream,duFreestreamdx)
    % Get previous time step velocities
    uprev = u(:,timeIndex);
    vprev = v(:,timeIndex);

    % Compute cell-centered velocities and their differences
    [uu, uHorAvg_uHorDiff, vv, vVerAvg_vVerDiff] = computeCellCenteredVelocities(M, N, dx, dy, uprev, vprev, uInlet, timeIndex);

    % Compute node-centered velocities and their differences
    [uv, uVerAvg_vHorDiff, vHorAvg_uVerDiff] = computeNodeCenteredVelocities(M, N, dx, dy, uprev, vprev, uFreestream, timeIndex);

    % Compute gamma
    gamma = min(1.2 * dt * max(max(abs(u(:,timeIndex))), max(abs(v(:,timeIndex)))), 1);
%     gamma = 1;

    % Compute duudx and duvdy
    [duudx, duvdy] = computeDuudxDuvdy(M, N, dx, dy, uprev, vprev, uu, uHorAvg_uHorDiff, uv, vHorAvg_uVerDiff, gamma, uInlet, timeIndex);

    % Compute dvudx and dvvdy
    [dvudx, dvvdy] = computeDvudxDvvdy(M, N, dx, dy, uprev, vprev, vv, vVerAvg_vVerDiff, uv, uVerAvg_vHorDiff, gamma, duFreestreamdx, timeIndex);

    % Construct the advection term
    ADV = [duudx + duvdy; dvudx + dvvdy];
end

function [uu, uHorAvg_uHorDiff, vv, vVerAvg_vVerDiff] = computeCellCenteredVelocities(M, N, dx, dy, uprev, vprev, uInlet, timeIndex)
    % Initialize variables
    nc = M * N;
    uu = zeros(nc,1);
    vv = zeros(nc,1);
    uHorAvg_uHorDiff = zeros(nc,1);
    vVerAvg_vVerDiff = zeros(nc,1);

    for iCell = 1:nc
        col = mod(iCell-1, M) + 1;
        row = ceil(iCell/M);
        ijleft = (row-1)*(M) + col - 1;
        ijright = (row-1)*(M) + col - 1 + 1;
        ijbot = (row - 2) * M + col;
        ijtop = (row - 2 + 1) * M + col;

        % Compute uu and uHorAvg_uHorDiff
        if col == 1
            uu(iCell) = 0.25 * (uInlet(row,timeIndex) + uprev(ijright))^2;
            uHorAvg_uHorDiff(iCell) = 0.25 * abs(uInlet(row,timeIndex) + uprev(ijright)) * (uprev(ijright) - uInlet(row,timeIndex));
        else
            uu(iCell) = 0.25 * (uprev(ijleft) + uprev(ijright))^2;
            uHorAvg_uHorDiff(iCell) = 0.25 * (uprev(ijleft) + uprev(ijright)) * (uprev(ijright) - uprev(ijleft));
        end

        % Compute vv and vVerAvg_vVerDiff
        if row == 1
            vv(iCell) = 0.25 * (vprev(ijtop) + 0)^2;
            vVerAvg_vVerDiff(iCell) = 0.25 * (vprev(ijtop) + 0) * (vprev(ijtop) - 0);
        elseif row == N
            vv(iCell) = 0.25 * (vprev(ijbot) + 0)^2;
            vVerAvg_vVerDiff(iCell) = 0.25 * (vprev(ijbot) + 0) * (vprev(ijbot) - 0);
        else
            vv(iCell) = 0.25 * (vprev(ijtop) + vprev(ijbot))^2;
            vVerAvg_vVerDiff(iCell) = 0.25 * (vprev(ijtop) + vprev(ijbot)) * (vprev(ijtop) - vprev(ijbot));
        end
    end
end

function [uv, uVerAvg_vHorDiff, vHorAvg_uVerDiff] = computeNodeCenteredVelocities(M, N, dx, dy, uprev, vprev, uFreestream, timeIndex)
    % Initialize variables
    nn = (M+1)*(N+1);
    uv = zeros(nn,1);
    uVerAvg_vHorDiff = zeros(nn,1);
    vHorAvg_uVerDiff = zeros(nn,1);
    dy(N+1) = dy(N);  % Extend dy to N+1

    for iuv = 1:nn
        row = floor((iuv-1)/(M+1)) + 1;
        col = mod(iuv-1, M+1) + 1;
        ijtop = (M-1) * (row-1) + col - 1;
        ijbot = ijtop - (M-1);
        ijright = (M) * (row-1) + col - M;
        ijleft = ijright - 1;

        % Handle boundaries and compute uv, uVerAvg_vHorDiff, vHorAvg_uVerDiff
        if (col == 1)
            uv(iuv) = 0;
            uVerAvg_vHorDiff(iuv) = 0;
            vHorAvg_uVerDiff(iuv) = 0;
        elseif (row == 1)
            uv(iuv) = 0;
            uVerAvg_vHorDiff(iuv) = 0;
            vHorAvg_uVerDiff(iuv) = 0;
        elseif (col == M+1) 
            hs = dy(row-1);
            hn = dy(row);
            he = dx(col-1);
            uv(iuv) = vprev(ijleft) * (hs * uprev(ijtop) + hn * uprev(ijbot)) / (hn + hs);
            uVerAvg_vHorDiff(iuv) = 0;
            vHorAvg_uVerDiff(iuv) = 0;
        elseif (row == N+1) && (col ~= M+1)
            hw = dx(col-1);
            he = dx(col);
            uv(iuv) = uFreestream(col-1, timeIndex) * (he * vprev(ijleft) + hw * vprev(ijright)) / (he + hw);
            uVerAvg_vHorDiff(iuv) = abs(uFreestream(col-1, timeIndex)) * (hw * vprev(ijright) - he * vprev(ijleft)) / (he + hw);
            vHorAvg_uVerDiff(iuv) = abs((hw * vprev(ijright) + he * vprev(ijleft)) / (he + hw)) * (uFreestream(col-1, timeIndex) - uprev(ijbot));
            uv(iuv) = 0;
            uVerAvg_vHorDiff(iuv) = 0;
            vHorAvg_uVerDiff(iuv) = 0;
%         elseif (col == M+1) && (row == N+1)
%             hs = dy(row-1);
%             hn = dy(row);
%             uv(iuv) = vprev(ijleft) * uFreestream(col-1, timeIndex);
%             uVerAvg_vHorDiff(iuv) = 0;
%             vHorAvg_uVerDiff(iuv) = abs(vprev(ijleft)) * (uFreestream(col-1, timeIndex) - uprev(ijbot));
        else
            hs = dy(row-1);
            hn = dy(row);
            hw = dx(col-1);
            he = dx(col);
            uv(iuv) = ((he * vprev(ijleft) + hw * vprev(ijright)) / (he + hw)) * ((hs * uprev(ijtop) + hn * uprev(ijbot)) / (hn + hs));
            uVerAvg_vHorDiff(iuv) = abs((hs * uprev(ijtop) + hn * uprev(ijbot)) / (hn + hs)) * (hw * vprev(ijright) - he * vprev(ijleft)) / (he + hw);
            vHorAvg_uVerDiff(iuv) = abs((hw * vprev(ijright) + he * vprev(ijleft)) / (he + hw)) * ((hs * uprev(ijtop) - hn * uprev(ijbot)) / (hn + hs));
        end
    end
end

function [duudx, duvdy] = computeDuudxDuvdy(M, N, dx, dy, uprev, vprev, uu, uHorAvg_uHorDiff, uv, vHorAvg_uVerDiff, gamma, uInlet, timeIndex)
    nu = M * N;
    duudx = zeros(nu,1);
    duvdy = zeros(nu,1);
    dx(M+1) = dx(M);  % Extend dx to M+1
    dy(N+1) = dy(N);  % Extend dy to N+1

    for iu = 1:nu
        row = floor((iu-1)/(M)) + 1;
        col = mod(iu-1, M) + 1;
        leftcell = (row-1) * (M) + col;
        rightcell = leftcell + 1;
        botnode = (row-1) * (M+1) + col + 1;
        topnode = botnode + M + 1;
        deltax = 0.5 * (dx(col) + dx(col+1));
        deltay = dy(row);

        if row == N % top bdry 
            duvdy(iu) = (0 - uv(botnode) - gamma * (0 - vHorAvg_uVerDiff(botnode))) / deltay;
        else
            duvdy(iu) = (uv(topnode) - uv(botnode) - gamma * (vHorAvg_uVerDiff(topnode) - vHorAvg_uVerDiff(botnode))) / deltay;
        end

        if col == M
            % Handle boundary at rightmost column
            vTprevprev = vprev(iu);
            if row > 1
                vBprevprev = vprev(iu - M);
            else
                vBprevprev = 0;
            end
            ug = -2 * deltax * (vTprevprev - vBprevprev) / deltay + uprev(iu - 1);
            uuright = 0.25 * (ug + uprev(iu))^2;
            uHorAvg_uHorDiff_GHOST = 0.25 * abs(ug + uprev(iu)) * (ug - uprev(iu));
            duudx(iu) = (uuright - uu(leftcell) - gamma * (uHorAvg_uHorDiff_GHOST - uHorAvg_uHorDiff(leftcell))) / deltax;
        else
            duudx(iu) = (uu(rightcell) - uu(leftcell) - gamma * (uHorAvg_uHorDiff(rightcell) - uHorAvg_uHorDiff(leftcell))) / deltax;
        end
    end
end

function [dvudx, dvvdy] = computeDvudxDvvdy(M, N, dx, dy, uprev, vprev, vv, vVerAvg_vVerDiff, uv, uVerAvg_vHorDiff, gamma, duFreestreamdx, timeIndex)
    nv = M * N;
    dvudx = zeros(nv,1);
    dvvdy = zeros(nv,1);
    dx(M+1) = dx(M);  % Extend dx to M+1
    dy(N+1) = dy(N);  % Extend dy to N+1

    for iv = 1:nv
        row = floor((iv-1)/(M)) + 1;
        col = mod(iv-1, M) + 1;
        botcell = (row-1) * (M) + col;
        topcell = botcell + M;
        leftnode = (row-1) * (M+1) + col + 1 + M;
        rightnode = leftnode + 1;
        deltax = dx(col);
        deltay = 0.5 * (dy(row) + dy(row+1));
        dvudx(iv) = (uv(rightnode) - uv(leftnode) - gamma * (uVerAvg_vHorDiff(rightnode) - uVerAvg_vHorDiff(leftnode))) / deltax;

        if row == N
            % Handle boundary at top row
            dvvdy(iv) = 0;
        else
            dvvdy(iv) = (vv(topcell) - vv(botcell) - gamma * (vVerAvg_vVerDiff(topcell) - vVerAvg_vVerDiff(botcell))) / deltay;
        end
    end
end

% function [L] = laplacianGen(M, N, dx, dy)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %          Generate Laplacian Matrix L using            %
%     %                 Sparse Matrices                       %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % Ensure dx and dy are column vectors
%     dx = dx(:);
%     dy = dy(:);
% 
%     % Initialize dimensions
%     nu = M * N;
%     nv = M * N;
%     nq = nu + nv;
% 
%     % Extend dx and dy to include ghost points
%     dx(M + 1) = dx(M);
%     dy(N + 1) = dy(N);
% 
%     % Initialize index and value arrays for Luxx components
%     I_Luxx = [];
%     J_Luxx = [];
%     S_Luxx = [];
% 
%     % Build Luxx
%     for iu = 1:nu
%         row = ceil(iu / M);
%         col = mod(iu - 1, M) + 1;
% 
%         if col == 1
%             hw = dx(1);
%             he = dx(2);
%             hc = 0.5 * (he + hw);
% 
%             % Coefficients
%             val_center = -2 / (hw * he);
%             val_east = 1 / (he * hc);
% 
%             % Append indices and values
%             I_Luxx = [I_Luxx; iu; iu];
%             J_Luxx = [J_Luxx; iu; iu + 1];
%             S_Luxx = [S_Luxx; val_center; val_east];
%         elseif col == M
%             hw = dx(M - 1);
%             he = dx(M);
%             hc = 0.5 * (he + hw);
% 
%             % Coefficients
%             val_center = -2 / (hw * he);
%             val_west = 1 / (hw * hc) + 1 / (he * hc);
% 
%             % Append indices and values
%             I_Luxx = [I_Luxx; iu; iu];
%             J_Luxx = [J_Luxx; iu; iu - 1];
%             S_Luxx = [S_Luxx; val_center; val_west];
%         else
%             he = dx(col + 1);
%             hw = dx(col);
%             hc = 0.5 * (he + hw);
% 
%             % Coefficients
%             val_center = -2 / (hw * he);
%             val_west = 1 / (hw * hc);
%             val_east = 1 / (he * hc);
% 
%             % Append indices and values
%             I_Luxx = [I_Luxx; iu; iu; iu];
%             J_Luxx = [J_Luxx; iu; iu - 1; iu + 1];
%             S_Luxx = [S_Luxx; val_center; val_west; val_east];
%         end
%     end
% 
%     % Build Luyy
%     I_Luyy = [];
%     J_Luyy = [];
%     S_Luyy = [];
% 
%     for iu = 1:nu
%         row = ceil(iu / M);
%         col = mod(iu - 1, M) + 1;
% 
%         if row == 1
%             hn = 0.5 * (dy(row) + dy(row + 1));
%             hs = dy(row);
%             hc = dy(row);
% 
%             % Coefficients
%             val_center = -2 / (hs * hn) - 1 / (hs * hc);
%             val_north = 1 / (hn * hc);
% 
%             % Append indices and values
%             I_Luyy = [I_Luyy; iu; iu];
%             J_Luyy = [J_Luyy; iu; iu + M];
%             S_Luyy = [S_Luyy; val_center; val_north];
%         elseif row == N
%             hn = dy(row);
%             hc = dy(row);
%             hs = 0.5 * (dy(row) + dy(row - 1));
% 
%             % Coefficients
%             val_center = -(2 * hs + hc) / (hc * hc * hs);
%             val_south = 1 / (hs * hc);
% 
%             % Append indices and values
%             I_Luyy = [I_Luyy; iu; iu];
%             J_Luyy = [J_Luyy; iu; iu - M];
%             S_Luyy = [S_Luyy; val_center; val_south];
%         else
%             hc = dy(row);
%             hs = 0.5 * (dy(row) + dy(row - 1));
%             hn = 0.5 * (dy(row) + dy(row + 1));
% 
%             % Coefficients
%             val_center = -2 / (hs * hn);
%             val_south = 1 / (hs * hc);
%             val_north = 1 / (hn * hc);
% 
%             % Append indices and values
%             I_Luyy = [I_Luyy; iu; iu; iu];
%             J_Luyy = [J_Luyy; iu; iu - M; iu + M];
%             S_Luyy = [S_Luyy; val_center; val_south; val_north];
%         end
%     end
% 
%     % Build Lu as sparse matrix
%     Luxx = sparse(I_Luxx, J_Luxx, S_Luxx, nu, nu);
%     Luyy = sparse(I_Luyy, J_Luyy, S_Luyy, nu, nu);
%     Lu = Luxx + Luyy;
% 
%     % Build Lvxx
%     I_Lvxx = [];
%     J_Lvxx = [];
%     S_Lvxx = [];
% 
%     for iv = 1:nv
%         row = ceil(iv / M);
%         col = mod(iv - 1, M) + 1;
% 
%         if col == 1
%             he = 0.5 * (dx(1) + dx(2));
%             hw = dx(1);
%             hc = dx(1);
% 
%             % Coefficients
%             val_center = -2 / (hw * he) - 1 / (hw * hc);
%             val_east = 1 / (he * hc);
% 
%             % Append indices and values
%             I_Lvxx = [I_Lvxx; iv; iv];
%             J_Lvxx = [J_Lvxx; iv; iv + 1];
%             S_Lvxx = [S_Lvxx; val_center; val_east];
%         elseif col == M
%             he = dx(M);
%             hc = dx(M);
%             hw = 0.5 * (dx(M) + dx(M - 1));
% 
%             % Coefficients
%             val_center = -1 / (hw * hc);
%             val_west = 1 / (hw * hc);
% 
%             % Append indices and values
%             I_Lvxx = [I_Lvxx; iv; iv];
%             J_Lvxx = [J_Lvxx; iv; iv - 1];
%             S_Lvxx = [S_Lvxx; val_center; val_west];
%         else
%             hw = 0.5 * (dx(col) + dx(col - 1));
%             he = 0.5 * (dx(col) + dx(col + 1));
%             hc = dx(col);
% 
%             % Coefficients
%             val_center = -2 / (hw * he);
%             val_west = 1 / (hw * hc);
%             val_east = 1 / (he * hc);
% 
%             % Append indices and values
%             I_Lvxx = [I_Lvxx; iv; iv; iv];
%             J_Lvxx = [J_Lvxx; iv; iv - 1; iv + 1];
%             S_Lvxx = [S_Lvxx; val_center; val_west; val_east];
%         end
%     end
% 
%     % Build Lvyy
%     I_Lvyy = [];
%     J_Lvyy = [];
%     S_Lvyy = [];
% 
%     for iv = 1:nv
%         row = ceil(iv / M);
%         col = mod(iv - 1, M) + 1;
% 
%         if row == 1
%             hn = dy(row + 1);
%             hs = dy(row);
%             hc = 0.5 * (hn + hs);
% 
%             % Coefficients
%             val_center = -2 / (hn * hs);
%             val_north = 1 / (hc * hn);
% 
%             % Append indices and values
%             I_Lvyy = [I_Lvyy; iv; iv];
%             J_Lvyy = [J_Lvyy; iv; iv + M];
%             S_Lvyy = [S_Lvyy; val_center; val_north];
%         elseif row == N
%             hs = dy(row);
%             hn = dy(row + 1);
%             hc = 0.5 * (hn + hs);
% 
%             % Coefficients
%             val_center = -2 / (hn * hs);
%             val_south = 1 / (hc * hs) + 1 / (hn * hc);
% 
%             % Append indices and values
%             I_Lvyy = [I_Lvyy; iv; iv];
%             J_Lvyy = [J_Lvyy; iv; iv - M];
%             S_Lvyy = [S_Lvyy; val_center; val_south];
%         else
%             hs = dy(row);
%             hn = dy(row + 1);
%             hc = 0.5 * (hn + hs);
% 
%             % Coefficients
%             val_center = -2 / (hn * hs);
%             val_south = 1 / (hc * hs);
%             val_north = 1 / (hc * hn);
% 
%             % Append indices and values
%             I_Lvyy = [I_Lvyy; iv; iv; iv];
%             J_Lvyy = [J_Lvyy; iv; iv - M; iv + M];
%             S_Lvyy = [S_Lvyy; val_center; val_south; val_north];
%         end
%     end
% 
%     % Build Lv as sparse matrix
%     Lvxx = sparse(I_Lvxx, J_Lvxx, S_Lvxx, nv, nv);
%     Lvyy = sparse(I_Lvyy, J_Lvyy, S_Lvyy, nv, nv);
%     Lv = Lvxx + Lvyy;
% 
%     % Shift indices for v components
%     [I_Lv, J_Lv, S_Lv] = find(Lv);
%     I_Lv_shifted = I_Lv + nu;
%     J_Lv_shifted = J_Lv + nu;
% 
%     % --------------------
%     % Assemble L
%     % --------------------
% 
%     % Get indices and values from Lu
%     [I_Lu, J_Lu, S_Lu] = find(Lu);
% 
%     % Combine indices and values
%     I_L = [I_Lu; I_Lv_shifted];
%     J_L = [J_Lu; J_Lv_shifted];
%     S_L = [S_Lu; S_Lv];
% 
%     % Build L as sparse matrix
%     L = sparse(I_L, J_L, S_L, nq, nq);
% end

% function [LBC] = laplacianRHS(M,N,dx,dy,u,v,q,dt,timeIndex,x,A,k,omega,uFreestream,duFreestreamdx)
%     nu = (M)*N;
%     nv = M*(N);
%     nq = nu+nv;
%     LBC = zeros(nq,1);
%     for iu = 1:M % bot bc
%         hs = dy(1);
%         hn = 0.5*(dy(1)+dy(2));
%         hc = dy(1);
%         LBC(iu) = LBC(iu) + 0;
% %         fprintf('LBC(%3d)=%18.15f\n',iu,LBC(iu));
%     end
%     x = 1.0;
%     for iu = nu:-1:nu-(M)+1 % TOP v component
%         % tangential component
%         idx = iu-nu+(M);
%         x = x - dx(idx);
% %         fprintf('%8.5f %3d %3d\n',x,idx,iu);
%         hn = dy(N);
%         hs = 0.5*(dy(N)+dy(N-1));
%         hc = dy(N);
%         LBC(iu) = LBC(iu) + 2*uFreestream(idx,timeIndex)/(hn*hc);
% 
%         % normal component
%         hc = dy(N);
%         LBC(iu+nu) = LBC(iu+nu) + duFreestreamdx(idx,timeIndex) *(hn+hs)/(hn*hc);
%     end
% 
%     for iu=M:M:nu % right bc u
%         vtprev = v(iu,timeIndex-1);
%         
%         if iu==M
%             vbprev=0;
%         else
%             vbprev = v(iu-M,timeIndex-1);
%         end
%         hc = dx(M);
%         he = dx(M);
%         hw = dx(M);
%         row = floor(iu/M);
%         deltay = dy(row);
%         LBC(iu) = LBC(iu) - 2*hw*(vtprev - vbprev)/deltay/(he*hc);
%     end
% 
% 
%     
% end

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

% function [pRHS] = pressure(M,N, dx, dy, q, u, v, dt, timeIndex, delta, U0, Re, A, k, w,t) % add uf freestream for left inlet
%     
%     nu = (M)*N;
%     nv = M*(N);
%     nq = nu+nv;
%     vprev = v(:,timeIndex-1);   % prev timestep v
%     uprev = u(:,timeIndex-1);   % prev timestep u
%     pRHS  = zeros(nq,1);
%     
%     iUstart = M;
%     iVNW = M; % indexation in vprev starts from 1
%     for iU = iUstart:M:nu
%         
%         iVSW = iVNW-M;
%         iVNWW = iVNW - 1;
%         iVSWW = iVSW - 1;
%         hc = dx(M);
%         hw = 0.5*(dx(M-1)+hc);
%         he = 0.5*hc; % v located at bdry
%         if iVSW==0
%             vInterpS = 0.0;
%         else
%             vInterpS = (hw+he)*(vprev(iVSW) - vprev(iVSWW))/hw + vprev(iVSWW);
%         end
%         vInterpN = (hw+he)*(vprev(iVNW) - vprev(iVNWW))/hw + vprev(iVNWW);
%         vInterp = 0.5*(vInterpS + vInterpN);
%         kineticEnergy = 0.5 * (uprev(iU)*uprev(iU)+vInterp*vInterp);
%         Redudx = (1/Re)*(uprev(iU)-uprev(iU-1))/dx(M);
%         stepfun = 0.5 * (1-tanh(uprev(iU)/(U0 * delta)));
%         pRHS(iU) = Redudx - (kineticEnergy * stepfun);
% %         fprintf(' pRHS(%3d)=%18.15f; VNW: %3d; VNWW: %3d; VSW: %3d; VSWW: %3d\n', iU,pRHS(iU), iVNW,iVNWW, iVSW,iVSWW);
%         iVNW = iVNW + M;
%     end
%     
%     % New code for the top boundary % take uf -> x-mom -> integrate over dx
%     x_coord = 0;
%     ctr = 1;
%     for iV = nv - M + 1 : nv  % Indices for the top boundary in v
%         % Determine grid indices
%         row = ceil(iV / M);   % Should be N
%         col = mod(iV - 1, M) + 1;
%         x_coord = x_coord + dx(ctr)*0.5;
%         pRHS(nu + iV) = -(2*A*(Re*k + Re*w) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2)^3 + 2*A*tan((k*x_coord)/2 + (t*w)/2)^2*(Re*k + Re*w - A*Re*k) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2))/(Re*k*(tan((k*x_coord)/2 + (t*w)/2)^2 + 1)^2);
%         x_coord = x_coord + dx(ctr)*0.5;
%         ctr = ctr+1;
%     end
% 
%     % extrapolate two pressure values from right boundary to top
%     deltay1=(dy(N)+dy(N-1))*0.5;
%     deltay2 = dy(N);
%     pCornerTopFromRightBC = deltay2/deltay1 * (pRHS(nu) - pRHS(nu-M)) + pRHS(nu);
%     pCornerRightFromTopBC = -(2*A*(Re*k + Re*w) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2)^3 + 2*A*tan((k*x_coord)/2 + (t*w)/2)^2*(Re*k + Re*w - A*Re*k) + 2*A*k^2*tan((k*x_coord)/2 + (t*w)/2))/(Re*k*(tan((k*x_coord)/2 + (t*w)/2)^2 + 1)^2);
%     deltap = pCornerTopFromRightBC - pCornerRightFromTopBC;
%     % shift the top values to match the right (solve up to a constant) 
%     % pRHS(nq-M+1:1:nq) = pRHS(nq-M+1:1:nq) + deltap;
%     pRHS(M:M:nu) = pRHS(M:M:nu) - deltap;
% end

function [uInlet] = inletProfile(y,A,w,time,N,ly,phi)
    const=((1) + 2e-18 * (0.5 - rand(1, N)));
    smooth = (A*cos(w*time+phi)*(0.5+0.5*tanh(30*(0.5*(y(2:end)+y(1:end-1))-(N-8)*ly/N)))');
%     uInlet=(1+A*cos(w*time)) + 2e-8 * (0.5 - rand(1, N)) + A*cos(w*time)*(0.5+0.5*tanh(1000*(y-9*ly/10)))';%e*[tanh(smootheningFactor*(0.5*(y(1:end-1)+y(2:end))))];
    uInlet = const + smooth;
    %     uInlet(1) = uInlet(1)/2;
end

function [uFreestream] = uFree(x,Amplitude,k,omega,timeI,phi)
    uFreestream = 1+Amplitude*cos(omega*timeI+k*x(2:end)+phi);
end

function [duFreestream] = duFree(Amplitude,k,omega,time,x,phi)
    duFreestream = -k*Amplitude*sin(omega*time + k*0.5*(x(2:end)+x(1:end-1))+phi); % du/dx
end

function [qp] = generateParticularSolutionFlatPlate(M, N, dxq, dyq, uInlet, timeNow)
    % M: Количество точек по оси x
    % N: Количество точек по оси y
    % dxq, dyq: сеточные шаги
    % uInlet: профиль скорости на входе
    % timeNow: текущий момент времени для временно-зависимых условий

    nu = M * N;  % Количество узлов для горизонтальной скорости (u)
    nq = 2 * nu; % Общее количество узлов (u + v)

    % Инициализация вектора для частного решения
    qp = zeros(nq, 1); 

    % Заполнение горизонтальной скорости (u)
    ctr = 1;
    uparticular = zeros(1, nu);  % массив для хранения горизонтальной скорости
    for i = 1:N
        % Горизонтальная скорость задается с профилем скорости uInlet
        uparticular(ctr:1:ctr+M-1) = uInlet(i, timeNow);
        ctr = ctr + M;  % Сдвигаем счетчик для заполнения следующего слоя
    end
    
    % Горизонтальная компонента скорости (u)
    qp(1:nu) = dyq .* uparticular;  % Умножаем на сеточные шаги по y для корректировки скорости
    
    % Вертикальная компонента скорости (v) равна нулю
    qp(nu+1:nq) = dxq .* 0.0;  % Вертикальная скорость v = 0 для потока вдоль пластины

end

function [px, py] = advectParticles(px, py, u, v, M, N, dt, X, Y)
    % Reshape velocity components for interpolation
    u_matrix = reshape(u(:, end), M, N)';
    v_matrix = reshape(v(:, end), M, N)';
    
    % Define a function to interpolate velocities
    interpVel = @(x, y) deal(interp2(X, Y, u_matrix, x, y, 'linear', 0), interp2(X, Y, v_matrix, x, y, 'linear', 0));
    
    % RK4 coefficients
    [k1x, k1y] = interpVel(px, py);
    [k2x, k2y] = interpVel(px + 0.5 * dt * k1x, py + 0.5 * dt * k1y);
    [k3x, k3y] = interpVel(px + 0.5 * dt * k2x, py + 0.5 * dt * k2y);
    [k4x, k4y] = interpVel(px + dt * k3x, py + dt * k3y);
    
    % Update positions
    px = px + (dt / 6) * (k1x + 2 * k2x + 2 * k3x + k4x);
    py = py + (dt / 6) * (k1y + 2 * k2y + 2 * k3y + k4y);
    
    % Remove particles that move out of the domain
    valid_idx = (px >= min(X(:)) & px <= max(X(:))) & (py >= min(Y(:)) & py <= max(Y(:)));
    px = px(valid_idx);
    py = py(valid_idx);
end

function plotParticles(px, py, X, Y)
    % Plot the particles on the streamlines plot
    scatter(px, py, 2, 'filled', 'MarkerFaceColor', [1, 1, 1]); % White particles
    axis tight;
end

function [max_val, scaled_coord] = findMaxVelocity(u, v, M, N, dx, dy)
    % Combine u and v components
    velocities = [u(:); v(:)];
    
    % Find maximum velocity by absolute value
    [max_val, idx] = max(abs(velocities));
    
    % Determine if the maximum is from u or v
    if idx <= numel(u)
        % The maximum value is from u
        [row, col] = ind2sub([N, M], idx);
    else
        % The maximum value is from v
        idx_v = idx - numel(u);
        [row, col] = ind2sub([N, M], idx_v);
    end
    
    % Calculate the x and y coordinates corresponding to the maximum value
    x_coord = sum(dx(1:col-1)) + dx(col) / 2;
    y_coord = sum(dy(1:row-1)) + dy(row) / 2;
    
    % Scale the coordinates to the [0, 1] domain
    lx = sum(dx);
    ly = sum(dy);
    scaled_coord = [x_coord / lx, y_coord / ly];
end

function Amplitude = smoothAmplitude(iter, targetValue, totalIterations)
    % linearAmplitude: Increases Amplitude linearly from 0 to the target value
    % over the specified number of iterations.
    %
    % Arguments:
    %   iter            - Current iteration number
    %   targetValue     - Target value for Amplitude
    %   totalIterations - Total number of iterations for the linear increase
    %
    % Returns:
    %   Amplitude       - Linearly increasing Amplitude value at iteration 'iter'

    if iter >= totalIterations
        Amplitude = targetValue; % Cap the amplitude at the target value
    else
        Amplitude = (targetValue / totalIterations) * iter;
    end
end




