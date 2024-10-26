clc; clear; close all

omega = 1;
k = (2*pi)/0.5;
Amplitude = 0.0;

stepPlot = 500;
M = 15;      % grid size in x-dir
N = 15;      % grid size in y-dir
nt =1000000;     % number of timesteps
dt = 1e-2;  % timestep size
Re = 3500;   % Reynolds number
U0 = 1;     % characteristic velocity scale
delta = 1e-5; % chosen non-dimensional positive constant that is suﬃciently small
filename = ['streamfunction_dim' num2str(sqrt(M*N)) '_it' num2str(nt) '_colored_Re' num2str(Re) '_Amplitude' num2str(Amplitude) '_k' num2str(k) '_omega' num2str(omega) '.gif'];
set(gcf, 'Position', get(0, 'Screensize'));
check_flag = 0;
kx = 1.15;  % grid ratio in x-dir
ky = 1.05;  % grid ratio in y-dir
kx = 1;
ky = 1;
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
[D,G,C] = divGradCurlGen(M,N);
pinvC = pinv(C);
pinvC = sparse(pinvC);
C = sparse(C);
[dxq, dyq] = delQ(M,N,dx,dy);
[u,v,q] = initialConditions(M,N,nt,dyq,dxq);
[MRinv,R] = scalingMatrices(M,N,dx,dy);
MRinv = sparse(MRinv);
R = sparse(R);

% [ADV] = advection(M,N,dx,dy,u,v,2); % test adv

[L] = laplacianGen(M,N,dx,dy);
% L = (L+L')/2;
L = sparse(L);
% L = zeros(nq,nq);
hatA = eye(nq,nq) - dt/Re * L;
A = C'*MRinv*hatA*C;
% A = (A+A')/2;
Rinv = R^(-1);

% qprev = C*psi;
qprev(1:nu) = dyq.*1;
qprev(nu+1:nq) = dxq.*(0.0);
psi = pinv(full(C))*qprev';
MM = MRinv*R;
qp = pinv(D)*boundary2(M, N, dx, dy, q, u, v,dt,1,x,Amplitude,k,omega);    % these two are different
% qp = qprev';
pinvD = sparse(pinv(D));
q(:,1) = qprev;
q(:,2) = qprev;
qprev = qprev';
bc2 = zeros(M*N,nt);
qpadd=qp;
for timeIndex = 4:nt
    [LBC] = laplacianRHS(M,N,dx,dy,u,v,q,dt,timeIndex,x,Amplitude,k,omega);
%     [ADV] = 1.5*advection(M,N,dx,dy,u,v,timeIndex-1)-0.5*advection(M,N,dx,dy,u,v,timeIndex-2);
    [ADV] = advection(M,N,dx,dy,u,v,timeIndex-1);
%     [ADV] = zeros(nq,1);

    % TEST BC ADD
%     qpadd(M:M:nu) = 0;
%     qpadd(nu+nv-M:1:nq) = 0;
%     qpadd(:) = 0;
%     qp(:)=0;

    [pbc] = pressure(M,N, dx, dy, q, u, v, dt, timeIndex, delta, U0, Re);
    B = dt*C'*(MM)*(1/Re*LBC-ADV-pbc) + C'*MRinv*qprev - C'*MM*(hatA)*Rinv * qpadd;
    psi = mldivide(A,B);
%     bc2(:,timeIndex) = outlet(M, N, dx, dy, q, u, v,dt,timeIndex-1,x,Amplitude,k,omega);
%     qp = pinvD*bc2(:,timeIndex);
%     qp(1:nu) = dyq.*1;
%     qp(nu+1:nq) = 0;
    
    
    
    q(:,timeIndex) = C*psi+qp;
    qnow = q(:,timeIndex);
    qprev = q(:,timeIndex);
    u(:,timeIndex) = qnow(1:nu)./dyq';
    v(:,timeIndex) = qnow(nu+1:nq)./dxq';
%     u(u > 1) = 1;
%     v(u > 1) = 1;
    
    if (mod(timeIndex,stepPlot)==0) || (timeIndex==3)
        fprintf('%3d\n',timeIndex);
        tiledlayout(4,6); 
        nexttile(1,[2 2]);
        qh = C*psi;
        uh = qh(1:nu);
        vh = qh(nu+1:nq);
        uploth = reshape(uh,M,N);
        vploth = reshape(vh,M,N);
        psiPlot = reshape(psi,M,N);
        
        contourf(X,flipud(Y),(psiPlot'), 20,'b--'); hold on;
        streamslice(X,flipud(Y),(uploth'),(vploth')); 
%         caxis([-0.1 0.1]);
        colorbar
        title (sprintf('psi contour at %3d''th time step (homogeneous solution)',timeIndex))
        nexttile(13,[2 2]);
%         contourf(reshape(q(1:nu,timeIndex),M-1,N)',  50,'LineColor','none');
        up = qp(1:nu)./dyq';
        vp = qp(nu+1:nq)./dxq';
        uplotp = reshape(up,N,M);
        vplotp = reshape(vp,N,M);
        uplotp = uplotp(1:N,:);
        vplotp = vplotp(:,1:M);
        contourf(X,flipud(Y),(uplotp.^2+vplotp.^2),50,'LineColor','none'); hold on; colorbar;
        sip = streamslice(X,flipud(Y),uplotp,vplotp);
%         sip = quiver(X,Y,uplotp,vplotp);
        axis tight
        set(sip,'Color','red');
%         caxis([0 2]);
        colorbar
        title(sprintf('particular velocity field at %3d''th time step',timeIndex));
        nexttile(10,[3 2]);
        psifinal = pinvC*qprev;
%         plot(q(:,timeIndex));
%         contourf(reshape(qp(1:nu),M-1,N)', 'LineColor','k', 'LevelStep',0.01);
%         title (sprintf('q(1:end) at %3d''th time step',timeIndex))
        uplot = reshape(u(:,timeIndex),M,N)';
        vplot = reshape(v(:,timeIndex),M,N)';
        uplot = uplot(1:N,:);
        vplot = vplot(:,1:M);
        contourf(X,flipud(Y),uplot.^2+vplot.^2,50,'LineColor','none'); hold on;
        VERT = plot(X,Y,'LineWidth',0.0625,"Color", [0.4, 0.4, 0.4, 1]); hold on; 
        HORIZ = plot(X.',Y.','LineWidth',0.0625,"Color", [0.4, 0.4, 0.4, 1]); hold on; 
        colorbar;
        si = streamslice(X,flipud(Y),uplot,vplot);
        set(si,'LineWidth',0.1)
        axis tight
        set(si,'Color','red');
        title(sprintf('total velocity field at %3d''th time step',timeIndex));

%         % plot bc top
        nexttile(3,[1 4]);
%         plot(x,pbc(nq-M:1:nq),'blue');
        hold on;
        plot(x,v(nv-M:1:nv,timeIndex),'red');
        title(sprintf('top velocity BC at %3d',timeIndex));
%         t = dt*timeIndex;
%         utop = 1+Amplitude*(cos(omega*t-k*x));
%         plot(x,utop);
%         ylim([0 2]);
%         title(sprintf('top tangential BC at %3d',timeIndex));
%         % plot bc left
        nexttile(9,[3,1]);
        plot([u(M:M:nu-1,timeIndex)]',y(1:N-1));
        title(sprintf('outlet u component BC at %3d',timeIndex));
%         row = 1;
%         cell = 1;
%         height = 0;
%         velocity1 = 1;
%         acos = Amplitude * cos(omega*t);
%         velocity2 = 1+acos;
%         deltav = (velocity2 - velocity1);
%         uinlet = zeros(1,N);
%         for i = 1:M-1:nu
%             height = height+dy(row);
%             uinlet(row) = 1+0.5*(deltav*(tanh(6*(height-0.5)))+deltav);
%             cell = cell+M;
%             row = row + 1;
%     %         fprintf('%18.15f\n',height);
%         end
%         plot(uinlet',y(2:N+1)); 
%         xlim([0 2]);
%         ylim([0 1]);
%         title(sprintf('inlet normal BC at %3d',timeIndex));
%         % plot bc right
        nexttile(12,[3,1]);
        plot(pbc(M:M:M*N)',y(1:N));
        title(sprintf('outlet pressure BC at %3d',timeIndex));
%         row = 1;
%         uoutlet = zeros(N,1);
%         for i = M:M:M*N
%             uoutlet(row) = -bc2(i,timeIndex)/dy(row);
%             row = row+1;
%         end
%         row = 1;
%         for i = M-1:M-1:(M-1)*N
%             uoutlet(row) = u(i-4,timeIndex);
%             row = row+1;
%         end
%         plot(uoutlet',y(1:N));
%         title(sprintf('outlet normal BC at %3d',timeIndex));
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

bc2 = boundary2(M, N, dx, dy, q, u, v,dt,timeIndex,x,A,k,omega);


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
function [D,G,C] = divGradCurlGen(M,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           DIVERGENCE          %
    %           & GRADIENT          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nu = (M)*N;
    nv = M*(N);
    nq = nu+nv;
    D = zeros(M*N,M*N+M*N);
    nodes = size(D,1);
    links = size(D,2);
    
    
    
    for i = 1:nodes
        D(i,i) = 1;
    end
    
    ctr = 0;
    for i = 1:N % rows
        ctr = ctr + 2;
        D(ctr,ctr-1) = -1;
        for j = 1:M-2 % columns
            ctr = ctr + 1;
            D(ctr,ctr-1) = -1;
        end
    end
    
    for i = 1:nodes
        D(i,M*N+i) = 1;
    end
    
    for i = M+1:nodes
        % fprintf("i=%d.mn+1=%d\n",i,m*(n-1)+i);
        D(i,M*(N-1)+i) = -1;
    end
    G = -D';
    
    C = zeros(links,(M)*(N));
    
    C(1:nu,:) = D(:,nu+1:nq);
    C(nu+1:nq,:) = -D(:,1:nu);
end
function [MRinv,R] = scalingMatrices(M,N,dx,dy)
    [dxq, dyq] = delQ(M,N,dx,dy);
    R = diag([dyq dxq]);
    m = M;
    n = N;
    
    nu = (m)*n;
    nv = m*(n);
    nq = nu + nv;
    
    M1 = zeros(nu, 1);
    M2 = zeros(nv, 1);
    
    MRinv = zeros(nq, nq); % Allocate a zero matrix of size nq x nq
    
    % Fill M1 array
    dx(M+1) = dx(M); % extend to ghost
    for I = 1:nu
        arg = floor((I-1)/(m));
        JJ = floor(arg) + 1;
        II = mod(I-1, m) + 1;
        
        M1(I) = (dx(II)+dx(II+1))/2 / dy(JJ);
    end    
    
    % Fill M2 array
    dy(N+1)=dy(N); % extend to ghost
    for I = 1:nv
        arg = floor((I-1)/m);
        jj = floor(arg) + 2;
        ii = mod((I-1), m) + 1;
%         fprintf('%4d. %4d. %18.15f %18.15f\n',jj-1,ii,dy(jj-1),dy(jj));
        M2(I) = (dy(jj)+dy(jj-1))/2 / dx(ii);
    end
    % Fill diagonal of MM with values from M1
    for I = 1:nu
        MRinv(I, I) = M1(I);
    end
    
    % Fill diagonal of MM with values from M2
    for I = nu+1:nq
        MRinv(I, I) = M2(I-nu);
    end
    
%     dyq = zeros((m-1)*n, 1);
%     for i = 1:n
%         dyq((i-1)*(m-1)+1:i*(m-1)) = repmat(dy(i), m-1, 1);
%     end
%     dxq = zeros((n-1)*m, 1); % Preallocate dxq with zeros
%     index = 1; % Initialize an index variable
%     for i = 1:n-1
%         dxq(index:index+length(dx)-1) = dx; % Fill dxq
%         index = index + length(dx); % Update the index
%     end
%     R = [diag([dyq; dxq])];
%     
%     
%     % M GEN
%     MM = zeros((m-1)*n+m*(n-1),1);
%     size(MM)
%     ctr = 0;
%     rowShift = 0;
%     for j = 1:n
%         for i = 2:m
%             MM((i-1)+rowShift) = 0.5*(dx(i)+dx(i-1));
%             ctr = ctr+1;
%         end
%         % ctr = ctr + 1; % SHIFT BY ONE ROW
%         rowShift = ctr;
%     end
%     rowShift = (m-1)*n;
%     ctr = 0;
%     for i = 2:n
%         for j = 1:m
%             MM(rowShift+(i-2)*(m)+j) = 0.5*(dy(i)+dy(i-1));
%         end
%     end
%     MM = diag(MM);
%     size(MM)
%     size(R)
%     MRinv = M*R^(-1);
 
end
% ^ DONE
function [ADV] = advection(M,N,dx,dy,u,v,timeIndex) % freestream at top hardcoded (uf)
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
            uu(iCell) = 0.25* (0+uprev(ijright))^2; 
        elseif col == M % rightmost cell
            uu(iCell) = 0.25* (uprev(ijleft)+uprev(ijright))^2;
        else
            uu(iCell) = 0.25* (uprev(ijleft)+uprev(ijright))^2; % every cell centre
        end
        
        % compute vv
        if row == 1 % BOT
            vv(iCell) = 0.25*(vprev(ijtop)+0)^2;
        elseif row == N % topmost cell
            vv(iCell) = 0.25*(vprev(ijtop)+vprev(ijbot))^2;
        else
            vv(iCell) = 0.25*(vprev(ijtop)+vprev(ijbot))^2;
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
        if (col==1 && row~=N+1) || (row==1)
            uv(iuv) = 0;
        elseif (col==M+1) && (row~=N+1) % use Dong bc dvtangential/dnormal = 0
            hs = dy(row-1);
            hn = dy(row);
            hw = dx(col-1);
            he = dx(col-1);
            uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijleft))/(he+hw)*(hs*uprev(ijtop) + hn*uprev(ijbot))/(hn+hs);
        elseif (row==N+1) && (col~=M+1) % (uf) freestream hardcoded
            hs = dy(row-1);
            hn = dy(row-1);
            if (col>1)
                hw = dx(col-1); % top bdry case
                uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijright))/(he+hw)*(hs*uprev(ijbot) + hn*uprev(ijbot))/(hn+hs);
            elseif col==1
                hw = dx(col); %  top left corner case
                uv(iuv) = (he*vprev(ijright) + hw*vprev(ijright))/(he+hw)*(hs*uprev(ijbot) + hn*uprev(ijbot))/(hn+hs);
            end
%             % used to be freestream based, 
%             % switched to dong case du/dy=0
%             % at top boundary
%             uf = 1;
%             ug = 2*uf - uprev(ijbot); % 0.5*(ug+ub)=uf; MIGHT NEED TO SET (fix) uTANGENTIAL TOP = UFreestream 
%                                   % since ug+ub = 2 from tangential and
%                                   % ug-ub=0 from dong
%                                   % or just use DONG,
%                                   % without tangential shear velocity
% %             uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijright))/(he+hw)*(hs*ug + hn*uprev(ijbot))/(hn+hs);
% %             uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijleft))/(he+hw)*(hs*uprev(ijtop) + hn*uprev(ijbot))/(hn+hs); % standard
        elseif (col==M+1) && (row==N+1) % mix of use:
            % RIGHT Dong bc dvtangential/dnormal = 0 at right corner 
            % TOP shear tangential for top
            hs = dy(row-1);
            hn = dy(row-1);
            hw = dx(col-1);
            he = dx(col-1);
%             uf = 1;
%             ug = 2*uf - uprev(ijbot); % switched to DONG BC
            uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijleft))/(he+hw)*(hs*uprev(ijbot) + hn*uprev(ijbot))/(hn+hs);
        else
            hs = dy(row-1);
            hn = dy(row);
            hw = dx(col-1);
            he = dx(col);
            uv(iuv) = (he*vprev(ijleft) + hw*vprev(ijright))/(he+hw)*(hs*uprev(ijtop) + hn*uprev(ijbot))/(hn+hs);
        end
%         fprintf('row:%3d, col:%3d, topU:%3d, botU:%3d\t leftV:%3d, rightV=%3d, uv(%3d)=%18.15f\n',row,col,ijtop,ijbot,ijleft,ijright,iuv,uv(iuv));
    end
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
        if col==M
            if row==1 % no v at bottom, v=0
                ug = 2*deltax/deltay*(-v(iu,timeIndex-1) + 0) + u(iu-1,timeIndex);
            else % use Dong du/dx=-dv/dy
                ug = 2*deltax/deltay*(-v(iu,timeIndex-1) + v(iu-M,timeIndex-1)) + u(iu-1,timeIndex);
            end
            uuright = 0.5*(ug+uprev(iu));
            duudx(iu) = (uuright - uu(leftcell))/(deltax);
        else
            duudx(iu) = (uu(rightcell) - uu(leftcell))/(deltax);
%             fprintf('iu=%3d, row:%3d, col:%3d, leftcell=%3d, rightcell=%3d, botnode=%3d, topnode=%3d, uuright=%8.4f, uuleft=%8.4f\n', ...
%             iu,row,col,leftcell,rightcell,botnode,topnode, uu(rightcell), uu(leftcell));
        end
        duvdy(iu) = (uv(topnode) - uv(botnode))/(deltay);
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
        dvudx(iv) = (uv(rightnode) - uv(leftnode))/deltax;
        if row==N
            if col==1
                vg = 2*deltay/deltax*(-u(iu,timeIndex-1) + 0) + v(iu-M,timeIndex);
            else
                vg = 2*deltay/deltax*(-u(iu,timeIndex-1) + u(iu-1,timeIndex-1)) + v(iu-M,timeIndex);
            end
            vvtop = 0.5*(vg+vprev(iv));
            dvvdy(iv) = (vvtop - vv(botcell))/deltay;
%             fprintf('\t %3d. vvtop = %18.15f, vvbot = %18.15f, dy=%18.15f, vg = %18.15f\n', iv, vvtop,vv(botcell),deltay,vg);
        else
            dvvdy(iv) = (vv(topcell) - vv(botcell))/deltay;
%             fprintf('\t vvtop = %18.15f, vvbot = %18.15f, dy=%18.15f\n', vv(topcell),vv(botcell),deltay);
        end
%         fprintf('\t dvvdy(%3d)=%18.15f\n',iv,dvvdy(iv));
    end

%     for i = 1:M
%         iV = nv - M + i;  % Index for v at top boundary
%         iU = nu - M + i;  % Index for u at top boundary
%         
%         % For v-component
%         if i < M
%             vN = vprev(iV);
%             vS = vprev(iV - M);
%             dudyV = (vN - vS) / dy(N);
%             ADV(iV) = -uprev(iU) * dudyV;
%         end
%         
%         % For u-component
%         if i > 1
%             uN = uprev(iU);
%             uS = uprev(iU - M);
%             dudyU = (uN - uS) / dy(N);
%             ADV(iU) = -uprev(iU) * dudyU;
%         end
%     end



    advU = duudx+duvdy;
    advV = dvudx+dvvdy;
    ADV = [advU; advV];

%     % Top boundary treatment
%     u0 = 1; % Freestream velocity
%     for i = 1:M
%         iV = nv - M + i;
%         iU = nu - M + i;
%         
%         % For u-component (tangential)
%         if i > 1 && i < M
%             uE = uprev(iU);
%             uW = uprev(iU-1);
%             dudx = (uE - uW) / (2*dx(i));
%             ADV(iU) = -u0 * dudx;
%         end
%         
%         % For v-component (normal)
%         if i < M
%             vN = vprev(iV);
%             vS = vprev(iV - M);
%             dvdy = (vN - vS) / dy(N);
%             ADV(iV) = -max(vN, 0) * dvdy; % Outflow-favoring upwind
%         end
%     end
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
            Luxx(iu,iu) = 1/(hw*hc); % new DONG BC % du/dx = -dv/dy ugnp1=dx/dy(vt-vb)n +ucnp1
            Luxx(iu,iu) = -2/(hw*he)+1/(hw*hc); %du/dx=-dv/dy
            Luxx(iu,iu-1) = 1/(hw*hc); % west coeff is 1/hw / hc
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
            Luyy(iu,iu) = -2/(hs*hn)+1/(hn*hc); % du/dy=0 ug=uc
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
            Lvyy(iv,iv) = -2/(hc*hs) + 1/(hc*hn); % dv/dy=0 vg = vc (at bdry)
            Lvyy(iv,iv-M) = 1/(hc*hs);


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
%     spy(L);
end
function [LBC] = laplacianRHS(M,N,dx,dy,u,v,q,dt,timeIndex,x,A,k,omega)
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
    t = dt*timeIndex;
    uf = zeros(1,M);
    for iu = nu:-1:nu-(M)+1 % TOP
        % tangential component
        idx = iu-nu+(M);
        x = x - dx(idx);
%         fprintf('%8.5f %3d %3d\n',x,idx,iu);
        hn = dy(N);
        hs = 0.5*(dy(N)+dy(N-1));
        hc = dy(N);
        acos = A * cos(omega*t - k*x);
        uf(idx) = 1.0+acos;

        % normal component
        hc = dy(N);
        deltax = dx(idx);
        urprev = u(iu,timeIndex-1);
        if iu == nu-M+1
            ulprev = 0;
        else
            ulprev = u(iu-1,timeIndex-1);
        end
        LBC(iu+nu) = LBC(iu+nu) - (urprev-ulprev)/deltax/hc;
    end

    for iu=M:M:nu % right bc u
        vtprev = v(iu,timeIndex-1);
        
        if iu==M
            vbprev=0;
        else
            vbprev = v(iu-M,timeIndex-1);
        end
        hc = dx(M);
        row = floor(iu/M);
        deltay = dy(row);
        LBC(iu) = LBC(iu) - (vtprev - vbprev)/deltay/hc;
    end


    
end
% Add this function to compute bc2
function [bc2] = boundary2(M, N, dx, dy, q, u, v,dt,timeIndex,x,A,k,omega)
    nu = (M-1)*N;
    nv = M*(N-1);
    nq = nu+nv;
    bc2 = zeros(M*N,1);
    
    row = 1;
    cell = 1;
    t = dt*timeIndex;
    height = 0;
    velocity1 = 1;
    acos = A * cos(omega*t);
    velocity2 = 1+acos;
    deltav = (velocity2 - velocity1);
    for i = 1:M:M*N
        height = height+dy(row);
        uinlet = 1+0.5*(deltav*(tanh(6*(height-0.5)))+deltav);
        uinlet = 1;
        bc2(cell) = uinlet*dy(row);
        cell = cell+M;
        row = row + 1;
%         fprintf('%4d\n',i);
    end
end
function [pRHS] = pressure(M,N, dx, dy, q, u, v, dt, timeIndex, delta, U0, Re) % add uf freestream for left inlet
    
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
    
    iVstart = nv-M+1;
    iUSE = iVstart;
    for iV = iVstart:1:nv
        iUSW = iUSE - 1;
        iUSSW = iUSW - M;
        iUSSE = iUSE - M;
        hc = dy(N);
        hs = 0.5*(dy(N-1)+hc);
        hn = 0.5*hc; % v located at bdry
        if iVSW==0
            uInterpW = 1.0;
        else
            uInterpW = (hs+hn)*(uprev(iUSW) - uprev(iUSSW))/hs + uprev(iUSSW);
        end
%         vInterpN = (hw+he)*(vprev(iVNW) - vprev(iVNWW))/hw + vprev(iVNWW);
        uInterpE = (hs+hn)*(uprev(iUSE) - uprev(iUSSE))/hs + uprev(iUSSE);
        uInterp = 0.5*(uInterpE + uInterpW);
        kineticEnergy = 0.5 * (vprev(iV)*vprev(iV)+uInterp*uInterp);
%         if iV == 781-400
%             uInterp
%         end
        Redudx = (1/Re)*(vprev(iV)-vprev(iV-M))/dy(N);
        stepfun = 0.5 * (1-tanh(vprev(iV)/(U0 * delta)));
        pRHS(iV+nu) = Redudx - (kineticEnergy * stepfun);
%         fprintf(' pRHS(%3d)=%18.15f; USW: %3d; USSW: %3d; USE: %3d; USSE: %3d\n', iV+nu,pRHS(iV+nu), iUSW,iUSSW, iUSE,iUSSE);
        iUSE = iUSE + 1;
    end
% 
    deltap = pRHS(nq) - pRHS(nu);
%     pRHS(M:M:nu) = pRHS(M:M:nu) + deltap;
    pRHS(nq-M:1:nq) = pRHS(nq-M:1:nq) + deltap;
end