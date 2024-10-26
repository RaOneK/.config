% Параметры задачи
Lx = 1;               % Длина области по x
Ly = 1;               % Длина области по y
Nx = 100;             % Количество пространственных шагов по x
Ny = 100;             % Количество пространственных шагов по y
x = linspace(0, Lx, Nx); % Координаты сетки по x
y = linspace(0, Ly, Ny); % Координаты сетки по y

T = 10 * pi;           % Время моделирования
Nt = 10000;            % Количество временных шагов
t = linspace(0, T, Nt); % Моменты времени

c = 1;                % Скорость адвекции
dx = Lx / (Nx - 1);   % Пространственный шаг по x
dy = Ly / (Ny - 1);   % Пространственный шаг по y
dt = T / (Nt - 1);    % Временной шаг

% Проверка условия устойчивости CFL
CFL_x = c * dt / dx;
CFL_y = c * dt / dy;
if CFL_x > 1 || CFL_y > 1
    error('Нарушено условие Куранта: CFL_x = %f, CFL_y = %f > 1', CFL_x, CFL_y);
end

% Начальные условия
u = zeros(Nx, Ny);

% Основной цикл по времени
for n = 1:Nt
    % Левая граница: v(0, t) = 1 + sin(t)
    u(1, :) = 1 + sin(t(n));
    
    % Явная схема вперёд по времени, назад по пространству
    for i = 2:Nx
        for j = 2:Ny
            u(i, j) = u(i, j) - CFL_x * (u(i, j) - u(i-1, j)) - CFL_y * (u(i, j) - u(i, j-1));
        end
    end
    
    % Правая граница: du/dx = 0 и du/dy = 0 (неявно сохраняем значения на границах)
    u(Nx, :) = u(Nx-1, :);
    u(:, Ny) = u(:, Ny-1);
    
%     % Визуализация
%     hold off;
%     [X, Y] = meshgrid(x, y);
%     streamslice(X, Y, u', zeros(size(u')));
%     axis([0 Lx 0 Ly]);
%     title(sprintf('Время t = %f', t(n)));
%     xlabel('x'); ylabel('y'); zlabel('u');
%     drawnow;

     % Визуализация
    if mod(n,100)==0
%         surf(x, y, u');
%         axis([0 Lx 0 Ly -0.5 2.5]);
%         title(sprintf('Время t = %f', t(n)));
%         xlabel('x'); ylabel('y'); zlabel('u');
%         drawnow;
        [X, Y] = meshgrid(x, y);
        streamslice(X, Y, u', zeros(size(u')));
        axis([0 Lx 0 Ly]);
        title(sprintf('Время t = %f', t(n)));
        xlabel('x'); ylabel('y'); zlabel('u');
        drawnow;
    end
end