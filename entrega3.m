% Entrega 3

clc
clear
close all

% Tipo de gráfica y tiempo simulación
tf = input("Tiempo de simulación: ");
dinamica = input("Plot dinámico [0 / 1]: ");
zoomed = input("Plot zoom [0 / 1]: ");
estatica = input("Plot estático [0 / 1]: ");

% Condiciones iniciales
r1 = [0, 0, 0];
r2 = [10.5990201138371, -5.98069207827745, 0];
r3 = [-3.96596970555168, -19.5396345015178, 0];

v1 = [0, 0, 0.1];
v2 = [-6.05028539819850, -11.2054399690566, 0.1];
v3 = [-8.31356419153732, 2.16208426107951, 0.1];

% Constantes del sistema
m1 = 2000;
m2 = 1.8;
m3 = 0.01;

G = 1;

% Ecuación de fuerza
Fg = @(r1, r2, m1, m2) (G * m1 * m2 / (sum((r2 - r1).^2))^(3/2)) * (r2 - r1);

% Ecuaciones energía
T = @(v1, m1) 1/2 * m1 * (sum(v1.^2));
U = @(r1, r2, m1, m2) -G * m1 * m2 / sqrt(sum((r2 - r1).^2));


%% Simulación
dt = 0.01;
iteraciones = floor(tf / dt);

r1_log = zeros(iteraciones, 3);
r2_log = zeros(iteraciones, 3);
r3_log = zeros(iteraciones, 3);

r_centro_log = zeros(iteraciones, 3);

E_t = zeros(iteraciones, 1);

% Desfase Leapfrog
F1 = Fg(r1, r2, m1, m2);
F2 = Fg(r2, r3, m2, m3);
F3 = Fg(r3, r1, m3, m1);
 
a1 = (F1 - F3)/m1;
a2 = (F2 - F1)/m2;
a3 = (F3 - F2)/m3;

v1 = v1 + a1*dt/2;
v2 = v2 + a2*dt/2;
v3 = v3 + a3*dt/2;

% Distancia colisión
dist_c = 0.02;

for i = 1:iteraciones
    t = i*dt;

    r1 = r1 + v1*dt;
    r2 = r2 + v2*dt;
    r3 = r3 + v3*dt;

    r1_log(i, :) = r1;
    r2_log(i, :) = r2;
    r3_log(i, :) = r3;

    F1 = Fg(r1, r2, m1, m2);
    F2 = Fg(r2, r3, m2, m3);
    F3 = Fg(r3, r1, m3, m1);
    
    a1 = (F1 - F3)/m1;
    a2 = (F2 - F1)/m2;
    a3 = (F3 - F2)/m3;

    v1 = v1 + a1*dt;
    v2 = v2 + a2*dt;
    v3 = v3 + a3*dt;

    r_centro = (m1*r1 + m2*r2 + m3*r3) / (m1 + m2 + m3);
    r_centro_log(i, :) = r_centro;

    T_i = T(v1, m1) + T(v2, m2) + T(v3, m3);
    U_i = U(r1, r2, m1, m2) + U(r2, r3, m2, m3) + U(r3, r1, m3, m1);
    E_t(i) = T_i + U_i;

    % Detección de colisión
    dist12 = sqrt(sum((r1 - r2).^2));
    dist23 = sqrt(sum((r2 - r3).^2));
    dist31 = sqrt(sum((r3 - r1).^2));

    if dist12 < dist_c
        fprintf("Colisión entre cuerpo 1 y 2 en tiempo %.4f\n", t)
        r1_log = r1_log(1:i, :);
        r2_log = r2_log(1:i, :);
        r3_log = r3_log(1:i, :);
        r_centro_log = r_centro_log(1:i, :);
        break
    elseif dist23 < dist_c
        fprintf("Colisión entre cuerpo 2 y 3 en tiempo %.4f\n", t)
        r1_log = r1_log(1:i, :);
        r2_log = r2_log(1:i, :);
        r3_log = r3_log(1:i, :);
        r_centro_log = r_centro_log(1:i, :);
        break
    elseif dist31 < dist_c
        fprintf("Colisión entre cuerpo 1 y 3 en tiempo %.4f\n", t)
        r1_log = r1_log(1:i, :);
        r2_log = r2_log(1:i, :);
        r3_log = r3_log(1:i, :);
        r_centro_log = r_centro_log(1:i, :);
        break
    end
end

figure(1)
X = (dt:dt:tf)';
plot(X, E_t)

Emax = max(E_t);
Emin = min(E_t);
Emed = (Emax + Emin) / 2;
ymax = Emax + 0.5;
ymin = Emin - 0.5;

hold on

err = dt^2 * E_t .* X;
plot(X, E_t + err, "r--")
plot(X, E_t - err, "r--")

hold off

axis([dt tf ymin ymax])
title("Energía mecánica total")
legend("Energía calculada", "Error")
xlabel("tiempo")
ylabel("Energía (u. naturales)")


%% Dibujar órbita estática

if estatica == 1
    figure (2)
    title("Trayectorias de las órbitas")
    plot3(r_centro_log(:, 1), r_centro_log(:, 2), r_centro_log(:, 3), "k:")
    hold on
    axis equal
    plot3(r1_log(:, 1), r1_log(:, 2), r1_log(:, 3), "Color", "red", "LineWidth", 0.5)
    plot3(r2_log(:, 1), r2_log(:, 2), r2_log(:, 3), "Color", "green")
    plot3(r3_log(:, 1), r3_log(:, 2), r3_log(:, 3), "Color", "blue")
    hold off
    legend("centro de masa", "Sol", "Júpiter", "asteroide")
end


%% Dibujar órbita dinámica

d0= 2;

if dinamica == 1
    figure(3)
    title("Simulación dinámica")
    
    xmin = min(r2_log(:,1)) - 0.1 * d0;
    xmax = max(r2_log(:,1)) + 0.1 * d0;
    ymin = min(r2_log(:,2)) - 0.1 * d0;
    ymax = max(r2_log(:,2)) + 0.1 * d0;
    zmin = min(r2_log(:,3)) - 0.1 * d0;
    zmax = max(r2_log(:,3)) + 0.1 * d0;
    
    axis equal    
    hold on
    
    axis([xmin xmax ymin ymax zmin zmax])

    points1 = animatedline("Color", "red","Marker","o");
    points2 = animatedline("Color", "green", "Marker", "o");
    points3 = animatedline("Color", "blue","Marker","o");
    orbit1 = animatedline('Color','red');
    orbit2 = animatedline('Color','green');
    orbit3 = animatedline('Color','blue');
    c_masa = animatedline("Color", "black", "LineStyle", ":");

    for i = 1:iteraciones
        addpoints(points1, r1_log(i, 1), r1_log(i, 2), r1_log(i, 3));
        addpoints(points2, r2_log(i, 1), r2_log(i, 2), r2_log(i, 3));
        addpoints(points3, r3_log(i, 1), r3_log(i, 2), r3_log(i, 3));

        addpoints(orbit1, r1_log(i, 1), r1_log(i, 2), r1_log(i, 3));
        addpoints(orbit2, r2_log(i, 1), r2_log(i, 2), r2_log(i, 3));
        addpoints(orbit3, r3_log(i, 1), r3_log(i, 2), r3_log(i, 3));
        addpoints(c_masa, r_centro_log(i, 1), r_centro_log(i, 2), r_centro_log(i, 3))

        drawnow limitrate

        clearpoints(points1)
        clearpoints(points2)
        clearpoints(points3)
    end
    hold off
end


%% Zoooooom

if zoomed == 1
    figure(4)
    title("Simulación con zoom")
    
    axis equal
    hold on

    points1 = animatedline("Color", "red","Marker","o");
    points2 = animatedline("Color", "green", "Marker", "o");
    points3 = animatedline("Color", "blue","Marker","o");
    orbit1 = animatedline('Color','red');
    orbit2 = animatedline('Color','green');
    orbit3 = animatedline('Color','blue');

    for i = 1:iteraciones
        addpoints(points1, r1_log(i, 1), r1_log(i, 2), r1_log(i, 3));
        addpoints(points2, r2_log(i, 1), r2_log(i, 2), r2_log(i, 3));
        addpoints(points3, r3_log(i, 1), r3_log(i, 2), r3_log(i, 3));

        addpoints(orbit1, r1_log(i, 1), r1_log(i, 2), r1_log(i, 3));
        addpoints(orbit2, r2_log(i, 1), r2_log(i, 2), r2_log(i, 3));
        addpoints(orbit3, r3_log(i, 1), r3_log(i, 2), r3_log(i, 3));

        drawnow

        clearpoints(points1)
        clearpoints(points2)
        clearpoints(points3)

        axis([r2_log(i,1)-d0, r2_log(i,1)+d0, r2_log(i,2)-d0,...
            r2_log(i,2)+d0, r2_log(i,3)-d0, r2_log(i,3)+d0])
    end
    hold off
end


