%% Animazione braccio 2 GDL da dati marker su CSV

clear; clc; close all

% Parametri (devono essere noti, non ricavati dai dati)
l1 = 0.3;   % Lunghezza omero [m]
l2 = 0.25;  % Lunghezza avambraccio [m]

% Parametri di animazione
anim_speed = 0.1;      % pausa tra i frame [s] (aumenta per rallentare)
stop_time = 10;       % istante di tempo (in secondi) a cui fermare l'animazione

% Caricamento dati marker (usa il percorso relativo corretto)
data = readtable('dati_marker_2.csv');
nframe = height(data);

% Preleva marker
x1 = data.x1; y1 = data.y1; % Spalla
x2 = data.x2; y2 = data.y2; % Gomito
x3 = data.x3; y3 = data.y3; % Mano
t = data.time;

figure('Name','Animazione 2 GDL: marker da CSV, ellisse forza e diagonale maggiore');
axis equal; grid on; hold on
xlabel('x [m]'); ylabel('y [m]');
title('Braccio planare 2 GDL: ellisse forza e diagonale maggiore');
xlim([-0.3 0.3]);
ylim([-0.05 l1+l2+0.1]);

for k = 1:nframe
    if t(k) > stop_time
        disp(['Animazione stoppata a t = ', num2str(t(k)), ' s (frame ', num2str(k), ')']);
        break
    end

    cla

    % Marker
    S = [x1(k); y1(k)];
    G = [x2(k); y2(k)];
    H = [x3(k); y3(k)];

    % Ricavo angoli articolari (cinematica inversa)
    theta1 = atan2(y2(k)-y1(k), x2(k)-x1(k));
    theta2 = atan2(y3(k)-y2(k), x3(k)-x2(k)) - theta1;

    % Jacobiano rispetto al punto H
    J = [-l1*sin(theta1) - l2*sin(theta1+theta2), -l2*sin(theta1+theta2);
          l1*cos(theta1) + l2*cos(theta1+theta2),  l2*cos(theta1+theta2)];

    % Ellisse della forza trasmissibile
    phi = linspace(0, 2*pi, 200);
    manipulability = sqrtm(inv(J*J'));
    ellisse = manipulability * [cos(phi); sin(phi)];
    ellisse = 0.08 * ellisse;

    % Calcolo degli assi principali (autovalori, autovettori)
    [V,D] = eig(manipulability*manipulability');
    [~, idx] = max(diag(D));
    v_max = V(:,idx);
    len_max = 0.08 * sqrt(D(idx,idx));
    pt1 = H + len_max * v_max;
    pt2 = H - len_max * v_max;

    % Omero
    plot([S(1), G(1)], [S(2), G(2)], 'k-', 'LineWidth', 4);
    % Avambraccio
    plot([G(1), H(1)], [G(2), H(2)], 'k-', 'LineWidth', 4);
    % Giunti
    plot(S(1), S(2), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    plot(G(1), G(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    plot(H(1), H(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

    % Ellisse trasmissibilità forza
    plot(H(1) + ellisse(1,:), H(2) + ellisse(2,:), 'm-', 'LineWidth', 2);

    % Diagonale maggiore (asse principale) in rosso
    plot([pt1(1) pt2(1)], [pt1(2) pt2(2)], 'r-', 'LineWidth', 3);

    % Etichette statiche
    text(S(1), S(2), '  Spalla', 'FontWeight', 'bold', 'FontSize', 10, 'VerticalAlignment', 'bottom');
    text(G(1), G(2), '  Gomito', 'FontWeight', 'bold', 'FontSize', 10, 'VerticalAlignment', 'top');
    text(H(1), H(2), '  Mano', 'FontWeight', 'bold', 'FontSize', 10, 'VerticalAlignment', 'bottom');

    % Mostra angoli attuali
    angolo1_gradi = rad2deg(theta1);
    angolo2_gradi = rad2deg(theta2);
    txt1 = sprintf('\\theta_1 = %.1f^\\circ', angolo1_gradi);
    txt2 = sprintf('\\theta_2 = %.1f^\\circ', angolo2_gradi);
    text(-0.28, l1+l2+0.10, txt1, 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'b');
    text(-0.28, l1+l2+0.05, txt2, 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'r');

    pause(anim_speed)
end