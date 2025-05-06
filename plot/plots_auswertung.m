% clear;

% load 'websocket_vicon_5x_v2500_100hz_abb_ereignisse.mat'

%% Farben definieren

% Für Bahnvergleich
blau = [0 0.4470 0.7410]; % Standard Blau
rot = [0.78 0 0]; % gewählter Rotton in Masterarbeit

% Für Plots Verfahren
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.9290 0.6940 0.1250];
c4 = [0.4940 0.1840 0.5560];
c5 = [0.4660 0.6740 0.1880];
c6 = [0.3010 0.7450 0.9330];
c7 = [0.6350 0.0780 0.1840];

tab1 = 0;
tab2 = 0;

plotiso = 0;
ploteckgerade = 0;
plotabweichungen = 0;

% if tab1
%% Plots für Abweichungen bei Segmenten und Trajektorien
%%%%% PLOTS FÜR WEBSOCKET UND VICON

for i = 1:size(trajectories_ist,2) % size(segments_ist,2)

    %%%%%%%%%%% TRAJEKTORIEN
    xt_dists_av_dtw(i) = struct_dtw{i}.dtw_average_distance*1000;
    xt_dists_av_frechet(i) = struct_frechet{i}.frechet_average_distance *1000;
    xt_dists_av_lcss(i) = struct_lcss{i}.lcss_average_distance*1000;
    xt_dists_av_sidtw(i) = struct_sidtw{i}.dtw_average_distance*1000;
    xt_dists_av_eucl(i) = struct_euclidean{i}.euclidean_average_distance*1000;

    xt_dists_max_dtw(i) = struct_dtw{i}.dtw_max_distance*1000;
    xt_dists_max_frechet(i) = struct_frechet{i}.frechet_max_distance *1000;
    xt_dists_max_lcss(i) = struct_lcss{i}.lcss_max_distance*1000;
    xt_dists_max_sidtw(i) = struct_sidtw{i}.dtw_max_distance*1000;
    xt_dists_max_eucl(i) = struct_euclidean{i}.euclidean_max_distance*1000;
end
%%
for i = 1:size(segments_ist,2) 
    %%%%%%%%%% SEGMETE
    xs_dists_av_dtw(i) = struct_dtw_segments{i}.dtw_average_distance*1000;
    xs_dists_av_frechet(i) = struct_frechet_segments{i}.frechet_average_distance *1000;
    xs_dists_av_lcss(i) = struct_lcss_segments{i}.lcss_average_distance*1000;
    xs_dists_av_sidtw(i) = struct_sidtw_segments{i}.dtw_average_distance*1000;
    xs_dists_av_eucl(i) = struct_euclidean_segments{i}.euclidean_average_distance*1000;

    xs_dists_max_dtw(i) = struct_dtw_segments{i}.dtw_max_distance*1000;
    xs_dists_max_frechet(i) = struct_frechet_segments{i}.frechet_max_distance *1000;
    xs_dists_max_lcss(i) = struct_lcss_segments{i}.lcss_max_distance*1000;
    xs_dists_max_sidtw(i) = struct_sidtw_segments{i}.dtw_max_distance*1000;
    xs_dists_max_eucl(i) = struct_euclidean_segments{i}.euclidean_max_distance*1000;
end

%% Mittel und Max Durchschnittswerte 
yy_max_dtw = mean(xt_dists_max_dtw);
yy_max_frechet = mean(xt_dists_max_frechet);
yy_max_lcss = mean(xt_dists_max_lcss);
yy_max_sidtw = mean(xt_dists_max_sidtw);
yy_max_eucl = mean(xt_dists_max_eucl);

% yy_mean_dtw = mean(y_dists_av_dtw);
% yy_mean_frechet = mean(y_dists_av_frechet);
% yy_mean_lcss = mean(y_dists_av_lcss);
% yy_mean_sidtw = mean(y_dists_av_sidtw);
% yy_mean_eucl = mean(y_dists_av_eucl);

% Für Plot falls Linie gewünscht
% a = 1:length(y_dists_max_dtw);
% b = yy_max_dtw*ones(length(a),1);

%% 
figure('Color','white');
hold on
% ylim([0 max(xt_dists_av_frechet)]);
% Segmente
plot(xs_dists_av_dtw,LineWidth=1)
plot(xs_dists_av_frechet,LineWidth=1)
plot(xs_dists_av_lcss,LineWidth=1)
plot(xs_dists_av_sidtw,LineWidth=1)
plot(xs_dists_av_eucl,LineWidth=1)
% Trajektorien
% plot(xt_dists_av_dtw,LineWidth=2.5,Color=c1)
% plot(xt_dists_av_frechet,LineWidth=2.5,Color=c2)
% plot(xt_dists_av_lcss,LineWidth=2.5,Color=c3)
% plot(xt_dists_av_sidtw,LineWidth=2.5,Color=c4)
% plot(xt_dists_av_eucl,LineWidth=2.5,Color=c5)
% % Linie die Trajektorien von Segmenten trennt

xlabel('Bahnsegmente (Trajektorien)');
ylabel('Abweichung in mm');
% legend('DTW','DFD','LCSS','SIDTW','Eukl. Dist.')
xline(size(trajectories_ist,2))
% axis equal
axis tight
hold off

% end

%% Berechnung der mitteleren und max Abweichung für Trajektorienvergleich

% if tab2
for i = 1:size(trajectories_ist,2)-1 % size(segments_ist,2)

    %%%%%%%%%%%TRAJEKTORIEN
    xt_dists_av_dtw(i) = struct_dtw{i}.dtw_average_distance*1000;
    xt_dists_av_sidtw(i) = struct_sidtw{i}.dtw_average_distance*1000;
    xt_dists_av_eucl(i) = struct_euclidean{i}.euclidean_average_distance*1000;
    % xt_dists_av_lcss(i) = struct_lcss{i}.lcss_average_distance*1000;

    xt_dists_max_dtw(i) = struct_dtw{i}.dtw_max_distance*1000;
    xt_dists_max_sidtw(i) = struct_sidtw{i}.dtw_max_distance*1000;
    xt_dists_max_eucl(i) = struct_euclidean{i}.euclidean_max_distance*1000;
    % xt_dists_max_lcss(i) = struct_lcss{i}.lcss_max_distance*1000;
end

% %%%%% Mittelwerte bilden
% 
yy_mean_dtw = mean(xt_dists_av_dtw)
yy_mean_sidtw = mean(xt_dists_av_sidtw)
yy_mean_eucl = mean(xt_dists_av_eucl)
% yy_mean_lcss = mean(xt_dists_av_lcss)

xx_max_dtw = mean(xt_dists_max_dtw)
xx_max_sidtw = mean(xt_dists_max_sidtw)
xx_max_eucl = mean(xt_dists_max_eucl)
% xx_max_lcss = mean(xt_dists_max_lcss)

% end


%% Plotten der Bewegungsroutine ISO-Würfel

if plotiso

figure('Color','white');
hold on
% plot3(trajectories_ist{1}(:,2),trajectories_ist{1}(:,3),trajectories_ist{1}(:,4),Color= c1,LineWidth=1.5)
% plot3(trajectories_ist{2}(:,2),trajectories_ist{2}(:,3),trajectories_ist{2}(:,4),Color= c1,LineWidth=1.5)
% plot3(trajectories_ist{3}(:,2),trajectories_ist{3}(:,3),trajectories_ist{3}(:,4),Color= c2,LineWidth=1.5)
% plot3(trajectories_ist{4}(:,2),trajectories_ist{4}(:,3),trajectories_ist{4}(:,4),Color= c3,LineWidth=1.5)
% plot3(trajectories_ist{5}(:,2),trajectories_ist{5}(:,3),trajectories_ist{5}(:,4),Color= c4,LineWidth=1.5)
plot3(trajectories_ist{6}(:,2),trajectories_ist{6}(:,3),trajectories_ist{6}(:,4),Color= c5,LineWidth=1.5)
% plot3(trajectories_ist{7}(:,2),trajectories_ist{7}(:,3),trajectories_ist{7}(:,4),Color= c6,LineWidth=1.5)
view(3)
grid on
% axis equal 
axis padded
xlabel('\textbf{x}'); ylabel('\textbf{y}'); zlabel('\textbf{z}');

end
%% Plotten der Eckpunkte mit Zuordnung

% figure('Color','white');
% hold on
% plot3(trajectories_ist{3}(:,2),trajectories_ist{3}(:,3),trajectories_ist{3}(:,4),Color= blau,LineWidth=1.5)
% plot3(trajectories_soll{3}(:,1),trajectories_soll{3}(:,2),trajectories_soll{3}(:,3),Color= rot,LineWidth=1.5)
% xlabel('\textbfx'); ylabel('\textbfy'); zlabel('\textbfz');
% view(90, 0);
% % view(2)
% axis padded

%% PLOTS DER ECK PUNKT ZUORDNUNG
% % Ist = trajectories_ist{6};
% % Soll = trajectories_soll{6};
% 
soll = struct_sidtw_segments{7}.dtw_X*1000;
ist = struct_sidtw_segments{7}.dtw_Y*1000;
N = length(struct_sidtw_segments{7}.dtw_path);

% % Indizes der zu löschenden Zeilen ermitteln
% rows_to_delete = 2:2:size(ist, 1);
% 
% % Die entsprechenden Zeilen löschen
% ist(rows_to_delete, :) = [];
% soll(rows_to_delete,:) = [];

figure('Color','white');
hold on
plot3(ist(:,1),ist(:,2),ist(:,3),Color= c1,LineWidth=1.5,Marker = "o",MarkerFaceColor= c1)
plot3(soll(:,1),soll(:,2),soll(:,3),Color= rot,LineWidth=1.5,Marker = "square",MarkerFaceColor=rot)

% % Darstellung der Punkte und Verbindungslinien
%     % Punkte
%     plot3(Ist(:,1), Ist(:,2), Ist(:,3), '-bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
%     plot3(Soll(:,1), Soll(:,2), Soll(:,3), '-ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% 
% 
for k = 1:N 
    % Verbindungslinien
    plot3([ist(k,1), soll(k,1)], [ist(k,2), soll(k,2)], [ist(k,3), soll(k,3)], 'k-');
end
% 
title('DTW');
% xlim([1376.5 1381]);
% ylim([-931 -926]);
% xlabel('\textbfx'); ylabel('\textbfy'); zlabel('\textbfz');
% % view(90, 0);
% view(2)
% % axis padded
%%%%%%% Vergleich DTW DFD 3D
% xlim([1372 1381]);
% ylim([-308 -299]);
xlabel('\textbfx'); ylabel('\textbfy'); zlabel('\textbfz');
grid on
% 
% % legend({'Istbahn', 'Sollbahn', 'Zuordnung'}, 'Location', 'best');

%% PLOTTEN DER ECKEN UND GERADEN UND ANWENDUNG DER VEFAHREN

% if ploteckgerade
%%%% ALLE VERFAHREN
% soll = struct_frechet{3}.dtw_X*1000;
% ist = struct_sidtw{3}.dtw_Y*1000;
N = length(struct_frechet{3}.frechet_path);

%%% FÜR FRECHET NOTWENIDIG
ix = struct_frechet{6}.frechet_path(1,:);
iy = struct_frechet{6}.frechet_path(2,:);

soll = trajectories_soll{6}(ix,:);
ist = trajectories_ist{6}(iy,2:4);

%%%%%% EUKLIDEAN
% soll = trajectories_soll{6};
% ist = trajectories_ist{6}(:,2:4);
% 
% intersections = struct_euclidean{6}.euclidean_intersections;
% numPoints = length(soll);
% 
% [xy,eucl_distances,~] = distance2curve(soll,ist,'linear');

figure('Color','white');
hold on
plot3(ist(:,1),ist(:,2),ist(:,3),Color= c1,LineWidth=1.5,Marker = "o",MarkerFaceColor= c1)
plot3(soll(:,1),soll(:,2),soll(:,3),Color= rot,LineWidth=1.5,Marker = "square",MarkerFaceColor=rot)
for k = 1:N 
    % Verbindungslinien
    plot3([ist(k,1), soll(k,1)], [ist(k,2), soll(k,2)], [ist(k,3), soll(k,3)], 'k-');
end
%%%%%%%%%% Nur bei Euklidean
% for i = 1:numPoints
%     % Zeichnen der Linie
%     line([ist(i,1), xy(i,1)], [ist(i,2), xy(i,2)], [ist(i,3), xy(i,3)], 'color', 'black');
% 
%     % Speichern der Koordinaten in der Matrix
%     lines(i, 1:3) = ist(i, :); % Startpunkt (x, y, z)
%     lines(i, 4:6) = xy(i, :);    % Endpunkt (x, y, z)
% end

%%%%%% ECKE

% xlim([1378 1381]);
% ylim([-931 -926]);
% set(gcf, 'Position', [100, 100, 300, 400]);

%%%%% GERADE
% xlim([1379.5 1381.5]);
% ylim([-642 -585]);
% set(gcf, 'Position', [100, 100, 200, 400]);


title('DFD');
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
xlabel('\textbfx'); ylabel('\textbfy'); zlabel('\textbfz');
grid on

%%%%%%% Vergleich DTW DFD 3D
xlim([1372 1381]);
ylim([-308 -299]);

% view(90, 0);
% pbaspect([1 2 1]); % Verhältnis 1:2 für x:y
view(2)
% axis padded

% end

%% ABB 15 x v1000 und v2000

if plotabweichungen
% Anzahl der Einträge im struct-Array
numEntries = length(struct_euclidean);
% numEntries = 3;

% Maximale Länge der dtw_distances-Vektoren herausfinden
maxLen = 0;
for i = 1:numEntries
    maxLen = max(maxLen, length(struct_euclidean{i}.euclidean_distances));
end

% Matrix initialisieren mit NaN-Werten für ungleiche Längen
xeucl_matrix = NaN(numEntries, maxLen);

% dtw_distances in die Matrix kopieren
for i = 1:numEntries
    xeucl = struct_euclidean{i}.euclidean_distances;
    xeucl_matrix(i, 1:length(xeucl)) = xeucl;
end

% Mittelwert über die Zeilen berechnen, unter Berücksichtigung von NaN-Werten
xmean_eucl = nanmean(xeucl_matrix, 1)*1000;


% Euklidische Distanzen plotten

figure('Color','white');
title('Eukl. Dist')
xlabel("Bahnpunkte"); ylabel("Abweichung in mm")
hold on
for i = 1:size(struct_euclidean,2)
    plot(struct_euclidean{i}.euclidean_distances*1000,Color=c1,LineWidth=0.5)
end
plot(xmean_eucl,Color=c2,LineWidth=1.5)
hold off

%%
% Anzahl der Einträge im struct-Array
numEntries = length(struct_dtw);
% numEntries = 3;

% Maximale Länge der dtw_distances-Vektoren herausfinden
maxLen = 0;
for i = 1:numEntries
    maxLen = max(maxLen, length(struct_dtw{i}.dtw_distances));
end

% Matrix initialisieren mit NaN-Werten für ungleiche Längen
xdtw_matrix = NaN(numEntries, maxLen);

% dtw_distances in die Matrix kopieren
for i = 1:numEntries
    xdtw = struct_dtw{i}.dtw_distances;
    xdtw_matrix(i, 1:length(xdtw)) = xdtw;
end

% Mittelwert über die Zeilen berechnen, unter Berücksichtigung von NaN-Werten
xmean_dtw = nanmean(xdtw_matrix, 1);

% Ergebnis anzeigen
xmean_dtw = xmean_dtw*1000';


figure('Color','white');
title('DTW')
xlabel("Bahnpunkte"); ylabel("Abweichung in mm")
hold on
for i = 1:numEntries
    plot(struct_dtw{i}.dtw_distances*1000,Color=c1,LineWidth=0.5)
end
plot(xmean_dtw,Color=c2,LineWidth=1.5)
%% SIDTW
% Anzahl der Einträge im struct-Array
numEntries = length(struct_sidtw);
% numEntries = 3;

% Maximale Länge der dtw_distances-Vektoren herausfinden
maxLen = 0;
for i = 1:numEntries
    maxLen = max(maxLen, length(struct_sidtw{i}.dtw_distances));
end

% Matrix initialisieren mit NaN-Werten für ungleiche Längen
xdtw_matrix = NaN(numEntries, maxLen);

% dtw_distances in die Matrix kopieren
for i = 1:numEntries
    xsidtw = struct_sidtw{i}.dtw_distances;
    xsidtw_matrix(i, 1:length(xsidtw)) = xsidtw;
end

% Mittelwert über die Zeilen berechnen, unter Berücksichtigung von NaN-Werten
xmean_sidtw = nanmean(xsidtw_matrix, 1)*1000;


figure('Color','white');
title('SIDTW')
xlabel("Bahnpunkte"); ylabel("Abweichung in mm")
hold on
for i = 1:numEntries
    plot(struct_sidtw{i}.dtw_distances*1000,Color=c1,LineWidth=0.5)
end
plot(xmean_sidtw,Color=c2,LineWidth=1.5)

%% LCSS

% Anzahl der Einträge im struct-Array
numEntries = length(struct_lcss);
% numEntries = 3;

% Maximale Länge der dtw_distances-Vektoren herausfinden
maxLen = 0;
for i = 1:numEntries
    maxLen = max(maxLen, length(struct_lcss{i}.lcss_distances));
end

% Matrix initialisieren mit NaN-Werten für ungleiche Längen
xlcss_matrix = NaN(numEntries, maxLen);

% dtw_distances in die Matrix kopieren
for i = 1:numEntries
    xlcss = struct_lcss{i}.lcss_distances;
    xlcss_matrix(i, 1:length(xlcss)) = xlcss;
end

% Mittelwert über die Zeilen berechnen, unter Berücksichtigung von NaN-Werten
xmean_lcss = nanmean(xlcss_matrix, 1)*1000;


figure('Color','white');
title('LCSS')
xlabel("Bahnpunkte"); ylabel("Abweichung in mm")
hold on
for i = 1:numEntries
    plot(struct_lcss{i}.lcss_distances*1000,Color=c1,LineWidth=0.5)
end
plot(xmean_lcss,Color=c2,LineWidth=1.5)

end
%% Standardabweichungen



% y_std_dtw = std(xmean_dtw);
% y_std_sidtw = std(xmean_sidtw);
% y_std_eucl = std(xmean_eucl);

av_dtw = mean(xt_dists_av_dtw)

av_frechet = mean(xt_dists_av_frechet)

av_lcss = mean(xt_dists_av_lcss)

av_sidtw = mean(xt_dists_av_sidtw)

av_eucl = mean(xt_dists_av_eucl)