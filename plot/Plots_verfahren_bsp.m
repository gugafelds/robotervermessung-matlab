clear

% Farben Für Bahnvergleich
blau = [0 0.4470 0.7410]; % Standard Blau
rot = [0.78 0 0];

% Für Plots Verfahren
c1 = [0 0.4470 0.7410];
c2 = [0.8500 0.3250 0.0980];
c3 = [0.9290 0.6940 0.1250];
c4 = [0.4940 0.1840 0.5560];
c5 = [0.4660 0.6740 0.1880];
c6 = [0.3010 0.7450 0.9330];
c7 = [0.6350 0.0780 0.1840];

load square1

% Umbenennen der Variablen
X = pathX'*1000;
Y = pathY'*1000;

% X = X(43:137,:);

pflag = 0;
% figure;
% plot3(X(:,1),X(:,2),X(:,3))
% hold on
% plot3(Y(:,1),Y(:,2),Y(:,3))

%% DTW
[distances, maxDistance, averageDistance, accdist, dtw_X, dtw_Y, path, ix, iy, localdist] = fkt_dtw3d(X,Y,pflag);
    
figure('Color','white','Name','DTW - Zuordnung der Bahnpunkte','NumberTitle','off')
    hold on;
    grid on;
    box on;
    plot3(dtw_X(:,1),dtw_X(:,2),dtw_X(:,3),Color= c1,LineWidth=1.5,Marker = "o",MarkerFaceColor= c1,MarkerSize=4);
    plot3(dtw_Y(:,1),dtw_Y(:,2),dtw_Y(:,3),Color= rot,LineWidth=1.5,Marker = "square",MarkerFaceColor=rot,MarkerSize=4);
    for i = 1:1:length(dtw_X)
        line([dtw_Y(i,1),dtw_X(i,1)],[dtw_Y(i,2),dtw_X(i,2)],[dtw_Y(i,3),dtw_X(i,3)],'Color','black')
    end
    view(220, 30)
    legend({'Istbahn','Sollbahn','Abweichung'},'Location','northeast',"FontWeight", "bold")
    xlabel("x","FontWeight","bold")
    ylabel("y","FontWeight","bold")
    zlabel("z","FontWeight","bold")
    hold off

%% DFD

% Plot ist in der Funktion umgesetzt
fkt_discreteFrechet(X,Y,pflag);


%% LCSS
[lcss_av, lcss_max, lcss_distances, lcss_accdist, lcss_path, lcss_X, lcss_Y,lcss_score,lcss_epsilon] = fkt_lcss(X,Y,pflag);
figure('Color','white','Name','LCSS - Zuordnung der Bahnpunkte','NumberTitle','off')
    hold on;
    grid on;
    box on;
    plot3(X(:,1),X(:,2),X(:,3),Color= c1,LineWidth=0.5,Marker = "o",MarkerFaceColor= c1,MarkerSize=4);
    plot3(Y(:,1),Y(:,2),Y(:,3),Color= rot,LineWidth=0.5,Marker = "square",MarkerFaceColor=rot,MarkerSize=4);
    for k = 1:length(lcss_path)    
        % Verbindungslinien
        plot3([lcss_X(k,1), lcss_Y(k,1)], [lcss_X(k,2), lcss_Y(k,2)], [lcss_X(k,3), lcss_Y(k,3)], 'k-');
    end
    plot3(lcss_X(:,1), lcss_X(:,2), lcss_X(:,3),Color= c1,LineWidth=2,Marker = "o",MarkerFaceColor= c1,MarkerSize=6);
    plot3(lcss_Y(:,1), lcss_Y(:,2), lcss_Y(:,3),Color= rot,LineWidth=2,Marker = "square",MarkerFaceColor=rot,MarkerSize=6);

    view(220, 30)
    legend({'Istbahn','Sollbahn','Abweichung'},'Location','northeast',"FontWeight", "bold")
    xlabel("x","FontWeight","bold")
    ylabel("y","FontWeight","bold")
    zlabel("z","FontWeight","bold")
    hold off
% LCSS Besser mit 10 mal Eukl Dist. 
figure('Color','white','Name','LCSS - Zuordnung der Bahnpunkte','NumberTitle','off')
    hold on;
    grid on;
    box on;
    plot3(lcss_X(:,1), lcss_X(:,2), lcss_X(:,3),Color= c1,LineWidth=1.5,Marker = "o",MarkerFaceColor= c1,MarkerSize=4);
    plot3(lcss_Y(:,1), lcss_Y(:,2), lcss_Y(:,3),Color= rot,LineWidth=1.5,Marker = "square",MarkerFaceColor=rot,MarkerSize=4);
    for k = 1:length(lcss_path)    
        % Verbindungslinien
        plot3([lcss_X(k,1), lcss_Y(k,1)], [lcss_X(k,2), lcss_Y(k,2)], [lcss_X(k,3), lcss_Y(k,3)], 'k-');
    end

    view(220, 30)
    legend({'Istbahn','Sollbahn','Abweichung'},'Location','northeast',"FontWeight", "bold")
    xlabel("x","FontWeight","bold")
    ylabel("y","FontWeight","bold")
    zlabel("z","FontWeight","bold")
    hold off


%% SIDTW
[distances, maxDistance, averageDistance, AccumulatedDistance, dtw_X, dtw_Y, MappingIndexes, ix, iy] = fkt_selintdtw3d(X,Y,pflag);

figure('Color','white','Name','SelectiveInterpolationDTW - Zuordnung der Bahnpunkte','NumberTitle','off')
    hold on;
    grid on;
    box on;
    plot3(dtw_X(:,1),dtw_X(:,2),dtw_X(:,3),Color= c1,LineWidth=1.5,Marker = "o",MarkerFaceColor= c1,MarkerSize=4);
    plot3(dtw_Y(:,1),dtw_Y(:,2),dtw_Y(:,3),Color= rot,LineWidth=1.5,Marker = "square",MarkerFaceColor=rot,MarkerSize=4);
    for i = 1:1:length(dtw_X)
        line([dtw_Y(i,1),dtw_X(i,1)],[dtw_Y(i,2),dtw_X(i,2)],[dtw_Y(i,3),dtw_X(i,3)],'Color','black')
    end
    view(220, 30)
    legend({'Istbahn','Sollbahn','Abweichung'},'Location','northeast',"FontWeight", "bold")
    xlabel("x","FontWeight","bold")
    ylabel("y","FontWeight","bold")
    zlabel("z","FontWeight","bold")
    hold off
