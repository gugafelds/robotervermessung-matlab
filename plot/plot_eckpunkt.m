figure
plot3(selection(:,1),selection(:,2),selection(:,3))
hold on 
mean_selection = mean(selection);
plot3(mean_selection(:,1),mean_selection(:,2),mean_selection(:,3),'ro')
axis equal

%%
mean_selection = mean(selection)

% 3D-Plot erstellen
figure('Color','white'); 
plot3(selection(:,1), selection(:,2), selection(:,3), 'LineWidth', 2) % Dickere Linie
hold on
plot3(mean_selection(:,1), mean_selection(:,2), mean_selection(:,3), 'rx', 'MarkerSize', 10, 'LineWidth', 2) % Mittelwert plotten
hold off

axis equal

% Achsbeschriftungen hinzufügen
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')

% Anzahl der Ticks reduzieren
xt = linspace(min(selection(:,1)), max(selection(:,1)), 5);
yt = linspace(min(selection(:,2)), max(selection(:,2)), 5);
zt = linspace(min(selection(:,3)), max(selection(:,3)), 5);

xticks(unique(round(xt, 1)))
yticks(unique(round(yt, 1)))
zticks(unique(round(zt, 1)))


% 3D-Ansicht festlegen
view(45, 30)

% Titel hinzufügen
% title('3D-Plot mit Mittelwert und angepassten Achsen')
%%
a = vecnorm(trajectory_ist(idx1)-selection_mean,2,2)
