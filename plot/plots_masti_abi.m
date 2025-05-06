%% PLOTTEN FÜR MASTI ABI

% % WENN NEUE DATEN ERZEUGEN
% t1 = 0:0.1:1;
% t2 = 0:0.035:1;
% X = ones(size(t1));
% X(3) = 1.3;
% X(8) = 0.3;
% Y = 2*rand(size(t2));
% 
% X = [X; t1]';
% Y = [Y; t2]';
% %%
% pflag = 1;
% 
% [distances, maxDistance, averageDistance, accdist, dtwX, dtwY, path] = fkt_dtw3d(X,Y,pflag);

% BEKANNTE DATEN LADEN
load Beispiel_Zeitreihe.mat

pflag = 1;

[distances, maxDistance, averageDistance, accdist, dtwX, dtwY, path, ix, iy] = fkt_dtw3d(X,Y,pflag);


% f = figure; 
% hold on;
% plot(X(:,2),X(:,1),'-bo','MarkerFaceColor','b','LineWidth',1.5);
% plot(Y(:,2),Y(:,1),'-ro','MarkerFaceColor','r','LineWidth',1.5);
% xlabel('Zeit');
% ylabel('Amplitude');
% legend('Zeitreihe X', 'Zeitreihe Y')
% f.Position(3) = 3*f.Position(3);
% axis([0 1 0 2.5])

%% PLOTTEN
% Signal mit Zuordnung
f = figure('Name','DTW - Zuordnung der Datenpunkte','NumberTitle','off'); 
hold on;
plot(dtwX(:,2),dtwX(:,1),'-bo','MarkerFaceColor','b','LineWidth',1.5);
plot(dtwY(:,2),dtwY(:,1),'-ro','MarkerFaceColor','r','LineWidth',1.5);
for i = 1:1:length(dtwX)
    line([dtwY(i,2),dtwX(i,2)],[dtwY(i,1),dtwX(i,1)],'Color','black')
end
xlabel('Zeit [Sekunden]');
ylabel('Datenwert [$X,Y$]');
legend('Zeitreihe $X$', 'Zeitreihe $Y$','Zuordnung')
f.Position(3) = 2*f.Position(3);
axis([0 1 0 2.3])

clear f

f = figure('Name','DTW - Kostenkarte und Mapping ','NumberTitle','off');
ax = gca;
hold on
% main=subplot('Position',[0.19 0.19 0.67 0.79]);           
imagesc(accdist)
colormap("turbo"); % "turbo" "jet" "hot" "cool"
colorb = colorbar;
colorb.Label.String = 'Akkumulierte Kosten';
% --------To-Do: Colorbar normen auf 1----------
% set(colorb,'FontSize',10,'YTickLabel','');
% set(colorb,'FontSize',10;
% hold on
plot(iy+0.5, ix+0.5,"-w","LineWidth",1.5)
grid on;
xlabel('Index Zeitreihe $Y$');
ylabel('Index Zeitreihe $X$');
axis([min(iy) max(iy) 1 max(ix)]);
% Grid funktioniert nicht!
ax.GridColor = "black";
ax.GridLineWidth = 1.0;
set(gca,'FontSize',10,'YDir', 'normal');


%%

x = [1 2 3];
y = [1 2 1];
plot(x,y, 'LineWidth',2);
hold on;
z = [1.5 2 1.5];
plot(x,z, 'LineWidth',2);