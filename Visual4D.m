topLevelFolder = pwd;
filename = strcat(topLevelFolder,'\Figures\Comparison.xlsx');
sheet = 2;
startPos=2;
endPos=617;
xlRange = strcat('A',int2str(startPos),':A',int2str(endPos));
k_merCol = xlsread(filename,sheet,xlRange);

xlRange = strcat('B',int2str(startPos),':B',int2str(endPos));
joinCol = xlsread(filename,sheet,xlRange);

xlRange = strcat('C',int2str(startPos),':C',int2str(endPos));
distCol = xlsread(filename,sheet,xlRange);

xlRange = strcat('D',int2str(startPos),':D',int2str(endPos));
rFCol = xlsread(filename,sheet,xlRange);

scatter3(joinCol,k_merCol,distCol,40,rFCol,'filled')    % draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('Joining Methods')
ylabel('Kmers')
zlabel('Distance Methods')

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'RF Distance';