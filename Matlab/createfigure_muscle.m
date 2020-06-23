function createfigure_muscle(X1, Y1, X2, Y2)
%CREATEFIGURE(X1,Y1,X2,Y2)
%  X1:  vector of x data
%  Y1:  vector of y data
%  X2:  vector of x data
%  Y2:  vector of y data

%  Auto-generated by MATLAB on 20-May-2010 11:47:48

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'LineWidth',2,'GridLineStyle','--',...
    'FontSize',18);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 20000]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0 20000]);
grid(axes1,'on');
hold(axes1,'all');

% Create title
title('Muscle Area (mm^2)','FontSize',24);

% Create xlabel
xlabel('Semi-manual CT-based method','FontSize',22);

% Create ylabel
ylabel('Automated MRI-based method','FontSize',24);

% Create plot
plot(X1,Y1,'Parent',axes1,'Marker','o','LineWidth',4,'LineStyle','none',...
    'Color',[1 0 0],...
    'DisplayName','muscle area scatterplot');

% Create plot
plot(X2,Y2,'Parent',axes1,'LineWidth',2,'Color',[1 0 0],...
    'DisplayName','linear model');

% Create text
text('Parent',axes1,'String','y=0.97*x+817','Position',[12000 8000 0],...
    'FontSize',24);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','SouthEast','LineWidth',2,'FontSize',24);

