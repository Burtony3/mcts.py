close all
figure
theta = linspace(0, 2*pi, 100);


% Plotting Orbits
plot(cos(theta), sin(theta), "k:", "HandleVisibility", "off")
hold on
plot(2*cos(theta), 2*sin(theta), "k:", "HandleVisibility", "off")

% Launch Dates
n = 5;
thetaL = linspace(0, 2*pi, n+1);
xL = cos(thetaL);
yL = sin(thetaL);
color = {'#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F'};

% First Flyby Dates
thetaP = linspace(0.05*2*pi, 0.2*2*pi, n+1);
for i = 1:n
    xL_ = xL(i);
    yL_ = yL(i);
    color_ = color{i};
    thetaL_ = thetaL(i);
    for j = 1:n
        xF = 2*cos(thetaL_ + thetaP(j));
        yF = 2*sin(thetaL_ + thetaP(j));
        scatter(xF, yF, 25, "filled", "MarkerEdgeColor", color_, "MarkerFaceColor", color_)
        plot([xL_ xF], [yL_ yF], 'color', color_)
    end
end

plot([0, xL(1)], [0, yL(1)], "k--")
plot([0, xL(2)], [0, yL(2)], "k--")

% Plotting Launch Dates
for i = 1:n
    scatter(xL(i), yL(i), 25, "k", "filled")% "filled", "MarkerEdgeColor", color{i}, "MarkerFaceColor", color{i})
end

% Plotting Sun
scatter(0, 0, 200, 'o', 'MarkerEdgeColor', [0.85,0.33,0.10], "MarkerFaceColor", [0.93,0.69,0.13], "LineWidth", 1.5, "HandleVisibility", "off")

% Setting Plot Properties
hold off
axis equal square
xlim([-2.25 2.25])
ylim([-2.25 2.25])