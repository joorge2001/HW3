function [] = orbitPlotter(COE,mu,name,title_graph,rB,flag,param)
% -------------------------------------------------------------------------
% Function to plot the orbit of an spacecraft and its position according to
% the problem.
%
% Inputs: 
%     - COE: vector containing:
%       - a: semimajor axis [km]
%       - e: eccentricity
%       - i: inclination [rad]
%       - raan: right ascension of the ascending node [rad]
%       - omega: argument of periapsis [rad]
%       - theta: true anomaly [rad]
%     - name: name of the celestial body and orbits in a column vector
%     - title: title of the plot
%     - rB: radius of the celestial body
%     - flag: vector of binary values:
%       - If flag (k) = 0 the position of the spacecraft is not shown. 
%       - If flag (k) = 1 the position of the spacecraft is shown.
%     - param: parameter obtained by trial and error to better plot
%     hyperbolas and parabolas. param = 0 if the orbit is elliptical. 
%     
% Outputs:
%     - Plot in 3D axes comprising one or more orbits, the celestial body for
%       reference and a highlighted point in the orbit marking its position if 
%       needed.
% -------------------------------------------------------------------------

% Initialize arrays for position vectors
% x = zeros(length(nu), height(COE));
% y = zeros(length(nu), height(COE));
% z = zeros(length(nu), height(COE));

figure
hold on  
Legend{1} = string(name(1,:));
[xBody,yBody,zBody] = sphere; % create a unit sphere
rBody = rB; % scale the sphere
surf(rBody*xBody,rBody*yBody,rBody*zBody,'EdgeColor','none','FaceColor','b'); % plot the sphere

% Convert COE to position vectors at each true anomaly
for k = 1:height(COE)

    a(k) = COE(k,1);
    e(k) = COE(k,2);
    i(k) = COE(k,3);
    raan(k) = COE(k,4);
    omega(k) = COE(k,5);

    % Define a range of true anomalies

    if e(k) < 1
        nu = linspace(0, 2*pi, 1000);
    else
        nu = linspace(-acosh(e(k))-param(k), acosh(e(k))+param(k), 1000);
    end

    for j = 1:length(nu)
    % Replace the true anomaly in the COE
        theta(k) = nu(j);

    % Convert COE to position vector
        [r,~] = COE2rv(a(k),e(k),i(k),raan(k),omega(k),theta(k),mu); % assuming COE2RV is a function that converts COE to position vector

        x(j,k) = r(1);
        y(j,k) = r(2);
        z(j,k) = r(3);
    end
    plot3(x(:,k), y(:,k), z(:,k))
    Legend{1+k} = string(name(1+k,:));    
    hold on
end

for k = 1:height(COE)
    if flag(k) == 1
        [pos,~] = COE2rv(a(k),e(k),i(k),raan(k),omega(k),COE(k,6),mu);
        p(k) = pos(1);
        q(k) = pos(2);
        r(k) = pos(3);
    plot3(p(k),q(k),r(k),'.r','MarkerSize',20)
    end
end    

view(3)
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title(title_graph)
legend(Legend,'Location','southoutside','Orientation','horizontal')
grid on
axis equal
hold off


end