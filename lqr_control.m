m0 = 1.5;
m1 = 0.5;
m2 = 0.75;
L1 = 0.75;
L2 = 1;
Q = diag([5 50 50 700 700 700]);
R = 1;
Nh = 40;
g = 10;


d0 = [m0 + m1 + m2, (m1/2 + m2)*L1, (m2*L2)/2; (m1/2 + m2)*L1, (m1/3 + m2)*L1^2, (m2*L1*L2)/2; (m2*L2)/2, (m2*L1*L2)/2, (m2*L2^2)/3];
g1 = [0, 0, 0; 0, -(m1/2 + m2)*L1*10, 0; 0, 0, -(m2*L2*g)/2];

a = [zeros(3), eye(3); -inv(d0)*g1, zeros(3)];
b = [zeros(3, 1); inv(d0) * [1; 0; 0]];
c = eye(6);
sys = ss(a, b, c, 0);

k = lqr(sys, Q, R);
sysnew = ss(a - b*k, b, c, 0);
curr1 = [0; deg2rad(30); -deg2rad(30); 0; 0; 0];
curr2 = [0; deg2rad(30); -deg2rad(30); 0; 0; 0];

dt = 0.001;
strt = 0.001;
en = 7;
tot = ((en - strt) / dt) + 1;

mu = [0 0 0 0 0 0];
sigmaVal = 0.01
sigma = sigmaVal * eye(6);
rng('default')  % For reproducibility
noise = mvnrnd(mu,sigma,ceil(tot));

y1 = [];
t1 = [];
x1 = [];
y2 = [];
t2 = [];
x2 = [];
i = 1;

for time = strt : dt : en
    [y1d, t1d, x1d] = initial(sysnew, curr1, 0 : dt : dt);
    y1 = vertcat(y1, y1d(2,:));
    t1 = vertcat(t1, time);
    x1 = vertcat(x1, x1d(2,:));
    curr1 = x1d(2,:);
    curr1 = curr1 + noise(i, :);
    [y2d, t2d, x2d] = initial(sysnew, curr2, 0 : dt : dt);
    y2 = vertcat(y2, y2d(2,:));
    t2 = vertcat(t2, time);
    x2 = vertcat(x2, x2d(2,:));
    curr2 = x2d(2,:);
    i = i + 1;
end


y1(:, 2 : 3) = rad2deg(y1(:, 2 : 3));
y1(:, 5 : 6) = rad2deg(y1(:, 5 : 6));
y2(:, 2 : 3) = rad2deg(y2(:, 2 : 3));
y2(:, 5 : 6) = rad2deg(y2(:, 5 : 6));

close all

for i = 1 : 1 : 6
    subplot(2, 3, i)
    plot(t1, y1(:, i), 'r', t2, y2(:, i), 'b')
    grid
end

figure();

i = 1;
video = VideoWriter('0.01.avi', 'Uncompressed AVI');
open(video);
ax = gca();
for k = 1 : 10 : length(t1)
    W = 4;  % cart width
    H = .5; % cart height
    mr = .3; % mass radius
   
    x11 = y1(k, 1);
    th1 = deg2rad(y1(k, 2));
    th2 = deg2rad(y1(k, 3));
    y11 = H/2;

    % positions
    % y = wr/2; % cart vertical position
    w1x = x11-.9*W/2;
    w1y = 0;
    w2x = x11+.9*W/2;
    w2y = 0;
    
    x22 = x11 + L1*sin(th1);
    y22 = y11 + L1*cos(th1);
    
    x33 = x22 + L2*sin(th2);
    y33 = y22 + L2*cos(th2);
    
    plot([-5*W 5*W], [0 0], 'k', 'LineWidth', 2)
    hold on
    
    rectangle('Position',[x11-W/2,0,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1]);
    
    plot([x11 x22], [y11 y22], 'k', 'LineWidth', 2)
    plot([x22 x33], [y22 y33], 'k', 'LineWidth', 2)
    rectangle('Position',[x22-0.1,y22-0.1,0.2,0.2],'Curvature',[1,1],'FaceColor',[1 0.1 0.1])
        
    xlim([-5*W 5*W]);
    ylim([-2 2.5]);
    
    drawnow
    
    if k == 1
        ax = gca();
    end
    
    F(i) = getframe(ax);
    i = i + 1;
    hold off
    
end

writeVideo(video, F);
close(video);

% figure();
% 
% for k = 1 : 10 : length(t2)
%     W = 4;  % cart width
%     H = .5; % cart height
%     mr = .3; % mass radius
%     
%     x11 = y2(k, 1);
%     th1 = deg2rad(y2(k, 2));
%     th2 = deg2rad(y2(k, 3));
%     y11 = H/2;
% 
%     % positions
%     % y = wr/2; % cart vertical position
%     w1x = x11-.9*W/2;
%     w1y = 0;
%     w2x = x11+.9*W/2;
%     w2y = 0;
%     
%     x22 = x11 + L1*sin(th1);
%     y22 = y11 + L1*cos(th1);
%     
%     x33 = x22 + L2*sin(th2);
%     y33 = y22 + L2*cos(th2);
%     
%     plot([-5*W 5*W], [0 0], 'k', 'LineWidth', 2)
%     hold on
%     
%     rectangle('Position',[x11-W/2,0,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1]);
%     rectangle('Position',[x22-0.1,y22-0.1,0.2,0.2],'Curvature',[1,1],'FaceColor',[1 0.1 0.1])
%     
%     plot([x11 x22], [y11 y22], 'k', 'LineWidth', 2)
%     plot([x22 x33], [y22 y33], 'k', 'LineWidth', 2)
%         
%     xlim([-5*W 5*W]);
%     ylim([-2 2.5]);
%     
%     drawnow
%     hold off
%     
% end
