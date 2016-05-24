clear all;
close all;
% rng(0);
for outer=1:1
hold on;
npoints = 15;

focal = [512 ; 512];
principle = [256 ; 256];
skew = 0;
K = [focal(1) skew principle(1) ; 0 focal(2) principle(2) ; 0 0 1];

C1 = K * [1 0 0 0; 0 1 0 0; 0 0 1 0]; % Trivial extrinsics.
C2Rod = [0 30 0] .* pi ./ 180;
C2R = expm([0 -C2Rod(3) C2Rod(2); C2Rod(3) 0 -C2Rod(1); -C2Rod(2) C2Rod(1) 0]);
C2t = [1 0 0]';
C2t = C2t ./ norm(C2t);
C2 = K * [C2R' -C2R'*C2t];

trueE = ([0 -C2t(3) C2t(2); C2t(3) 0 -C2t(1); -C2t(2) C2t(1) 0] * C2R)';

worldPoints = zeros(4,npoints);
for i=1:npoints
    worldPoints(1:2,i) = rand(2,1) - 0.5;
    worldPoints(3,i) = rand() + 3;
    worldPoints(4,i) = 1;
end

points1 = C1*worldPoints;
points1(1,:) = points1(1,:) ./ points1(3,:);
points1(2,:) = points1(2,:) ./ points1(3,:);
points1(3,:) = ones(1,npoints);

wpoints2 = worldPoints;
for i=1:npoints
    wpoints2(:,i) = wpoints2(:,i) - [C2t; 0];
    wpoints2(:,i) = [[C2R' [0;0;0]]; 0 0 0 1] * wpoints2(:,i);
end

points2 = C1 * wpoints2; %hrmmm
points2(1,:) = points2(1,:) ./ points2(3,:);
points2(2,:) = points2(2,:) ./ points2(3,:);
points2(3,:) = ones(1,npoints);
   
Q1 = points1;
Q2 = points2;
for i=1:npoints
    Q1(1,i) = (Q1(1,i) - principle(1)) ./ focal(1);
    Q1(2,i) = (Q1(2,i) - principle(2)) ./ focal(2);
    Q2(1,i) = (Q2(1,i) - principle(1)) ./ focal(1);
    Q2(2,i) = (Q2(2,i) - principle(2)) ./ focal(2);
end
key = randsample(npoints, 5, false);
q1 = zeros(3,5);
q2 = zeros(3,5);
for i=1:5
    q1(:,i) = Q2(:,key(i));
    q2(:,i) = Q1(:,key(i));
end

% figure; plot(q1(1,:), q1(2,:));
% figure; plot(q2(1,:), q2(2,:));

% uh1 = q1;
% uh2 = q2;
% for i=1:npoints
%     uh1(:,i) = uh1(:,i) ./ norm(uh1(:,i));
%     uh2(:,i) = uh2(:,i) ./ norm(uh2(:,i));    
% end

G1 = [0  0 0; 0 0 -1;  0 1 0];
G2 = [0  0 1; 0 0  0; -1 0 0];
G3 = [0 -1 0; 1 0  0;  0 0 0];
P = [1 0 0; 0 1 0];
a = [0 pi/6 0 0 pi/2]';
a(1:3) = sph2cart(2*pi*rand(), asin(2*rand()-1), pi*rand().^(1/3));
[a(4:5) ~] = sph2cart(2*pi*rand(), asin(2*rand()-1), 1);

N = 100;
rn = zeros(1,N);
r = zeros(5,1);
J = zeros(5,5);
tdef = [0;0;1];
for n=1:N
    R = expm(a(1)*G1 + a(2)*G2 + a(3)*G3);
    Rt = expm(a(4)*G1 + a(5)*G2);
    th = Rt*tdef;
    thx = [0 -th(3) th(2); th(3) 0 -th(1); -th(2) th(1) 0];
    E = thx*R;
    
    for i=1:5
        r(i) = q2(:,i)'*E*q1(:,i);
    end
    rn(n) = norm(r);
    
    % residuals...
    tmp = P*E*q1;
    g1 = 1 ./ sqrt(tmp(1,:).^2 + tmp(2,:).^2)';
    tmp = P*E'*q2;
    g2 = 1 ./ sqrt(tmp(1,:).^2 + tmp(2,:).^2)';
    
    Rq1 = R*q1;
    q2thx = q2'*thx;
    for i=1:5
        J(i,1) = q2thx(i,:) * G1 * Rq1(:,i);
        J(i,2) = q2thx(i,:) * G2 * Rq1(:,i);
        J(i,3) = q2thx(i,:) * G3 * Rq1(:,i);

        tmp = Rt * G1 * tdef;
        J(i,4) = q2(:,i)' * [0 -tmp(3) tmp(2); tmp(3) 0 -tmp(1); -tmp(2) tmp(1) 0] * Rq1(:,i);
        tmp = Rt * G2 * tdef;
        J(i,5) = q2(:,i)' * [0 -tmp(3) tmp(2); tmp(3) 0 -tmp(1); -tmp(2) tmp(1) 0] * Rq1(:,i);
    end
    
    a = a - J \ r;
end
plot(log10(rn));
drawnow;
end
hold off;

disp(a' .* 180 ./ pi);
disp(g1');
disp(g2');
% 
% 
% for i=1:1000
% hold on;
% [x,y,z] = sph2cart(2*pi*rand(), asin(2*rand()-1), rand().^(1/3));
% plot3(x,y,z,'.');
% drawnow;
% end

allr = zeros(1,npoints);
for i=1:npoints
    allr(i) = Q2(:,i)'*E*Q1(:,i);
end
figure; histogram(allr);