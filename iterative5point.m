clear all;
close all;
rng(1);

focal = [512 ; 512];
principle = [256 ; 256];
skew = 0;
K = [focal(1) skew principle(1) ; 0 focal(2) principle(2) ; 0 0 1];

C1 = K * [1 0 0 0; 0 1 0 0; 0 0 1 0]; % Trivial extrinsics.
C2Rod = [0 -30 0] .* pi ./ 180;
C2R = expm([0 -C2Rod(3) C2Rod(2); C2Rod(3) 0 -C2Rod(1); -C2Rod(2) -C2Rod(1) 0]);
C2t = [1 0 0]';
C2t = C2t ./ norm(C2t);
C2 = K * [C2R C2t];

trueE = C2R * [0 -C2t(3) C2t(2); C2t(3) 0 -C2t(1); -C2t(2) C2t(1) 0];

worldPoints = zeros(4,5);
for i=1:5
    worldPoints(1:2,i) = rand(2,1) - 0.5;
    worldPoints(3,i) = rand() + 3;
    worldPoints(4,i) = 1;
end

points1 = C1*worldPoints;
points1(1,:) = points1(1,:) ./ points1(3,:);
points1(2,:) = points1(2,:) ./ points1(3,:);
points1(3,:) = ones(1,5);

wpoints2 = worldPoints;
for i=1:5
    wpoints2(:,i) = wpoints2(:,i) - [C2t; 0];
    wpoints2(:,i) = [[C2R' [0;0;0]]; 0 0 0 1] * wpoints2(:,i);
end

points2 = C1 * wpoints2; %hrmmm
points2(1,:) = points2(1,:) ./ points2(3,:);
points2(2,:) = points2(2,:) ./ points2(3,:);
points2(3,:) = ones(1,5);
   
upoints1 = points1;
upoints2 = points2;
for i=1:5
    upoints1(1,i) = (upoints1(1,i) - principle(1)) ./ focal(1);
    upoints1(2,i) = (upoints1(2,i) - principle(2)) ./ focal(2);
    upoints2(1,i) = (upoints2(1,i) - principle(1)) ./ focal(1);
    upoints2(2,i) = (upoints2(2,i) - principle(2)) ./ focal(2);
end

% figure; plot(upoints1(1,:), upoints1(2,:));
% figure; plot(upoints2(1,:), upoints2(2,:));

uh1 = upoints1;
uh2 = upoints2;
for i=1:5
    uh1(:,i) = uh1(:,i) ./ norm(uh1(:,i));
    uh2(:,i) = uh2(:,i) ./ norm(uh2(:,i));    
end

a = [0 0 0 0 0]';
G1 = [0  0 0; 0 0 -1;  0 1 0];
G2 = [0  0 1; 0 0  0; -1 0 0];
G3 = [0 -1 0; 1 0  0;  0 0 0];
% Why the above pattern appears is left as an exercise to the reader.
P = [1 0 0; 0 1 0];
R1old = expm(a(1).*G1 + a(2).*G2 + a(3).*G3);
R2old = expm(a(4).*G1 + a(5).*G2);

rrmin = 9999;
tikhonov = true;
lambda = 0.3;
plotFrequency = 1000000;
h = waitbar(0, '...');
N = 1000;
rrlog = zeros(1,N);
alog = zeros(1,N);
for n=1:N
    waitbar(n./N, h, 'Thinking...');
    R1 = expm(a(1).*G1 + a(2).*G2 + a(3).*G3) * R1old;
    R2 = expm(a(4).*G1 + a(5).*G2) * R2old;

    w1 = P * R1 * uh1;
    w2 = P * R2 * uh2;

    wh1 = w1 ./ repmat(arrayfun(@(idx) norm(w1(:,idx)), 1:size(w1,2)), [2,1]); % Norm along one dimension only
    wh2 = w2 ./ repmat(arrayfun(@(idx) norm(w2(:,idx)), 1:size(w2,2)), [2,1]); % https://www.mathworks.com/matlabcentral/newsreader/view_thread/172116

    if mod(n,plotFrequency) == 0
        figure; plot(0,0,'.'); hold on;
        for i=1:5
            plot([wh1(1,i) wh2(1,i)], [wh1(2,i) wh2(2,i)]);            
        end
        hold off; xlim([-1.1 1.1]); ylim([-1.1 1.1]);
    end

    r = acos(sum(wh1.*wh2, 1))';
    s = cross([wh1;zeros(1,5)], [wh2;zeros(1,5)]);
    s = sign(s(3,:))';
    r = r .* s;
    rr = sum(r.^2);
    if (rr < rrmin)
        rrmin = rr;
        best1 = R1;
        best2 = R2;
    end
    rrlog(n) = rr;

    J = zeros(5,5);
    for i=1:5
        W1 = w1(:,i);
        W2 = w2(:,i);
        UH1 = R1old * uh1(:,i);
        UH2 = R2old * uh2(:,i);

        q = sum(W1.*W1) + sum(W2.*W2);
        p = W1(1).*W2(1) - W1(2).*W2(2);
        scale = -1 / sqrt(q*q * (q - p*p));

        dtdW = zeros(1,4);
        dtdW(1) = W1(2)^2*W2(1) - W1(1)*W1(2)*W2(2) + W2(1)^3 + W2(1)*W2(2)^2;
        dtdW(2) = W1(1)^2*W2(2) - W1(2)*W1(1)*W2(1) + W2(1)^2*W2(2) + W2(2)^3;
        dtdW(3) = W1(1)^3 + W1(1)*W1(2)^2 + W1(1)*W2(2)^2 - W2(1)*W1(2)*W2(2);
        dtdW(4) = W1(1)^2*W1(2) - W2(2)*W1(1)*W2(1) + W1(2)^3 + W1(2)*W2(1)^2;
        dtdW = dtdW .* scale;

        dWda = zeros(5,4);
        dWda(1,:) = [   0     UH1(3)  0      0    ];
        dWda(2,:) = [-UH1(3)    0     0      0    ];
        dWda(3,:) = [ UH1(2) -UH1(1)  0      0    ];
        dWda(4,:) = [   0       0     0    UH2(3) ];
        dWda(5,:) = [   0       0  -UH2(3)   0    ];

        dtda = sum(repmat(dtdW, [5,1]) .* dWda, 2);
        J(:,i) = dtda;
    end

    JTJ = J'*J;
    if tikhonov
        regularization = diag(lambda .* diag(JTJ));
    else
        regularization = lambda .* eye(5,5);
    end
    delta = (JTJ+regularization) \ (J'*r);
    a = delta;
    alog(n) = norm(a);

    R1old = R1;
    R2old = R2;
end
close(h);
figure; plot(log10(rrlog)); hold on; plot(log10(alog)); hold off; legend('residual', 'gradient');

E = best2' * [0 -1 0; 1 0 0; 0 0 0] * best1;
err1 = sum(sum((trueE - E).^2));
E = -E;
err2 = sum(sum((trueE - E).^2));
min(err1, err2)


