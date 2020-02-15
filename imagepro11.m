function f = imagepro11()
Z_img=imread('Z.jpg');
B_img = imread('B_test.jpg');
X_img = imread('X.jpg');
Z1 = imresize(Z_img,[16 16]); 
% disp(I2)
Z1=rgb2gray(Z1); % use if the image containing RGB value 3
% newin = im2double(In);
% disp(newin)
Zout  = imbinarize(Z1);
Z_oned = Zout(:);
Z_onedd = Z_oned - mean(Z_oned);
[n m] = size(Z_oned)
H = hadamard(n);    % Hadamard matrix
Z_walsh = H*Z_oned;
Z_walsh1 = H*Z_onedd;


B1 = imresize(B_img,[16 16]); 
B1=rgb2gray(B1); % use if the image containing RGB value 3
Bout  = imbinarize(B1);
B_oned = Bout(:);
B_onedd = B_oned - mean(B_oned);
[g h] = size(B_oned)
H1 = hadamard(g);    % Hadamard matrix
B_walsh = H1*B_oned;
B_walsh1 = H1*B_onedd;

X1 = imresize(X_img,[16 16]); 
X1=rgb2gray(X1); % use if the image containing RGB value 3
Xout  = imbinarize(X1);
X_oned = Xout(:);
X_onedd = X_oned - mean(X_oned);
[g h] = size(X_oned)
H = hadamard(g);    % Hadamard matrix
X_walsh = H*X_oned;
X_walsh1 = H*X_onedd;

% plotting figure comparison
figure
subplot(2,1,1)
plot(Z_walsh,'r')
hold on
plot(B_walsh,'g')
hold on
plot(X_walsh,'b')
legend('Z','B','X')
title('frequency plot for Walsh Of Z,B,X')

subplot(2,1,2)
plot(Z_walsh1,'r')
hold on
plot(B_walsh1,'g')
hold on
plot(X_walsh1,'b')
legend('Z','B','X')
title('frequency plot for Walsh Z,B,X - Mean')


% code for fast walsh hadamard needs package installation for that
% Z_fwht = fwht(Z_walsh);
% Z_fwht1 = fwht(Z_walsh1);
% 
% subplot(2,2,[3 4]);
% plot(Z_fwht)
% title('frequency plot for fast Walsh Of Z')
% subplot(2,3,8);
% plot(Z_fwht)
% title('frequency plot for fast Walsh Of Z-Mean(Z)')



% y = fwht(newin,N,'sequency'); %Perform Fast-walsh-hadamard-transform with order 128
% disp(H)
% plot(H)
% imshow(H); %Display image transform
% figure;
% imshow(In);
% disp(I)
% whos In
% out  = imbinarize(In);




