% OCR (Optical Character Recognition).

warning off 
% Clear all
clc, close all, clear all
% Read image
imagen=imread('TEST_3.jpg');
% Show image
imshow(imagen);
title('Inpute Image')
% Convert to gray scale
if size(imagen,3)==3 %RGB image
    imagen=rgb2gray(imagen);
end
% Convert to BW
threshold = graythresh(imagen);
imagen =~im2bw(imagen,threshold);
% Remove all object containing fewer than 30 pixels
imagen = bwareaopen(imagen,30);
%Storage matrix word from image
word=[ ];
re=imagen;
%Opens text.txt as file for write
fid = fopen('text.txt', 'wt');
% Load templates
load templates
global templates
% Compute the number of letters in template file
num_letras=size(templates,2);
while 1
    %Fcn 'lines' separate lines in text
    [fl re]=lines(re);
    imgn=fl;
    %Uncomment line below to see lines one by one
    imshow(fl);pause(0.5)    
    %-----------------------------------------------------------------     
    % Label and count connected components
    [Ls Ne] = bwlabel(imgn);    
    for n=1:Ne
        [r,c] = find(Ls==n);
        % Extract letter
        n1=imgn(min(r):max(r),min(c):max(c));  
        % Resize letter (same size of template)
        img_r=imresize(n1,[42 24]);
        %Uncomment line below to see letters one by one
         imshow(img_r);pause(0.5)
        %-------------------------------------------------------------------
        % Call fcn to convert image to text
        letter=read_letter(img_r,num_letras);
        % Letter concatenation
        word=[word letter];
    end
    %fprintf(fid,'%s\n',lower(word));%Write 'word' in text file (lower)
    fprintf(fid,'%s\n',word);%Write 'word' in text file (upper)
    % Clear 'word' variable
    word=[ ];
    %*When the sentences finish, breaks the loop
    if isempty(re)  %See variable 're' in Fcn 'lines'
        break
    end    
end
fclose(fid);
%Open 'text.txt' file
winopen('text.txt')

%*_*_*_*_*_*__*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_**_
% HADAMARD TRANSFORM
%*_*_*_*_*_*_*_*_

Z_img=imread('Z.jpg');
B_img = imread('B_test.jpg');
X_img = imread('X.jpg');
Z1 = imresize(Z_img,[32 32]); 
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


B1 = imresize(B_img,[32 32]); 
B1=rgb2gray(B1); % use if the image containing RGB value 3
Bout  = imbinarize(B1);
B_oned = Bout(:);
B_onedd = B_oned - mean(B_oned);
[g h] = size(B_oned)
H1 = hadamard(g);    % Hadamard matrix
B_walsh = H1*B_oned;
B_walsh1 = H1*B_onedd;

X1 = imresize(X_img,[32 32]); 
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
title('Frequency plot for Walsh Of Z ,B ,X  Letters')

subplot(2,1,2)
plot(Z_walsh1,'r')
hold on
plot(B_walsh1,'g')
hold on
plot(X_walsh1,'b')
legend('Z','B','X')
title('Frequency plot for Walsh Z ,B ,X Letters - Mean/Normalized Value')


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

%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*


%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
% Forier transform of alphabets
%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*

Fs = 1000;            % Sampling frequency                    
Time = 1/Fs;             % Sampling period       
Li = 1500;             % Length of signal
ti = (0:Li-1)*Time; 

Zfft = fft(Z_onedd);
Bfft = fft(B_onedd);
Xfft = fft(X_onedd); 

Zfft1 = abs(fft(Z_oned));
Bfft1 = abs(fft(B_oned));
Xfft1 = abs(fft(X_oned)); 

% plot(Zfft)
% Zfft = abs(fft(Z_oned));
% Bfft = abs(fft(B_oned));
% Xfft = abs(fft(X_oned)); 
% P2 = abs(Zfft/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% plot(Fs,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
 
figure
subplot(2,1,1)
plot(1000*ti(1:512),Zfft(1:512),'r') 
hold on
plot(1000*ti(1:512),Bfft(1:512),'b') 
hold on  
plot(1000*ti(1:512),Xfft(1:512),'g') 
legend('Z','B','X')
title('Fourier Transform for Z ,B ,X Letters ')    

subplot(2,1,2)
plot(1000*ti(1:512),Zfft1(1:512),'r') 
hold on
plot(1000*ti(1:512),Bfft1(1:512),'b') 
hold on  
plot(1000*ti(1:512),Xfft1(1:512),'g') 
legend('Z','B','X')
title('Fourier Transform for Z ,B ,X Letters - With Absolute Values')

% y = fwht(newin,N,'sequency'); %Perform Fast-walsh-hadamard-transform with order 128
% disp(H)
% plot(H)
% imshow(H); %Display image transform
% figure;
% imshow(In);
% disp(I)
% whos In
% out  = imbinarize(In);
%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*

%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_



%__*_*_*_*_*_*__*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_
%print CERTRIOD
%__*_*_*_*_*_*__*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_

Im = imread('TEST_3.jpg');
Ibw = imbinarize(Im);
Ibw = imfill(Ibw,'holes');
%Ilabel = bwlabel(Ibw);
stat = regionprops(~Ibw,'centroid');

figure
imshow(Im); hold on;

for lp = 1: numel(stat)
      plot(stat(lp).Centroid(1),stat(lp).Centroid(2),'*');
end
%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_

All_letters = ["A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" ... 
    "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z" "1" "2" "3" "4" "5" "6" "7" "8" "9" "0"];

for pi = 1:length(templates)
    
    disp('Featurs of Letter : '+All_letters(pi));
    stats=regionprops(cell2mat(templates(pi)),'Area','Centroid','MajorAxisLength','MinorAxisLength','FilledArea','Perimeter');
    %stats1=regionprops(cell2mat(templates(pi)),'Centroid');
    %stats2=regionprops(cell2mat(templates(pi)),'MinorAxisLength');
    %Ma_A_Length=stats1.MajorAxisLength;
    %Mi_A_Length=stats2.MinorAxisLength;
    %disp("MajorAxisLength of : ");
    %disp(Ma_A_Length)
    %disp("MinorAxisLength of T : ");
    %disp(Mi_A_Length)
    t0=struct2table(stats);
    %t1=struct2table(stats1);
    %t2=struct2table(stats2);
    disp(t0);
    %writetable(t0,'test_1.txt');
    %disp(t1);
    %disp(t2);
    
    %[~, numberOfClosedRegions] = bwlabel(1-(cell2mat(templates(pi))));
        %disp('no of regions')
        %disp(numberOfClosedRegions)
        
     centsx(pi)=stats.Centroid(1);
     centsy(pi)=stats.Centroid(2);
     area(pi)=stats.Area;
     filled_area(pi)=stats.FilledArea;
end


%path = 'H:\Desktop\OCR';
%writetable(t0,'test_1.txt');

figure
plot(centsx,centsy,'.r')
text(centsx,centsy,All_letters,'VerticalAlignment','bottom','HorizontalAlignment','right')
xlim([5 20])
ylim([15 30])
title('Centroid of all Characters')
%*_*_*_*_*_*_*_*_*_*_*_*_*_
%AREA
%*_*_*_*_*_*_*_*_*_*_*_*_
figure
x=(1:36);
plot(x,area,'*r')
text(x,area,All_letters,'VerticalAlignment','bottom','HorizontalAlignment','right')
ylim([min(area)-50 max(area)+50])
title('Area value of all Characters')

%*_*_*_*_*_*_*_*_*_*_*_*_
%FilledArea
%*_*_*_*_*_*_*_*_*_*_*_*_

figure
x=(1:36);
plot(x,area,'.b')
text(x,filled_area,All_letters,'VerticalAlignment','bottom','HorizontalAlignment','right')
ylim([min(filled_area)-50 max(filled_area)+50])
title('Filled Area value of all Characters')

%*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
fprintf('For more details contact Tejas B, Sudhanshu C, Arjun S.')

clear all