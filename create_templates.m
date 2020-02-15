%CREATE TEMPLATES
%Letter

clc;
A=imread('letters_numbers\A.bmp');B=imread('letters_numbers\B.bmp');
C=imread('letters_numbers\C.bmp');D=imread('letters_numbers\D.bmp');
E=imread('letters_numbers\E.bmp');F=imread('letters_numbers\F.bmp');
G=imread('letters_numbers\G.bmp');H=imread('letters_numbers\H.bmp');
I=imread('letters_numbers\I.bmp');J=imread('letters_numbers\J.bmp');
K=imread('letters_numbers\K.bmp');L=imread('letters_numbers\L.bmp');
M=imread('letters_numbers\M.bmp');N=imread('letters_numbers\N.bmp');
O=imread('letters_numbers\O.bmp');P=imread('letters_numbers\P.bmp');
Q=imread('letters_numbers\Q.bmp');R=imread('letters_numbers\R.bmp');
S=imread('letters_numbers\S.bmp');T=imread('letters_numbers\T.bmp');
U=imread('letters_numbers\U.bmp');V=imread('letters_numbers\V.bmp');
W=imread('letters_numbers\W.bmp');X=imread('letters_numbers\X.bmp');
Y=imread('letters_numbers\Y.bmp');Z=imread('letters_numbers\Z.bmp');
%Number
one=imread('letters_numbers\1.bmp');  two=imread('letters_numbers\2.bmp');
three=imread('letters_numbers\3.bmp');four=imread('letters_numbers\4.bmp');
five=imread('letters_numbers\5.bmp'); six=imread('letters_numbers\6.bmp');
seven=imread('letters_numbers\7.bmp');eight=imread('letters_numbers\8.bmp');
nine=imread('letters_numbers\9.bmp'); zero=imread('letters_numbers\0.bmp');
%*-*-*-*-*-*-*-*-*-*-*-
letter=[A B C D E F G H I J K L M...
    N O P Q R S T U V W X Y Z];

%All_letters = ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

number=[one two three four five...
    six seven eight nine zero];
character=[letter number];
templates=mat2cell(character,42,[24 24 24 24 24 24 24 ...
    24 24 24 24 24 24 24 ...
    24 24 24 24 24 24 24 ...
    24 24 24 24 24 24 24 ...
    24 24 24 24 24 24 24 24]);
save ('templates','templates')


for i = 1:length(templates)
    stats=regionprops(cell2mat(templates(i)));
    stats1=regionprops(cell2mat(templates(i)),'MajorAxisLength');
    stats2=regionprops(cell2mat(templates(i)),'MinorAxisLength');
    Ma_A_Length=stats1.MajorAxisLength;
    Mi_A_Length=stats2.MinorAxisLength;
    %disp("MajorAxisLength of : ");
    %disp(Ma_A_Length)
    %disp("MinorAxisLength of T : ");
    %disp(Mi_A_Length)
    t=struct2table(stats);
    t1=struct2table(stats1);
    t2=struct2table(stats2);
    disp(t);
    disp(t1);
    disp(t2);
    
  
end

figure
disp(stats1);

%for i = 1:length(letter)
    %stats=regionprops(cell2mat(templates(i)),'Centroid');
    %centr=stats.Centroid;
    %s = regionprops(letter(i));
    %t=struct2table(s);
    %disp(t)
%end


%*******************************
%t=struct2table(letter);
%disp(t)
%*********************************
%celldisp(templates(20));

      
    %stats=regionprops(cell2mat(templates(20)),'Centroid');
    %centr=stats.Centroid;
    %disp("Centroid of T : ");
    %disp(centr)
    
   %for x = 1: numel(stats)
      %plot(stats(x).Centroid(1),stats(x).Centroid(2),'ro');
   %end
 
 %data=load('templates');
%  disp(data)
 
%  n = numel(templates);
%  disp(n)
 
% p = table(templates,'RowNames');
 %disp(p)

%stats=regionprops(cell2mat(templates(20)),'MajorAxisLength');
%Ma_A_Length=stats.MajorAxisLength;
%disp("MajorAxisLength of T : ");
%disp(Ma_A_Length)

%stats=regionprops(cell2mat(templates(20)),'MinorAxisLength');
%Mi_A_Length=stats.MinorAxisLength;
%disp("MinorAxisLength of T : ");
%disp(Mi_A_Length)

%stats=regionprops(cell2mat(templates(20)),'PixelValues');
%pix_val_t=stats.PixelValues;
%disp("PixelValues of T : ");
%disp(pix_val_t)

%celldisp(templates(21));
%stats=regionprops(cell2mat(templates(21)),'Centroid');
%centr=stats.Centroid;
%disp("Centroid of  : ");
%disp(centr)

%stats=regionprops(cell2mat(templates(21)),'MajorAxisLength');
%Ma_A_Length=stats.MajorAxisLength;
%disp("MajorAxisLength of U : ");
%disp(Ma_A_Length)

%stats=regionprops(cell2mat(templates(21)),'MinorAxisLength');
%Mi_A_Length=stats.MinorAxisLength;
%disp("MinorAxisLength of U : ");
%disp(Mi_A_Length)

%stats=regionprops(cell2mat(templates(21)),'PixelValues');
%pix_val_U=stats.PixelValues;
%disp("PixelValues of U : ");
%disp(pix_val_U)

%MajorAxisLength
%MinorAxisLength
%PixelValues

clear all

