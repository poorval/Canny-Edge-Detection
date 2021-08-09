tic;
clc();
C=imread('cameraman.tif');
%Gaussian Filter
%C=rgb2gray(C);
m=input('Enter size of filter=');
n=input('Sigma=');
G=gauss2d(m,n);
S=uint8(round(conv2(G,C),0));
%SobelX and SobelY
Gx=[-1 0 1;-2 0 2;-1 0 1];
Gy=[1 2 1;0 0 0;-1 -2 -1];
SGx=round(conv2(Gx,S),0);
SGy=round(conv2(Gy,S),0);
Theta=180*(atan2(SGy,SGx))/pi;
[w,h]=size(Theta);
%Adjusting to nearest 0, 45, 90, and 135 degree
for i=1:w
    for j=1:h
        if((Theta(i,j)<22.5 && Theta(i,j)>=0) || (Theta(i,j)<=180 && Theta(i,j)>=157.5) || (Theta(i,j)<-157.5 && Theta(i,j)>-180) || (Theta(i,j)<0 && Theta(i,j)>=-22.5))
            Theta(i,j)=0;
        elseif((Theta(i,j)<67.5 && Theta(i,j)>=22.5) || (Theta(i,j)<-22.5 && Theta(i,j)>=-67.5))
            Theta(i,j)=45;
        elseif((Theta(i,j)<=112.5 && Theta(i,j)>=67.5) || (Theta(i,j)<-67.5 && Theta(i,j)>=-112.5) )
            Theta(i,j)=90;
        elseif((Theta(i,j)<157.5 && Theta(i,j)>=112.5) || (Theta(i,j)<-112.5 && Theta(i,j)>=-157.5))
            Theta(i,j)=135;
        end
    end
end
%Non-Maximum Suppression
mag=sqrt((SGx.^2)+(SGy.^2));
NM = zeros (w, h);
for i=2:w-1
    for j=2:h-1
        if (Theta(i,j)==0)
            if(mag(i,j) == max([mag(i,j), mag(i,j+1), mag(i,j-1)]))
                NM(i,j) = mag(i,j);
            end
        elseif (Theta(i,j)==45)
            if(mag(i,j) == max([mag(i,j), mag(i+1,j-1), mag(i-1,j+1)]))
                NM(i,j) = mag(i,j);
            end
        elseif (Theta(i,j)==90)
            if(mag(i,j) == max([mag(i,j), mag(i+1,j), mag(i-1,j)]))
                NM(i,j) = mag(i,j);
            end
        elseif (Theta(i,j)==135)
            if(mag(i,j) == max([mag(i,j), mag(i+1,j+1), mag(i-1,j-1)]))
                NM(i,j) = mag(i,j);
            end
        end
    end
end
disp("Enter between 0 to 1 (Threshold 1 < Threshold 2)");
Th1=input("Threshold 1=");
Th2=input("Threshold 2=");
Th1=Th1*max(max(NM));
Th2=Th2*max(max(NM));
Result=zeros(w,h);
for i=1:w
    for j=1:h
        if(NM(i,j)<Th1)
            Result(i,j)=0;
        end
        if(NM(i,j)>Th2)
            Result(i,j)=255;
        end
    end
end
subplot(3,2,1),imshow(C);
title("Original Image");
subplot(3,2,2),imshow(S);
title("Gaussian Filtered Image");
subplot(3,2,3),imshow(uint8(SGx));
title("Sobel-X Filtered Image");
subplot(3,2,4),imshow(uint8(SGy));
title("Sobel-Y Filtered Image");
subplot(3,2,5),imshow(NM);
title("Non-Max Suppresed Image");
subplot(3,2,6),imshow(uint8(Result));
title("Resultant Image");
toc;