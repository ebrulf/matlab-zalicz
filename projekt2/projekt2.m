%% Przygotowanie
% Tr�jk�ty
W =[0,0,10,1];
A=[-22,0,20,1];
B=[11,19,20,1];
C=[11,-19,20,1];
WABeps = [1, 0.5, 0.25];
WACeps = [0.25, 1, 0.5];
WBCeps = [0.5, 0.25, 1];
Iamb = [30, 30, 30];
Imax = [255, 255, 0];
% do cz�ci 1
eye = [65, 0, 0, 1];
plight1 = [10, 0, 0, 1];
plight2 = [-7, 0, 0, 1];
m = 30;
alpha = 0.001;
brzegkolor = [255,255,255]; % legacy
% p�aszczyzna obrazu - OXY
%% Cz�� 1 - generacja czworo�cianu
X = 640; Y = 480;
dx = 0.1; dy = 0.1;%wcze�niej by�a 1/100 i si� dziwi�, �e wykracza poza
tlo = [0,0,0]; % kolor

%tr�jk�ty rysowa� umiemy
%tym razem rzutowanie r�wnoleg�e
Im = uint8(zeros(Y,X));
%eksperyment
%[~, plight1, eye, plight2] = trojkat(Im, plight1, eye, plight2, brzegkolor, dx, dy);
% nieudany, wykrzacza si� na lambercie
[Im1, WW, AA, BB] = trojkat(Im, W, A, B, brzegkolor, dx, dy);
[Im2, ~,~,CC] = trojkat(Im, W, A, C, brzegkolor, dx, dy);
[Im3, ~,~,~] = trojkat(Im, W, B, C, brzegkolor, dx, dy);
[Im4,~,~,~] = trojkat(Im, A, B, C, brzegkolor, dx, dy);% to jest zas�oni�te i tak
% p�aszczyzny

Imm = zmiesc(Im1, Im2, Im3, dx, dy, WW, AA, BB, CC, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight1, m, alpha, 0);
% najpierw o�wietli�, potem z-buforowa�?
figure(1); imshow(Im1);
figure(2);imshow(Im2);
figure(3);imshow(Im3);
figure(4);imshow(Im4);
figure(5);imshow(Imm);
imwrite(Imm, 'lambert1.png');
Imm = zmiesc(Im1, Im2, Im3, dx, dy, W, A, B, C, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight2, m, alpha, 0);
figure(6);imshow(Imm);%teraz po poprawce mamy dwa r�ne obrazki
imwrite(Imm, 'lambert2.png');
Imm = zmiesc(Im1, Im2, Im3, dx, dy, W, A, B, C, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight1, m, alpha, 1);
figure(7);imshow(Imm);%teraz po poprawce mamy dwa r�ne obrazki
imwrite(Imm, 'phong1.png');
Imm = zmiesc(Im1, Im2, Im3, dx, dy, W, A, B, C, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight2, m, alpha, 1);
figure(8);imshow(Imm);%teraz po poprawce mamy dwa r�ne obrazki
imwrite(Imm, 'phong2.png');
%% Cz�� 2
% animacja
dt = 0.01; %od 0 do 1
%to b�dzie mno�nik
liczba_klatek = 200;
for L = 0:1
    video = VideoWriter(['czesc',int2str(L+2),'.avi']); % patrzcie, konkatenacja tablic!
    open(video); % tu wypada ostrze�enie o potencjalnie nieograniczonej d�ugo�ci filmu
    for K = 0:liczba_klatek
        if L==1
            eye = [60+5*cos(2*pi*K*dt), 60+10*sin(2*pi*K*dt), 0, 1]; % w poprzednim by� zwyk�y sin
        end
        plight = [10*cos(2*pi*K*dt), 10*sin(2*pi*K*dt), 0, 1];
        [Im1, WW, AA, BB] = trojkat(Im, W, A, B, brzegkolor, dx, dy);
        [Im2, ~,~,CC] = trojkat(Im, W, A, C, brzegkolor, dx, dy);
        [Im3, ~,~,~] = trojkat(Im, W, B, C, brzegkolor, dx, dy);
        Imm = zmiesc(Im1, Im2, Im3, dx, dy, W, A, B, C, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight2, m, alpha, L);
        writeVideo(video, Imm);
        %imshow(Imm); % przyda�o si�
    end
    close(video);
end
%% Cz�� dla ambitnych
%czym jest HSV w Matlabie?
for L = 0:1
    video = VideoWriter(['czesc4',int2str(L+2),'.avi']); % patrzcie, konkatenacja tablic!
    open(video); % tu wypada ostrze�enie o potencjalnie nieograniczonej d�ugo�ci filmu
    for K = 0:liczba_klatek
        if L==1
            eye = [60+5*cos(2*pi*K*dt), 60+10*sin(2*pi*K*dt), 0, 1]; % w poprzednim by� zwyk�y sin
        end
        plight = [10*cos(2*pi*K*dt), 10*sin(2*pi*K*dt), 0, 1];
        if K<=100
        Imax = round(hsv2rgb([1-K*dt, 1, 1])*255); % to za�atwi konwersj�
        else
        Imax = round(hsv2rgb([1-(K-100)*dt, 1, 1])*255); % to za�atwi konwersj�
        end
        [Im1, WW, AA, BB] = trojkat(Im, W, A, B, brzegkolor, dx, dy);
        [Im2, ~,~,CC] = trojkat(Im, W, A, C, brzegkolor, dx, dy);
        [Im3, ~,~,~] = trojkat(Im, W, B, C, brzegkolor, dx, dy);
        Imm = zmiesc(Im1, Im2, Im3, dx, dy, W, A, B, C, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight2, m, alpha, L);
        writeVideo(video, Imm);
        %imshow(Imm); % przyda�o si�
    end
    close(video);
end
%% Pozosta�e
function [Im, AA1, AA2, AA3] = trojkat(Im1, A1, A2, A3, brzegkolor, dx, dy)
[Y, X] = size(Im1);
%a co je�li najpierw przekszta�c� wsp�rz�dne, a potem zrzutuj�?
AA1 = rzutR3(A1);
AA2 = rzutR3(A2);
AA3 = rzutR3(A3);
AA1 = mat2pix(X, Y, dx, dy, AA1);
AA2 = mat2pix(X, Y, dx, dy, AA2);% tutaj to eksploduje
AA3 = mat2pix(X, Y, dx, dy, AA3);

% potem przekszta�ci� wsp�rz�dne (przesun�� punkt zerowy na koniec)
%albo zrobimy na chama, a co si� stanie
% AA1 = [AA1(1), AA1(2), AA1(4)];
% AA2 = [AA2(1), AA2(2), AA2(4)];
% AA3 = [AA3(1), AA3(2), AA3(4)];
% niech brzegi b�d� bia�e
% kroki: rysujemy brzegi
Im = linia(Im1,AA1(2),AA1(1),AA2(2),AA2(1), brzegkolor(1));
Im = linia(Im,AA2(2),AA2(1),AA3(2),AA3(1), brzegkolor(1));
Im = linia(Im,AA3(2),AA3(1),AA1(2),AA1(1), brzegkolor(1));
% potem wype�niamy
Im = floodfill(Im, round((AA1(2)+AA2(2)+AA3(2))/3.0), round((AA1(1)+AA2(1)+AA3(1))/3.0), 0, 255);
end
function A = rzutR(X)
[a,b] = size(X);
if a>b
    X = X';
end
T = [1.0,0,0;
    0,1.0,0;
    0,0,1.0];
A = (T*X');
A = A/A(3);
end
function A = rzutR3(X)
[a,b] = size(X);
if a>b
    X = X';
end
T = [1.0,0,0,0;
    0,1.0,0,0;
    0,0,0,0;
    0,0,0,1.0];
A = (T*X');
A = A/A(4);
A = [A(1) A(2) A(4)];
end
function A = mat2pix(m, n, dx, dy, punkt)
vec = punkt;
T = [1/dx, 0, m/2;
    0 -1/dy,  n/2;
    0,0,1];
A = (T*vec');
A = round(A/A(3));
end
function A = mat2pix3(m, n, dx, dy, punkt)
vec = punkt;
T = [1/dx, 0, 0, m/2;
    0 -1/dy, 0, n/2;
    0,0,0,0;
    0,0,0,1];
A = (T*vec');
A = round(A/A(4));% chwilowa zamiana
end
function [x, y] = pix2mat(m, n, dx, dy, punkt)%to jest najcz�ciej wywo�ywana funkcja
vec = punkt;
T = [dx, 0, -m*dx/2;
    0, -dy, n*dy/2;
    0,0,1];
A = (T*vec');
% A = (A/A(3)); % uprosci�
x = A(1)/A(3);
y = A(2)/A(3);
end
function [x, y] = pix2mat3(m, n, dx, dy, punkt)
vec = punkt;
T = [dx, 0,0, -m*dx/2;
    0, -dy,0, n*dy/2;
    0,0,0,0;
    0,0,0,1];
A = (T*vec');
A = (A/A(4));
x = A(1);
y = A(2);
end
function Im2 = linia(Im1,yp,xp,yk,xk, kolor)
[Y, X] = size(Im1);
Im2 = Im1;
kolor = kolor(1);
if (abs(yk-yp) > abs(xk-xp))                   % stroma linia
    x0 = yp;y0 = xp; x1 = yk;y1 = xk;          % zamiana wsp�rzednych
    zamiana_XzY = 1;
else
    x0 = xp;y0 = yp; x1 = xk;y1 = yk;
    zamiana_XzY = 0;
end
if(x0 >x1)                                      % te� zamiana wsp�rz�dnych
    temp1 = x0; x0 = x1; x1 = temp1;
    temp2 = y0; y0 = y1; y1 = temp2;
end
dx = abs(x1 - x0);                              % odleg�o�� x
dy = abs(y1 - y0);                               % odleg�o�� y
sx = sign(x1 - x0);                              % kierunek x
sy = sign(y1 - y0);                              % dodatnie nachylenie
x = x0; y = y0;
% inicjalizacja
if x>0 && x<=X && y>0 && y<=Y                    % rysowanie punktu
    if (zamiana_XzY ==0)
        Im2(y,x) = kolor;           % poprawka na kolor cofni�ta
    else
        Im2(x,y) = kolor;
    end
end
error = 2*dy - dx;                              % inicjalizacja b��du (z wiki)
for i = 0:dx-1
    error = error + 2*dy;                        % modyfikacja b��du
    if (error >0)
        y = y +1*sy;
        error = error - 2*dx;
    end
    x = x + 1*sx;                                % zwi�kszamy x
    % rysowanie punktu
    if (zamiana_XzY ==0)
        if x>0 && x<=X && y>0 && y<=Y
            Im2(y,x) =  kolor;
        end
    else
        if x>0 && x<=Y && y>0 && y<=X
            Im2(x,y) =  kolor;
        end
    end
end
end
function Im2 = floodfill(Im,y,x,T,C)
[Y, X] = size(Im);
%co je�li x i y s� negatywne?
if x<=0
    x=1;
end
if y<=0
    y=1;
end
if x>X
    x=X;
end
if y>Y
    y=Y;
end
Im2 = Im;
que_init = 1; que_end = 1;
yy(que_init)=y; xx(que_init)=x;
while que_init <= que_end
    %zamiana, bo wykrzaczy�o si� na pierwszem
    if xx(que_init)-1>0
        if Im2(yy(que_init),xx(que_init)-1) == T % Sprawdzanie, czy lewy sasiad nie jest tlem
            Im2(yy(que_init),xx(que_init)-1) = C; que_end = que_end+1; % Wypelnienie
            yy(que_end)=yy(que_init); xx(que_end) = xx(que_init)-1;
        end
    end
    
    if yy(que_init)-1 > 0
        if Im2(yy(que_init)-1,xx(que_init)) == T % Sprawdzanie, czy dolny sasiad nie jest tlem
            Im2(yy(que_init)-1,xx(que_init)) = C; que_end = que_end+1; % Wypelnienie
            yy(que_end)=yy(que_init)-1; xx(que_end) = xx(que_init);
        end
    end
    if xx(que_init)+1 <= X
        if Im2(yy(que_init),xx(que_init)+1) == T % Sprawdzanie, czy prawy s�siad nie jest t�em
            Im2(yy(que_init),xx(que_init)+1) = C; que_end = que_end+1; % Wype�nienie
            yy(que_end)=yy(que_init); xx(que_end) = xx(que_init)+1;
        end
    end
    if yy(que_init)+1 <= Y
        if Im2(yy(que_init)+1,xx(que_init)) == T % Sprawdzanie, czy g�rny sasiad nie jest t�em
            Im2(yy(que_init)+1,xx(que_init)) = C; que_end = que_end+1; % Wype�nienie
            yy(que_end)=yy(que_init)+1; xx(que_end) = xx(que_init);
        end
    end
    que_init = que_init+1;
end
end
function [A,B,C,D] = plane(p1,p2,p3)
A = [p1(2), p1(3), 1;
    p2(2), p2(3), 1;
    p3(2), p3(3), 1];
B = [p1(1), p1(3), 1;
    p2(1), p2(3), 1;
    p3(1), p3(3), 1];
C = [p1(1), p1(2), 1;
    p2(1), p2(2), 1;
    p3(1), p3(2), 1];
D = [p1(1), p1(2), p1(3);
    p2(1), p2(2), p2(3);
    p3(1), p3(2), p3(3)];
A = det(A);
B = -det(B);
C = det(C);
D = -det(D);
end
function d = dist2planeR(A,B,C,D,x,y)
%Rzutowanie r�noleg�e
ro = (A*x+B*y+D)/(C);
d = -ro;
end
function s = cosvec(v, w)
s = dot(v, w)/(norm(v)*norm(w)); %lekkie wyszukanie w guglu
end % nie pomog�o, wi�c zaimplementowa�em sam
function Iamb = ambience(amb, eps)
Iamb = eps.*amb;
end
function Ilam = lambert(p, Imax, eps, plight, alpha, n)
Ilam = eps;
Ilam = (Ilam.*Imax)/(1+alpha*sum((p-plight).^2));
Ilam = Ilam * cosvec(n, plight-p);%a jak dam norm�?
end % tak! o to chodzi�o!
function Iph = phong(p, Imax, eps, plight, alpha, eye, n, m)
r = 2 * dot(plight-p,n).*n-(plight-p);
Iph = eps;
Iph = (Iph.*Imax)/(1+alpha*sum((p-plight).^2));
Iph = Iph.*1; % tu by�oby co� wi�cej
Iph = Iph * cosvec(r, eye-p).^m;
end
function I = oswietlenie(eps, amb, p, Imax, plight, alpha, n, eye, m, phon)
I = ambience(amb, eps) + lambert(p, Imax, eps, plight, alpha,n);
if phon == 1 % upsik
    I = I + phong(p, Imax, eps, plight, alpha, eye, n, m);
end
end
function Imm = zmiesc(Im1, Im2, Im3, dx, dy, W, A, B, C, WABeps,WACeps,WBCeps, Iamb, Imax, eye, plight, m, alpha, phong)
%tym si� r�ni od poprzedniego, �e s� cztery tr�jk�ty i o�wietlenie
[Y, X] = size(Im1);
Imm = uint8(zeros(Y,X,3)); % nowy kolorowy obrazek
Zbuf = 537467*ones(Y,X); % maksymalna g��bia piksela
%przeniesione tutaj, aby zmniejszy� ilo�� wymaganych zmiennych
[PA1,PB1,PC1,PD1] = plane(W,A,B);%to pewnie b��d ze zwyk�ym W A B
[PA2,PB2,PC2,PD2] = plane(W,A,C);
[PA3,PB3,PC3,PD3] = plane(W,B,C);
for i = 1:X
    for j = 1:Y
        [x,y] = pix2mat(X,Y,dx,dy,[i,j,1]); % Rzutowanie na wsp�rz�dne matematyczne
        if Im1(j,i) ~= 0 % Sprawdzanie, czy piksel nie jest t�em
            d = dist2planeR(PA1,PB1,PC1,PD1,x,y); % Inne parametry
            if d<Zbuf(j,i)
                Zbuf(j,i) = d;
                %tu o�wietlenie (drugi wektor to normalna p�aszczyzny
                % brzegi pomijam
                Imm(j,i,:) = oswietlenie(WABeps, Iamb, [x,y,d,1], Imax, plight, alpha, [PA1,PB1,PC1, 1], eye, m, phong);
                % Wype�nienie kolorem, je�eli jest blizej umiejscowiony
            end
        end
        if Im2(j,i) ~= 0 % Sprawdzanie, czy piksel nie jest tlem
            d = dist2planeR(PA2,PB2,PC2,PD2,x,y);
            if d<Zbuf(j,i)
                Zbuf(j,i) = d;
                Imm(j,i,:) = oswietlenie(WACeps, Iamb, [x,y,d,1], Imax, plight, alpha, [PA2,PB2,PC2, 1],eye, m, phong);
            end
        end
        if Im3(j,i) ~= 0 % Sprawdzanie, czy piksel nie jest t�em
            d = dist2planeR(PA3,PB3,PC3,PD3,x,y);
            if d<Zbuf(j,i)
                Zbuf(j,i) = d;
                Imm(j,i,:) = oswietlenie(WBCeps, Iamb, [x,y,d,1], Imax, plight, alpha, [PA3,PB3,PC3, 1],eye, m, phong);
            end
        end
    end
end
end