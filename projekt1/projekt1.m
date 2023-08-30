%% Przygotowanie
% Tr�jk�ty
% import java.util.LinkedList
A1 = [2,0,10, 1];
A2 = [-1, -1,12, 1];
A3 = [-1,1,8, 1];
Akolor = [255,0,0]; % trzeba b�dzie potem podzieli� przez 255
B1 = [0,-2,10, 1];
B2 = [-1,2,11, 1];
B3 = [1,2,9, 1];
Bkolor = [0,255,0];
f = 5;
srodek = [0,0,0];
% o� kamery - pokryta z OZ
% p�aszczyzna obrazu - OXY
%% Cz�� 1 - generacja obrazu
X = 640; Y = 480;
dx = 0.01; dy = 0.01;
tlo = [0,0,0]; % kolor
brzegkolor = [255,255,255];
Im = uint8(zeros(Y,X)); % p�achta ju� w kolorze wybranym
% tr�jk�ty generujemy
Im1 = trojkat(Im, A1, A2, A3, brzegkolor, f, dx, dy);
[PA1,PB1,PC1,PD1] = plane(A1,A2,A3);

Im2 = trojkat(Im, B1, B2, B3, brzegkolor, f, dx, dy);
[PA2,PB2,PC2,PD2] = plane(B1,B2,B3);
%figure(1); imshow(Im); % nic nie wida�
% trzeba dokona� rzutowania
% z-bufor jeszcze ugh
Imm = zmiesc(Im1, Im2, dx, dy, f, PA1, PB1, PC1, PD1, Akolor, PA2, PB2, PC2, PD2, Bkolor);
%figure(2); imshow(Imm); % a tu wida�
%% Cz�� 2 - animacja
% obroty robi� wok� osi OZ
a = 6; % 2 + a/10;
b = 6; % -2 -b/10;
Aroll = (2+a/10)*pi/180.0; % od razu na radiany
Broll = (-2-b/10)*pi/180.0; % b�d� si� kumulowa�y b��dy przybli�enia
Apitch = 0;
Bpitch = 0;
Ayaw = 0;
Byaw = 0;
Aprzesun = [0 0 0];
Bprzesun = [0 0 0];
liczba_klatek = 200; % co najmniej, liczymy od zera
% Po prawdzie wystarczy klatek 139, ale ta liczba jest daleka od ok. 200.
% Wyk�adowca podpowiada ROLL, czyli RPY model obrotu.
% Dobra nasza! Nie musimy obraca� ka�dego punktu tr�jk�ta, tylko rogi
% i od nowa pokolorowa�.
% kopia zapasowa
AA1 = A1;
AA2 = A2;
AA3 = A3;
BB1 = B1;
BB2 = B2;
BB3 = B3;
for L = 0:1
    video = VideoWriter(['czesc',int2str(L+2),'.avi']); % patrzcie, konkatenacja tablic!
    open(video); % tu wypada ostrze�enie o potencjalnie nieograniczonej d�ugo�ci filmu
    for K = 0:liczba_klatek
        
%         if K ==17
%             AA1 = AA1;
%         end
        %if K>0 % moje sklecone na �ywio�
            if L==1
                Bprzesun(3) = sin(K*pi/50);
            end % przesuwamy stopniowo, liniowo - zmieniamy, bo 3 nie dzia�a
            AA1 = ObrPrz(A1, K*Aroll, K*Apitch, K*Ayaw, Aprzesun);% wcze�niej by�o Obr(AA1, ..)
            AA2 = ObrPrz(A2, K*Aroll, K*Apitch, K*Ayaw, Aprzesun);% (.., Aroll, ...)
            AA3 = ObrPrz(A3, K*Aroll, K*Apitch, K*Ayaw, Aprzesun);
            BB1 = ObrPrz(B1, K*Broll, K*Bpitch, K*Byaw, Bprzesun);
            BB2 = ObrPrz(B2, K*Broll, K*Bpitch, K*Byaw, Bprzesun);
            BB3 = ObrPrz(B3, K*Broll, K*Bpitch, K*Byaw, Bprzesun);
        %end
        %przekopiuj
        % Im = uint8(zeros(Y,X)); % analiza wykaza�a, �e Im nie ulega zmianie
        % tr�jk�ty generujemy
        Im1 = trojkat(Im, AA1, AA2, AA3, brzegkolor, f, dx, dy); %problem jest w �rodku tego
        Im2 = trojkat(Im, BB1, BB2, BB3, brzegkolor, f, dx, dy);%by�o Im, ale przy klatce 18 si� wykrzacza
        [PA1,PB1,PC1,PD1] = plane(AA1,AA2,AA3); %nie naprawi�o si�
        [PA2,PB2,PC2,PD2] = plane(BB1,BB2,BB3);
        %imshow(Im1); imshow(Im2); % zast�pione faktycznym debugowaniem
        Imm = zmiesc(Im1, Im2, dx, dy, f, PA1, PB1, PC1, PD1, Akolor, PA2, PB2, PC2, PD2, Bkolor);
        writeVideo(video, Imm);
        if K==0 && L==0 % sklecone na �ywio�
            imwrite(Imm, 'czesc1.png');
        end
        %imshow(Imm); % przyda�o si�
    end
    close(video);
end
%% Cz�� 3 - dodatek
% drugi tr�jk�t oscyluje wok� OZ
% [0,0,sin(K*pi/50)] - wektor przesuni�cia
%% Pozosta�e
% punkty przekszta�camy
function A = rzut(X, f)
if f==0
    print("Ogniskowa r�wna zero");
    return;
end
[a,b] = size(X);
if a>b % dlatemu si� wykrzacza
    X = X';
end
T = [1.0,0,0,0;
    0,1.0,0,0;
    0,0,0,0;
    0,0,-1/f,1.0];
A = (T*X');
% to si� wykrzacza przy zmiennych rzeczywistych, dlaczemu
%A = round(A/A(4));
% Dlaczego bez zaokr�gle�?
% Bo nie musimy zaokr�gla� (plus na slajdzie nie ma o tym wzmianki).
% Patrz ni�ej uwaga.
A = A/A(4);
A = [A(1) A(2) A(4)];

end
% poni�sza funkcja przekopiowana z lab�w 4, bo napisa�em to ju�
% ze wsp�rz�dnych rzeczywistych (zera w centrum)
% do obrazkowych (zera w rogu)
function A = mat2pix(m, n, dx, dy, punkt)
%dlaczego punkt jst nan?
vec = punkt;
T = [1/dx, 0, m/2;
    0 -1/dy, n/2;
    0,0,1];
A = (T*vec');
A = round(A/A(3)); % bez dzielenia przez A(3)
% Tylko w zamianie wsp�rz�dnych
% z wyra�anych liczbami rzeczywistymi
% na te wyra�ane liczbami naturalnymi
% potrzeba zaokr�gla� liczby do ca�o�ci.
% if A(1)>m
%     A(1)=m;
% end
% if A(2)>n
%     A(2)=n;
% end
% if A(1)<0
%     A(1)=0;
% end
% if A(2)<0
%     A(2)=0;
% end
% x = A(1);
% y = A(2);
end
function [x, y] = pix2mat(m, n, dx, dy, punkt)
vec = punkt;
T = [dx, 0, -m*dx/2;
    0, -dy, n*dy/2;
    0,0,1];
A = (T*vec');
A = (A/A(3)); % bez dzielenia przez A(3)
% Wsp�rz�dne matematyczne s� wyra�ane
% liczbami rzeczywistymi.
% Nie trzeba zaokr�gla�.
% if A(1)>m
%     A(1)=m;
% end
% if A(2)>n
%     A(2)=n;
% end
% if A(1)<0
%     A(1)=0;
% end
% if A(2)<0
%     A(2)=0;
% end
x = A(1);
y = A(2);
end
% ukradni�te z zaj�� funkcje nr 1
function Im2 = linia(Im1,yp,xp,yk,xk, kolor)
% kre�li BIA�� LINI� na obrazie kolorowym Im1
% gdzie kolor to tak naprawd� odcie� szaro�ci
% od punktu (yp,xp) do punktu (yk, xk)
% zwraca obraz Im2, na kt�rym ta linia jest ju� narysowana
% algorytm Bresenhama
[Y, X] = size(Im1); % zamieniam
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
% rekurencyjnie? zabrak�o pami�ci. nie idziemy w t� stron�
% a zatem kolejka. Matlab niby ma kolejk�,
% ale nie umiem z niej korzysta�
% ani zaimportowa� z javy, jak radzi stackoverflow,
% wi�c trzeba b�dzie samemu j� zaimplementowa�.
% ech
function Im2 = floodfill(Im,y,x,T,C)
[Y, X] = size(Im);
Im2 = Im;
que_init = 1; que_end = 1;
yy(que_init)=y; xx(que_init)=x;
% if Im2(x,y)~=T % a co je�li?
%     return;
% end % w ko�cu to by�a pozosta�o�� z rekurencyjnej implementacji
while que_init <= que_end % kolejno�� s�siad�w wg wyk�adu 4 str 8
    if xx(que_init)+1 <= X
        % floodfill(Im, y, x+1, T, C);
        if Im2(yy(que_init),xx(que_init)+1) == T % Sprawdzanie, czy prawy s�siad nie jest t�em
            Im2(yy(que_init),xx(que_init)+1) = C; que_end = que_end+1; % Wype�nienie
            yy(que_end)=yy(que_init); xx(que_end) = xx(que_init)+1;
        end
    end
    if xx(que_init)-1>0
        % floodfill(Im, y, x-1, T, C);
        if Im2(yy(que_init),xx(que_init)-1) == T % Sprawdzanie, czy lewy sasiad nie jest tlem
            Im2(yy(que_init),xx(que_init)-1) = C; que_end = que_end+1; % Wypelnienie
            yy(que_end)=yy(que_init); xx(que_end) = xx(que_init)-1;
        end
    end
    
    if yy(que_init)-1 > 0
        % floodfill(Im, y-1, x, T, C);
        if Im2(yy(que_init)-1,xx(que_init)) == T % Sprawdzanie, czy dolny sasiad nie jest tlem
            Im2(yy(que_init)-1,xx(que_init)) = C; que_end = que_end+1; % Wypelnienie
            yy(que_end)=yy(que_init)-1; xx(que_end) = xx(que_init);
        end
    end
    if yy(que_init)+1 <= Y
        % floodfill(Im, y+1, x, T, C);
        if Im2(yy(que_init)+1,xx(que_init)) == T % Sprawdzanie, czy g�rny sasiad nie jest t�em
            Im2(yy(que_init)+1,xx(que_init)) = C; que_end = que_end+1; % Wype�nienie
            yy(que_end)=yy(que_init)+1; xx(que_end) = xx(que_init);
        end
    end
    que_init = que_init+1;
end
end
function [A,B,C,D] = plane(p1,p2,p3)
% Danymi wej�ciowymi sa punkty po�o�one na p�aszczy�nie
% Danymi wyj�ciowymi s� warto�ci A B C D r�wnania
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

% Funkcja obliczaj�ca odleg�o�� od p�aszczyzny
function d = dist2plane(A,B,C,D,x,y,z,l,m,n)
% Danymi wej�ciowymi s� warto�ci A B C D p�aszczyzny, polozenie punktu itd.
% Danymi wyj�ciowymi jest odleg�o�� od p�aszczyzny
% Wed�ug prezki z wyk�adu 4 slajd 29
% Rzutowanie perspektywistyczne

ro = (A*x+B*y+C*z+D)/(A*l+B*m+C*n); % Obliczanie ro wg slajdu 29 after all
d = abs(ro)*sqrt(l^2+m^2+n^2); % Obliczanie odleg�o�ci

end
function Imm = zmiesc(Im1, Im2, dx, dy, f, PA1, PB1, PC1, PD1, Akolor, PA2, PB2, PC2, PD2, Bkolor)
[Y, X] = size(Im1);
Imm = uint8(zeros(Y,X,3)); % nowy kolorowy obrazek
Zbuf = 537467*ones(Y,X); % maksymalna g��bia piksela
for i = 1:X
    for j = 1:Y
        [x,y] = pix2mat(X,Y,dx,dy,[i,j, 1]); % Rzutowanie na wsp�rz�dne matematyczne
        if Im1(j,i) ~= 0 % Sprawdzanie, czy piksel nie jest t�em
            d = dist2plane(PA1,PB1,PC1,PD1,x,y,0,-x,-y,f); % Obliczanie odleg�o�ci od kamery dla pierwszego tr�jkata
            if d<Zbuf(j,i)
                Zbuf(j,i) = d; Imm(j,i,:) = Akolor'; % Wype�nienie kolorem, je�eli jest blizej umiejscowiony
            end
        end
        if Im2(j,i) ~= 0 % Sprawdzanie, czy piksel nie jest tlem
            d = dist2plane(PA2,PB2,PC2,PD2,x,y,0,-x,-y,f); % Obliczanie odleglosci od kamery dla drugiego trojkata
            if d<Zbuf(j,i)
                Zbuf(j,i) = d; Imm(j,i,:) = Bkolor'; % Wypelnienie kolorem, jezeli jest blizej umiejscowiony
            end
        end
    end
end
end
function w = ObrPrz(punkt, roll, pitch, yaw, p)
[a,b] = size(punkt);
if a>b % dlatemu si� wykrzacza
    punkt = punkt';
end
S1 = sin(roll); C1 = cos(roll);
S2 = sin(pitch); C2 = cos(pitch);
S3 = sin(yaw); C3 = cos(yaw);
R = [C1*C2, C1*S2*S3-S1*C3, C1*S2*C3+S1*S3, p(1);
    S1*C2, S1*S2*S3+C1*C3, S1*S2*C3-C1*S3, p(2);
    -S2, C2*S3, C2*C3, p(3);
    0, 0, 0, 1];
w = (R*punkt');
% bez rzutowania, tym zajmie si� tr�jk�t.
end
% dla uproszczenia dajmy gotow� funkcj�
function Im = trojkat(Im1, A1, A2, A3, brzegkolor, f, dx, dy)
[Y, X] = size(Im1);
AA1 = rzut(A1, f);
AA2 = rzut(A2, f);
AA3 = rzut(A3, f);
% potem przekszta�ci� wsp�rz�dne (przesun�� punkt zerowy na koniec)
AA1 = mat2pix(X, Y, dx, dy, AA1);
AA2 = mat2pix(X, Y, dx, dy, AA2);
AA3 = mat2pix(X, Y, dx, dy, AA3);
% a teraz rysujemy tr�jk�t
% niech brzegi b�d� bia�e
% kroki: rysujemy brzegi (algorytmem bresenhama na przyklad)
Im = linia(Im1,AA1(2),AA1(1),AA2(2),AA2(1), brzegkolor(1));
Im = linia(Im,AA2(2),AA2(1),AA3(2),AA3(1), brzegkolor(1));
Im = linia(Im,AA3(2),AA3(1),AA1(2),AA1(1), brzegkolor(1));
% potem wype�niamy
Im = floodfill(Im, round((AA1(2)+AA2(2)+AA3(2))/3.0), round((AA1(1)+AA2(1)+AA3(1))/3.0), 0, 255);%na razie monochrom
end % w trzeciej zmiennej wy�ej by�y jedynki :facepalm: i nie w tej kolejno�ci (to akurat nie)