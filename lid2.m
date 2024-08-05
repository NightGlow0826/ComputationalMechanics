%% Initial
clear; clc;
%40, 80, 160
dim = 40;
h = 1/dim; % 根据格子数计算h值
x = 0:h:1; y = 0:h:1;
% 100, 400, 1000
Re = 1000;
rho = 0.05;
tol = 1e-5;
ite = 0;
len = size(x, 2);
psi = zeros(len, len);
w = zeros(len, len);
psi_new=psi;
w_new=w;

while 1
    for i=2:len-1
        for j=2:len-1
            psi_new(i,j)=0.25*rho*(w(i,j)*h^2+psi(i+1,j)+psi_new(i-1,j)+psi(i,j+1)+psi_new(i,j-1))...
                +(1-rho)*psi(i,j);
            w_new(i,j)=0.25*rho*(w(i+1,j)+w_new(i-1,j)+w(i,j+1)+w_new(i,j-1)...
                +0.25*Re*((psi(i,j+1)-psi_new(i,j-1))*(w(i+1,j)-w_new(i-1,j))...
                -(psi(i+1,j)-psi_new(i-1,j))*(w(i,j+1)-w_new(i,j-1))))...
                +(1-rho)*w(i,j);
        end
    end
    w_new(:,1)=-2*(psi(:,2)-psi(:,1))/h^2;
    w_new(:,len)=-2*(psi(:,len-1)-psi(:,len))/h^2;
    w_new(len,:)=-2*(psi(len-1,:)-psi(len,:))/h^2;
    w_new(1,:)=-2*(psi(2,:)-psi(1,:)+h)/h^2;
    ite=ite+1;
    cur_err = max(max(max(abs(w_new-w))),max(max(abs(psi_new-psi))));
    if cur_err <=tol
        break
    else
        psi=psi_new; w=w_new;
    end
end
%% PLot
er=max(max(psi_new)); contourf(x,y,flipud(psi_new))
[t,s]=title(['Re = ',num2str(Re),', ','err = ',num2str(tol),', ',num2str(dim),' X ',num2str(dim)]); t.FontSize=14; s.FontSize=14; colorbar
v=-0.005:0.00025:er;
hold on
contour(x,y,flipud(psi_new),v,'k')

disp('done')
u = zeros(len);
v = zeros(len);
u(1) = -1;
mid = floor(len/2);
for i = 2:len-1
    u(i) = (psi(i+1, mid) - psi(i-1, mid)) / (2*h);
    v(i) = (psi(mid, i+1) - psi(mid, i-1)) / (2*h);
end


figure
plot(u(:, 1), y)
% plot(x, v(:, 1))
%
[t,s]=title(['Re = ',num2str(Re),',     ',  num2str(dim),' X ',num2str(dim)]); t.FontSize=14; s.FontSize=14;
xlabel('u')
ylabel('y')
%
% xlabel('x')
% ylabel('u')



