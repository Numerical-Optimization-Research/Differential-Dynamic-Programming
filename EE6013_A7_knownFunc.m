%% EE6013_A6.m


%% Clean Up

clear all;
close all;
clc


%% Variables

N = 50;
delta = 1;


%% Initialize

x(1,:) = 1/5*pi:(1.5*pi-1/5*pi)/(N-1):1.5*pi;
x(2,:) = 3/8*pi:(1.8*pi-3/8*pi)/(N-1):1.8*pi;

% Where u contains the control signals for N-1 timesteps
u = 1/36*pi:(1/41*pi-1/36*pi)/(N-2):1/4*pi;

timeMax = 10;
deltaT = 0.01;
maxIteration = timeMax/deltaT;

Xrecord1 = zeros(maxIteration,N);
Xrecord2 = zeros(maxIteration,N);

threshold = 0.1;

%% Main

% [x,u] = forwardPass(x,u);
[fx,fu] = f(x,u);
[lx,lu] = costCalcDeriv(x,u);
[lxx,lux,luu] = costCalcSecondDeriv(x);
[k,K] = backwardPass(lx,lu,lxx,lux,luu,fx,fu,N);
[uOut,xNew] = updateU(x,u,k,K,N);


for i = 1:maxIteration
    
    if abs(cost(x(:,1),uOut,N) - cost(x(:,1),u,N)) < threshold
        disp(i);
        return;
        
    elseif cost(x(:,1),uOut,N) < cost(x(:,1),u,N)
        
        u = uOut;
        [fx,fu] = f(x,u);
        [lx,lu] = costCalcDeriv(x,u);
        [lxx,lux,luu] = costCalcSecondDeriv(x);
        [k,K] = backwardPass(lx,lu,lxx,lux,luu,fx,fu,N);
        [uOut,xNew] = updateU(x,u,k,K,N);
        
    elseif cost(x(:,1),uOut,N) >= cost(x(:,1),u,N)
        
%         [fx,fu] = f(x,u);
%         [lx,lu] = costCalcDeriv(x,u);
%         [lxx,lux,luu] = costCalcSecondDeriv(x);
        [k,K] = backwardPass(lx,lu,lxx,lux,luu,fx,fu,N);
        [uOut,xNew] = updateU(x,u,k,K,N);
        
    end
    
    Xrecord1(i,:) = xNew(1,:);    
    Xrecord2(i,:) = xNew(2,:);
    
   
end

% Calculate the difference between the starting and final u.




% plot the space

[X1,X2] = meshgrid(-1.5*pi:0.1:2*pi,-1.5*pi:0.1:2*pi);

fx = sin(X1) + cos(X2) + pi/2;

figure('name','Mesh Plot - P1 - Zoomed Out')
surf(X1,X2,fx)
xlabel('X1');
ylabel('X2');
zlabel('F');

figure('name','Contour Plot - P1 - Zoomed Out')
contour(X1,X2,fx);
xlabel('X1');
ylabel('X2');
title('Differential Dynamic Programming');
grid on
hold on 
plot(x(1,:),x(2,:),'-*');
hold off
% for i = 2:eq
%     if val(i) ~= 0
%         plot(outputX1(i,:),outputX2(i,:));
%     end
% end
%% Functions

function [X,U] = forwardPass(x,u)

    lx(1,:) = cos(x(1,:));
    lx(2,:) = -sin(x(2,:));
    lu(:) = 2*u(:);
    
    lxx(1,:) = -sin(x(1,:));
    lxx(2,:) = -cos(x(2,:));
    lux(1,:) = zeros(length(x(1,:)),1);
    lux(2,:) = zeros(length(x(1,:)),1);
    luu = 2*ones(1,length(u));
    
    X = [lx;lxx;lux];
    U = [lu,luu];

end


function [k,K] = backwardPass(lx,lu,lxx,lux,luu,fx,fu,N)
    Vx = zeros(1,N);
    Vxx = zeros(1,N);
    fxx = zeros(1,N);
    fux = zeros(1,N);
    fuu = zeros(1,N);
    k = zeros(1,N);
    K = zeros(1,N);
    
    for i = N-1:-1:1
        
        Qx(i) = lx(i) + fx(i)'*Vx(i+1);
        Qu(i) = lu(i) + fu(i)'*Vx(i+1);
        
        Qxx(i,:) = lxx(:,i) + fx(i)'*Vxx(i+1)*fx(i) + Vx(i+1)*fxx(i);
        Qux(i) = lux + fu(i)'*Vxx(i+1)*fx(i) + Vx(i+1)*fux(i);
        Quu(i) = luu + fu(i)'*Vxx(i+1)*fu(i) + Vx(i+1)*fuu(i);
        
        k(i) = -inv(Quu(i))*Qu(i);
        K(i) = -inv(Quu(i))*Qux(i);
        Vx(i) = Qx(i) - K(i)'*Quu(i)*k(i);
        Vxx(i) = Qxx(i) - K(i)'*Quu(i)*K(i);
        
    end

end


function [fx,fu] = f(x,u)

    out(1,:) = ones(1,length(u)+1) ;
    out(2,:) = ones(1,length(u)+1) ;
    out_u(3,:) = 2*ones(1,length(u));
    fx = [out(1,:);out(2,:)];
    fu = out_u(3,:);

end


function [Qxx,Qxu,Quu] = hessianCalc(x)


    Qxx = [-sin(x(1,:));-cos(x(2,:))];
    Qxu = 0;
    Quu = 2;

end


function cost = costCalc(x,u,N)

    cost(1) = sin(x(1,1)) + cos(x(2,1)) + u(1)^2;
    for i = 2:N-1
        cost(i) = cost(i-1) + sin(x(1,i)) + cos(x(2,i)) + u(i)^2;
    end
    cost(N) = cost(N-1) + sin(x(1,N)) + cos(x(2,N));

end


function [cost,costU] = costCalcDeriv(x,u)

    cost(1,:) = cos(x(1,:));
    cost(2,:) = -sin(x(2,:)) ;
    costU = 2*u(:);

end



function c = cost(x0,u,N)
    xnext = [0;0];
    c = sin(x0(1)) + cos(x0(2)) + u(1)^2;
    for j = 2:N-1
        
        c = c + sin(xnext(1)) + cos(xnext(2)) + u(j)^2;
        xnext(1) = x0(1) + u(j);
        xnext(2) = x0(2) + u(j);
        
    end
    c = c + sin(xnext(1)) + cos(xnext(2));

end

function [cost,costXU,costU] = costCalcSecondDeriv(x)

    cost(1,:) = -sin(x(1,:));
    cost(2,:) = -cos(x(2,:)) ;
    costXU = 0;
    costU = 2;

end


function [uOut,xNew] = updateU(x,u,k,K,N)
    
    xNew(:,1) = x(:,1);
    for i = 1:N-1
        uOut(i) = u(i) + k(i) + K(i)*max(xNew(:,i)-x(:,i));
        xNew(1,i+1) = xNew(1,i) + uOut(i)^2;
        xNew(2,i+1) = xNew(2,i) + uOut(i)^2;
    end
    
        
end


function xOut = newX(x,u)

    xOut(1,:) = x(1,:) + 2*u(:);
    xOut(2,:) = x(2,:) + 2*u(:);

end

% old function


%% Edits
%
%
%
%
%