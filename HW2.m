
clear;
tic
%===============Numerically finding the interpolation threshold N==========

for N = 90:102

x = (0:1/N:1);
y = exp(-400.*(x-0.5).^2);
xi = 0:1/2048:1;
yi = interp1(x,y,xi);
y = exp(-400.*(xi-0.5).^2);

norm(y-yi,inf)
end




%==============Implement the numerical solver for 2D wave equation=========
f = @(x)...
        exp(-400*(x-0.5)^2);
    

   
    
    
    

    
    

%------------------------------fine grid-----------------------------------    
N = 1024;
deltat = 1/(2*N);
deltax = 1/N;
U = zeros(N+1,N+1,N+1);
%--------------Incoporate initial conditions-------------------------------
for i = 1:N+1
    for j = 1:N+1
        U(i,j,2) = deltat*f((i-1)*deltax)*f((j-1)*deltax);
    end
end


for n = 3:N+1
    for i = 2:N
        for j = 2:N
            U(i,j,n) = (deltat^2/deltax^2)*...
                (U(i+1,j,n-1)+U(i,j+1,n-1)+U(i-1,j,n-1)+U(i,j-1,n-1)-4*U(i,j,n-1))...
                -U(i,j,n-2)+2*U(i,j,n-1);
        end
    end
end

%------------------------comparing with fine grid--------------------------
error_vec = zeros(3,1);
for k = 7:9
    N = 2^k;
    deltat = 1/(2*N);
    deltax = 1/N;
    Unew = zeros(N+1,N+1,N+1);

    for i = 1:N+1
        for j = 1:N+1
            Unew(i,j,2) = deltat*f((i-1)*deltax)*f((j-1)*deltax);
        end
    end


    for n = 3:N+1
        for i = 2:N
            for j = 2:N
                Unew(i,j,n) = (deltat^2/deltax^2)*...
                    (Unew(i+1,j,n-1)+Unew(i,j+1,n-1)+Unew(i-1,j,n-1)+Unew(i,j-1,n-1)-4*Unew(i,j,n-1))...
                    -Unew(i,j,n-2)+2*Unew(i,j,n-1);
            end
        end
    end
    error_vec(k-6) = norm(reshape(Unew-U(1:2^(10-k):1025,...
        1:2^(10-k):1025,1:2^(10-k):1025),[],1),inf);
end
plot((7:9),-log2(error_vec),'*')
line((7:9),-log2(error_vec),'LineStyle','-')
line([7,9],[14,18]);

toc








