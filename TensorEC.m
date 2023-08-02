function [S0,Z,obj]=TensorEC(H,c,alpha,beta)
V = size(H,3);
N = size(H,2); %sample number

for v=1:V
    S{v} = eye(N);
    Z{v} = eye(N);
    Y1{v} = zeros(N,N);
    Y2{v} = zeros(N,N);
    A{v} = eye(N,N);
end
w = 1/V*ones(1,V); %1*v
Isconverg = 0;epson = 1e-8;
iter = 0;
mu = 10e-6; max_mu = 10e10; pho_mu = 3;
rho = 0.01; max_rho = 10e12; pho_rho = 3;
tic;

sX = [N, N, V];

while(Isconverg == 0)
    fprintf('----processing iter %d--------\n', iter+1);
    for v=1:V
        %1 update Z^k and S^k
        B1 = Z{v} - Y1{v}/mu;
        B2 = A{v} - Y2{v}/rho;
        B12 = (2 * w(1,v) * H(:,:,v) + mu * B1 + rho * B2)/(2 * w(1,v) + rho + mu);
        B = (B12 + B12')/2;
        S{v} = project_fantope(B,c);
        Z{v} = prox_l1((S{v}+Y1{v}/mu),alpha/mu);
        Y1{v} =  Y1{v} + mu*(S{v} - Z{v});
    end
    
    
    %update G
    S_tensor = cat(3, S{:,:});
    Y2_tensor = cat(3, Y2{:,:});
    s = S_tensor(:);
    y2 = Y2_tensor(:);
    [a, ~] = wshrinkObj(s+1/rho*y2,beta/rho,sX,0,3);
    A_tensor = reshape(a, sX); 
    
    % update omega
    for i = 1:V
        w(1,i) = 0.5/norm(H(:,:,i)-S{i},'fro');
    end
    
    %update variables
    Y2_tensor = Y2_tensor + rho*(S_tensor - A_tensor);
    
    %record the iteration information
    %history.objval(iter+1)   =  objV;
    
    %% coverge condition
    Isconverg = 1;
    maxvalue = 0;
    for v=1:V
        A{v} = A_tensor(:,:,v);
        Y2{v} = Y2_tensor(:,:,v);
        if (norm(S{v}-A{v},inf)>epson)
            history.norm_S_A = norm(S{v}-A{v},inf);
            maxvalue = max(maxvalue,history.norm_S_A);
            %%fprintf('@[%d], history.norm_S_A %7.10f    \n', v, history.norm_S_A);
            Isconverg = 0;
        end
        if (norm(S{v}-Z{v},inf)>epson)
            history.norm_S_Z = norm(S{v}-Z{v},inf);
            maxvalue = max(maxvalue,history.norm_S_Z);
            %%fprintf('#history.norm_S_Z %7.10f    \n', history.norm_S_Z);
            Isconverg = 0;
        end
    end
    history.objval(iter+1)=maxvalue;
    if (iter>30)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
    rho = min(rho*pho_rho, max_rho);
end
S0 = zeros(N);
for v=1:V
    S0 = S0 + abs(Z{v})+abs(Z{v}');
end
S0 = S0 -diag(diag(S0));
obj =history.objval;
