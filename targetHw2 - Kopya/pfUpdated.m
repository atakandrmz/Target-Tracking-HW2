% hold on;

%----------system definition------------------
    T = 1;
    stepSize = 5
    vx = 5
    vy = 10

    A = [1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1]
    G = [T^2/2 0;T 0;0 T^2/2;0 T];
    C = [1 0 0 0;0 0 1 0]
    
    
    Qtilda = [0.3 0;0 0.1];
    zeroMeanProcessNoise = [0 0]
    zeroMeanMeasNoise = [0 0]
    
    Q = G*Qtilda*G'
    R = [0.1 0;0 0.1];
    
    rng('default')  % For reproducibility
    wk = mvnrnd(zeroMeanProcessNoise,Qtilda,1000)';
    vk = mvnrnd(zeroMeanMeasNoise,R,1000)';

%------------------system model------------------
    %acc = [wx;wy]
    %measNoise = [vx;vy]
    %xState = [x;vx;y;vy]
    %xStateNew = A*xState+G*acc
    %y = C*xState + measNoise

%----------True Position----------------------
    
    x = zeros(4,stepSize); %initialization
    x(2,:) = vx             % constant velocity assumption
    x(4,:) = vy
    y = zeros(2,stepSize);
    
    for i=1:stepSize
    
        x(:,i+1) = A*x(:,i) +G*wk(:,i)
        y(:,i) = C*x(:,i) + vk(:,i)
    end



%------------PF-----------------------
rng('default')  % For reproducibility
    pzgz = eye(4);
    xzgz = x(:,1)
    mu = (xzgz)'
    Sigma = pzgz

xi = mvnrnd(mu,Sigma,1000);
xi = xi'
    
sigmaPtPos = xi
    
    xi = [sigmaPtPos]
    wi = 1/(length(xi));
    wi = ones(length(xi),1)*wi
    
for i =1:length(xi)  
        yi(:,i) = C*(A*xi(:,i)+G*wk(:,i)) + vk(:,i)
end

% for j = 1:stepSize  
        clf;
        hold on;

        
            muBar = [0 0 0 0]

    for i =1:length(xi) 
        
        %yi(:,i) = C*(A*xi(:,i)+G*wk(:,i)) + vk(:,i)
        muBar = muBar + wi(i)*(A*xi(:,i)+G*wk(:,i))'
    end
    
        SigmaBar = zeros(4,4);

    for i =1:length(xi) 
        SigmaBar = SigmaBar + wi(i).*((A*xi(:,i)+G*wk(:,i))'-muBar)'*((A*xi(:,i)+G*wk(:,i))'-muBar)
    end
        mmse1 = [0 0]

    for i =1:length(xi)

        p = mvnpdf((C*xi(:,i)-[xi(1,i) xi(3,i)]),[],R)
        pw(:,i) = p(1)*(wi(i))'
        pw(:,i) = pw(:,i)./(sum(pw,"all"))

        mmse1 = mmse1 + pw(:,i)*[xi(1,i) ; xi(3,i)]'

     end
%     for i =1:length(xi)

    [xiRe,wiRe] = resample(xi,wi)

    xi = xiRe
    pw = wiRe
%     end
plot(xi(1,:),xi(3,:),'r*')
hold on
xi = [yi(1,:) ; vx*ones(1,length(yi));yi(2,:);vy*ones(1,length(yi))]
muBarRe = [0 0 0 0] 
for i =1:length(R) 
    muBarRe = muBarRe + pw(i)*(A*xi(:,i)+G*wk(:,i))'
end

SigmaBarRe = zeros(4,4);

for i =1:length(R) 
    SigmaBarRe = SigmaBarRe + pw(i).*((A*xi(:,i)+G*wk(:,i))'-muBarRe)'*((A*xi(:,i)+G*wk(:,i))'-muBarRe)
end
mmse2 = [0 0]
for i =1:length(xi)  
            mmse2 = mmse2 + pw(:,i)*[xi(1,i) ; xi(3,i)]'

        yi(:,i) = C*(A*xi(:,i)+G*wk(:,i)) + vk(:,i)
end
plot(xi(1,:),xi(3,:),'c*')
hold on;
plot(yi(1,:),yi(2,:),'g*')

plot(mmse1(1,1),mmse1(1,2),'b*')
plot(mmse2(1,1),mmse2(1,2),'b*')

%plot(mmse3(1,1),mmse3(1,2),'b*')

error_ellipse(C*Sigma*C' +R, C*mu',0.99)
error_ellipse(C*SigmaBar*C' +R, C*muBar',0.99)
error_ellipse(C*SigmaBarRe*C'+R, C*(A*muBar'),0.99)
title('Particle Filter')
legend('Initial samples', 'Time Update Result', 'Measurement Updated Time Update Result','MMSE','','Resampling Result')
grid minor
%  for i=1:length(xi)-stepSize
%     
%         xi(:,i+1) = A*xi(:,i) +G*wk(:,i)
%         yi(:,i) = C*xi(:,i) + vk(:,i)
%     end

% xi = [yi(1,:) ; vx*ones(1,length(yi));yi(2,:);vy*ones(1,length(yi))]
% wi = pw(1,:)'
% %
% end


