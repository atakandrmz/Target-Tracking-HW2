hold on;

nx = 2;
% mu = [2 3];
% Sigma = [1 1.5; 1.5 3];
mu = [-1 2];
muMeas =[0 0];

Sigma = [0.3 0.05; 0.05 0.1];
SigmaMeas = [0.1 0.05; 0.05 0.1];
rng('default')  % For reproducibility
acc = normrnd(0,0.1)
R = mvnrnd(mu,Sigma,1000);
% plot(R(:,1),R(:,2),'+')
v = 10;
w0 = 0.3;
sigmaPt0 = mu;
% ------------------UT sqrtm---------------------
% sigmaPtPos = mu + sqrtm((nx/(1-w0).*Sigma)); 
% sigmaPtNeg = mu - sqrtm((nx/(1-w0).*Sigma)); 
% % ------------------UT cholcov---------------------
% sigmaPtPos = mu + cholcov((nx/(1-w0).*Sigma));
% sigmaPtNeg = mu - cholcov((nx/(1-w0).*Sigma));
% 
% xi = [sigmaPt0;sigmaPtPos;sigmaPtNeg]
% wi = [w0;w0/2;w0/2;w0/2;w0/2];%nx/(1-w0);nx/(1-w0);nx/(1-w0);nx/(1-w0)];
% 
% wi = [w0;(1-w0)/(2*nx);(1-w0)/(2*nx);(1-w0)/(2*nx);(1-w0)/(2*nx)];
% 
% muBar = [0 0 0 0]
% for i =1:2*nx+1 %length(R)+1 
%     muBar = muBar + wi(i).*g(xi(i,1),v,xi(i,2),v,acc)'
% end
% 
% SigmaBar = zeros(4,4);
% 
% for i =1:2*nx+1 %length(R)+1 
%     SigmaBar = SigmaBar + wi(i).*(g(xi(i,1),v,xi(i,2),v,acc)'-muBar)'*(g(xi(i,1),v,xi(i,2),v,acc)'-muBar)
% end

%----------------Monte Carlo--------------------

sigmaPtPos = R

xi = [sigmaPtPos]
wi = 1/(length(R));

muBar = [0 0 0 0]
for i =1:length(R) 
    muBar = muBar + wi.*g(xi(i,1),v,xi(i,2),v,acc)'
end

SigmaBar = zeros(4,4);

for i =1:length(R) 
    SigmaBar = SigmaBar + wi.*(g(xi(i,1),v,xi(i,2),v,acc)'-muBar)'*(g(xi(i,1),v,xi(i,2),v,acc)'-muBar)
end

% ------------------EKF TT1---------------------
% muBar = g(mu(1),mu(2))
% SigmaBar = gBar(mu(1),mu(2))*Sigma*(gBar(mu(1),mu(2)))'
% 
% error_ellipse(Sigma, mu,0.99)
% hold on
% error_ellipse(SigmaBar, muBar,0.99,'style', '--')

%---------------UT----------------
% hold on;
% plot(sigmaPtPos(:,1),sigmaPtPos(:,2),'r*')
% plot(sigmaPtNeg(:,1),sigmaPtNeg(:,2),'g*')
% plot(sigmaPt0(:,1),sigmaPt0(:,2),'m*')
% error_ellipse(Sigma, mu,0.99)
% error_ellipse([SigmaBar(1,1) 0;0 SigmaBar(3,3)], [muBar(1) muBar(3)],0.99,'style', '--')
% error_ellipse(Sigma + SigmaMeas, h(muBar(1),muBar(3),SigmaMeas),0.99)
% 
% utPos = g(sigmaPt0(1,1),v,sigmaPt0(1,2),v,acc)'
% utPos = [utPos(1) utPos(3)]
%     for i =1:nx %length(R)+1 
%         utPosNew = g(sigmaPtPos(i,1),v,sigmaPtPos(i,2),v,acc)'
%         utPosNew = [utPosNew(1) utPosNew(3)]
%         utNegNew = g(sigmaPtNeg(i,1),v,sigmaPtNeg(i,2),v,acc)'
%         utNegNew = [utNegNew(1) utNegNew(3)]
%         utPos = [utPos;utPosNew ;utNegNew]
%     end
% 
% plot(utPos(:,1),utPos(:,2),'r+')
% hold on
% % plot(utNeg(:,1),utNeg(:,1),'g+')
% % plot(ut0(1),ut0(2),'m+')
% grid minor
% % title('Unscented Transform with sqrtm function')
% title('Unscented Transform with cholcov function')

%--------------Monte Carlo------------------
hold on;
plot(sigmaPtPos(:,1),sigmaPtPos(:,2),'r*')
plot(sigmaPt0(:,1),sigmaPt0(:,2),'m*')
error_ellipse(Sigma, mu,0.99)
error_ellipse([SigmaBar(1,1) 0;0 SigmaBar(3,3)], [muBar(1) muBar(3)],0.99,'style', '--')

utPos = g(sigmaPt0(1,1),v,sigmaPt0(1,2),v,acc)'
utPos = [utPos(1) utPos(3)]
    for i =1:length(R)
        utPosNew = g(sigmaPtPos(i,1),v,sigmaPtPos(i,2),v,acc)'
        utPosNew = [utPosNew(1) utPosNew(3)]
        utPos = [utPos; utPosNew]

    end

plot(utPos(:,1),utPos(:,2),'g+')
hold on
grid minor

legend('Initial Random Data','-','Transformed Data','--')

