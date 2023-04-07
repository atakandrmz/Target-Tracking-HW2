function [x,w]=resample(x,w)
    % Multinomial sampling with sort
    N = 4
   % Multinomial sampling with sort
u=rand(N, 1);
wc=cumsum(w) ;
wc=wc/wc(N) ;
[dum,ind1]=sort([u;wc]);
ind2=find(ind1<=N);
ind=ind2-(0:N-1)';
x=x(ind,:) ;
w=ones (1, N) ./N;
end