close all
clear all

n = 50;
h = 1./(n+1);
A = 2*eye(n) - diag( ones(n-1,1), 1) - diag(ones(n-1,1),-1);
A = -A/h/h;


x0 = ones(n,1);
% d/dt x = -Laplace x
[T,X] = ode23s( @(t,x) A*x, linspace(0,1), x0 );

figure(1);
for i = 1:length(T)
    plot(X(i,:));
    drawnow;
    pause
end

%%
[V,S,~] = svd(X', 'econ');
s = diag(S)
sum(s)
figure(2);
plot(cumsum(s)/sum(s));
xlabel('l');
ylabel('ratio');
title('"POD energy"');

figure(3);
semilogy( diag(S) );
title('Singular value decay');
hold on
%%

figure(4);
U = V(:, 1:5);
for k = 1:size(U,2)
    subplot(1,size(U,2),k);
    plot(U(:,k));
    title(sprintf('%d  basis mode', k));
end

%%
AN = U'*A*U;
x0 = rand(n,1);
%x0 = ones(n,1);
[TN, XN] = ode23s( @(t,x) AN*x, linspace(0,1), U'*x0);

%%
figure(5);
for k = 1:length(TN)
    plot([X(k,:)' U*XN(k,:)'])
    pause
end