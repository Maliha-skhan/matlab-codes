clear all; close all;
pause on;
Nmax = 60;
delta = 1;
M = 1000;
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)*0.01 scrsz(4)*0.15 scrsz(3)*0.5 scrsz(4)*0.75])
for N=1:delta:Nmax
%figure('Position',[scrsz(3)*0.01 scrsz(4)*0.15 scrsz(3)*0.5 scrsz(4)*0.75])
    clf;
%    x = rand(N,M);
    x = normrnd(0.5,0.3,N,M);
%    x = (rand(N,M)-0.5)/10+0.05+(randi(2,N,M)-1)*0.9;
%    x = randi(2,N,M)-1;
    for m=1:M
        y(m) = mean(x(:,m));
    end
    [h xbin]=hist(y,20);
    hist(y,20);
    axis([0 1 0 1.1*max(h)]);
    sig(N) = std(y);
    xlabel('Mean of N random numbers');
    ylabel('Frequency');
    hold on;
    gsigma = sig(1)/sqrt(N);
    xx = xbin(1)+[0:0.01:100]*(xbin(end)-xbin(1));
%    yy = M*(xbin(2)-xbin(1))/(sqrt(2*pi)*gsigma).*gaussmf(xbin,[gsigma 0.5]);
    yy = M*(xbin(2)-xbin(1))/(sqrt(2*pi)*gsigma).*gaussmf(xx,[gsigma 0.5]);
    plot(xx,yy,'-r','linewidth',4);
    txt = sprintf('N=%d, std=%5.3f',N,sig(N));
    legend(txt,'gauss with sigma = sigma(1)/sqrt(N)');
    pause(3/N);
    hold off;
end

figure('Position',[scrsz(3)*0.6 scrsz(4)*0.3 scrsz(3)*0.4 scrsz(4)*0.6])
xx = (1:Nmax);
plot(xx,sig,'linewidth',4);
hold on;
plot(xx,sig(1)./sqrt(xx),'-r','linewidth',2);
legend('st.dev. of mean of N random measurements', 'sigma(1)/sqrt(N)');
xlabel('Number to average, N');
ylabel('Standard deviation of distribution');
