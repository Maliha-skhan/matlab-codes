function [U,dt] = diffusion(Lx,Ly,J,K,Dup,Ddown,Dright,Dleft,U0,tf)
dx = 2 * Lx / J;
dy = 2 * Ly / K;
[x,y] = ndgrid(-Lx:dx:Lx, -Ly:dy:Ly);
U = U0;
maxD = max(max(Dup,Ddown),max(Dleft,Dright));
dt0 = 0.25 * min(dx,dy)^2 / max(max(maxD));
N = ceil(tf / dt0);
dt = tf / N;
mu_x = dt / (dx*dx);
mu_y = dt / (dy*dy);
for n=1:N
U(2:J,2:K) = U(2:J,2:K) + ...
mu_y * Dup .* ( U(2:J,3:K+1) - U(2:J,2:K) ) - ...
mu_y * Ddown .* ( U(2:J,2:K) - U(2:J,1:K-1) ) + ...
mu_x * Dright .* ( U(3:J+1,2:K) - U(2:J,2:K) ) - ...
mu_x * Dleft .* ( U(2:J,2:K) - U(1:J-1,2:K) );
end