

function I = cross_field_re(I0, G, lambda, beta, eps, eta_sqr, phi_alpha, phi_eps, iternum)

show = 0;

S = ones(size(I0));
I = I0;

[m, n] = size(I);

Cx = get_Cx(m, n); Cy = get_Cy(m, n);
Cxt = get_Cxt(m, n); Cyt = get_Cyt(m, n);

Gx = cal(Cx, G); Gy = cal(Cy, G);

Px = 1 ./ (msign(Gx) .* max(abs(Gx), eps));
Py = 1 ./ (msign(Gy) .* max(abs(Gy), eps));

G2 = Gx .* Gx + Gy .* Gy;
S1 = eta_sqr ./ (G2 + 2 * eta_sqr);
S2 = (G2 + eta_sqr) ./ (G2 + 2 * eta_sqr);

Gs = sqrt(G2);
Vx = Gx ./ max(abs(Gs), eps);
Vy = Gy ./ max(abs(Gs), eps);

if show
    figure;
    subplot(221); mesh(S1); view(2); colorbar;
    subplot(222); mesh(S2); view(2); colorbar;
    subplot(223); mesh(Vx); view(2); colorbar;
    subplot(224); mesh(Vy); view(2); colorbar;
end

ls1 = zeros(1, iternum);
ls2 = zeros(1, iternum);

li1 = zeros(1, iternum);
li2 = zeros(1, iternum);

for i = 1:iternum

    Ix = cal(Cx, I); Iy = cal(Cy, I);

    Ax = phi(S - Px .* Ix, phi_alpha, phi_eps);
    Ay = phi(S - Py .* Iy, phi_alpha, phi_eps);

    lhs = S_func(Cx, Cxt, Cy, Cyt, S1, S2, Vx, Vy, beta, Ax, Ay);
    rhs = Ax .* Px .* Ix + Ay .* Py .* Iy;
    rhs = reshape(rhs, [m*n, 1]);
    
    tol = 0; max_iter = 100; display_id = 10;
    Sinit = reshape(S, [m*n, 1]);
    S = cg(lhs, rhs, Sinit, tol, max_iter, show, num2str(i), display_id);
    S = reshape(S, [m, n]);
    
    % TODO:
    S(Gs<0.005) = 0;

    ls1(i) = sum(sum((S - Px .* Ix).^2 .* Ax))...
               + sum(sum((S - Py .* Iy).^2 .* Ay));
           
    L = L_func(Cx, Cxt, Cy, Cyt, S1, S2, Vx, Vy);
    St = reshape(S, [m*n, 1]);
    ls2(i) = beta * sum(sum(St' * L * St));
    

    Ax = phi(S - Px .* Ix, phi_alpha, phi_eps);
    Ay = phi(S - Py .* Iy, phi_alpha, phi_eps);

    B = phi(I - I0, phi_alpha, phi_eps);

    lhs = I_func(Cx, Cxt, Cy, Cyt, Px, Py, Ax, Ay, lambda, B);
    rhs = cal(Cxt, Px .* Ax .* S) + cal(Cyt, Py .* Ay .* S) + lambda * B .* I0;
    rhs = reshape(rhs, [m*n, 1]);
    
    tol = 0; max_iter = 100; display_id = 11;
    Iinit = reshape(I, [m*n, 1]);
    I = cg(lhs, rhs, Iinit, tol, max_iter, show, num2str(i), display_id);
    I = reshape(I, [m, n]);
    
    I(I<0) = 0;
    I(I>1) = 1;
    
    li1(i) = sum(sum((S - Px .* Ix).^2 .* Ax))...
               + sum(sum((S - Py .* Iy).^2 .* Ay));
    li2(i) = lambda * sum(sum((I-I0).^2 .* B));
    
    
    if show
        figure(1); clf;
        plot(1:iternum, ls1, '--o', 'DisplayName', 'ls1');
        hold on;
        plot(1:iternum, ls2, '--*', 'DisplayName', 'ls2');
        hold on;
        plot(1:iternum, li1, 'r', 'DisplayName', 'li1');
        hold on;
        plot(1:iternum, li2, 'b', 'DisplayName', 'li2');
        legend('show');
        
        figure(2); clf;
        Sa = abs(S);
        Sn = (Sa - min(Sa(:))) ./ (max(Sa(:)) - min(Sa(:)));
        subplot(211); imshow(Sn, []); title(num2str(i));
        subplot(212); imshow(I); title(num2str(i));
        
        drawnow;
    end
    
end

end


function y = S_func(Cx, Cxt, Cy, Cyt, S1, S2, Vx, Vy, beta, Ax, Ay)

num = numel(S1);

y = sparse(1:num, 1:num, Ax(:), num, num) + ...
    sparse(1:num, 1:num, Ay(:), num, num) + ...
    beta * L_func(Cx, Cxt, Cy, Cyt, S1, S2, Vx, Vy);

end

function y = L_func(Cx, Cxt, Cy, Cyt, S1, S2, Vx, Vy)

num = numel(S1);
t = S1 .* Vx.^2 + S2 .* Vy.^2;
t = sparse(1:num, 1:num, t(:), num, num);
first = Cxt * t * Cx;

t = S2 .* Vx.^2 + S1 .* Vy.^2;
t = sparse(1:num, 1:num, t(:), num, num);
second = Cyt * t * Cy;

t = (S1 - S2) .* Vx .* Vy;
t = sparse(1:num, 1:num, t(:), num, num);
third = 2 * Cyt * t * Cx;

y = first + second + third;
end


function y = I_func(Cx, Cxt, Cy, Cyt, Px, Py, Ax, Ay, lambda, B)

num = numel(Px);
t = Px.^2 .* Ax;
t = sparse(1:num, 1:num, t(:), num, num);
first = Cxt * t * Cx;

t = Py.^2 .* Ay;
t = sparse(1:num, 1:num, t(:), num, num);
second = Cyt * t * Cy;

t = lambda * B;
third = sparse(1:num, 1:num, t(:), num, num);
y = first + second + third;

end


function y = phi(x, alpha, eps)
y = 1 ./ (abs(x).^(2-alpha) + eps);
end


function y = msign(x)
y = sign(x);
y(y==0) = 1;
end

function Tc = get_Cx(m, n)

dx = [1, -1];
T = convmtx2(dx, m, n);
Tc = T(m+1:end, :);
Tc(end-m:end, :) = 0;

end

function Tc = get_Cxt(m, n)

dx = [1, -1];
T = convmtx2(dx, m, n);
Tc = T(m+1:end, :);
Tc = Tc';
Tc(1:m, :) = 0;

end

function Tc = get_Cy(m, n)

dy = [1; -1];
T = convmtx2(dy, m, n);
Tc = T;
Tc(1:m+1:end, :) = [];
Tc(m:m:end, :) = 0;

end

function Tc = get_Cyt(m, n)

dy = [1; -1];
T = convmtx2(dy, m, n);
Tc = T;
Tc(1:m+1:end, :) = [];
Tc = Tc';
Tc(1:m:end, :) = 0;

end

function y = cal(Cm, x)
y = reshape(Cm*x(:), size(x));
end

