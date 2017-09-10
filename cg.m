function x = cg(A, b, x, tol, max_iter, show, tstr, display_id)

if isa(A,'numeric')
  explicitA = true;
elseif isa(A,'function_handle')
  explicitA = false;
else
  error('cg:Atype','%s','A must be numeric or a function handle');
end

if explicitA
    r = b - A * x;
else
    r = b - A(x);
end

p = r;

loss = zeros(1, max_iter);

for k = 1:max_iter

    if explicitA
        Ap = A * p;
    else
        Ap = A(p);
    end
    
    alpha = sum(sum(r .* r)) / sum(sum(p .* Ap));

    x = x + alpha * p;
    r_new = r - alpha * Ap;
    
    beta = sum(sum(r_new .* r_new)) / sum(sum(r .* r));

    p = r_new + beta * p;
    r = r_new;
    
    loss(k) = sum(sum(r .* r));
    if loss(k) < tol
        break;
    end
    
    if show && k == max_iter
        figure(display_id); clf;
        plot(1:max_iter, loss, '-', 'DisplayName', '||Ax-b||^2');
        title(tstr);
        drawnow;
    end

end
    
end