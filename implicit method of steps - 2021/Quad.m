function I_Riemann = Quad(y,DOM)
% Quadrature of data y on DOM. Assumes equal grid of the form
% [DOM(1),(DOM(2)-DOM(1))/N,2(DOM(2)-DOM(1))/N,...,DOM(2)]. Simple code,
% included only for readability of other code.
if isintval(y)
    I_Riemann = (DOM(2)-DOM(1))/intval(size(y,1))*sum(y,1);
else
    I_Riemann = (DOM(2)-DOM(1))/size(y,1)*sum(y,1);
end
end