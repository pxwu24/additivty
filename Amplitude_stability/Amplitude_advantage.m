% generate r dimensional pure state from 2r real vectors. params is a 2r
% real vector
function f = pure_state_density(params)
    r = length(params);
    if (-1)^r == 1
        r = r/2;
    else
        disp('Not a valid even length real vector!')
    end
    % Extract real and imaginary parts
    realPart = params(1:r);
    imagPart = params(r + 1:2*r);
    % Construct the complex state vector
    psi = complex(realPart, imagPart);
    % Normalize the state vector
    f = psi / norm(psi);   
end

% compute the advantange of amplitude damping with d-dimensional ancillary
% system given input rho; rho is 2d by 2d density matrix.
function f = compute_advantage(gamma,d,rho)
    V = amp_damp(gamma,'isom');
    V_amp = kron(eye(d), V);
    sigma_VBE = V_amp*rho*V_amp';
    sigma_VB = TrX(sigma_VBE,2,[2*d,2]);
    sigma_VE = TrX(sigma_VBE,2,[d,2,2]);
    sigma_BE = TrX(sigma_VBE,1,[d,4]); 
    % f1 is S(B)- S(E); f2 is S(BV) - S(EV)
    f1 = real(VNent(TrX(sigma_BE,2,[2,2]))-VNent(TrX(sigma_BE,1,[2,2])));
    f2 = real(VNent(sigma_VB)-VNent(sigma_VE));
    f = f1 - f2;
end
