"Historical finite-order running expression, ported from alphas.m (not a modern threshold-matched solver)."
function alpha_s(s::Real)
    nf = s < (2*1.85)^2 ? 3 : s < (2*5.6)^2 ? 4 : 5
    a0 = 0.327/pi
    zeta3 = 1.202056903
    b0=(11-2nf/3)/4
    b1=(102-38nf/3)/16/b0
    b2=(2857/2-5033nf/18+325nf^2/54)/64/b0
    b3=(149753/6+3564zeta3-(1078361/162+6508zeta3/27)*nf+
         (50065/162+6472zeta3/81)*nf^2+1093nf^3/729)/256/b0
    K=1+b0*a0*log(s/1.7769^2)
    K>0 || throw(DomainError(s,"outside historical perturbative running domain"))
    L=log(K)
    a=a0/(K+b1*a0*L+b1^2*a0^2*(1-K+L)/K)
    c2=b2*(1-K)
    c3=b3/2*(1-K^2)+b1*b2*K*(K-1-L)+b1^3/2*(L^2-(1-K)^2)
    return pi*a*(1+c2*a^2+c3*a^3)
end

function r_pqcd(s::Real)
    nf,R,eta = s < (2*1.85)^2 ? (3,2.0,0.0) : s < (2*5.62)^2 ? (4,10/3,2/15) : (5,11/3,1/33)
    a=alpha_s(s)/pi
    c2=1.9857-0.1152nf
    c3=-6.63694-1.20013nf-0.00518nf^2-1.24eta
    c4=-156.61+18.775nf-0.7974nf^2+0.0215nf^3+(17.828-0.575nf)*eta
    return R*(1+a+c2*a^2+c3*a^3+c4*a^4)
end
