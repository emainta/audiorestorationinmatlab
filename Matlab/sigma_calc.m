function [sigma] = sigma_calc(N_,p_,m_,a_,sig_)

c = zeros(p_ ,1);
c00 = zeros(1);

%% Definition of c(x)
i=0;
for k = (p_+1):N_
    for j = 1:p_
        c(j,i+1) = c(j,i+1) + sig_(k-i)*sig_(k-j);
    end
end


for k = (p_+1):N_
        c00 =  c00 + sig_(k)*sig_(k);
end

sigma = (1/(N_- p_- m_))*(c00+dot(a_(2:end),c));


end

