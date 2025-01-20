function h = filtre_canal(m,A,tau,Ts,L)
h = 0;
for ii=1:length(tau)
    h = h + A(ii)*nyquist(m*Ts-tau(ii)*Ts,Ts,0.5);
end
end