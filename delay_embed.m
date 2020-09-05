function [X_out]=delay_embed(X,E,tau)
% delay_embed embeds a scalar time series X into an E-dimensional time series

X_out = zeros((length(X)-(E-1)*tau),E);
for t=1:(length(X)-(E-1)*tau)
    X_out(t,:) = X((E-1)*tau+t:-tau:(E-1)*tau+t-(E-1)*tau);
end
