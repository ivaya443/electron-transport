function [y_i] = slice_size(x_i,Nx,Ny)
%SLICE_SIZE Summary of this function goes here
%   Detailed explanation goes here

y_i = real(Ny + 2.*round(Nx - sqrt(Nx^2 - x_i.^2)));

end

