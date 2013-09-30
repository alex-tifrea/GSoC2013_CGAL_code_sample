%Tifrea Alexandru 311CA
function [x y z] = readFile(filename, D)
	fid = fopen(filename);
	aux = fscanf(fid, '%f');
	n = length(aux) / D;
	x = zeros(1, n);
	y = zeros(1, n);
	z = zeros(1, n);

	x(:) = aux(1:D:D*n);
	y(:) = aux(2:D:D*n);
	z(:) = aux(3:D:D*n);

	x_trig = [x(1:D+1),x(1)];
	y_trig = [y(1:D+1),y(1)];
	z_trig = [z(1:D+1),z(1)];

	x_edges = x(D+2:end);
	y_edges = y(D+2:end);
	z_edges = z(D+2:end);

%	clf;
%	plot(x_trig, y_trig);
%	hold on;
%	scatter(x_edges, y_edges, 'r');
endfunction
