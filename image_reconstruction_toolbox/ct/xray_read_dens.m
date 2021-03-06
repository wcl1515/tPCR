 function [density, mtype] = xray_read_dens(mtype, varargin)
%function [density, mtype] = xray_read_dens(mtype, [options])
%
% Read densities for  given material type(s).
% 
% in
%	mtype			'aluminum', 'copper', ...
%				(optionally can be a cell array of several)
% option
%	'units'		cm | mm	default: cm
% out
%	density			g/cm^3, scalar, or vector if mtype is cell array
%
% Copyright 2004-05-1, Jeff Fessler, The University of Michigan

% default is to show example
if nargin < 1, help(mfilename), error(mfilename), end
if streq(mtype, 'test'), xray_read_dens_test, return, end

arg.units = 'cm';
arg = vararg_pair(arg, varargin);

if iscell(mtype)
	for ll=1:length(mtype)
		density(ll) = xray_read_dens(mtype{ll}, 'units', arg.units);
	end
return
end

% densities, all in g/cm^3
% http://www.flukebiomedical.com/RMS/productManuals/76-430-1.pdf

switch mtype
case 'aluminum'
	density = 2.7;
case 'cesium'
	density = 1.9; % 1.879 ?
case 'copper'
	density = 0.92;
case 'csi'
	density = 4.51;
case 'gadolinium'
	density = 7.9;
case 'iodine'
	density = 4.93;
case 'lead'
	density = 11.34;
case 'water'
	density = 1;
case 'polyethylene'
	density = 0.950;
case {'lexan', 'polycarbonate'} % C16 H14 O3
	density = 1.190;
case 'polystyrene'
	density = 1.050;
case 'acrylic'
	density = 1.180;
case 'teflon'
	density = 2.214;
otherwise
	error('unknown material density "%s"', mtype)
end

switch arg.units
case 'mm'
	density = density / 1000;
case 'cm'
	;
otherwise
	error('bad units "%s"', arg.units)
end


%
% example usage
%
function xray_read_dens_test
mtype = {'aluminum', 'copper', 'lead', 'lexan'}
density = xray_read_dens(mtype)
