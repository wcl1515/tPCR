 function display(ob)
%function display(ob)
% "display" method for this class

name = inputname(1);
printm('"%s" is an object of class "%s":', name, class(ob))
ob = struct(ob);
disp(ob)

fnames = fieldnames(ob);
for ii=1:length(fnames)
	fname = fnames{ii};
	if isstruct(ob.(fname))
		printf('%s.%s :', inputname(1), fname)
		if streq(fname, 'meth') && ~isempty(ob.docs)
			mnames = fieldnames(ob.meth);

			nmax = 1; % find max length method name
			for im=1:length(mnames)
				nmax = max(nmax, length(mnames{im}));
			end
			nmax = min(nmax, 40);

			hmax = 1; % find max length handle name
			for im=1:length(mnames)
				hmax = max(hmax, ...
				length( func2str(ob.meth.(mnames{im})) ));
			end
			hmax = min(hmax, 40);

			format = ['\t%' sprintf('%d', nmax) 's: ' ...
				 '@%' sprintf('%d', hmax) 's %s'];
			for im=1:length(mnames)
				printf(format, mnames{im}, ...
					func2str(ob.meth.(mnames{im})), ...
					ob.docs{im})
			end
		else
			disp(ob.(fname))
		end
	end
end
