% Gtomo2_moj_test
% Test the Mojette case of the Gtomo2_table object
% Copyright 2005-12-14, Jeff Fessler, The University of Michigan

list = {'Gtomo2_moj_test1', 'Gtomo2_moj_test3'};
if 3 == exist('Gtomo2_dd')
	list{end} = 'Gtomo2_moj_test2';
end

run_mfile_local(list);

if 0 % put it through the paces...
	Gtab2 = Gtab;
	Gtab2.arg.nthread = 2;
	Gtomo2_test(Gtab, mask, Gtab2)
return
end
