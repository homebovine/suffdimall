for i = 1:4
		for itr = 1:100
		cuvx =  mx0(:, :, itr) * otherresult(itr).betahat(:, 1, i);
		cuvy = my0(:, itr);
		cuvnewx = mnewx(:, :, itr) * otherresult(itr).betahat(:, 1, i)  ;
		%fity = csaps(cuvx, cuvy,'xx',  cuvnewx)%fit([cuvx(:, 1)],  cuvy, 'lowess');
		[m, m0] = reg_local_linear_smoothing(cuvx, cuvy, 2, 0.3, cuvnewx); 
		%cross = median((m0 - newy).^2); 
		cross(itr, i) = mean((m0 - mnewy(:, itr)).^2);

		%cuvx =  x0 * otherresult.betahat(:, 2:3, 3);
		%cuvy = y0;
		%cuvnewx = newx * otherresult.betahat(:, 1:2, i);
		%fity = fit([cuvx(:, 1) cuvx(:, 2)],  cuvy, 'poly23');
		%fity = feval(fity, [cuvnewx(:, 1), cuvnewx(:, 2)]);
		%cross(itr, i) = median((fity - newy).^2);
end
end
%save(sprintf('bootresult.mat'), 'mx0', 'my0', 'mnewx', 'mnewy', 'otherresult'); 
