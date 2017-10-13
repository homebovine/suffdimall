datasrd = load('datasrdreal.tab');
nsim = 100
for itr = 1:nsim
n = size(datasrd, 1)/100;
p = size(datasrd, 2)/2 - 1;
y0 = datasrd(((itr - 1) * n + 1): (itr * n), 1);
x0 = datasrd(((itr - 1) * n + 1): (itr * n), 2:7);
newy = datasrd(((itr - 1) * n + 1): (itr * n), 8);
newx = datasrd(((itr - 1) * n + 1): (itr * n), 9:end);

itr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                         %
%                                                                         %
%                             DIRECTIONAL  REGRESSION                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x0; y = y0; dim = 4;
ss = cov(x,0); mu = mean(x);
z = (x - ones(n,1)*mu)*(inv(ss))^(1/2);
beta0 =[1,2.211599e-01,1.867076e+00,1.215989e+00,-1.273452e-01,-2.537102e-01];% [1,2.198969e-01,1.255190e+00,1.014425e+00,-2.727293e-01,-1.231710e-01];% [1,5.443134e-01,5.776732e-01,7.712757e-01,5.739989e-01,-5.448373e-01];
%beta0 = 0.25 * ones(p, p); %; %0.65 * ones(p, p); %
beta0 = repmat((beta0)', 1, p);%0.65 * ones(p, p); %;

for i=1:4
    if i==1
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_SIR(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    if i==2
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_SAVE(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    if i==3
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_DR(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    if i==4
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i), ...
        div1(:, :, i)]=select_main_PHD(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    i
end
div
otherresult = struct('LIC', LIC, 'ee', ee, 'eflag', eflag, 'dis', ...
                     dis, 'betahat', betahat, 'div1', div1, 'div', div);       
for i = 1:4
cuvx =  x0 * otherresult.betahat(:, 1, i) ;
cuvy = y0;
cuvnewx = newx * otherresult.betahat(:, 1, i)  ;
fity = csaps(cuvx, cuvy, 'xx',  cuvnewx);
cross(itr, i) = median((fity - newy).^2);
end
end

    
save(sprintf('subapplication_result.mat'), 'cross'); 



