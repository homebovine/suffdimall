function [rcorr, pdist, xcorr, init, mave, sr, best] = Example44
datasrd = load('datasrdrealfull.tab');
times = 1

n = 278; %(size(datasrd, 1)+ 1)/2;
p = size(datasrd, 2) - 1;


d = 1;


for ii = 1:times
    ii
    [ init(:,:,ii),mave(:,:,ii)] = SIRMIM(n,p, ii, datasrd, d);
end



      function [rcorr, pdist, xcorr, init, mave, sr, best, varbest] = SIRMIM(n,p, ii, datasrd,  d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           THE SIR METHOD FOR MULTIPLE INDICES MODEL                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            DATA GENERATION                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = d
%  [datasrd0, ind] = datasample(datasrd, n, 'Replace', false);
ind = 1:278;
y0 = datasrd(ind, 1);
x0 = datasrd(ind, 2: (p + 1));
datasrd1 =datasrd(setdiff(1:size(datasrd, 1),ind), :);
newy = datasrd1(:, 1);
newx = datasrd1(:, 2: (p+ 1));
n = n;
x = x0;
y = y0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    SLICED INVERSE REGRESSION                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss = cov(x,0); mu = mean(x);
z = (x - ones(n,1)*mu)*(inv(ss))^(1/2);
invfest = kernelregression(y,z,y, 'kernel','epanechinikov');
bLam = invfest'*invfest/n;
[U,S,V] = svd(bLam); 
    
mave = dmave(z,y,dim);
for j=1:dim
    i = 1;
    while (mave(i,j)==0)
        i = i+1;
    end
    if (mave(i,j)<0)
        mave(:,j) = -mave(:,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       FMINCON SEARCH ALGORITHM                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           BACK-TRANSFORMATION                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mave = mave * inv(mave(1:dim, 1:dim))
allOneString = sprintf('%d,' ,  mave)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        RESULTS SUMMARY                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% The objective function for fmincon
function ee = objcon(best,invfest,z,p,n)
zfest = kernelregression(z*best,z,z*best, 'kernel','epanechinikov');
zres = z - zfest;
zfest = kernelregression(z*best,invfest,z*best, 'kernel','epanechinikov');
resinv = invfest - zfest;
eematrix = resinv'*zres/n;
ee = (eematrix(:))'*(eematrix(:));


function ee = gobjcon(best,invfest,z,p,n)
zfest = kernelregression(z*best,z,z*best, 'kernel','epanechinikov');
zres = z - zfest;
zfest = kernelregression(z*best,invfest,z*best, 'kernel','epanechinikov');
resinv = invfest - zfest;
innersum = zeros(p * p, p*p, n);
for i = 1:n
inner = (resinv(i, :)' * zres(i, :));
innersum(:, :, i) =  inner(:) * inner(:)'; 
end
						  ee = mean(innersum, 3);
				  


	 function ee = hobjcon(best2,invfest,z,p,d, n)
n = size(z,1);
best2 = reshape(best2, (p-2), d);
best = [eye(d); best2];
zfest = kernelregression(z*best,z,z*best, 'kernel','epanechinikov');
zres = z - zfest;
zfest = kernelregression(z*best,invfest,z*best, 'kernel','epanechinikov');
resinv = invfest - zfest;
						    innersum = zeros(p, p, n);
for i = 1:n
						    innersum(:, :, i) = (resinv(i, :)' * zres(i, :)); 
end
						  ee = mean(innersum, 3);
ee = ee(:);

function J=numjacobian(x, invfest, z, p, d, n)
n = size(z, 1)

	     nn=length(x); % n=2 in your case, but can be higher 
	     fx=hobjcon(x, invfest, z, p, d, n); % evaluate the function at point x
	     step=1e-6; % difference step you choose, can be 1e-10 if you like
size(fx)
	     for i=1:nn
	     xstep = x;
   % compute the i-th component partial derivative 
   % numerically using first order forward difference approximation
	     xstep(i)=x(i)+step;
	     J(:,i)=(hobjcon(xstep, invfest, z, p, d, n)-fx)/step;
	     end




% the constraint function for fmincon
function [c,ceq]=constraint(best)
[p,d] = size(best);
ceq = best'*best-eye(d);
c = zeros(d,1);
for j = 1:d
    i = 1;
    while (sign(best(i,j))==0)
        i = i+1;
    end
    c(j) = -sign(best(i));
end


% the constraint function for fmincon
function [c,ceq]=constraint1(best)
[p,d] = size(best);
beq = zeros(p, d);
beq(1:d, 1:d) = eye(d);
Aeq = zeros(p, p);
Aeq(1:d, 1:d) = eye(d);
ceq = Aeq * best - beq;
c = zeros(d,1);

    
% the kernel regression
function [f,u] = kernelregression(x,y,x0,varargin)

%   f = kernelregression(X,Y) computes the estimate of regression function
%   evaluated at x0. The default estimate is based on a normal kernel
%   function using a window parameter (bandwidth) that is a function of the
%   number of points in X.

%   [f,u] = kernelregression(...) also returns the bandwidth of the kernel
%   smoothing window.

%
%   [f,u] = kernelregression(...,'PARAM1',val1,'PARAM2',val2,...) specifies
%   parameter name/value pairs to control the regression estimation.  Valid
%   parameters are the following:
%
%      Parameter    Value
%      'kernel'     The type of kernel smoother to use, chosen from among
%                   'normal' (default), 'box', 'triangle', and
%                   'epanechinikov'.
%      'width'      The bandwidth of the kernel smoothing window.  The default
%                   is optimal for estimating normal densities, but you
%                   may want to choose a smaller value to reveal features
%                   such as multiple modes.
%
%
%   Example:
%      x = randn(30,2); y = exp(x(:,1));
%      f = kernelregression(x,y,x, 'kernel','epanechinikov','width',1);
%      plot(x(:,1),f,'go');


   
[n, p] = size(x);
xmin = min(x);
xmax = max(x);



% Process additional name/value pair arguments
okargs = {'width'  'kernel'};
defaults = {[] 'normal'};
[emsg,u,m,kernel] = statgetargs(okargs, defaults, varargin{:});
error(emsg);

% Default window parameter is optimal for normal distribution
if (isempty(u)),
   med = median(x);
   sig = median(abs(x-repmat(med,n,1))) / 0.6745;
   if sig>0
      u = sig * (4/(3*n))^(1/6);
   else
      u = ones(1,p);
   end
end


xsize = size(x0);
m = size(x0,1);

okkernels = {'normal' 'epanechinikov' 'box' 'triangle'};
if isempty(kernel)
   kernel = okkernels{1};
elseif ~(isa(kernel,'function_handle') | isa(kernel,'inline'))
   if ~ischar(kernel)
      error('Smoothing kernel must be a function.');
   end
   knum = strmatch(lower(kernel), okkernels);
   if (length(knum) == 1)
      kernel = okkernels{knum};
   end
end



blocksize = 1e8;
if n*m<=blocksize
   for ii = 1:p
       z(:,:,ii) = (repmat(x0(:,ii)',n,1)-repmat(x(:,ii),1,m))/u(ii);
       fff(:,:,ii) = feval(kernel, z(:,:,ii))/u(ii);
   end
   ff = prod(fff,3);
   f1 = y'*ff;
   f2 =  sum(ff,1);
   f = f1./(ones(size(y,2),1)*f2+1E-6);
else   
    display('error: the dimensionality of data set is too large')
end

f = f';

% -----------------------------
% The following are functions that define smoothing kernels k(z).
% Each function takes a single input Z and returns the value of
% the smoothing kernel.  These sample kernels are designed to
% produce outputs that are somewhat comparable (differences due
% to shape rather than scale), so they are all probability
% density functions with unit variance.
%
% The density estimate has the form
%    f(x;k,h) = mean over i=1:n of k((x-y(i))/h) / h

function f = normal(z)
%NORMAL Normal density kernel.
%f = normpdf(z);
f = exp(-0.5 * z .^2) ./ sqrt(2*pi);

function f = epanechinikov(z)
%EPANECHINIKOV Epanechinikov's asymptotically optimal kernel.
a = sqrt(5);
z = max(-a, min(z,a));
f = .75 * (1 - .2*z.^2) / a;

function f = box(z)
%BOX    Box-shaped kernel
a = sqrt(3);
f = (abs(z)<=a) ./ (2 * a);

function f = triangle(z)
%TRIANGLE Triangular kernel.
a = sqrt(6);
z = abs(z);
f = (z<=a) .* (1 - z/a) / a; 




function [eid,emsg,varargout]=statgetargs(pnames,dflts,varargin)
emsg = '';
eid = '';
nparams = length(pnames);
varargout = dflts;
unrecog = {};
nargs = length(varargin);

% Must have name/value pairs
if mod(nargs,2)~=0
    eid = 'WrongNumberArgs';
    emsg = 'Wrong number of arguments.';
else
    % Process name/value pairs
    for j=1:2:nargs
        pname = varargin{j};
        if ~ischar(pname)
            eid = 'BadParamName';
            emsg = 'Parameter name must be text.';
            break;
        end
        i = strmatch(lower(pname),pnames);
        if isempty(i)
            % if they've asked to get back unrecognized names/values, add this
            % one to the list
            if nargout > nparams+2
                unrecog((end+1):(end+2)) = {varargin{j} varargin{j+1}};
                
                % otherwise, it's an error
            else
                eid = 'BadParamName';
                emsg = sprintf('Invalid parameter name:  %s.',pname);
                break;
            end
        elseif length(i)>1
            eid = 'BadParamName';
            emsg = sprintf('Ambiguous parameter name:  %s.',pname);
            break;
        else
            varargout{i} = varargin{j+1};
        end
    end
end

varargout{nparams+1} = unrecog;


function [B, cv] = rMAVE(xx, y, h, nd)
%Searching the effective dimension reduction subspace of model
%  			y = g(B^TX) + e
%Useage: directions(x, y, h, d)
%Input: 
%     x --- expanaltory variables
%     y --- response
%     h --- bandwidth
%     d --- working dimension of the space
%Output:
%     Directions B
%Reference: Y. Xia, H. Tong, W.K. Li and L.X. Zhu, 
% "adaptive estimation of the effictive dimension space",(2001)  
%-------------------------------------------------------------
% Example
%	n = 100;
%	x = randn(n,4);
%	beta1 = [1 2 3 0]';
%	beta1 = beta1/sqrt(beta1'*beta1);
%	beta2 = [-2 1 0 1]';
%	beta2 = beta2/sqrt(beta2'*beta2);
%	y = (x*beta1).^2 + x*beta2 + 0.2*randn(n,1);
%	B = rMAVE(x, y, 0.5, 2)
%  % estimation errors
%	B0 = [beta1 beta2];
%	error = B'*(eye(4)-B0*B0')*B

[n,p] = size(xx);
[ny, py] = size(y); 
mm = mean(xx);
xx = xx - repmat(mm,n,1);
ss = inv(xx'*xx/n)^0.5;
xx = xx*ss;


if (ny ~= n)
   disp('Error: matrix of x and y dont match')
   return
end
if (p <= 1)
   disp('Error: please check your matrix x')
   return
end
if (py > 1)
   disp('Error: please check your matrix y')
   return
end
if (sum(abs(mean(xx,1))) > 0) | (abs(mean(std(xx,1))-1)>0) 
%   disp('Warning: covariates must be standardized')
end;   
b = (1:p)';             
b = b/sqrt((b')*b);
c = b;
Jzf1 = ones(n, 1);
txt = Jzf1;
J1 = ones(n,1);
zfdd = zeros(1,p);
zfa0 = zeros(n, p*p);
zfa1 = zeros(n, p);
zfb0 = zeros(n, p);
zfb1 = ones(n, 1);
zfb2 = ones(n, 1);
nref = 10;
niter = 10;
zfh0 = mean(std(xx))/n^(1/(p+4));	
for zfia = 1:n;
   	zfxy = xx - repmat(xx(zfia,:), n,1);
      zfx = zfxy.*zfxy*ones(p,1);
      zfK = exp(- zfx /(2*zfh0*zfh0*p) );
      K = zfK/sum(zfK);
      Kv = repmat(K,1,p).*zfxy;
      zfa0(zfia,:) = reshape((Kv')*zfxy, 1, p*p);
      zfa1(zfia,:) = ones(1,n,1)*Kv;
      zfb0(zfia,:) = (y')*Kv;
      zfb1(zfia) = (K')*y;
   end;
   
	zfbb0 = [];
	zfcc0 = zeros(p, nd);
	for ik = 1 : min(nd+1,p);
      if ik > 1; 
         b = (eye(p) - zfbb0*(zfbb0'))*ones(p,1);
   	end;
   	for iter = 1 : niter;
        b0 = b;
        zfbb = [zfbb0 b];
        zfAs = 0.0;
        zfBs = 0.0;
     	for zfia = 1 : n;
        zfaa = reshape(zfa0(zfia,:), p, p);
        zfat12 = zfa1(zfia,:)*zfbb ;
        zfat22 = (zfbb')*zfaa*zfbb;
        zfat1 = [1
           (zfat12')];
        zfat2 = [zfa1(zfia,:)*zfbb 
           (zfbb')*zfaa*zfbb ];
        zfat = [zfat1 zfat2]  + (1.0E-5)*eye(1+ik);
        zfbt = [zfb1(zfia)   zfb0(zfia,:) * zfbb ]';
        zfbt = inv(zfat)*zfbt;
        a = zfbt(1);

        zfbe1 = 0;
        if ik >1 ;
        		zfbe1 = zfbt(2:ik);
        end;
        d = zfbt(1+ik);

		  if ik == 1 
        		zfbb0 = zeros(p,1);
        end
        zfBs = zfBs + d*( (zfb0(zfia,:)- zfa1(zfia,:)*a)' - zfaa*zfbb0*zfbe1  );
        zfAs = zfAs + d*d*zfaa;
        if ik == 1 
           zfbb0=[];
        end
     	end;

      zfAs = zfAs + 1.0e-10*eye(p);
      b = inv(zfAs)*zfBs;
      if ik > 1 
      		b = (eye(p) - zfbb0*(zfbb0'))*b;
      end;
      b = b/sqrt((b')*b);

      if  ( (b0')*(eye(p) - b*(b'))*b0 < 1.0E-5 )
                iter = niter + 1;
      end;
   	end;
		if ik == 1 
         zfbb0 = b;
      end
      if ik > 1 
         zfbb0 = [zfbb0 b];
      end
   end;
   error1 = 0.0;
   error2 = 0.0;
   for kk = 1 : nref;
        	zfbbi = zfbb0(:,1:nd);
        	xxk = xx*zfbb0;
        	[n1,nd1] = size(xxk);  
 		  	zfh = h;
        	for zfia = 1 : n;
         	zfxyk = xxk - repmat(xxk(zfia,:), n,1);
         	zfxy = xx - repmat(xx(zfia,:), n, 1);
         	zfx = zfxyk.*zfxyk*ones(nd1,1);
         	zfK = exp(- zfx /(2*zfh*zfh) );
            K = zfK/sum(zfK);
         	Kv = repmat(K,1,p).*zfxy;
         	zfa0(zfia,:) = reshape((Kv')*zfxy, 1, p*p);
         	zfa1(zfia,:) = ones(1,n)*Kv;
         	zfb0(zfia,:) = (y')*Kv;
         	zfb1(zfia) = (K')*y;
        	end;
    	  	zfbb0 = zfbb0(:,1:nd);
    	  	for ik = 1 : nd;
    	  		b = zfbb0(:,ik);
    			for iter = 1 : niter;
        			b0 = b;
        			zfAs = 0.0;
        			zfBs = 0.0;
        			for zfia = 1 : n;
        				zfaa = reshape(zfa0(zfia,:), p, p);
        				zfat12 = zfa1(zfia,:)*zfbb0 ;
        				zfat22 = (zfbb0')*zfaa*zfbb0;
        
        				zfat1 = [1 zfat12]';
        				zfat2 = [zfa1(zfia,:)*zfbb0 
           				(zfbb0')*zfaa*zfbb0] ;
        				zfat = [zfat1 zfat2]  + (1.0E-5)*eye(1+nd);
        				zfbt = [zfb1(zfia)   zfb0(zfia,:)*zfbb0]';

        				zfbt = inv(zfat)*zfbt;
        				a = zfbt(1);
        				zfbe1 = zfbt(2:(nd+1));
        				d = zfbe1(ik);
                        DD(zfia, ik) = d;
        				zfbe1(ik) = 0.0;
        				zfBs = zfBs + d*( (zfb0(zfia,:)- zfa1(zfia,:)*a)' - zfaa*zfbb0*zfbe1  );
        				zfAs = zfAs + d*d*zfaa;
        			end;
        			zfAs = zfAs + 1.0e-10*eye(p);
        			b = inv(zfAs)*zfBs;
        			zfbb0(:,ik) = zeros(p,1);
        			b = (eye(p) - zfbb0*(zfbb0'))*b;
        			b = b/sqrt((b')*b);
        			zfbb0(:,ik) = b;
        		end;
       end;
       error1 =  sum( diag( zfbb0'*(eye(p)-zfbbi*zfbbi')*zfbb0 ) );   
       if error2 + error1 < 0.0001 
          break;
       end;
       error2 = error1;
   end;
   B = zfbb0;

xB = xx*B;
ye = y;
for i = 1:n;
    xi = xB - repmat(xB(i,:), n,1);
    kx = exp(-sum(xi.^2,2)/(2*h*h));
    kx(i) = 0;
    xi1 = [ones(n,1) xi];
    kxi1 = xi1.*repmat(kx, 1, nd+1);
    beta = inv( kxi1'*xi1 )*kxi1'*y;
    ye(i) = beta(1);
end
cv = (y-ye)'*(y-ye)/n;

DD = DD-repmat(mean(DD), n, 1);
[v, d] = eig(DD'*DD);
d = diag(d);
[d, I] = sort(d);
v = v(:,I);
B = B*v;

B = ss*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end

function B = dmave(x, y, m0);

% Input : x (n x p); y (n x 1); m0 (integer)
% Output: B (p x m0) 


[n,p] = size(x);
onen = ones(n,1);

mm = mean(x);
x = x - repmat(mm,n,1);
ss = inv(x'*x/n)^0.5;
x = x*ss;

m = p;
B = eye(p);
Ba = B;
BI = Ba;
B = B(:,1:m);
noip0 = 0;
noip1 = 1;
iterstop = 0;
Btmp = B;
rige = std(y)*mean(std(x,[],1));
y = (y-mean(y))/std(y);
x = x./repmat(std(x,[],1),n,1);

yc = y;
x0 = x;
xc = x;

if 1 == 1
[a, yt] = sort(y);
[a, yt] = sort(yt);
ycc = yt/max(yt);
for i = 1:size(x,2);
    a = x(:,i);   
    [b, a] = sort(a);
    [b, a] = sort(a);
    xc(:,i) = a/std(a);
end
end

niter = floor(p*3/2);
%ch = (sqrt(p)/n^(1/(p+4))*n^(2/(m0+4)))^(1/niter);
for iter = 1:niter;
    x = xc;
    if iter >= p;
        x = x0;
    end
  	adj = p^(1/iter)*n^(1/(m0+4)-1/(m+4));
  	hy = std(yc)/n^(1/5)*p^(1/(iter+1)); %*adj;
	ky = repmat(yc, 1, n);
    ky = [ycc/sqrt(2*pi)/hy*n^(1/iter) exp(-(ky-ky').^2/(2*hy*hy))/sqrt(2*pi)/hy];
    n1 = size(ky,2);
    h = p^(1/(iter+1))*mean(std(x*Ba,[],1))/n^(1/(m+4));
    h2 = 2*h*h*adj;
   
   ABI = zeros(m, n*n); 
   for iiter = 1:max(1, (m < p)*floor(m/2));
	dd = zeros(m*p, m*p);
   	dc = zeros(m*p,1);
   	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        sxij = sum((xij*Ba).^2,2) + sum(xij.^2,2)/1.5^iter;
      	ker = exp(-sxij/h2);  
        rker = repmat(ker, 1, p+1);
   		onexi = [xij*B onen];
   		xk = (onexi.*rker(:, 1:m+1))';
   		abi = inv(xk*onexi+eye(m+1)/n)*xk*ky;
      
      	kxij = (xij.*rker(:,1:p))';
      	kxijy = kxij*(ky - repmat(abi(m+1,:),n,1));
      	ddx = kxij*xij;
      	for k1 = 1:m
         	ka = (k1-1)*p+1:k1*p;
         	dc(ka) = dc(ka) + kxijy*abi(k1,:)';
      	end
        tmp = abi(1:m, :)*abi(1:m, :)';
        dd = dd + kron(tmp, ddx);
      
        ABI(:, (j-1)*n1+1:j*n1) = abi(1:m,:); 
   	end
    
	   B = pinv(dd+ rige*eye(length(dc))/n)*dc;
       
   	   B0 = reshape(B, p, m)*ABI;
       [B, R] = eig(B0*B0');
       B = B(:,p-m+1:p);
       Ba = B;
       if (max(svd(B*B'-BI*BI')) < 0.001);
           break
       end
       BI = B;
       
   end

   mb = m;
   ma = max(m-1,m0);
   m = ma; 
   B = B(:,mb-m+1:mb);
   Ba = B;
   if (max(svd(B*B'-BI*BI')) < 0.001)*(iter>p+3);
       break
   end
   BI = Ba;   
   
end

cv = 0;
kye = ky;
	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        sxij = sum((xij*Ba).^2,2) + 0*sum(xij.^2,2)/1.5^iter;
        ker = exp(-sxij/h2);  
        ker(j) = 0;
        if mean(ker) > 1/n/n;
           kye(j,:) = ker'*ky/sum(ker);
        end      
  end   
cv = (ky-kye).^2;
cv = sum(sum(cv))/sum(sum(cv>0));

B = ss*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end


function [B, cv, Binitial] = SR(x, y, nslice, m0, initial, kappa);
% INPUT
% nslice: number of slices
%    m0 : specified efficient dimension
%    initial=1, if OPG is used for initial value
%            2, if OLS is used for initial value
%            0, others (a gradually reducing procedure)
%  
% OUTPUT
% B     : estimated dimension
% cv    : CV value
% Binitial: initial estimator of B.


[n,p] = size(x);
onen = ones(n,1);

mm = mean(x);
x = x - repmat(mm,n,1);
ss = inv(x'*x/n)^0.5;
x = x*ss;


m = p;
B = eye(p);
Ba = B;
BI = Ba;
B = B(:,1:m);
noip0 = 0;
noip1 = 1;
iterstop = 0;
Btmp = B;
rige = std(y)*mean(std(x,[],1));

yc = y;
x0 = x;
xc = x;

[yc, I] = sort(y);
y = y(I);
x = x(I,:);
x0 = x0(I,:);
xc = xc(I,:);
ky = [];
for i = 1:nslice
    ends = i*floor(n/nslice);
    if i == nslice;  
        ends = n;
    end
    I = (i-1)*floor(n/nslice)+1:ends;
    yi = zeros(n,1);
    yi(I) = y(I)*0+1;
    ky = [ky yi];
end


if initial == 0
    m = p;
end
if initial == 1
    h = 1/n^(1/(p+4))*kappa;
    h2 = 2*h*h;
   	BB0 = 0;
	onexi = ones(n,p+1);
   	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        dis = sum((xij).^2,2);
        disd = sort(dis);
        h2j = max(h2, disd(p+1));
        ker0 = (sum((xij).^2,2))/h2j;
      	ker = exp(-ker0);
   		onexi(:,1:p) = xij;
   		xk = (onexi.*repmat(ker, 1, p+1))';
   		abi = inv(xk*onexi+eye(p+1)/n)*xk*ky;        
        BB0 = BB0 + abi(1:p,:)*abi(1:p,:)';
    end   
    [B, R] = eig(BB0/n);
    Binitial = B(:,p-m0+1:p);  
    B0 = Binitial;
    m = m0;
    if p > 10
        B0 = B(:,p-9:p);
        m = 10;
    end
    Ba = B0;
    B = B0;
end

BB0 = 0;
if initial == 2
    m = m0;
    for i = 1:nslice
        ends = i*floor(n/nslice);
        if i == nslice;  
            ends = n;
        end
        I = (i-1)*floor(n/nslice)+1:ends;
        xi = [x(I,:) ones(length(I),1)];
        ncol=length(xi(1,:));D=1.0E-6*eye(ncol);    %ridge effect added by Hanson
        bi = inv(xi'*xi+D)*(xi'*y(I));              %ridge effect added by Hanson
        BB0 = BB0 + bi(1:p)*bi(1:p)';
    end   
    [B, R] = eig(BB0/n);
    B0 = Binitial;
    m = m0;
    if p > 10
        B0 = B(:,p-9:p);
        m = 10;
    end
    Ba = B0;
    B = B0;
end

cv = 0;
Binitial = [];
if (initial ==1)+(initial ==2)
B0 = ss*B0;
for i = 1:size(Binitial,2);
    Binitial(:,i) = Binitial(:,i)/sqrt(Binitial(:,i)'*Binitial(:,i));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

niter = floor(p*3/2);
for iter = 1:niter;
    x = xc;
    if iter >= p;
        x = x0;
    end
  	adj = p^(1/iter)*n^(1/(m0+4)-1/(m+4));
    n1 = size(ky,2);
    h = p^(1/(iter+1))*mean(std(x*Ba,[],1))/n^(1/(m+4));
    h2 = 2*h*h*adj;
   
   ABI = zeros(m, n*n); 
   for iiter = 1:max(1, (m < p)*floor(m/2));
	dd = zeros(m*p, m*p);
   	dc = zeros(m*p,1);
   	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        dis = sum((xij*Ba).^2,2);
        disd = sort(dis);
        h2j = max(h2, disd(p+1));
        sxij = dis + sum(xij.^2,2)/1.5^iter;
      	ker = exp(-sxij/h2j);  
        rker = repmat(ker, 1, p+1);
   		onexi = [xij*B onen];
   		xk = (onexi.*rker(:, 1:m+1))';
   		abi = inv(xk*onexi+eye(m+1)/n)*(xk*ky);
      
      	kxij = (xij.*rker(:,1:p))';
      	kxijy = kxij*(ky - repmat(abi(m+1,:),n,1));
      	ddx = kxij*xij;
      	for k1 = 1:m
         	ka = (k1-1)*p+1:k1*p;
         	dc(ka) = dc(ka) + kxijy*abi(k1,:)';
      	end
        tmp = abi(1:m, :)*abi(1:m, :)';
        dd = dd + kron(tmp, ddx);
      
        ABI(:, (j-1)*n1+1:j*n1) = abi(1:m,:); 
   	end
    
	   B = inv(dd+ rige*eye(length(dc))/n)*dc;
       
   	   B0 = reshape(B, p, m)*ABI;
       [B, R] = eig(B0*B0');
       B = B(:,p-m+1:p);
       Ba = B;
       if (max(svd(B*B'-BI*BI')) < 0.001);
           break
       end
       BI = B;
       
   end

   mb = m;
   ma = min(max(m-1,m0), 10);
   m = ma; 
   B = B(:,mb-m+1:mb);
   Ba = B;
   if (max(svd(B*B'-BI*BI')) < 0.001)*(iter>p+3);
       break
   end
   BI = Ba;   
%   [iter,max(svd(B*B'-BI*BI'))]
end

cv = 0;
kye = ky;
	for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        sxij = sum((xij*Ba).^2,2) + 0*sum(xij.^2,2)/1.5^iter;
        ker = exp(-sxij/h2);  
        ker(j) = 0;
        if mean(ker) > 1/n/n;
            kye(j,:) = ker'*ky/sum(ker);
        end      
  end   
cv = (ky-kye).^2;
cv = sum(sum(cv))/sum(sum(cv>0));


B = ss*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end


