function [LICa,LICb,ee,eflag,dis,betahat,div, div1]=select_main_SAVE(x,y,dim,beta,ti, pini, betapre)
    warning off;
    [n,p] = size(x);
    z  = x;
    zz = [];
    for ii = 1:p
        zz = [zz, repmat(z(:,ii),1,p).*z];
    end
    for ii = 1:p
        invfest(:,ii) = kernelregression(y,z(:,ii),y, 'kernel',...
            'epanechinikov');
    end
    for ii = 1:p
        for jj = 1:ii
            invsqfest = kernelregression(y,z(:,ii).*z(:,jj),y, 'kernel',...
                'epanechinikov');
            blam(ii,jj,:) = invsqfest - invfest(:,ii).*invfest(:,jj);
            blam(jj,ii,:) = blam(ii,jj,:);
        end
    end
    for ii = 1:n
        blam11(ii,:) = ((eye(p) - blam(:,:,ii))*(x(ii,:) - invfest(ii,:))')';
        blam22(:,:,ii) = eye(p) - blam(:,:,ii);
    end 
    eflag  = 6*ones(1,dim);
    ii=1;
    betahat = [];    
iiend = 0;
    for ii = 1:dim
        Large = 4000; Small = 10^(-8);
        options = optimset('display','off','MaxIter',Large,'MaxFunEvals',Large,...
            'TolX',Small,'TolCon',Small,'TolFun',Small, 'algorithm','active-set');%,,'LargeScale','off',
        %betainit = beta(ii+1:p,1:ii);
iiinit = iiend + 1;
iiend = iiinit + ii - 1;
betainit = betapre(ii + 1:p, iiinit:iiend )
if pini ~= 1
        [best,ee(1,ii),eflag(1,ii)] = fminsearch(@(best)objcon(best,blam11,blam22,z,p,zz),...
                betainit,options);
        [best,ee(1,ii),eflag(1,ii),output,grad,hessian] = fminunc(@(best)objcon(best,blam11,blam22,z,p,zz),...
                best * 1.01,options);
end
if pini == 1
        [best,ee(1,ii),eflag(1,ii),output,grad,hessian] = fminunc(@(best)objcon(best,blam11,blam22,z,p,zz),...
                betainit,options);
end

%output
%         [best,ee(1,ii),eflag(1,ii)] = fminunc(@(best)objcon(best,blam11,blam22,z,p,zz),...
%                 betainit,options);
%         [best,ee(1,ii),eflag(1,ii)] = fminunc(@(best)objcon(best,blam11,blam22,z,p,zz),...
%                 best*1.0,options);    
        if ii==ti 
            [best betainit]
            dis = sum(sum(abs(best - betainit)));
	    zBtemp = zeros((p - ii)*ii, 1 );
            sma   = 0.001;
            for bi=1:((p-ii) * ii)
	      Btemp = zBtemp;
	      Btemp(bi) = sma;
	      Btemp = reshape(Btemp, p-ii, ii);
                b  = bfunction(best-Btemp,blam11,blam22,z,p,zz);
                b2 = bfunction(best+Btemp,blam11,blam22,z,p,zz);
                deriv(:,bi) = 2*b'*((mean(b2,2)-mean(b,2))/sma/2);
                A(:,bi) = (mean(b2,2)-mean(b,2))/sma/2;
            end
%             Ch=(2*A'*A - hessian/n)./abs(hessian/n)
%             [tt1,tt2] = eig(hessian/n);
%             hessian_inv = tt1*tt2^(-1)*tt1';
%             hessian_inv = inv(hessian/n);
            hessian_inv = inv(2*A'*A);
            COV=hessian_inv*cov(deriv)*hessian_inv/n;
            div=diag(COV).^0.5;
            div1 = COV;
	    size(div)
        end
        betahatU = best(1,:);
        betahatL = best(2:p-ii,:);
        gamma1  = [betahatL-zeros(p-ii-1,1)*betahatU zeros(p-ii-1,1)]; 
        gamma2  = [betahatL-ones(p-ii-1,1)*betahatU ones(p-ii-1,1)];
        b1        = bfunction(gamma1,blam11,blam22,z,p,zz);
        b2        = bfunction(gamma2,blam11,blam22,z,p,zz);     
        [W1a,W1b] = Wb(b1,(p-ii-1)*(ii+1));  
        [W2a,W2b] = Wb(b2,(p-ii-1)*(ii+1)); 
        L1a        = objconV(gamma1,blam11,blam22,z,p,zz,W1a);
        L1b        = objconV(gamma1,blam11,blam22,z,p,zz,W1b);
        L2a        = objconV(gamma2,blam11,blam22,z,p,zz,W2a);
        L2b        = objconV(gamma2,blam11,blam22,z,p,zz,W2b);
        LICa(1,ii)    = L1a+L2a;
        LICb(1,ii)    = L1b+L2b;
        betahat      = [betahat [eye(ii);best]];        
        ii = ii+1;
    end
%%
% The objective function for fmincon
function ee = objcon(best,blam11,blam22,z,p,zz)
n = size(z,1);
b = bfunction(best,blam11,blam22,z,p,zz);
ee = mean(b,2)'*mean(b,2)*n;%'

function b =bfunction(best,blam11,blam22,z,p,zz)
best = [eye(length(best(1,:))); best];
n = size(z,1);
zfest = kernelregression(z*best,z,z*best,'kernel','epanechinikov');
zres = z - zfest;
invsqfesttemp=kernelregression(z*best,zz,z*best,'kernel','epanechinikov');
for ii = 1:p
    for jj = 1:ii
        invsqfest = invsqfesttemp(:,(ii-1)*p+jj);
        blam(ii,jj,:) = invsqfest - zfest(:,ii).*zfest(:,jj);
        blam(jj,ii,:) = blam(ii,jj,:);
    end
end
for ii = 1:n
    blam33(:,:,ii) = blam22(:,:,ii) * blam(:,:,ii);
end
b=[];
for i=1:n
    %a=zres(i,:)'*resinv(i,:);%'
    a=blam11(i,:)'*zres(i,:)-blam33(:,:,i);%'
    b=[b a(:)];
end


function ee = objconV(best,blam11,blam22,z,p,zz,W)
n = size(z,1);
b = bfunction(best,blam11,blam22,z,p,zz);
ee = mean(b,2)'*W*mean(b,2)*n;%'

function [W0a,W0b]=Wb(b,v)
bb=[];
for i=1:length(b(:,1))
    W0a(i,i) = var(b(i,:))^(-1);
    bb=[bb b(i,:)];
end
W0b = var(bb)^(-1)*eye(length(b(:,1)));
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
