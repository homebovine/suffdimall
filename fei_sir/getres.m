bcover = []
for i = 1:times
  temp = best(3:6, :, i);
temp1 = beta(3:6, :);
normbestm(i, :) = (temp(:) - temp1(:))./sqrt(diag(varbest(:, :, i)/n));
end 

for i = 1:times
  temp = best(3:6, :, i);
bestm(i, :) = temp(:);
end 

Q3 = quantile((normbestm), 0.75);
Q1 = quantile((normbestm), 0.25);
IQ = Q3- Q1;
lower1 = quantile((normbestm(1:times, :)), 0.05);%Q1 -  1.5 * IQ;%%
upper1 = quantile((normbestm(1:times, :)), 1);

lowerb1 = quantile((bestm(1:times, :)), 0.025);%Q1 -  1.5 * IQ;%%
upperb1 = quantile((bestm(1:times, :)), 0.975);

lowerv = diag(quantile(sqrt(varbest/n), 0.025, 3));%Q1 -  1.5 * IQ;%%
upperv = diag(quantile(sqrt(varbest/n), 0.975, 3));%Q1 -  1.5 * IQ;%%
for i = 1 : times
	     revar = reshape(sqrt(diag(varbest(:, :, i) / n)), [p-d, 2]); 
lower = best(3:6, :, i) - 1.96 * revar;
upper = best(3:6, :, i) + 1.96 * revar;
	     %lowerupper(:, 1, i) = lower(:, 1); 
	     %lowerupper(:, 2, i) = upper(:, 1); 
	     %lowerupper(:, 3, i) = lower(:, 2); 
	     %lowerupper(:, 4, i) = upper(:, 2);
cover(:, :, i) = ((beta(3:6, :) <= upper) & (beta(3:6, :) >= lower));
	     end 



for i = 1:times
	
temp1 = beta(3:6, :);
normbest = normbestm(i, :);%(temp(:) - temp1(:))./sqrt(diag(varbest(:, :, i)/n));
incm(i, :)= normbest(:) <= upper1(:) & normbest(:) >= lower1(:); %all(all(normbest > -1 & normbest < 1));
	     end


for i = 1:times
	
temp1 = best(3:6, :, i);
incb(i, :)= temp1(:) <= upperb1(:) & temp1(:) >= lowerb1(:); %all(all(normbest > -1 & normbest < 1));
	     end

for j = 1:8

temp = incm(1:times,j );
temp1= mean(cover(:, :, temp(:)), 3);
temp1 = temp1(:);
mcover(j)= temp1(j);
end
reshape(diag(sqrt(median((varbest(:, :,: )/n), 3))), p-2, 2)
abs(median(best(:, :, :), 3) - beta)
1.4826 * mad(best(:, :, :), 1, 3)
round(mean((cover), 3), 3)
round(reshape((mcover), p-2, 2), 3)
for i = 1:times


temp = sqrt(diag(varbest(:, :, i)/n))
incv(i, :)= temp < upperv(:) & temp > lowerv(:); %all(all(normbest > -1 & normbest < 1));
	     end


for j = 1:8

temp = incv(:,j );
temp1= diag(sqrt(median((varbest(:, :, temp(:) )/n), 3)));
mvarbest(j)= temp1(j);
end

for j = 1:8
temp = incb(:,j );
temp1= 	   	      std(best(:, :, temp(:)), 0, 3);
temp1 = temp1(3:6, :);
temp1 = temp1(:);
mbest(j)= temp1(j);
end

for i = 1:times
		       temp = best(3:6, :, i);
temp1 = beta(3:6, :);
normbest = (temp(:) - temp1(:))./sqrt(diag(varbest(:, :, i)/n));
inc(i)= all(all(normbest < upper(:) & normbest > lower(:))); %all(all(normbest > -1 & normbest < 1));
	     end



diag(sqrt(mean((varbest(:, :,: )/n), 3)))
diag(sqrt(trimmean((varbest(:, :,: )/n),20,  3)))
reshape(diag(sqrt(median((varbest(:, :,: )/n), 3))), p-2, 2)
abs(mean(best(:, :, :), 3) - beta)
1.4826 * mad(best(:, :, :), 1, 3)
1.253 * mad(best(:, :, :), 0, 3)
std(best(:, :, :), 1, 3)
	     
	     lowerupper = zeros(p, 4 , times); 

	     

	     
	     
 for i = 1 : times
	     revar = reshape(sqrt(diag(varbest(:, :, i) / n)), [p-d, 2]); 
lower = best(3:6, :, i) - 1.96 * revar;
upper = best(3:6, :, i) + 1.96 * revar;
	     %lowerupper(:, 1, i) = lower(:, 1); 
	     %lowerupper(:, 2, i) = upper(:, 1); 
	     %lowerupper(:, 3, i) = lower(:, 2); 
	     %lowerupper(:, 4, i) = upper(:, 2);
cover(:, :, i) = ((beta(3:6, :) <= upper) & (beta(3:6, :) >= lower));
	     end 



	     lower = chi2inv(0.025, 36)
	     upper = chi2inv(0.95, 36- 8)
	     for i = 1 : times
		       rebest = reshape(best(3:6, :, i) - beta(3:6, :), (p-2)*2, 1); 
	     invbest = inv(varbest(:, :, i)) * n;
	     chitest = rebest' * invbest * rebest
	     chicover(i) = (chitest <= upper);
	     end 

	     cover = cover(:, :, ix1)
for i = 1: times
minit(:, :, i)= init(:, :, i) * inv(init(1:dim, 1:dim, i)) * beta(1:dim, 1:dim);
end

for i = 1: times
resbest(:, :, i)= best(:, :, i)  * inv(best(1:dim, 1:dim, i));
end

Eudisini = [];
Eudisbest = [];
Group = cell(2000, 1);
for i = 1:times
tini = init(3:6, :, i) - beta(3:6, :);
tbest = best(3:6, :, i) - beta(3:6, :);
Eudisini(i) = tini(:)' * tini(:);
Eudisbest(i)= tbest(:)'*tbest(:);
end
Eudis = [Eudisini; Eudisbest]';
Eudis = Eudis(Eudis(: ,1)<= 2, :);
Group(1:times) = {'DR'};
Group(1001:2000) ={'S-DR'}; 
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'all');
h = boxplot(Eudis, 'labels', {'DR', 'S-DR'}, 'colors', [0, 0, 0], 'symbol', 'b', 'whisker', 1.5);
%set(h(7,:),'Visible','off');
ylim([-0.1 2.2])
set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')
saveas(figure1,'DRnorm.pdf') 

	       tempsd = diag(sqrt(median((varbest(3, 3,: )/n), 3)))
figure1 = figure;
axes1 = axes('Parent',figure1)
hold(axes1,'all');
h = histogram((best(3, 1, :) - beta(3, 1))./sqrt(varbest(3, 3, :)/n));
hold on;
h1 = histogram((randn(times, 1)))
%set(h(7,:),'Visible','off');
saveas(figure1,'hist.pdf')

j = 6; k = 2;
figure1 = figure;
axes1 = axes('Parent',figure1)
hold(axes1,'all');
temp = (best(j, k, :) - beta(j, k))./sqrt(varbest(j * k - 2 * k, j * k - 2 * k, :)/n);
h = qqplot(temp(:));
saveas(figure1,'qq.pdf')
%%%save srd66, DRsrd77 noncenter

allOneString = sprintf('%d,' , otherresult.betahat(:, 1, 2))
allOneString = sprintf('%d,' , otherresult.div1(:, :, 2))
allOneString = allOneString(1:end-1)
load('app_result_3_10.mat')
round([otherresult.betahat(2:end, 1, 1), otherresult.div(:, 1), otherresult.betahat(2:end, 1, 2), otherresult.div(:, 2), otherresult.betahat(2:end, 1, 3), otherresult.div(:, 3), otherresult.betahat(2:end, 1, 4), otherresult.div(:, 4)], 3)



figure1 = figure;
axes1 = axes('Parent',figure1)
hold(axes1,'all');
plot(fity, [cuvx(:, 1), cuvx(:, 2)], cuvy)
saveas(gcf,'qq.pdf')
  
