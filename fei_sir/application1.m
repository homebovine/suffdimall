% cd /gpfs/home/luz15/YanyuanSemiDR/JASA1/Study10
clc;clear;
data = [3	1	3	26	1	1	0	32
1	1	14	38	1	0	0	39.1
1	1	12	35	0	0	0	33.2
2	1	8	40	7	0	0	30.6
3	1	3	28	0	1	0	29
3	1	3	24	0	0	0	30.5
3	1	4	27	0	0	0	30
3	1	8	33	2	1	0	27
1	1	4	62	0	0	0	34
3	1	9	31	0	0	0	29.5
3	1	9	34	2	0	0	26.8
2	1	8	37	8	0	0	31.3
2	1	9	37	0	0	0	31.2
2	1	10	58	6	0	0	34.7
3	1	4	33	0	0	0	30
3	1	3	27	0	0	0	31
3	1	6	30	0	0	0	27
2	1	8	37	9	0	0	29.6
3	1	5	44	6	0	0	32.6
2	1	4	29	3	0	0	29.6
3	1	4	36	2	0	0	29.5
2	1	3	28	3	1	0	31
1	1	5	45	0	0	0	28.5
2	1	3	33	4	1	0	26.7
3	1	3	24	1	1	0	30.75
3	1	3	27	1	1	0	29.5
2	1	16	60	6	0	0	42.2
1	1	13	48	0	0	0	37.6
1	1	12	40	6	0	0	34
2	1	4	33	7	0	0	33
1	1	7	35	4	0	0	28.76
1	1	11	44	0	0	0	35.4
3	1	3	43	8	1	0	31
2	1	18	46	2	0	0	38.8
2	1	14	42	0	0	0	34.3
1	1	19	47	0	0	0	35
3	1	3	25	2	0	1	34.6
2	1	2	30	4	0	0	28.5
1	1	11	40	0	0	0	29.5
3	1	3	26	2	1	0	30.5
3	1	5	32	1	0	0	34.2
1	1	15	51	0	0	0	43.6
5	1	7	35	0	0	1	33.5
3	1	12	37	1	0	0	33
1	1	18	44	0	0	0	45.3
1	1	17	53	3	1	0	38.8
1	1	10	40	0	0	0	29.9
3	1	5	51	10	1	0	31.2
1	1	15	42	0	0	0	34
2	1	2	53	0	0	0	30.45
1	1	3	58	3	1	0	35.5
1	1	4	44	10	0	1	34
2	1	7	31	0	0	0	29.1
1	1	8	64	0	0	0	29.65
3	1	15	47	1	0	0	29.2
3	1	9	37	0	0	1	29.8
2	1	16	46	0	0	0	33.5
1	1	8	55	0	0	0	34
1	1	9	39	0	0	0	29.6
3	1	18	51	0	0	0	34
2	2	3	37	8	0	0	37.25
2	2	6	30	3	1	0	33
3	2	4	26	0	0	0	28.6
5	2	5	41	1	0	1	36
3	2	4	34	4	0	1	37.3
2	2	7	57	4	1	0	29.9
1	2	11	53	8	0	0	31.5
3	2	5	32	4	0	1	41.4
1	2	17	44	5	0	0	32.74
3	2	3	25	1	1	0	33.5
1	2	5	31	9	0	0	32
1	2	9	50	0	0	0	30.8
5	2	3	47	3	0	1	42
3	2	4	35	0	1	0	34
2	2	16	43	0	0	0	32.5
2	2	9	46	10	0	0	31.7
5	2	3	35	0	1	0	36.5
3	2	4	22	0	1	0	33
2	2	8	58	0	0	0	31.2
5	2	8	40	0	0	0	34
3	2	6	30	0	0	0	33
5	2	4	29	4	0	0	33.9
1	2	3	31	9	0	1	39
2	2	12	52	18	1	0	34.92
5	2	3	33	5	1	0	39
1	2	8	49	0	0	0	34
2	2	6	34	7	0	0	31.9
5	2	3	26	1	1	0	37
5	2	4	28	0	1	0	34
5	2	3	35	2	0	0	36.4
1	2	15	47	1	0	1	38.2
1	2	15	51	0	0	0	35.3
3	2	3	26	2	1	0	34.5
3	2	12	33	0	0	0	30.5
4	2	2	27	2	1	0	30
5	2	8	34	0	0	1	37.3
4	2	5	29	0	0	0	40.2
3	2	5	27	0	1	0	35.5
1	2	11	43	0	0	0	35
3	2	4	36	3	0	0	38
1	2	9	38	0	0	0	35.3
2	2	14	60	0	0	0	34.1
3	3	4	43	5	0	1	43.2
2	3	15	48	5	0	0	36.1
5	3	7	32	3	0	0	34.6
3	3	5	31	0	1	0	36
5	3	7	29	2	0	0	36.2
3	3	7	35	0	0	0	37.5
3	3	4	37	12	0	0	41
2	3	10	43	0	0	0	35.6
3	3	5	33	5	0	0	39.8
4	3	11	58	4	0	1	41.3
3	3	9	44	7	0	0	42.5
3	3	4	37	8	0	1	45.8
5	3	5	48	6	0	0	34.9
5	3	4	26	0	1	0	41.5
3	3	5	25	0	0	0	38
4	3	6	38	0	0	0	35
3	3	6	41	0	0	0	40
3	3	5	29	0	1	0	36
2	3	9	59	0	0	0	33.7
2	3	5	29	4	1	0	36.3
3	3	3	27	2	0	1	38
5	3	4	30	0	0	0	39.5
2	3	7	34	5	0	0	36.3
3	3	8	35	2	0	0	32.5
2	3	12	50	6	0	0	37
5	3	3	33	1	0	0	32.6
3	3	4	26	0	0	0	36
5	3	3	36	0	0	0	35
5	3	3	33	5	0	1	43.6
3	3	8	47	0	0	0	33.8
1	3	21	51	0	0	0	35.3
1	3	16	42	6	0	0	42.4
5	3	5	31	0	1	0	39.5
2	3	25	62	10	0	0	43.5
5	3	6	46	1	1	0	42
3	3	21	60	9	0	0	40.3
4	3	6	43	5	1	0	44
1	3	25	53	2	0	0	40.66
3	3	13	38	1	0	0	39.7
5	3	6	39	5	0	0	45
5	3	7	35	0	0	0	43.9
4	3	8	40	3	0	0	38
5	3	5	32	3	0	0	39.02
5	4	5	33	3	1	0	44.5
5	4	4	30	1	1	0	41
5	4	6	37	3	1	0	44
5	4	6	30	0	1	0	44
5	4	5	32	4	0	0	42.5
5	4	7	37	3	0	0	40.26
5	4	5	29	1	1	0	44.5
1	4	13	50	9	0	0	35.5
5	4	6	29	0	1	0	42.5
5	4	7	32	0	0	0	44
5	4	6	31	2	1	0	45
2	4	15	47	4	0	0	44.4
3	4	17	44	0	0	0	38
5	4	4	27	0	1	0	41.8
1	4	23	55	0	1	0	45.5
3	4	5	52	4	1	0	42.5
5	4	3	50	12	0	0	44
3	4	19	59	8	0	1	54.3
3	4	26	47	0	0	0	44.8
3	4	6	43	4	1	0	47
5	4	15	41	0	0	0	43.8
1	4	12	39	4	0	1	48
5	4	9	39	0	0	0	42.7
3	4	14	40	1	0	1	48.5
3	4	16	49	0	0	0	42
2	4	16	53	1	0	0	45.5
3	4	11	37	0	0	0	44.5
2	4	13	40	2	0	0	51.2
5	5	7	34	0	1	0	47.5
5	5	8	37	0	0	0	44.5
5	5	8	31	0	1	0	47
5	5	6	41	10	1	0	47
3	5	17	46	4	0	0	43.1
5	5	8	37	2	1	0	49
5	5	8	33	0	1	0	48.5
3	5	8	35	5	0	0	45
5	5	16	49	5	0	0	52.5
5	5	6	33	2	1	0	47.5
5	5	7	31	0	1	0	48
5	5	8	49	4	1	0	46.5
5	5	12	40	2	0	0	61.5
5	5	9	37	2	0	0	50
5	5	12	46	2	0	0	61.8
4	5	16	43	0	0	0	43
5	5	11	36	1	1	0	47
5	5	9	37	6	0	0	58.5
5	5	16	40	7	1	0	55
3	5	24	54	3	1	0	57
5	5	17	57	1	1	0	57
5	6	14	49	0	1	0	60
3	6	13	41	0	1	0	60
5	6	19	59	4	1	0	59
5	6	12	51	0	1	0	60
5	6	20	45	0	1	0	65
5	6	20	56	1	1	0	52
5	6	22	57	0	1	0	58
4	6	21	53	0	1	0	60
5	6	39	65	0	1	0	74
3	6	34	60	0	1	0	95
5	6	36	61	0	1	0	97
5	6	32	62	0	1	0	88
5	6	35	59	0	1	0	94];

x = data(:,2:7); [n,p] = size(x);
SS = cov(x); 
x0 = (x - ones(n,1) * mean(x)) * (inv(cov(x)))^(1/2);
y0 = data(:,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                         %
%                                                                         %
%                             DIRECTIONAL  REGRESSION                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = x0; y = y0; dim = 4;
ss = cov(x,0); mu = mean(x);
z = (x - ones(n,1)*mu)*(inv(ss))^(1/2);
beta0 = 0.65*ones(p,p); 

for i=1:4
    if i==1
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i)]=select_main_SIR(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n)
    end
    if i==2
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i)]=select_main_SAVE(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    if i==3
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i)]=select_main_DR(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    if i==4
        [LICa,LICb,ee(i,:),eflag(i,:),dis(i,:),betahat(:,:,i),div(:,i)]=select_main_PHD(z,y,dim,beta0,1);
        LIC(i,:)=(LICa+LICb)/2+([1:dim]*p)*log(n);
    end
    i
end
div
otherresult = struct('LIC', LIC, 'ee', ee, 'eflag', eflag, 'dis', dis, 'betahat', betahat);           
save(sprintf('application_result.mat'), 'otherresult'); 



