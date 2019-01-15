%########## THIS CODE IS EXACTLY SAME AS FILE GA_tersoff.m EXCEPT FOR CALLING GA_tersoff_objfun_RandomlySelectNconfig.m WHICH RANDOMLY SELECTS SOME CONFIGS 
clc();
clear();
Nind=25;
MAXGEN=1;
NVAR=13;%11;
PRECI=50;
GGAP=0.619;

%making A and B as a linear function of r

%solution from neural net
% lowerBound = [1.0288,-113.5513,-1.1884,0.2502, 0.00, 0.0,33.0331,10.1063,100000.0085,998.2959,41.9552,0.0006,-122.3017];
% upperBound = [1.0288,-113.5513,-1.1884,0.2502, 0.00, 0.0,33.0331,10.1063,100000.0085,998.2959,41.9552,0.0006,-122.3017];
%Trial 1
% lowerBound = [  -1e0,-0.5e3,-0.5e1,  -1e0, 0.00, 0.0,  1e1,      0,1e4,1e2,-1e0,-1e0,-0.5e3];
% upperBound = [ 0.5e1,   1e0,   1e0, 0.1e1, 0.00, 0.0,0.5e2, 0.25e2,1e6,1e3, 1e2, 1e1,   1e0];
%Trial 2
% lowerBound = [     0,-252.9782,  -0.5e1,     0, 0.00, 0.0,22.4409, 7.8711,        1e4,887.577,45.4815,  0e0,-142.7092];
% upperBound = [3.1151,      1e0, -0.0667, 0.685, 0.00, 0.0,  0.5e2, 0.25e2,100998.7203,    1e3,    1e2,0.162,      1e0];

% lowerBound = [     0,-252.9782,  -0.5e1,      0, 0.00, 0.0,22.4409, 7.8711, 30438.0581,891.5652,    40,   0e0,-142.7092];
% upperBound = [1.2886, -65.7572, -0.8371, 0.5673, 0.00, 0.0,39.8311,13.7491,100998.7203,    1e3,45.4815,0.1515,-129.4048];

lowerBound = [0.9835, -114.153,    -1.2, 0.243, 0.00,  0.0, 28.1118, 8.9059, 97332.0533,996.1527,    40,    0e0,-128.0627];
upperBound = [1.2124,-112.3126, -1.0258, 0.2597, 0.00, 0.0,35.6218, 10.4845,100998.7203,     1e3,42.7962,0.0247,-120.4048];

% [1.016,-112.7209,-1.1801,0.2575,0,0,35.0602,9.4722,98714.8491,999.9641,41.5017,0.0198,-125.9528;]



FieldD = [rep([PRECI],[1,NVAR]); lowerBound; upperBound; rep([1;0;1;1],[1,NVAR])];

% FieldD = [rep([PRECI],[1,NVAR]); rep([-512;512],[1,NVAR]); rep([1;0;1;1],[1,NVAR])];
Chrom = crtbp(Nind,NVAR*PRECI);
gen=0;
%%%ABOVE is GA Initialization


clust = 3;%GA
dat= xlsread('O3_UMP4_form_corrected.xls');%GA
%dat= xlsread('Co2_MP4SDQ_ConfigEnergy_sorted.xls');
%dat= xlsread('test_C_C.xls');%GA
temp=size(dat);%GA
totalConfig=temp(1);%GA

target = dat(:,4);
sumSQerr = zeros(MAXGEN,Nind);
[ObjV,sumSQerr] = GA_tersoff_objfun_RandomlySelectNconfig_NoGLOBALvariable(bs2rv(Chrom,FieldD),target,Nind,totalConfig,dat,MAXGEN,1,sumSQerr);
gen=1;
while gen < MAXGEN
    FitnV = ranking(ObjV);
    SelCh = select('sus',Chrom,FitnV,GGAP);
%     SelCh = recombin('xovdp',SelCh,exp(-gen/MAXGEN));
    
    SelCh = recombin('xovsp',SelCh,0.9);
    SelCh=mut(SelCh,0.8);
    
    [ObjVSel,sumSQerr] = GA_tersoff_objfun_RandomlySelectNconfig_NoGLOBALvariable(bs2rv(SelCh,FieldD),target,Nind,totalConfig,dat,MAXGEN,gen,sumSQerr);
    
    [Chrom ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
    
    gen = gen+1
    min(ObjV)
end

phen=bs2rv(Chrom,FieldD);


