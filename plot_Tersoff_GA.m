function plot_Tersoff_GA(dat,Q,clust_size,type)
err=0.0;
R=1.815;
D=0.335;
% variables=[1.0288,-113.5513,-1.1884,0.2502, 0.00,0.0,33.0331,10.1063,100000.0085,998.2959,41.9552,0.0006,-122.3017];%-NN Solution
% variables=[3.1151,-252.9782,-0.0667,0.685,0,0,22.4409,7.8711,100998.7203,887.577,45.4815,0.162,-142.7092;];
% % variables=[1.2886,-65.7572,-0.8371,0.5673,0,0,39.8311,13.7491,30438.0581,891.5652,55.77,0.1515,-129.4048;];
% variables=[1.2287,-112.3126,-0.8812,0.3936,0,0,28.1118,11.278,89561.5597,959.3034,44.2672,0.1148,-135.9967;];
% variables=[1.1024,-124.0561,-0.9165,0.3761,0,0,34.1069,10.7918,94299.8204,976.9135,42.2781,0.1038,-125.6333;];
variables=[1.0172,-113.1257,-1.1918,0.2443,0,0,35.4668,9.3547,100384.8445,999.5204,41.2344,0.0081,-125.6201;];
%%O-O-O
param(2,2,2,1)=variables(1);% type_1, type_2, type_3, A
param(2,2,2,2)=variables(2);% type_1, type_2, type_3, B
param(2,2,2,3)=variables(3);% type_1, type_2, type_3, lamda1
param(2,2,2,4)=variables(4);% type_1, type_2, type_3, lamda2
param(2,2,2,5)=0.0;% type_1, type_2, type_3, lamda3
param(2,2,2,6)=0.0;% type_1, type_2, type_3, alpha
param(2,2,2,7)=variables(7);% type_1, type_2, type_3, beta
param(2,2,2,8)=variables(8);% type_1, type_2, type_3, eta
param(2,2,2,9)=variables(9);% type_1, type_2, type_3, c
param(2,2,2,10)=variables(10);% type_1, type_2, type_3, d
param(2,2,2,11)=variables(11);% type_1, type_2, type_3, h
param(2,2,2,12)=R;% type_1, type_2, type_3, R
param(2,2,2,13)=D;% type_1, type_2, type_3, D

for iQ=1:1:Q
    rs(1,2,iQ)=dat(iQ,1);
    rs(2,1,iQ)=dat(iQ,1);
    rs(1,3,iQ)=dat(iQ,2);
    rs(3,1,iQ)=dat(iQ,2);
    rs(2,3,iQ)=dat(iQ,3);
    rs(3,2,iQ)=dat(iQ,3);
end

%looping over no. of data pts. Q
for iQ=1:1:Q
    %converting rs to my format
    for i=1:1:clust_size
         for j=1:1:clust_size
             if i~=j
                 r(i,j)=rs(i,j,iQ);
             end
         end
    end
    
    %calculating Tersoff energy for iQ configuration
    Vhat_iQ=0.0;
    for i=1:1:clust_size
        for j=1:1:clust_size
            if i==j
                continue;
            end
            
            A=variables(12)*r(i,j)+variables(1);
            B=variables(13)*r(i,j)+variables(2);%param(type(i),type(j),type(j),2);
            lambda1=param(type(i),type(j),type(j),3);
            lambda2=param(type(i),type(j),type(j),4);
    %         lambda3=param(type(i),type(j),type(j),5);
    %         alpha=param(type(i),type(j),type(j),6);
            beta=param(type(i),type(j),type(j),7);
            eta=param(type(i),type(j),type(j),8);
    %         c=param(type(i),type(j),type(j),9);
    %         d=param(type(i),type(j),type(j),10);
    %         h=param(type(i),type(j),type(j),11);
            R=param(type(i),type(j),type(j),12);
            D=param(type(i),type(j),type(j),13);        
            [fC]=fc(i,j,r,R,D);
            [fR]=fr(i,j,r,A,lambda1);
            [fB]=fb(i,j,r,B,lambda2);
            zetaij=0.0;
            for k=1:1:clust_size
                if i==j || i==k
                    continue;
                end
    %             A=param(type(i),type(j),type(k),1);
    %             B=param(type(i),type(j),type(k),2);
    %             lambda1=param(type(i),type(j),type(k),3);
    %             lambda2=param(type(i),type(j),type(k),4);
    %             lambda3=param(type(i),type(j),type(k),5);
    %             alpha=param(type(i),type(j),type(k),6);
    %             beta=param(type(i),type(j),type(k),7);
    %             eta=param(type(i),type(j),type(k),8);
                c=param(type(i),type(j),type(k),9);
                d=param(type(i),type(j),type(k),10);
                h=param(type(i),type(j),type(k),11);
                R=param(type(i),type(j),type(k),12);
                D=param(type(i),type(j),type(k),13); 
                [zetaij] = zetaij+zeta(i, j, k, r, c, d, h, R,D);		
            end

            beta_eta=beta^eta;
            zeta_eta=zetaij^eta;
            [bij]=b(zetaij,beta_eta, zeta_eta, eta);
            Vhat_iQ=Vhat_iQ+fC*(fR+bij*fB);
        end
    end
    Vhat(iQ)=0.5*Vhat_iQ;
    err=err+ (Vhat(iQ)-dat(iQ,4))^2.0;
       
end %iQ
for i=1:1:Q
    rplot(i)=rs(2,3,i);
end
err
plot(rplot,Vhat);
hold on
plot(rplot,dat(:,4));