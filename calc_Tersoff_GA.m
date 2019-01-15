function [Vhat]=calc_Tersoff_GA(variables,dat)
R=1.815;
D=0.335;
clust_size=3;
type=[2,2,2];
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

for iQ=1:1:1
    rs(1,2,iQ)=dat(iQ,1);
    rs(2,1,iQ)=dat(iQ,1);
    rs(1,3,iQ)=dat(iQ,2);
    rs(3,1,iQ)=dat(iQ,2);
    rs(2,3,iQ)=dat(iQ,3);
    rs(3,2,iQ)=dat(iQ,3);
end

%looping over no. of data pts. Q
for iQ=1:1:1
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
       
end %iQ
