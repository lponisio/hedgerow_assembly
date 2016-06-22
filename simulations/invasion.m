function [In indx]=invasion(M,PoA,nInv,k_level,opW)

% Function that performs the introduction of nInv species of plants (PoA=0)
% or animals (PoA=1) to network M. Introduced species can be generalists 
% (k_level=1) or specialists (k_level=0), and they can be linked more likely to
% generalist natives (opW=1), to specialists (opW=2) or randomly to any native
% (opW=0)

[m n]=size(M);
    
% Defining if the introduced species are plants or animals to know 
% the degree of its native counterparts and to define the alien's degree
    
    if PoA==0
        knative=sum(M,2);% degree of natives of the same guild (plants)
        knative2=sum(M);% degree of natives of the other guild (animals)
        InvM=zeros(nInv,n); % preparing the sub-matrix for the alien
        fk=round(length(knative)*0.3);% number of 30% species of the same guild
    elseif PoA==1
        knative=sum(M);% same guild (animals)
        knative2=sum(M,2)';% the other guild (plants)
        InvM=zeros(m,nInv); % preparing the sub-matrix for the alien
        fk=round(length(knative)*0.3);
    end
    
% Defining the degree of the introduced species    

    ascend_knative=sort(knative);
    
    if k_level==0
        k=round(mean(ascend_knative(1:fk)));% degree of alien is equal to
                                            % the mean of 30% most
                                            % specialist natives
    elseif k_level==1
        k=round(mean(ascend_knative(end-(fk-1):end)));% degree of alien is equal
                                            % to the mean of 30% most
                                            % generalist natives
    end
    
% Defining the species to which the introduced one is linked to    
    native_label2=1:length(knative2);
    knative2_label=[native_label2; knative2];    

for i_nInv=1:nInv
    
    if opW==0
         R=randperm(length(knative2));
         Links_Inv(i_nInv,:)=R(1:k);
    
     elseif opW==1
         
         Links_Inv=zeros(nInv,k);
         knative2_labelB=knative2_label;
         
         for l=1:k
             Pm=knative2_labelB(2,:)./sum(knative2_labelB(2,:));% p porportional to the degree, normalized
             Pm2=suma_one(Pm);

             R = mnrnd(1,Pm2');   % rand numbers from multinomial distribution
             fW=find(R);          % find the chosen index ofthe linked native 
            
             Links_Inv(i_nInv,l)=knative2_labelB(1,fW);
             knative2_labelB(:,fW)=[];
                        
         end

    elseif opW==2
        [degree sp_label]=sort(knative2,'descend');
        kmin=min(degree(1:k));
        generalists=sp_label(degree>=kmin);
        perm_generalists=generalists(randperm(length(generalists)));
        Links_Inv=perm_generalists(1:k);
    
    end
            
% Making the introduction
    if PoA==0
        InvM(i_nInv,Links_Inv(i_nInv,:))=1;
    elseif PoA==1
        InvM(Links_Inv(i_nInv,:),i_nInv)=1;
    end

end

if PoA==0
    In=[ InvM ; M ];
elseif PoA==1
    In=[InvM M];
end

indx=1:nInv;
end

    
    

