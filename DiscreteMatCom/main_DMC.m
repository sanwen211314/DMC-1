% The main function for reading the data matrix, obtaining the results and then evaluating them.
% Bipartite matching and Minimum Error Correction(MEC) are used as performance matrix. 

function main_DMC(folder)

    cd ..
    cd simData;
    cd(folder);     
    fprintf('%s :',folder);
    
    if(exist('Data.mat','file'))
        load('Data.mat','X','B_true','r_label');
        [N,V]=size(X);
        S = size(B_true,1);          
    else
        filename='assembly.consensus.fragments.snv.mat.categorized';

        fileID=fopen(filename,'r');     


        l_all=0;
        C_all=[];
        dup_list={};
        C=[];
        l=[];
        S=0;
        r_label=[];

        while ~feof(fileID)

            line=fgetl(fileID);
            line=strsplit(line,'\t');
            read=line{1};
            dup=line{3}(end-4:end);
            if(all(read=='n')); continue; end;

            l_all=l_all+1;
            C_all(l_all,:)=read;


            find_dup=find(strcmp(dup,dup_list));
            if(isempty(find_dup))
                S=S+1;
                dup_list{S}=dup;            
                C{S}(1,:)=read;
                l(S)=1;
                r_label(l_all,:)=[l_all S];
            else
                l(find_dup)=l(find_dup)+1;
                C{find_dup}(l(find_dup),:)=read;
                r_label(l_all,:)=[l_all find_dup];
            end  

        end
        
    
    
        fclose(fileID);
        
        [N,V]=size(C_all); 
        X=zeros(N,V);
        X(C_all=='1')=-1;
        X(C_all=='.')=1;
        X=sparse(X);
        clear C_all;
        
        %% Ground Truth Extraction  

        B_true=zeros(S,V);

        for s=1:S
            buffer=zeros(size(C{s}));
            buffer(C{s}=='1')=1;
            buffer(C{s}=='.')=-1;

            B_true(s,:)=sum(buffer);

        end

        B_true(B_true>0)=1;
        B_true(B_true<0)=-1;
        save('Data.mat','X','B_true','r_label');
    
    end
    
    cd ..
    cd ..
    cd DiscreteMatCom;
   
    
    [indx.j, indx.i]=find(X');
    prev_i=1;
    pos_table=zeros(N,2);
    pos_table(1,1)=indx.j(1);
    for ele=1:length(indx.i)
        if(prev_i==indx.i(ele))
            continue;
        else
            pos_table(prev_i,2)=indx.j(ele-1);
            prev_i=indx.i(ele);
            pos_table(prev_i,1)=indx.j(ele);            
        end
    end
    pos_table(N,2)=indx.j(ele);   
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Running of Discrete Matrix Completion
    
   
    fprintf(', Running DMC :\n');   
    fname=strcat('HapDMC',folder,'.txt');
    Hap_recon_DMC(S,X,folder);
   
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    fprintf('Input of Segments from Output File\n');    
    
    cd HapDMC
    fileID=fopen(fname,'r'); 
    fgetl(fileID);
    Time=fgetl(fileID);
    fgetl(fileID);
    Temp=zeros(V,S);
    for v=1:V
        Temp(v,:)=str2num(fgetl(fileID));
    end  
    
    B_est=Temp';
    B_est=2*B_est-1;
    B_est(B_est==-3)=0;

    fclose(fileID);
    cd ..


    Hamm =@(x,y) sum(x~=y);

    %% Clustering of Reads to Duplications
    cluster_mat=zeros(S,S);    
    
    MEC=0;
    for row=1:N        

        temp=pos_table(row,1):pos_table(row,2);       
        dH_min=Inf;
        min_s=0;
            for s=1:S
                dH=Hamm(X(row,temp),B_est(s,temp));   
                    if(dH<dH_min)
                        dH_min=dH;
                        min_s=s;
                    end
            end  
            MEC=MEC+dH_min;

             cluster_mat(min_s,r_label(row,2))=cluster_mat(min_s,r_label(row,2))+1;

    end    

    fprintf('MEC = %d \n',MEC);
     
    [val,~,~]=bipartite_matching(cluster_mat);
    bip_score=val/sum(cluster_mat(:));



%% Hamming Distance Matrix   

    Hamm_mat=zeros(S,S);

    for s1=1:S
        dH_min=Inf;
        min_s=0;        
        for s2=1:S
            temp_est=B_est(s1,:);
            temp_est(B_true(s2,:)==0)=0;
            Hd=Hamm(temp_est,B_true(s2,:));
            Hamm_mat(s1,s2)=Hd;
            if(Hd<dH_min)                
                dH_min=Hd;
                min_s=s2;
            end
        end
        fprintf('  (%d,%d)  ',dH_min,min_s);
    end
   fprintf('\n'); 

    %% Output to File
    cd ResDMC;
    fout=strcat('ResDMC',folder,'.txt');
    Fmt=[repmat('%d\t\t',1,S),'\n'];
    Fmt_tab=[repmat('%d\t\t\t',1,S),'\n'];
    fileID=fopen(fout,'w');
    fprintf(fileID,'#Duplcations = %d\n',S);
    fprintf(fileID,'N = %d , V = %d\n',N,V);
    fprintf(fileID,'MEC=%d\n',MEC);
    fprintf(fileID,'%s\n',Time);
    fprintf(fileID,'Max. Bipartite Match Score = %.4f',bip_score);

    fprintf(fileID,'\n\nClustering Matrix -->\n');
    fprintf(fileID,Fmt,cluster_mat');

    fprintf(fileID,'\n\nHamming Distance Matrix -->\n');
    fprintf(fileID,Fmt_tab,Hamm_mat');

    fclose(fileID);   
    cd ..


end




