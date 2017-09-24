%% Discrete Stage-wise Alternating Minimization
% Function to implement an EM-type discrete refinement scheme on top of real valued optimization.

function DiscreteStageAltMinMC(S, X, fout)  
    
    Hamm =@(x,y) sum(x~=y);
    e=12.5;   
    tic    
   
    [N,V]=size(X);
    
    A_est=zeros(N,S);  
    
    [A_tilde, ~, B_tilde] = svds(0.75*X,1);
    
    I=2*(sum(A_tilde)>0)-1;
    A_initial=A_tilde*diag(I);
    B_initial=(B_tilde*diag(I))';
    
    MEC_Best=Inf;  
    
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
    Width_table=pos_table(:,2)-pos_table(:,1)+1;
       

    B_est=opt_StageAltMin(A_initial, B_initial,S,X);                            

                        
    %% Post Processing            

    max_iter=50;
    MEC_inner=Inf;
    for outer=1:max_iter
            % Estimating A    
            MEC=0;

            for row=1:N
                temp=pos_table(row,1):pos_table(row,2);
                W=Width_table(row);   
                dH_min=Inf;
                for s=1:S                      
                    dH=Hamm(X(row,temp),B_est(s,temp));
                    if(dH<dH_min)
                        dH_min=dH;                        
                    end
                    A_est(row,s)=(1-e*0.01)^(W-dH)*(e*0.01)^dH;
                end   
                MEC=MEC+dH_min;              

                if(sum(A_est(row,:))==0)
                    A_est(row,:)=zeros(1,S);
                else
                    A_est(row,:)=A_est(row,:)/sum(A_est(row,:));
                end

            end 

            if(MEC_inner<=MEC); break; end

            MEC_inner=MEC;


            % Estimating B

            weight=zeros(S,V);    

            for row=1:N            
                temp=pos_table(row,1):pos_table(row,2);        
                weight(:,temp)= weight(:,temp)+A_est(row,:)'*X(row,temp);          
            end 

            B_est=sign(weight);

    end   
    
    
    Time=double(toc);
    fprintf('MEC Best = %d \n',MEC_Best);
    cd HapStageAlt;
    fname=strcat('HapStageAlt',fout,'.txt');    
    fid=fopen(fname,'w');
    fprintf(fid,'MEC: %d\n',MEC);
    fprintf(fid,'Run Time: %f sec\n', Time);    
    fprintf(fid,'Recovered Haplotype : \n');
    B_est=B_est';
    for m = 1:V
       for n = 1:S
        fprintf(fid, '%d ', B_est(m, n));
       end
    fprintf(fid,'\n');
    end
    fclose(fid);

    cd ..

end
