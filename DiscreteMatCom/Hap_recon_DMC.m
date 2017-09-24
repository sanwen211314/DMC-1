%% Iterative Discrete Matrix Completion
% The EM-like framework to refine the estimates from real-valed matrix completion.
% The archetypes and latent factors are discovered by taking consensus from the reads. 

function Hap_recon_DMC(S, X, fout)  
    
    Hamm =@(x,y) sum(x~=y);
    e=12.5;   
    tic    
   
    [N,V]=size(X);
    
    A_est=zeros(N,S);  
    
    [A_tilde, ~, B_tilde] = svds(X,S);
    
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
    
    
    for k=S:-1:2                    
                  
            fprintf('Rank = %d \n',k);
            B_est=opt_alt_min(A_initial, B_initial,k,X);                            
                              
                        
            %% Post Processing            
            
            max_iter=50;
            MEC_inner=Inf;
            for outer=1:max_iter
                    % Estimating A    
                    MEC=0;
                    score=zeros(1,S);                   
                    for row=1:N
                        temp=pos_table(row,1):pos_table(row,2);
                        W=Width_table(row);   
                        dH_min=Inf;
                        for s=1:S                      
                            dH=Hamm(X(row,temp),B_est(s,temp));
                            if(dH<dH_min)
                                dH_min=dH;
                                min_s=s;
                            end                            
                            A_est(row,s)=(W-dH)*log((1-e*0.01))+ dH*log(e*0.01);
                        end   
                        MEC=MEC+dH_min;
                        score(min_s)=score(min_s)+dH_min;


                       % Using Logs to Prevent Underflow
                       shift = max(A_est(row,:));
                       A_est(row,:) = exp(A_est(row,:) - shift);
                       A_est(row,:)=A_est(row,:)/sum(A_est(row,:));
                       

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

            
            if(MEC<MEC_Best)
                    B_star=B_est;
                    MEC_Best=MEC;
            else
                break;
            end
            

          
            if(k==2); break; end          
            


            [~,best]=min(score(1:k));


            swapB=B_est(k,:); 
            B_est(k,:)=B_est(best,:);
            B_est(best,:)=swapB;
            
            swapA=A_est(k,:);
            A_est(k,:)=A_est(best,:);
            A_est(best,:)=swapA;         
           
            A_initial=A_est;
            B_initial=B_est;
                         
    end
    
    
    
    
    B_est=B_star;
    B_est(B_est==0)=-3;
    Time=double(toc);
    fprintf('MEC Best = %d \n',MEC_Best);
    cd HapDMC;
    fname=strcat('HapDMC',fout,'.txt');
    fid=fopen(fname,'w');
    fprintf(fid,'MEC: %d\n',MEC);
    fprintf(fid,'Run Time: %f sec\n', Time);    
    fprintf(fid,'Recovered Haplotype : \n');
    B_est=B_est';
    for m = 1:V
       for n = 1:S
        fprintf(fid, '%d ', (B_est(m, n)+1)/2);
       end
    fprintf(fid,'\n');
    end
    fclose(fid);

    cd ..

end
