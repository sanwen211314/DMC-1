% Implementation of Stage-wise optimization Steps. At each iteration, a low rank archetype is recovered. 
% These non-convex algorithms are sensitive to proper initialization. Here, SVD of the orthogonal subspace is used. 
function [B]=opt_StageAltMin(A_initial,B_initial,S,X)

    tol_obj=1e-4;
    tol_normB=1e-2;
    tol_normA=1e-2;
    iter_alt=100;
    max_step=50;

    Index=double(X~=0);
     A=A_initial;
     B=B_initial;

    for k=1:S
 
        fprintf('Rank = %d \n',k);
        P=(A*B-X).*Index;
        opt_new=Inf;   

        if(k>1)    
            tao= A*B - 0.75*P;
            [A_initial, ~, B_initial] = svds(tao,k);
            A=A_initial;
            B=B_initial';
        end


        for alternate=1:iter_alt


            %Alternating Minimization      

            norm_chngB=Inf;      
            step=1;
            while(norm_chngB>tol_normB && step<max_step) 

                %Pojected Gradient Descent on B
                B_old=B;
                grad_B=A'*P;
                w_B=0.25*norm(grad_B,'fro')^2/norm((A*grad_B).*Index,'fro')^2;
                B_tilde=B-w_B*grad_B;     


                B_tilde(B_tilde<-1)=-1;
                B_tilde(B_tilde>1)=1;
                B = B_tilde;


                P=(A*B-X).*Index;
                norm_chngB=norm(B-B_old,'fro');
                step=step+1;      
                
            end

            norm_chngA=Inf;
            step=1;            
            while(norm_chngA>tol_normA  && step<max_step)

                %Pojected Gradient Descent on A
                A_old=A;  
                grad_A=P*B';
                w_A=0.25*norm(grad_A,'fro')^2/norm((grad_A*B).*Index,'fro')^2;
                A_tilde=A-w_A*grad_A;


                A_tilde(A_tilde<0)=0;
                A_tilde(A_tilde>1)=1;
                A=A_tilde;

                P=(A*B-X).*Index;
                norm_chngA=norm(A-A_old,'fro');          
                step=step+1;
                
            end


             opt_old=opt_new;
             opt_new=0.5*norm(P,'fro');             

             if((opt_old-opt_new) < tol_obj);  break;     end
        end  

    end
    
   B=sign(B);
end
