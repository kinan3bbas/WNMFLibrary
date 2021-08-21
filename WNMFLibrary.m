classdef WNMFLibrary
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
   
    
    methods (Static) 
        %%
        function [X , G , F , t , initf , f] = EM_WMU_NMF( W ,X , G , F ,Iter_max_M_step , Iter_max_E_step , Nesterov_Max_Iterations , Nesterov_Min_Iterations )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            %E Step
            tic;
            for i=1:Iter_max_E_step
              WX = W2.*X;
              Xnew=WX+ (I-W).*(G*F); 
              for j=1: Iter_max_M_step
              %M Step
              [ G , F ]=MU_NMF( Xnew , G , F , Iter_max_M_step );
              end
            end

            t = toc;
            f = norm( W.*(X-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );

        end
        %%
        function [X , G, F , t , initf , f] = EM_WNE_NMF( W ,X , G , F ,  Iter_max_M_step , Iter_max_E_step , Nesterov_Max_Iterations , Nesterov_Min_Iterations  )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            WX = W2.*X;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            reducedDim=size(G,2);
            tic;
            %E Step
            for i=1:Iter_max_E_step 
              Result=WX+ (I-W).*(G*F); 
                for j=1:Iter_max_M_step
                  %M Step
                  [G , F ,it,ela,HIS]=NeNMF(Result,reducedDim,'W_INIT',G,'H_INIT',F,'MAX_ITER',Nesterov_Max_Iterations,'MIN_ITER',Nesterov_Min_Iterations);
                end
            end
            t = toc;
            f = norm( W.*(Result-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );
            
        end

            %%
        function [X , G, F , t , initf , f] = EM_L1R_WNE_NMF( W ,X , G , F , Iter_max_M_step , Iter_max_E_step , Nesterov_Max_Iterations , Nesterov_Min_Iterations  )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            reducedDim=size(G,2);
            tic;
            %E Step
            for i=1:Iter_max_E_step
              WX = W2.*X;
              Xnew=WX+ (I-W).*(G*F); 
               
              %M Step
               for j=1:Iter_max_M_step
                [G , F ,it,ela,HIS]=NeNMF(Xnew,reducedDim,'TYPE','L1R','W_INIT',G,'H_INIT',F,'MAX_ITER',Nesterov_Max_Iterations,'MIN_ITER',Nesterov_Min_Iterations);
               end
            end
            t = toc;
            f = norm( W.*(X-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );
        end

          %%
        function [X , G, F , t , initf , f] = EM_L2R_WNE_NMF( W ,X , G , F  ,Iter_max_M_step , Iter_max_E_step , Nesterov_Max_Iterations , Nesterov_Min_Iterations )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            reducedDim=size(G,2);
            tic;
            %E Step
            for i=1:Iter_max_E_step
              WX = W2.*X;
              Xnew=WX+ (I-W).*(G*F); 
              for j=1:Iter_max_M_step
              %M Step
              [G , F ,it,ela,HIS]=NeNMF(Xnew,reducedDim,'TYPE','L2R','W_INIT',G,'H_INIT',F,'MAX_ITER',Nesterov_Max_Iterations,'MIN_ITER',Nesterov_Min_Iterations);
              end 
            end
            t = toc;
            f = norm( W.*(X-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );
        end

    %%
        function [X , G, F , t , initf , f] = EM_GWNE_NMF( W ,X , G , F , Iter_max_M_step , Iter_max_E_step , Nesterov_Max_Iterations , Nesterov_Min_Iterations  )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            reducedDim=size(G,2);
            tic;
            S=constructW(X);
            %E Step
            for i=1:Iter_max_E_step
              WX = W2.*X;
              Xnew=WX+ (I-W).*(G*F); 
              for j=1:Iter_max_M_step
              %M Step
              [G , F ,it,ela,HIS]=NeNMF(Xnew,reducedDim,'TYPE','MR','S_MTX',S,'W_INIT',G,'H_INIT',F,'MAX_ITER',Nesterov_Max_Iterations,'MIN_ITER',Nesterov_Min_Iterations);
              end
            end
            t = toc;
            f = norm( W.*(X-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );
        end
    end
end

