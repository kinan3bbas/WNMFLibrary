classdef WNMFLibrary
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
   
    
    methods (Static) 
        %%
        function [ G , F ] = EM_WMU_NMF( W ,X , G , F , Iter_max_M_step, Iter_max_E_step )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            %E Step
            tic;
            for i=1:Iter_max_E_step
              WX = W2.*X;
              Xnew=WX+ (I-W).*(G*F); 

              %M Step
              [ G , F ]=MU_NMF( Xnew , G , F , Iter_max_M_step );
              X=Xnew;
            end

            t = toc;
            f = norm( W.*(X-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );

        end
        %%
        function [ G , F ] = EM_WNE_NMF( W ,X , G , F , Iter_max_E_step )
            initf = norm( W.*(X-G*F) , 'fro' )^2;
            W2 = W.^2;
            [n1,n2]=size(W);
            I=ones(n1,n2);
            reducedDim=50;
            tic;
            %E Step
            for i=1:Iter_max_E_step
              WX = W2.*X;
              Xnew=WX+ (I-W).*(G*F); 

              %M Step
              [G , F ,it,ela,HIS]=NeNMF(Xnew,reducedDim,'TYPE','L1R');
              X=Xnew;
            end
            t = toc;
            f = norm( W.*(X-G*F) , 'fro' )^2;
            fprintf('\n### Elapse time: %d sec.\n###   Initial objective value: %d\n###   Objective value: %d \n' , t , initf , f );
        end

    end
end
