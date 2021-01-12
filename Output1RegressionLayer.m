classdef Output1RegressionLayer < nnet.layer.RegressionLayer
    % Example custom regression layer with mean-absolute-error loss
    
    methods
        function layer = Output1RegressionLayer(name)
            % layer = maeRegressionLayer(name) creates a
            % mean-absolute-error regression layer and specifies the layer
            % name.
			
            % Set layer name.
            layer.Name = name;
           

            % Set layer description.
            layer.Description = 'Mean square error';
        end
        	
        
        function loss = forwardLoss(layer, Y, T)
            % loss = forwardLoss(layer, Y, T) returns the loss between
            % the predictions Y and the training targets T.
            
           global x;
           setGlobalx(x+1);
           
           inputs = linspace(1e-4, 1e-2,50);

            M=Y;
            for i=1:size(Y,2)
                
                mu = inputs(int8(T(i)*10^6-100));
                tau = double(extractdata(Y(i)));
                tau = abs(tau);
                
                if(x<=5)
                    tau=0;
                end
                 
                [errors,~,~,~]=C_main2D('Test7',3,tau,mu);
                M(i)=dlarray(errors.Time_error);
                
                
            end
            loss= mse(M,T,'DataFormat','T')/size(Y,2);            
           
            
        end
    end
end
