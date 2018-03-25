classdef Sigmoid
    %implements the normalized sigmoid. Pick derivative at beta, zero at
    %zero, and one at one.
    properties 
        alpha
        beta
        gamma
        delta
    end
    properties(Static=true)
        
    end
    
    methods
        function obj = Sigmoid(alpha,beta)
            obj.alpha = alpha;
            obj.beta  = beta;
            obj.gamma  = obj.baseFunction(0);
            obj.delta = obj.baseFunction(1)-obj.gamma;
        end
        
        function sig = baseFunction(obj,x)
            sig = obj.xMinusBeta(x)./sqrt(obj.onePlusXminusBeta2(x));
        end
        
        function xMb = xMinusBeta(obj,x)
            xMb = x-obj.beta;
        end
        
        function oPaXmb2 = onePlusXminusBeta2(obj,x)
            oPaXmb2 = 1+obj.alpha*obj.xMinusBeta(x).^2;
        end
        
        function x = backFunction(obj,sig)
            x = sig./sqrt(1-obj.alpha.*sig.*sig)+obj.beta;
        end
        
        function sig = sigmoid(obj,x)
            sig = (obj.baseFunction(x)-obj.gamma)/obj.delta;
        end
        
        function x = backSigmoid(obj, sig)
            sig(sig>1)=1; 
            sig(sig<0)=0;
            x = obj.backFunction(sig*obj.delta + obj.gamma);
        end
        
        function dSigDx = derivative(obj,x)
            div = (obj.onePlusXminusBeta2(x)-1)./obj.onePlusXminusBeta2(x);
            dSigDx = (1./obj.delta) * (1./sqrt(obj.onePlusXminusBeta2(x))).*(1-div);
        end
        
        
        function weightStatistics(obj,range, a, e, color)
            x = 0:0.01:1;
            s = zeros(size(x));
            m = max(a*(x-0.5).^e+0.5);
            
            %b = zeros(size(x));
%            fit.type = 'polynomial degree 4';
%            fit.coeff =  [-3.8783 7.8107 -6.3799 2.4700 0.7147];
%            f = (x.^4*fit.coeff(1)+x.^3*fit.coeff(2)+x.^2*fit.coeff(3)+x*fit.coeff(2)+fit.coeff(1));
%            f = max(f)-f+1;


            j = 1;
            for i = 1:50000
                obj.beta  = rand(1,1);
                if (m*rand(1,1)>a*(obj.beta-0.5).^e+0.5)
                    continue;
                end 
                obj.alpha = (range(2) - range(1))*rand(1,1)+range(1);
                obj.gamma  = obj.baseFunction(0);
                obj.delta = obj.baseFunction(1)-obj.gamma;
                s = s+obj.derivative(x);                
                j = j+1;
            end
            plot(x,s/j,color);
        end
        
       function test(obj)
            ah = subplot(1,3,1);
            bh = subplot(1,3,2);
            ch = subplot(1,3,3);
            x = 0:0.01:1;
            a = zeros(size(x));
            i = 1;
            while(1)
                s = obj.sigmoid(x);
                plot(ah,x,s,'-');
                hold(ah) 
                obj.backSigmoid(s)
                 plot(ah,x,obj.backSigmoid(s),'r-');
                 hold(ah)
                plot(bh,x,obj.derivative(x));
                a = a+obj.derivative(x);
                plot(ch,x,a/i,'-');
                i = i + 1;
                temp = input('alpha beta ');
                obj.alpha = temp(1);
                obj.beta  = temp(2);
                obj.gamma  = obj.baseFunction(0);
                obj.delta = obj.baseFunction(1)-obj.gamma;
            end
        end
    end
    
    methods(Static=true)
        
         function sigmoid = getSigmoid(randomStream)
             range = [0 10];
             e = 8;
             a = 1000;
             
             x = 0:0.01:1;
             m = max(a*(x-0.5).^e+0.5);
             
             beta1  = randomStream.rand(1,1);
             while (m*randomStream.rand(1,1)>a*(beta1-0.5).^e+0.5)
                    beta1  = randomStream.rand(1,1);
             end 
             alpha1 = (range(2) - range(1))*randomStream.rand(1,1)+range(1);
             sigmoid = Sigmoid(alpha1,beta1);
         end
         
        % function sigmoid = getSigmoid(numOfSigmoids,number)
        %     alpha1 = 100;
        %     
        %     beta1 = (1/numOfSigmoids)*number;
        %     
        %     sigmoid = Sigmoid(alpha1,beta1);
        % end
    end
end

        