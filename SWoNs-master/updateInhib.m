function [dA1,dA2] = updateInhib (a1,a2,dR,B,K,eps,dt)

dA1 =  (- B*a2-K*a1...
                + dR)*dt +  normrnd(dR,eps); %Update for R1
dA2 =   (B*a1 +K*a2...
                 - dR)*dt + normrnd(-dR,eps); %Update for R1

end

