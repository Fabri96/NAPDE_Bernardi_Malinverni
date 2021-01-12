mu_vec=10.^-[2:0.2:4];
h=0.0442;
beta=[1 1];
    
for i=1:length(mu_vec)
    
    
    tau_in(i)=1 / sqrt(30*mu_vec(i) /h^4 + norm(beta)^2 / h^2);
    
    [errors,solutions,femregion,Dati]=C_main2D(TestName,nRef,tau_in(i),mu_vec(i));
    
    for j=1:10
        
        tau_new(i)=tau_in(i)+ (-1)^j*0.05;
        [errors_new,solutions,femregion,Dati]=C_main2D(TestName,nRef,tau_new(i),mu_vec(i));
        
        if(errors_new.Error_L2 < errors.Error_L2)
            tau_in(i)=tau_new(i);
        end
        
    end
    
    
  
    
end