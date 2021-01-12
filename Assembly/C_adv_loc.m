function [ADV_loc]=C_adv_loc(Grad,dphiq,b_x,b_y,w_2D,nln,BJ,pphys_2D)

% It takes as input also the functions b_x, b_y, and pphys_2D

ADV_loc=sparse(nln,nln);

for i=1:nln
    for j=1:nln
        for k=1:length(w_2D)
            x = pphys_2D(k,1);                         % x quadrature node
            y = pphys_2D(k,2); 
            %x = Region.coord(:,1);
            %y = Region.coord(:,2);
            Binv=inv(BJ(:,:,k));                       % inverse
            Jdet=det(BJ(:,:,k));                       % determinant 
            ADV_loc(i,j)=ADV_loc(i,j) + (Jdet.*w_2D(k)) .*  dphiq(1,k,i) *( [b_x(x,y),b_y(x,y)]*(Grad(k,:,j) * Binv )');
        end
    end
end



                                              
                                              

