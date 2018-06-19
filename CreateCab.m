function C = CreateCab(J)

C = sym('C%d%d',length(J)) ;
for i = 1:length(J) 
    for j = 1:length(J)
        if(J(i,j)~=0)
            C(i,j) = J(i,j)/J(i,j) ;
        else
            C(i,j) = 0 ;
        end     
    end 
end

