%
% Copyright (c) 2015, Flavio Prattico
%

function E=lzentropy(rd)

n=length(rd);
L=zeros(1,n);
L(1)=1;

for i=2:n
    
    sub=rd(i);              
    
    match=rd(1:i-1)==sub;   
    
    if all(match==0)==1     
        L(i)=1;
    else                    
        k=1;
        
        while k<i  
            
            if i+k>n      
                L(i)=0;
                break
            end
            
            sub=rd(i:i+k);  
            
            for j=1:i-1      
                
                match=rd(j:j+length(sub)-1)==sub;
                
                if all(match==1)==1
                    break;
                end
            end
            
            L(i)=length(sub);
            if all(match==1)==0
                k=i;
            end
            k=k+1;
        end
    end
end

E=1/(1/n * sum(L))*log(n);

end