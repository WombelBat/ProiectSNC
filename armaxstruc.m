function [V] = armaxstruc(eData,vData,m)
    
    [len,nr] = size(m);

    V=zeros(nr,len+1);
    for i =1:len
        sys_m = armax(eData, m(i,:) );
        p= pe(sys_m,vData);
        mse_val = mean( p.OutputData .^2);
        V(:,i) = [mse_val ;m(i,1:end-1)' ];
    
    end

end