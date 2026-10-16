function delta=iteration(newenergy,sigma,l,m,n)
delta=ones(l,m,n);
s=newenergy.^2;
xmax=zeros(1,m);

for i=2:m
    for h=1:i
        xmax(1,h)=1-s(i+1-h)./s(i);
    end
    intvalue0=zeros(l,1);
    if i~=2
        for k=2:i-1
            b=s(i);
            %y=A+B*x;
            A=(s(i)-s(i-k))*(sigma(1:l,i-k+1,n)-sigma(1:l,i-k,n))/(s(i-k+1)-s(i-k))+sigma(1:l,i-k,n);
            B=-s(i)*(sigma(1:l,i-k+1,n)-sigma(1:l,i-k,n))/(s(i-k+1)-s(i-k));
            intvalue0=intvalue0+g1(xmax(k+1),b,A,B)-g1(xmax(k),b,A,B);
        end
    end
    for j=2:n;
        k=1;
        b=s(i);
        A=(s(i)-s(i-k))*(sigma(1:l,i-k+1,j-1)-sigma(1:l,i-k,j-1))/(s(i-k+1)-s(i-k))+sigma(1:l,i-k,j-1);
        B=-s(i)*(sigma(1:l,i-k+1,j-1)-sigma(1:l,i-k,j-1))/(s(i-k+1)-s(i-k));
        intvalue=intvalue0+g1(xmax(k+1),b,A,B)-g1(xmax(k),b,A,B);
        for h=1:l
            delta(h,i,j)=intvalue(h)/sigma(h,i,j-1);
            sigma(h,i,j)=sigma(h,i,1)/delta(h,i,j);
        end
    end
end

end