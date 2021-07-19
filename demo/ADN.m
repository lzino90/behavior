function A = ADN(N,m,alpha)
%TVn(N,m,alpha) constructs a Time-varying network with N nodes, activation probabilities of the
%nodes given by vector alpha. Each node construcs m links.
A=zeros(N);
for i=1:N
    if rand<alpha(i)
        temp=randperm(N);
        for j=1:m
            if i==temp(j)
                A(i,temp(N))=1;
                A(temp(N),i)=1;
            else
                A(i,temp(j))=1;
                A(temp(j),i)=1;
            end
        end
    end       
end
end