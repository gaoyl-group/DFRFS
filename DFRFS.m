function [Index] = DFRFS(X,beta,m,c,k,T)%,
%% 初始化
t=1;
[d,n]=size(X);
zr=1e-6;
M=2;

Xe=repmat(X,1,c);

Vs=eye(m,c);

Ws=eye(d,m);

Us=rand(n,c);
UUs=Us(:).^(M/2);
UUs=diag(UUs);

Es=[];
for ii=1:c
    Es = blkdiag(Es,Us(:,ii));
end
Es=Es.^(M/2);

Ps=ones(1,n)./n;
% Ps=rand(1,n);
PPs=repmat(Ps,1,c).^(0.5);
PPs=diag(PPs);

Ds=eye(d);

%% 迭代求解
while t<=T && (t==1 ||  sum(sum((Us-U).^2))>=zr || sum(sum((Vs-V).^2))>=zr || sum(sum((Ws-W).^2))>=zr || sum(sum((Ps-P).^2))>=zr)
% while t<=T
    W=Ws;
    UU=UUs;
    PP=PPs;
    E=Es;
    D=Ds;
    V=Vs;
    U=Us;
    P=Ps;

    Ws=(Xe*UU*(PP*PP')*UU'*Xe'+beta*D)\(Xe*UU*(PP*PP')*E*V');
    Wt=zeros(1,d);
    Ds=zeros(d);
    for ii=1:d
        Wt(ii)=Ws(ii,:)*Ws(ii,:)';
        Ds(ii,ii)=1/(2*sqrt(Wt(ii)));
    end

    DD = SquaredEuDist(V,Ws'*X)+eps;
    DD = DD.^(1/(M-1))+eps;
    Us = 1./DD./sum(1./DD,1);
    Us=Us';
    UUs=Us(:).^(M/2);
    UUs=diag(UUs);
    Es=[];
    for ii=1:c
        Es = blkdiag(Es,Us(:,ii));
    end
    Es=Es.^(M/2);

    G1=sum((Ws'*Xe*UUs-V*Es').^2);
    G2=zeros(c,n);
    for ii=1:c
        G2(ii,:)=G1(n*(ii-1)+1:n*ii);
    end
    G=sum(G2);
    [Gs,~]=sort(G,'ascend');
    G_sum=sum(Gs(1:k));
    for ii=1:n
        Ps(ii)=(Gs(k+1)-G(ii))/(k*Gs(k+1)-G_sum);
    end
    Ps(Ps<0)=0;
    PPs=repmat(Ps,1,c).^(0.5);
    PPs=diag(PPs);

    [Ub,~,Vb]=svd(Es'*(PPs*PPs')*UUs'*Xe'*Ws);%%%
    Imc=eye(size(Vb,1),size(Ub,1));
    Vs=Vb*Imc*Ub';

    t=t+1;

end

%features selection
[~,Index]=sort(Wt,'descend');

end

