function [L,lamF,lamB]=getlambda(E,L,dt,dstore,xn,embd_lib,embd_test,npt)

N        = size(embd_test,1);      % num points in the test attractor
DF_store = zeros(length(L)+1,npt); % store forward separation distances
DB_store = zeros(length(L)+1,npt); % store backward separation distances

[dsort I]=sort(dstore,'ascend'); % sort the nearest neighbors (NN) by phase space distance
                                 % if testing less points than are available i.e. npt < n_test, 
                                 % give preference first to points with closest near neighbors

% flag points on test attractor that are okay to use 
% i.e. must have length(L) to calculate distance spreading 
flag = zeros(length(I),1);
ind=1;
for ii=1:length(I) 
    if I(ii)+length(L)<size(embd_test,1) & I(ii)-length(L)>0 &...
            xn(I(ii))+length(L)<size(embd_lib,1) & xn(I(ii))-length(L)>0
        if ind<=npt
            flag(I(ii))=1;
            ind=ind+1;
        end
        
    end
end


ind = 1;
for k=1:length(I) % march through test attractor and measure separation distance forwards and backwards
    if flag(k)==1
        dsqF=(embd_test(k:k+length(L),:)-embd_lib(xn(k):xn(k)+length(L),:)).^2; % square of distances forward and backward from point k (along each dimension)    
        dsqB=(embd_test(k:-1:k-length(L),:)-embd_lib(xn(k):-1:xn(k)-length(L),:)).^2; 
        DF=sqrt(sum(dsqF,2))';    % total euclidean distance forward time 
        DB=sqrt(sum(dsqB,2))';    % total euclidean distance backward time 
        DF_store(:,ind)=DF/DF(1); % normalize by initial separation distance
        DB_store(:,ind)=DB/DB(1); 
        ind=ind+1;
    end
end

% if npt points could not be queried, then remove zeros from DB_store and DF_store, 
f_zero = find(DB_store(1,:)==0);
if numel(f_zero)>0
    f_zero=f_zero(1);
    DB_store(:,f_zero:end)=[];
    DF_store(:,f_zero:end)=[];
end

% get lambda^+ and lambda^- as in Equation 2 in manuscript
lamF=(1./L).*log(mean(DF_store(2:end,:),2)); 
lamB=(1./L).*log(mean(DB_store(2:end,:),2));
