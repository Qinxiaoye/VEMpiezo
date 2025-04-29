function GK = globalK(node,elem,d,e,rmt)

sumElem = size(elem,1); % the number of element
sumNode = size(node,1);

elemLen = cellfun('length',elem); 

nnz = sum((3*elemLen).^2);

ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 

ia = 0;
for n = 1:sumElem
    index = elem{n};

    coor = node(index,:);
    AK = elemKVEM(d,e,rmt,coor);
    AB = reshape(AK,1,[]);


    indexDof = [index, index+sumNode, index+2*sumNode];
    Ndof = length(indexDof);

    ii(ia+1:ia+Ndof^2) = repmat(indexDof, Ndof, 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;

end

GK = sparse(ii,jj,ss,sumNode*3,sumNode*3);


