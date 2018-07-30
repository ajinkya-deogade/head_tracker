function [head, tail] = extractHeadTail(contour)

display('Inside Extract Head Tail .........')
dataLength = length(contour);
curVec = nan(1, dataLength);
xOutline = (contour(:,1));
yOutline = (contour(:,2));
distCurv = 100;

parfor i = 1:dataLength
    index1 = i;
    index_2 = index1 + distCurv;
    index_3 = index1 - distCurv;
    
    if index_2 <= dataLength
        point_index_2 = contour(index_2,:);
    elseif index_2 > dataLength
        indexMove = dataLength - index_2;
        indexNo = distCurv - indexMove;
        point_index_2 = contour(indexNo,:);
    end
    
    if index_3 >= 1
        point_index_3 = contour(index_3,:);
    elseif index_3 < 1
        indexNo = dataLength + index_3;
        point_index_3 = contour(indexNo,:);
    end
    
    p1 = contour(i,:);
    p2 = point_index_2;
    p3 = point_index_3;
    
    angle1 = atan2((p2(1,1)-p1(1,1)),(p2(1,2) -p1(1,2)))-atan2((p3(1,1)-p1(1,1)),(p3(1,2) -p1(1,2)));
    curVec(i) = angle1;
end

idxMax0 = find(curVec == max(curVec));
idxMax = idxMax0(1);

curVec2 = curVec;
curVec2(idxMax0) = [0];

index1 = idxMax;

% some points forwards
index2 = index1 + distCurv;
if index2 <= dataLength
    curVec2(index1:index2) = [0];
elseif index2 > dataLength
    indexMove = dataLength - index1;
    indexNo = distCurv - indexMove;
    curVec2(index1:dataLength) = [0];
    curVec2(1:indexNo) = [0];
end

% some points backwards
index3 = index1 - distCurv;
if index3 >= 1
    curVec2(index3:index1) = [0];
elseif index3 < 1
    indexNo = dataLength+index3;
    curVec2(1:index1) = [0];
    curVec2(indexNo:dataLength) = [0];
end
%
% Finally, we find the second maximum
idxMax2=find(curVec2==max(curVec2));
% and pick any, in case there are two or more values together
idxMax2 = idxMax2(1);

head = contour(idxMax,:);
tail = contour(idxMax2,:);
