function [X_test] = calO(Data)
    X = Data.X{1};
    Y = Data.X{2};
    Z = Data.X{3};

    RID_Y = Y(:,1);
    RID_Z = Z(:,1);

    [common_RID, idx_X, idx_Y] = intersect(RID_Z, RID_Y, 'stable');

   
    a = zeros(size(X, 1), 1);
    a(idx_X) = 1;
    O_X = diag(a);

    b = zeros(size(Y, 1), 1);
    b(idx_Y) = 1;
    O_Y = diag(b);

   
    X1 = O_X * X;
    idx_X = find(X1(:, 1) ~= 0);
    X2 = X1(idx_X, :); 
    X = X2(:, 2:end);   

    Z1 = O_X * Z;
    idx_Z = find(Z1(:, 1) ~= 0);
    Z2 = Z1(idx_Z, :);  
    Z = Z2(:, 2:end);    

    Y1 = O_Y * Y;
    idx_Y = find(Y1(:, 1) ~= 0);
    Y2 = Y1(idx_Y, :);  
    Y = Y2(:, 2:end);   
    
    X_test = {X, Y, Z}; 


end