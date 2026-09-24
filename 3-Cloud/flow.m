function C = flow(W)
% Input:  W - 4 x N matrix, each column is one problem instance.
%             Components of each column must sum to zero.
% Output: C - 7 x N matrix of non-negative integer coefficients.

    c4 = max(W(4,:), 0);
    c7 = max(-W(4,:), 0);
    S  = W(2,:) + W(3,:) + c4;
    a  = max(W(2,:), min(S, 0));
    b  = S - a;
    C  = [max(a,0); max(b,0); a - W(2,:); c4; max(-a,0); max(-b,0); c7];
end
