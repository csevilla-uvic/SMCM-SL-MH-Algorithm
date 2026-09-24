function C = flow(W)

    w2 = W(2,:); w3 = W(3,:); w4 = W(4,:); w5 = W(5,:);

    % Step 1: Handle node 5 (edges 4->5 and 5->1)
    c7  = max(w5, 0);
    c11 = max(-w5, 0);

    % Step 2: Reduced 4-node subproblem parameter
    T = w2 + w3 + w4 + c7;

    % Step 3: Evaluate objective at 10 candidate (a, b) points
    z  = zeros(size(w2));
    ca = [z;  z;  z;     z;      T;  w2; w2+w3; T-w3; w2; w2  ];
    cb = [z;  w3; T;     w2+w3;  z;  z;  z;     w3;   w3; T-w2];

    % Feasibility
    feas = (ca >= w2) & (ca + cb >= w2 + w3);

    % Objective (set infeasible to inf)
    F = 2.*max(ca, 0) + abs(cb) + max(cb, w3) + abs(T - ca - cb);
    F(~feas) = inf;

    % Select uniformly at random among tied winners
    N   = size(W, 2);
    Fmin = min(F, [], 1);
    tied = (F == Fmin);
    tied_counts = sum(tied, 1);
    r = ceil(rand(1, N) .* tied_counts);
    cs = cumsum(tied, 1);
    chosen = (cs == r) & tied;
    [idx, ~] = find(chosen);
    % In case of exact floating-point ties across >1 row hitting the same
    % cumsum bucket (impossible here since tied_counts is exact), take first.
    [~, first] = unique((1:N), 'first');
    idx = idx(first);

    lin = sub2ind(size(ca), idx(:)', 1:N);
    a   = ca(lin);
    b   = cb(lin);

    % Step 4: Recover all 11 coefficients
    c4_v = max(0, w3 - b);
    c5_v = a - w2 - c4_v;
    c6_v = b + c4_v - w3;
    d    = T - a - b;

    C = [max(a, 0);      % c1  (1->2)
         max(b, 0);      % c2  (1->3)
         max(d, 0);      % c3  (1->4)
         c4_v;            % c4  (2->3)
         c5_v;            % c5  (2->4)
         c6_v;            % c6  (3->4)
         c7;              % c7  (4->5)
         max(-a, 0);      % c8  (2->1)
         max(-b, 0);      % c9  (3->1)
         max(-d, 0);      % c10 (4->1)
         c11];            % c11 (5->1)
end
