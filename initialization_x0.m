function x0 = initialization_x0(q0, eta0, qs0_c, width1, width2, slope1, slope2)
    total_weight = width1 * slope1 + width2 * slope2;
    w1 = (width1 * slope1) / total_weight;
    w2 = (width2 * slope2) / total_weight;

    q1_init = q0 * w1;
    q2_init = q0 * w2;

    qs1_c_init = qs0_c * w1;
    qs2_c_init = qs0_c * w2;

    eta1_init = eta0 - 0.2;
    eta2_init = eta0 - 0.2;

    x0 = [q1_init, q2_init, eta1_init, eta2_init, qs1_c_init, qs2_c_init];
end