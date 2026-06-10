function A_approx = rangefinder(A, k, sketchType)
    Omega = scegliSketch(sketchType, k, size(A, 1));
    Y = Omega * A;
    [Q, ~] = qr(Y', 0);
    A_approx = (A * Q) * Q';
end