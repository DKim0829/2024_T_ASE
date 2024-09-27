function LUT = Compute_LUT(vec, last, Primes)
    j = length(vec);
    LUT = Primes(last);
    for i = 2:j
        LUT = LUT * Primes(vec(i));
    end
end