use itertools::Itertools;

const PRIMES_NUM: usize = 169;
const PRIMES: [i32; PRIMES_NUM] = [
    -1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
    271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
    433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
    601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673,
    677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
    769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857,
    859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
    953, 967, 971, 977, 983, 991, 997,
];

fn main() {
    kraitchik::<false>(7);
}

fn trial_division(number: i32) -> [u8; PRIMES_NUM] {
    // This stores the prime (up to 1000) factors
    let mut exponent_vector = [0; PRIMES_NUM];

    // If the number is negative then make it positive but note its
    // first 'prime factor' is -1
    if number < 0 {
        exponent_vector[0] = 1;
    }

    // We know number is positive
    let mut number: u64 = number.abs().try_into().unwrap();

    for (i, p) in PRIMES.into_iter().skip(1).enumerate() {
        // Since we have skiped -1, p is postive
        let p: u64 = p.try_into().unwrap();

        // Try to keep dividing the number until it is not possible
        // to anymore
        loop {
            if number % p != 0 {
                break;
            }
            number /= p;

            // Add 1 to the (i+1)th position
            if let Some(n) = exponent_vector.get_mut(i + 1) {
                *n += 1
            }
        }
    }

    exponent_vector
}

// Calculate x - y, allowing it to be negative
fn sub_u32(a: u32, b: u32) -> i32 {
    if a > b {
        i32::try_from(a - b).unwrap()
    } else {
        -i32::try_from(b - a).unwrap()
    }
}

fn gcd(mut m: i32, mut n: i32) -> i32 {
    while m != 0 {
        let temp = m;
        m = n % m;
        n = temp;
    }
    n.abs()
}

fn kraitchik<const PRINT: bool>(to_be_factorized: u32) -> (i32, i32) {
    // This is Q(x) in Kraitchik's method. Due to avalible number
    // operations we have to manually compute the posiblly negative
    // division. A closure is used instead of a inline function
    // because it can capture its environment (i.e. to_be_factorised)
    let q = |x: u32| {
        let x_p = x.pow(2);
        sub_u32(x_p, to_be_factorized)
    };

    // Compute ⌈sqrt(n)⌉
    let start_point = (to_be_factorized as f64).sqrt().ceil() as u32;

    // These store the current number
    let mut t_up = start_point;
    let mut t_down = start_point - 1;

    // Preinitialise u and v
    let (u_actual, v_actual);

    // This contains computed Q(x) values and their prime
    // factorisations
    let mut bad = vec![];

    loop {
        // Preinitialise variables
        let mut uv_vec = vec![];
        let mut factor_vec = vec![];

        loop {
            let mut iterate = |t: u32| {
                // Compute Q(x) values and its prime factorisation
                let qt = q(t);
                let qt_pf = trial_division(qt);

                factor_vec.push((t, qt_pf));

                // Find a subset where all powers are square
                'exp_loop: for exponents in
                    factor_vec.iter().powerset()
                {
                    let mut base = [0; PRIMES_NUM];
                    // This does colounm-wise summation of the arrays
                    for precise_exp in &exponents {
                        for (index, base_n) in
                            base.iter_mut().enumerate()
                        {
                            *base_n += precise_exp.1[index];
                        }

                        // If we can find a subset where all indices
                        // are multiples of 2, we have found a square
                        if base.iter().all(|b| b % 2 == 0) {
                            uv_vec = exponents
                                .into_iter()
                                .copied()
                                .collect();
                            if bad.contains(&uv_vec) {
                                continue 'exp_loop;
                            }
                            if PRINT {
                                println!(
                                    "base = {:?},\nuv_vec = {:?}",
                                    base, uv_vec
                                );
                            }
                            return true;
                        }
                    }
                }
                false
            };

            if iterate(t_up) || iterate(t_down) {
                // factor_vec contains a subset that is a square
                break;
            }

            t_up += 1;
            t_down -= 1;
        }

        // This calculates u and v by multipling u's together and
        // colounm-wise summation of the v's.
        let (u, v) = uv_vec
            .clone()
            .into_iter()
            .reduce(|(accum_u, mut accum_v), (a, b)| {
                for (index, base_n) in accum_v.iter_mut().enumerate()
                {
                    *base_n += b[index];
                }
                ((accum_u * a) % to_be_factorized, accum_v)
            })
            .unwrap_or((0, [0; PRIMES_NUM]));

        // Convert v to a number
        let v: u32 = v
            .into_iter()
            .zip(PRIMES)
            .skip(1)
            .map::<u32, _>(|(exp, prime)| {
                // Since we skip -1, prime is always positive
                prime.pow(u32::from(exp) / 2).try_into().unwrap()
            })
            .product::<u32>()
            % to_be_factorized;

        // Now u, v are mod n

        // If u mod n != ±v mod n then, we have found u and v
        if u != v && u + v != to_be_factorized {
            u_actual = u;
            v_actual = v;
            break;
        }

        bad.push(uv_vec);
    }

    (
        gcd(
            sub_u32(u_actual, v_actual),
            to_be_factorized.try_into().unwrap(),
        ),
        gcd(
            (u_actual + v_actual).try_into().unwrap(),
            to_be_factorized.try_into().unwrap(),
        ),
    )
}


#[test]
fn tests() {
    let nums = [2, 20, 89, 42, 77, 1000, -9, 9979];
    for i in nums {
        println!("{:?}", trial_division(i))
    }

    let _nums_1 = [481, 793, 949, 1079, 1261, 9913, 9979];
    let nums_2 = [9913];
    for i in nums_2 {
        let (a, b) = kraitchik::<true>(i);
        println!("{:?}", (a, b, a * b))
    }
}
