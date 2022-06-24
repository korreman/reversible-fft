main = do
    let sines = [round $ sin (x/2**4 * pi * 2) * 2**16 | x <- [0..255]]
    print sines
