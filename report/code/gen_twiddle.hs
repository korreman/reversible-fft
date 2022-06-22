import System.Environment (getArgs)

main :: IO ()
main = do
    args <- getArgs
    case args of
        [log2N, fixpoint] -> do
            let (twiddles_m, twiddles_i) = compute (read log2N) (read fixpoint)
            print twiddles_m
            print twiddles_i
        _ -> putStrLn usage

compute :: Int -> Int -> ([Int], [Int])
compute log2N fixpoint =
    let n = 2 ^ log2N
        angles = [2.0 * pi * fromIntegral x / fromIntegral n | x <- [1..(n `div` 4) - 1]]
        twiddles_r = map (\t -> cos t) angles
        twiddles_i = map (\t -> sin t) angles
        twiddles_m = zipWith (\r i -> -(r - 1.0) / i) twiddles_r twiddles_i
        fix = map (\x -> floor $ x * fromIntegral (2 ^ fixpoint)) in
    (fix twiddles_m, fix twiddles_i)

usage = "usage: gen_twiddle gamma fixpoint\nwhere N = 2^gamma"
