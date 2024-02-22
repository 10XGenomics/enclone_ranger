// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Print the ANSI 8-bit colors.

fn main() {
    for i in 0..256 {
        println!("[38;5;{i}m█[0m {i}");
    }
}
