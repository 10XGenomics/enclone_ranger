// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Unoptimized functions for packing and unpacking some data structures.

use zstd::bulk::{compress, decompress};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Compression and decompression.  We use zstd rather than gzip because when tested it yielded
// slightly smaller compression size and much lower compression time.  Note that using zstd
// appears to add about 4 MB to the executable size.  If this is really true, it's not obvious
// that it's a good tradeoff.

pub fn compress_bytes(x: &[u8]) -> Vec<u8> {
    compress(x, 0).unwrap()
}

pub fn uncompress_bytes(x: &[u8], uncompressed_size: usize) -> Vec<u8> {
    decompress(x, uncompressed_size).unwrap()
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn u32_bytes(x: usize) -> [u8; 4] {
    (x as u32).to_le_bytes()
}

pub fn u32_from_bytes(x: &[u8]) -> u32 {
    u32::from_le_bytes([x[0], x[1], x[2], x[3]])
}

pub fn f32_bytes(x: usize) -> [u8; 4] {
    (x as f32).to_le_bytes()
}

pub fn f32_from_bytes(x: &[u8]) -> f32 {
    f32::from_le_bytes([x[0], x[1], x[2], x[3]])
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_string(x: &String) -> Vec<u8> {
    let mut bytes = Vec::<u8>::new();
    bytes.extend(u32_bytes(x.len()));
    bytes.extend(x.as_bytes());
    bytes
}

pub fn restore_string(x: &[u8], pos: &mut usize) -> Result<String, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + k > x.len() {
        return Err(());
    }
    let s = String::from_utf8(x[*pos..*pos + k].to_vec());
    if s.is_err() {
        return Err(());
    }
    *pos += k;
    Ok(s.unwrap())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_string(x: &[String]) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 + 5 * x.len());
    bytes.extend(u32_bytes(x.len()));
    for xi in x {
        bytes.extend(u32_bytes(xi.len()));
        bytes.extend(xi.as_bytes());
    }
    bytes
}

pub fn restore_vec_string(x: &[u8], pos: &mut usize) -> Result<Vec<String>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![String::new(); n];
    for yj in &mut y {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + k > x.len() {
            return Err(());
        }
        let s = String::from_utf8(x[*pos..*pos + k].to_vec());
        match s {
            Err(_) => return Err(()),
            Ok(s) => {
                *pos += k;
                *yj = s;
            }
        }
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_string_comp(x: &[String]) -> Vec<u8> {
    let z = save_vec_string(x);
    let mut y = compress_bytes(&z);
    let mut bytes = Vec::<u8>::with_capacity(8 + y.len());
    bytes.extend(u32_bytes(y.len()));
    bytes.extend(u32_bytes(z.len()));
    bytes.append(&mut y);
    bytes
}

pub fn restore_vec_string_comp(x: &[u8], pos: &mut usize) -> Result<Vec<String>, ()> {
    if *pos + 8 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let uncompressed_size = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + n > x.len() {
        return Err(());
    }
    let uncomp = uncompress_bytes(&x[*pos..*pos + n], uncompressed_size);
    *pos += n;
    let mut posx = 0;
    restore_vec_string(&uncomp, &mut posx)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_string(x: &[Vec<String>]) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 + 5 * x.len());
    bytes.extend(u32_bytes(x.len()));
    for xi in x {
        bytes.append(&mut save_vec_string(xi));
    }
    bytes
}

pub fn restore_vec_vec_string(x: &[u8], pos: &mut usize) -> Result<Vec<Vec<String>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<String>::new(); n];
    for yj in &mut y {
        *yj = restore_vec_string(x, pos)?;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_u8(x: &Vec<Vec<u8>>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 + 5 * x.len());
    bytes.extend(u32_bytes(x.len()));
    for vi in x {
        bytes.extend(u32_bytes(vi.len()));
        bytes.extend(vi);
    }
    bytes
}

pub fn restore_vec_vec_u8(x: &[u8], pos: &mut usize) -> Result<Vec<Vec<u8>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<u8>::new(); n];
    for yj in &mut y {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + k > x.len() {
            return Err(());
        }
        *yj = x[*pos..*pos + k].to_vec();
        *pos += k;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_vec_u32(x: &[Vec<u32>]) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 + 8 * x.len());
    bytes.extend(u32_bytes(x.len()));
    for vi in x {
        bytes.extend(u32_bytes(vi.len()));
        for xj in vi {
            bytes.extend(xj.to_le_bytes());
        }
    }
    bytes
}

pub fn restore_vec_vec_u32(x: &[u8], pos: &mut usize) -> Result<Vec<Vec<u32>>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    let mut y = vec![Vec::<u32>::new(); n];
    for yj in &mut y {
        if *pos + 4 > x.len() {
            return Err(());
        }
        let k = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
        *pos += 4;
        if *pos + 4 * k > x.len() {
            return Err(());
        }
        yj.reserve(4 * k);
        for _ in 0..k {
            yj.push(u32_from_bytes(&x[*pos..*pos + 4]));
            *pos += 4;
        }
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_bool(x: &Vec<bool>) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 + x.len());
    bytes.extend(u32_bytes(x.len()));
    for &xi in x {
        bytes.push(if xi { 1 } else { 0 });
    }
    bytes
}

pub fn restore_vec_bool(x: &[u8], pos: &mut usize) -> Result<Vec<bool>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + n > x.len() {
        return Err(());
    }
    let mut y = vec![false; n];
    for (yj, &xj) in y[..n].iter_mut().zip(x[*pos..].iter()) {
        *yj = xj == 1;
        *pos += 1;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_u32(x: &[u32]) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 * (x.len() + 1));
    bytes.extend(u32_bytes(x.len()));
    for &xi in x {
        bytes.extend(xi.to_le_bytes());
    }
    bytes
}

pub fn restore_vec_u32(x: &[u8], pos: &mut usize) -> Result<Vec<u32>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = u32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + 4 * n > x.len() {
        return Err(());
    }
    let mut y = vec![0; n];
    for yj in &mut y {
        *yj = u32_from_bytes(&x[*pos..*pos + 4]);
        *pos += 4;
    }
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_u32(x: u32) -> [u8; 4] {
    x.to_le_bytes()
}

pub fn restore_u32(x: &[u8], pos: &mut usize) -> Result<u32, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let y = u32_from_bytes(&x[*pos..*pos + 4]);
    *pos += 4;
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_bool(x: bool) -> [u8; 1] {
    [x as u8]
}

pub fn restore_bool(x: &[u8], pos: &mut usize) -> Result<bool, ()> {
    if *pos + 1 > x.len() {
        return Err(());
    }
    let y = x[*pos] != 0;
    *pos += 1;
    Ok(y)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn save_vec_f32(x: &[f32]) -> Vec<u8> {
    let mut bytes = Vec::<u8>::with_capacity(4 * (x.len() + 1));
    bytes.extend(f32_bytes(x.len()));
    for &xi in x {
        bytes.extend(xi.to_le_bytes());
    }
    bytes
}

pub fn restore_vec_f32(x: &[u8], pos: &mut usize) -> Result<Vec<f32>, ()> {
    if *pos + 4 > x.len() {
        return Err(());
    }
    let n = f32_from_bytes(&x[*pos..*pos + 4]) as usize;
    *pos += 4;
    if *pos + 4 * n > x.len() {
        return Err(());
    }
    let mut y = vec![0.0; n];
    for yj in &mut y[..n] {
        *yj = f32_from_bytes(&x[*pos..*pos + 4]);
        *pos += 4;
    }
    Ok(y)
}
